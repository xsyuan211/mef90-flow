// ---------------------------------------------------------------------------
// PhaseFieldFlow.h — Phase Field Fracture + Darcy + Biot Coupling
//
// deal.II port of mef90 vDef + m_Flow (Darcy-Biot poroelastic coupling)
//
// Design:
//   - Staggered AltMin:  StepV(damage) → StepP(pressure) → StepU(displacement)
//   - Macaulay strain decomposition (Ambati et al. hybrid model)
//   - Fixed-stress Biot iteration (IPACS pattern, Wheeler+Wick+Lee CMAME 2020)
//   - Permeability degradation k(v) = k_min + (k_max-k_min)*v_cellAvg
//   - History variable H_old = max(H_current, H_previous) for irreversibility
//   - Adaptive AMR + parallel (p4est)
//
// Reference implementations:
//   - deal.II Gallery Phase_field_fracture_model_in_3D (Wasim Munshi 2024)
//   - RACCOON MOOSE phase field framework (https://hugary1995.github.io/raccoon)
//   - IPACS fixed-stress Biot (Wheeler, Wick, Lee, CMAME 2020)
// ---------------------------------------------------------------------------

#ifndef PHASE_FIELD_FLOW_H
#define PHASE_FIELD_FLOW_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_point_data.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_newton.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/parameter/handler.h>

#include <mpi.h>
#include <cmath>
#include <iostream>
#include <fstream>

// deal.II 9.7+ 用 Trilinos 或 PETSc 后端
#define FORCE_USE_OF_TRILINOS
namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
    !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error "deal.II must be built with PETSC or TRILINOS"
#endif
} // namespace LA

namespace PhaseFieldFlow
{
  using namespace dealii;

  /**
   * 历史变量 — Macaulay 应变能量密度 H_old
   * deal.II 原生 CellDataStorage 在 AMR 后自动 pack/unpack
   */
  class HistoryData : public TransferableQuadraturePointData
  {
  public:
    HistoryData() = default;
    virtual ~HistoryData() = default;

    unsigned int number_of_values() const override { return 1; }

    virtual void pack_values(std::vector<double> &scalars) const override
    {
      Assert(scalars.size() == 1, ExcInternalError());
      scalars[0] = H_old;
    }

    virtual void unpack_values(const std::vector<double> &scalars) override
    {
      Assert(scalars.size() == 1, ExcInternalError());
      H_old = scalars[0];
    }

    double H_old = 0.0;
  };


  /**
   * 主类 — 三个子问题 + AltMin 循环 + Biot 耦合
   *
   * 对照 vDef F90 模块:
   *   PhaseFieldFlow      ≈ vDef.F90 主程序 + AltMin 控制
   *   elastic_subproblem  ≈ m_DefMech (力学平衡)
   *   damage_subproblem   ≈ m_DefMech (相场演化)
   *   flow_subproblem     ≈ m_Flow (Darcy 压力)
   *   biot_coupling       ≈ m_Flow::FlowApplyBiotForce
   */
  template <int dim>
  class PhaseFieldFlow
  {
  public:
    PhaseFieldFlow();
    void run();

  private:
    // ===== 1. 参数 =====
    void declare_parameters();
    void parse_parameters(const std::string &filename);

    // ===== 2. 网格 + BC =====
    void setup_mesh();
    void setup_boundary_conditions();

    // ===== 3. 系统分配 =====
    void setup_system_elastic();
    void setup_system_damage();
    void setup_system_flow();

    // ===== 4. 组装 =====
    void assemble_system_elastic();
    void assemble_system_damage();
    void assemble_system_flow();

    // ===== 5. 求解 =====
    void solve_elastic_subproblem();
    void solve_damage_subproblem();
    void solve_flow_subproblem();

    // ===== 6. Biot 耦合 =====
    void compute_biot_body_force();          // f_p = -α·∇p
    void apply_biot_force_to_elastic_rhs();  // u-方程 RHS += f_p
    void update_permeability();              // k(v) = k_min + (k_max-k_min)*v_cellAvg

    // ===== 7. AltMin 控制 =====
    bool check_altmin_convergence();
    void update_history_field();             // H_old = max(H_current, H_old)

    // ===== 8. AMR + 输出 =====
    void refine_grid();
    void output_results(const unsigned int load_step) const;

    // ===== 9. 辅助 =====
    double H_plus(const SymmetricTensor<2, dim> &strain) const;  // Macaulay 正部
    double lambda() const { return (E * nu) / ((1 + nu) * (1 - 2 * nu)); }
    double mu()     const { return E / (2 * (1 + nu)); }
    double gc_at(const Point<dim> &p) const;                      // 可能的 Gc 空间变化

    // ===== MPI + Timer =====
    MPI_Comm mpi_communicator;
    ConditionalOStream pcout;
    TimerOutput computing_timer;

    // ===== 并行网格 =====
    parallel::distributed::Triangulation<dim> triangulation;

    // ===== 弹性子问题 (u) =====
    const FESystem<dim>      fe_elastic;       // 3 (x,y,z) 位移分量
    DoFHandler<dim>          dof_handler_elastic;
    IndexSet                 locally_owned_dofs_elastic;
    IndexSet                 locally_relevant_dofs_elastic;
    AffineConstraints<double> constraints_elastic;
    LA::MPI::SparseMatrix    system_matrix_elastic;
    LA::MPI::Vector          system_rhs_elastic;
    LA::MPI::Vector          locally_relevant_solution_elastic;
    LA::MPI::Vector          completely_distributed_solution_elastic_old;
    LA::MPI::Vector          completely_distributed_solution_elastic;
    const QGauss<dim>        quadrature_formula_elastic;

    // ===== 损伤子问题 (d) =====
    const FE_Q<dim>          fe_damage;          // 1 个标量（相场损伤变量）
    DoFHandler<dim>          dof_handler_damage;
    IndexSet                 locally_owned_dofs_damage;
    IndexSet                 locally_relevant_dofs_damage;
    AffineConstraints<double> constraints_damage;
    LA::MPI::SparseMatrix    system_matrix_damage;
    LA::MPI::Vector          system_rhs_damage;
    LA::MPI::Vector          locally_relevant_solution_damage;
    LA::MPI::Vector          completely_distributed_solution_damage_old;
    LA::MPI::Vector          completely_distributed_solution_damage;
    const QGauss<dim>        quadrature_formula_damage;

    // ===== 流动子问题 (p) =====  ← 新增！对应 m_Flow
    const FE_Q<dim>          fe_flow;            // 1 个标量（压力）
    DoFHandler<dim>          dof_handler_flow;
    IndexSet                 locally_owned_dofs_flow;
    IndexSet                 locally_relevant_dofs_flow;
    AffineConstraints<double> constraints_flow;
    LA::MPI::SparseMatrix    system_matrix_flow;   // Darcy 扩散 + 压缩项
    LA::MPI::Vector          system_rhs_flow;       // 源项 + 边界压力
    LA::MPI::Vector          locally_relevant_solution_flow;
    LA::MPI::Vector          completely_distributed_solution_flow;
    const QGauss<dim>        quadrature_formula_flow;

    // ===== 历史变量 (AMR 后自动传递) =====
    CellDataStorage<typename Triangulation<dim>::cell_iterator, HistoryData>
      quadrature_point_history;

    // ===== Biot 耦合存储 =====
    LA::MPI::Vector          biot_body_force;      // f_p = -α·∇p (cell vector 组装)

    // ===== 物理参数 =====
    double E;              // 杨氏模量 [Pa]
    double nu;             // 泊松比
    double GC;             // 能量释放率 [N/m]
    double length_scale;   // 相场长度尺度 l [m]

    // ===== Darcy/Biot 参数 =====
    double k_min;          // 最小渗透率 [m²] (完整损伤区)
    double k_max;          // 最大渗透率 [m²] (完整介质区)
    double viscosity;      // 流体粘度 [Pa·s]
    double storativity;    // 储水系数 [Pa⁻¹]
    double biot_coeff;     // Biot-Willis 系数 α
    double source_term;    // 体积源项 [1/s]

    // ===== 求解控制 =====
    unsigned int num_load_steps;
    unsigned int max_altmin_iters;
    double       altmin_tol;
    double       load_increment;

    // ===== 参数文件 =====
    ParameterHandler prm;
  };

} // namespace PhaseFieldFlow

#endif // PHASE_FIELD_FLOW_H
