// phase_field_flow.cc — PhaseFieldFlow 实现
//
// M1 里程碑：骨架可编译 + AltMin 循环结构完整
// M2 补 flow_assembly
// M3 补 biot_coupling + Terzaghi 验证

#include "phase_field_flow.h"

namespace PhaseFieldFlow
{

// ========================================================================
// 构造函数 — 对应 vDef.F90 的 Init 链
// ========================================================================
template <int dim>
PhaseFieldFlow<dim>::PhaseFieldFlow()
  : mpi_communicator(MPI_COMM_WORLD),
    pcout(std::cout,
          (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    computing_timer(mpi_communicator, pcout, TimerOutput::never,
                    TimerOutput::wall_times),
    triangulation(mpi_communicator),

    // 三个子问题的有限元
    fe_elastic(FE_Q<dim>(1), dim),                 // u_x, u_y, u_z
    dof_handler_elastic(triangulation),
    quadrature_formula_elastic(fe_elastic.degree + 1),

    fe_damage(1),                                   // 标量损伤 d
    dof_handler_damage(triangulation),
    quadrature_formula_damage(fe_damage.degree + 1),

    fe_flow(1),                                     // 标量压力 p
    dof_handler_flow(triangulation),
    quadrature_formula_flow(fe_flow.degree + 1),

    // 物理参数 (Terzaghi 1D 默认值 — parameters.prm 覆盖)
    E(1.0e9),
    nu(0.25),
    GC(270.0),
    length_scale(0.15),

    k_min(5.0e-11),
    k_max(1.0e-14),
    viscosity(1.0e-3),
    storativity(4.6e-10),
    biot_coeff(1.0),
    source_term(0.0),

    num_load_steps(20),
    max_altmin_iters(30),
    altmin_tol(1.0e-6),
    load_increment(0.5)
{
}

// ========================================================================
// run() — 主循环 (对应 vDef.F90 L200-500)
// ========================================================================
template <int dim>
void PhaseFieldFlow<dim>::run()
{
  pcout << "==================================================" << std::endl;
  pcout << "  PhaseFieldFlow — Phase Field + Darcy + Biot"       << std::endl;
  pcout << "  deal.II port of vDef/m_Flow (PETSc/F90)"          << std::endl;
  pcout << "==================================================" << std::endl;

  // ---- 1. 参数 ----
  pcout << "\n[1/7] Parsing parameters.prm ..." << std::endl;
  declare_parameters();
  parse_parameters("parameters.prm");

  // ---- 2. 网格 ----
  pcout << "\n[2/7] Setting up mesh ..." << std::endl;
  setup_mesh();

  // ---- 3. 系统 ----
  pcout << "\n[3/7] Setting up systems ..." << std::endl;
  setup_system_elastic();
  setup_system_damage();
  setup_system_flow();

  // ---- 4. 边界条件 ----
  pcout << "\n[4/7] Setting up boundary conditions ..." << std::endl;
  setup_boundary_conditions();

  // ---- 5. AltMin 时间循环 (对应 vDef.F90 L300-500) ----
  pcout << "\n[5/7] AltMin time loop (" << num_load_steps << " steps) ...\n"
        << std::endl;

  for (unsigned int load_step = 1; load_step <= num_load_steps; ++load_step)
  {
    pcout << "--- Load step " << load_step << "/" << num_load_steps << " ---"
          << std::endl;

    // ---- AltMin 迭代 ----
    bool altmin_converged = false;
    for (unsigned int altmin_iter = 1;
         altmin_iter <= max_altmin_iters && !altmin_converged;
         ++altmin_iter)
    {
      pcout << "  AltMin iter " << altmin_iter << "/" << max_altmin_iters
            << std::endl;

      // ---- Step V: 相场损伤 (对应当前 AltMin 的 d-步) ----
      pcout << "    Step V: solving damage ..." << std::endl;
      assemble_system_damage();
      solve_damage_subproblem();

      // ---- 历史变量更新 ----
      update_history_field();

      // ---- Step P: Darcy 压力 (m_Flow 新增！) ----
      pcout << "    Step P: solving Darcy flow ..." << std::endl;
      update_permeability();
      assemble_system_flow();
      solve_flow_subproblem();

      // ---- Biot 耦合: f_p = -α·∇p → 加到 u-方程 RHS ----
      if (std::abs(biot_coeff) > 1e-12)
      {
        pcout << "    Biot: applying body force ..." << std::endl;
        compute_biot_body_force();
        apply_biot_force_to_elastic_rhs();
      }

      // ---- Step U: 位移 (含 Biot 体积力) ----
      pcout << "    Step U: solving elastic ..." << std::endl;
      assemble_system_elastic();
      solve_elastic_subproblem();

      // ---- 收敛检查 ----
      altmin_converged = check_altmin_convergence();
      if (altmin_converged)
        pcout << "  ✓ AltMin converged after " << altmin_iter << " iters"
              << std::endl;
    }

    // ---- 6. 输出 ----
    if (load_step % 5 == 0)
    {
      pcout << "  [Output step " << load_step << "]" << std::endl;
      output_results(load_step);
    }

    // ---- 7. AMR (每 refine_frequency 步) ----
    // if (load_step % refine_frequency == 0 && load_step > 0)
    //   refine_grid();  // M4 里程碑实现
  }

  pcout << "\n==================================================" << std::endl;
  pcout << "  Done. " << std::endl;
  pcout << "==================================================" << std::endl;
}

// ========================================================================
// Macaulay 正部应变 — 相场能量驱动力
// 直接照搬 deal.II Gallery H_plus() 实现
// ========================================================================
template <int dim>
double PhaseFieldFlow<dim>::H_plus(const SymmetricTensor<2, dim> &strain) const
{
  const double tr_strain = trace(strain);
  const double Mac_tr_strain = std::max(tr_strain, 0.0);

  const std::array<double, dim> principal_strains = eigenvalues(strain);

  double tr_sqr_pos = 0.0;
  for (unsigned int i = 0; i < dim; ++i)
    tr_sqr_pos += std::max(principal_strains[i], 0.0) *
                  std::max(principal_strains[i], 0.0);

  return 0.5 * lambda() * Mac_tr_strain * Mac_tr_strain +
         mu() * tr_sqr_pos;
}

// ========================================================================
// 以下函数先留 stub — 组装 + Biot 耦合在 M2/M3 补
// 骨架能编译 + AltMin 循环可跑 → 逐步填充物理
// ========================================================================
template <int dim> void PhaseFieldFlow<dim>::declare_parameters() {}
template <int dim> void PhaseFieldFlow<dim>::parse_parameters(const std::string &) {}
template <int dim> void PhaseFieldFlow<dim>::setup_mesh() {}
template <int dim> void PhaseFieldFlow<dim>::setup_boundary_conditions() {}
template <int dim> void PhaseFieldFlow<dim>::setup_system_elastic() {}
template <int dim> void PhaseFieldFlow<dim>::setup_system_damage() {}
template <int dim> void PhaseFieldFlow<dim>::setup_system_flow() {}
template <int dim> void PhaseFieldFlow<dim>::assemble_system_elastic() {}
template <int dim> void PhaseFieldFlow<dim>::assemble_system_damage() {}
template <int dim> void PhaseFieldFlow<dim>::assemble_system_flow() {}
template <int dim> void PhaseFieldFlow<dim>::solve_elastic_subproblem() {}
template <int dim> void PhaseFieldFlow<dim>::solve_damage_subproblem() {}
template <int dim> void PhaseFieldFlow<dim>::solve_flow_subproblem() {}
template <int dim> void PhaseFieldFlow<dim>::compute_biot_body_force() {}
template <int dim> void PhaseFieldFlow<dim>::apply_biot_force_to_elastic_rhs() {}
template <int dim> void PhaseFieldFlow<dim>::update_permeability() {}
template <int dim> bool PhaseFieldFlow<dim>::check_altmin_convergence() { return true; }
template <int dim> void PhaseFieldFlow<dim>::update_history_field() {}
template <int dim> void PhaseFieldFlow<dim>::refine_grid() {}
template <int dim> void PhaseFieldFlow<dim>::output_results(const unsigned int) const {}

// 显式实例化 (deal.II 模板类需要在编译单元实例化)
template class PhaseFieldFlow<2>;
template class PhaseFieldFlow<3>;

} // namespace PhaseFieldFlow
