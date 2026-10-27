# vDef (F90/PETSc) → deal.II 迁移映射表

> 目标：把 vDef + m_Flow (Darcy + Biot + 相场断裂) 迁移到 deal.II C++
> 起点代码：deal.II Gallery `Phase_field_fracture_model_in_3D/phase_field.cc`
> 关键参考：IPACS (Wheeler+Wick+Lee, CMAME 2020) — Biot + 固定应力迭代 + 自适应 AMR

---

## 1. 架构级对照

| 模块 | vDef F90 (MEF90/PETSc) | deal.II C++ | 说明 |
|------|------------------------|-------------|------|
| 网格 | EXO + DMplex + SectionReal | `parallel::distributed::Triangulation` + p4est | deal.II 用 GridGenerator 或 GMSH 导入 |
| 自由度 | DMScal + DMVect (stratum 分层) | `DoFHandler<dim>` × N（每个 PDE 独立） | deal.II 可以分开或 FESystem 合体 |
| 边界条件 | DMMeshGetStratumIS + MEF90ISCreateCelltoVertex | `VectorTools::interpolate_boundary_values` + `AffineConstraints` | deal.II 自动处理 hanging nodes |
| YAML 解析 | PetscBagRegister* + PetscBagSetFromOptions | `ParameterHandler` + `InputReader` | deal.II 原生支持 .prm 参数文件 |
| 组装 | DiffusionInterface (预封装) + 手写 Gauss 循环 | FEValues + 手写积分循环 | deal.II 无 Darcy 封装，自己写 |
| 求解器 | SNES/F90 + KSP/CG | `SolverCG` + `SolverNewton` 或 Trilinos/PETSc | deal.II 可以用 PETSc 后端 |
| 历史变量 | VecAXPy + max(d, d_old) | `CellDataStorage` 或 `QuadraturePointHistory` | deal.II 原生 AMR 后自动传递 |
| 并行 | MPI + MEF90ISAllGatherMerge | `parallel::distributed::*` 自动处理 | deal.II 内置 |
| AMR | ❌ (统一网格) | `SolutionTransfer` + `parallel::distributed::GridRefinement` | deal.II Gallery 已实现 |
| 后处理 | EXODUS II (`MEF90ViewEXO`) | `DataOut` (VTK + HDF5) | deal.II 原生 |

---

## 2. 子例程级映射

### 2.1 组装核心

| vDef F90 子例程 | 位置 | deal.II C++ 对应 |
|-----------------|------|-----------------|
| `FlowOperator` (Darcy div·k∇p) | m_Flow/Assembly L180 | **手写** `assemble_system_flow()` — `fe_values.JxW()` + 组装扩散矩阵 + 质量 |
| `FlowBilinearForm` (Jacobian) | m_Flow/Assembly L210 | 同上，刚度矩阵 + 压缩项 |
| `FlowIFunction` (非线性残差) | m_Flow/Assembly L260 | 残差 = 组装 RHS |
| `SetSourceCst` (体积源项) | m_Flow/Assembly L620 | `fe_values.get_quadrature_points()` → 插值源项 |
| `ApplyDirichletResidual` (强制 BC) | m_Flow/Assembly L376 | `AffineConstraints::distribute()` — deal.II 自动处理 |
| `FlowComputePressureGradient` (∇p 体积加权) | m_Flow/Assembly L849 | `fe_values[x].gradient(p, q)` 直接取 |
| `FlowApplyBiotForce` (f_p = -α∇p) | m_Flow/Assembly L920 | **关键**：Biot 体积力是 u-方程的 RHS |

### 2.2 vDef 初始化链

| vDef F90 | deal.II C++ |
|----------|-------------|
| `MEF90CtxCreate` (petsc.binary) | `ParameterHandler` 读 `.prm` 文件 |
| `MEF90DefMechCtxCreate` + SetSections | PhaseField 类构造函数 |
| `MEF90DefMechCreateVectors` + SetFromOptions | `setup_system_elastic()` + `setup_system_damage()` + **新增 `setup_system_flow()`** |
| `MEF90FlowCtxCreate` + baseForce 分配 | 同上 `setup_system_flow()` |
| `MEF90DefMechSetTransients` + `VecCopy(force→baseForce)` | AltMin 循环内：flow RHS = baseForce + (-α)·∇p |

### 2.3 AltMin 耦合循环（核心算法一致）

```
vDef F90 (我们现在的):                deal.II C++ (PhaseField gallery + 扩展):
┌────────────────────────────────┐    ┌───────────────────────────────────────┐
│ DO time_step = 1, N_steps       │    │ DO load_step = 1, num_load_steps      │
│   Call SetTransients (重力等)   │    │   setup_constraints_elastic()        │
│   VecCopy(force → baseForce)   │    │   Call SetBoundaryValues(流固耦合)   │
│   DO AltMin iter = 1, N_AltMin │    │   DO AltMin iter = 1, N_AltMin       │
│     StepV: SNES(damage)        │    │     StepV: solve_damage_subproblem   │
│     UpdatePerm: k(v)           │    │     update_history_field()           │
│     StepP: SNES(p)             │    │     StepP: solve_flow_subproblem ← NEW│
│     ApplyBiot: force += -α·∇p  │    │     compute_biot_body_force() ← NEW  │
│     StepU: SNES(u)             │    │     StepU: solve_elastic_subproblem  │
│     Check converged?           │    │     check_convergence()             │
│   END AltMin                   │    │   END AltMin                          │
│   UpdateDamageHistory          │    │   output_results()                   │
│ END time_step                  │    │   refine_grid() ← deal.II 原生 AMR   │
└────────────────────────────────┘    └───────────────────────────────────────┘
```

### 2.4 关键公式对照

| 物理量 | vDef F90 | deal.II C++ |
|--------|----------|-------------|
| 应力 | `σ = (1-d)·(λ tr(ε) I + 2μ ε)` | 相同，deal.II 用 `SymmetricTensor<2,dim>` |
| Biot 有效应力 | `σ_total = σ_eff - αp I` | 相同 |
| Biot 体积力 | `f_p = -α·∇p` | 加到 u-方程 RHS |
| Darcy 通量 | `q = -k(v)/μ · ∇p` | deal.II 组装扩散项 |
| k(v) 渗透率 | `k_min + (k_max-k_min)·v_cellAvg` | `damage→permeability()` 函数 |
| 历史变量 | `d_history = max(d_history, d_current)` | `CellDataStorage` 存积分点历史 |
| 应变分解 | MACAULAY(tr_strain) + 主应变 | deal.II Gallery `H_plus()` 已实现 ✅ |

---

## 3. deal.II 原型文件结构

```
dealii-port/
├── CMakeLists.txt              # deal.II 9.7+ + Trilinos 或 PETSc
├── README.md                   # 编译运行说明
├── parameters.prm              # deal.II ParameterHandler 输入
├── include/
│   ├── phase_field_flow.h      # 主类：AltMin + 三个子问题
│   ├── elastic_assembly.h      # u-方程组装（含 Biot 体积力）
│   ├── damage_assembly.h       # d-方程组装（MACAULAY + 历史变量）
│   ├── flow_assembly.h         # p-方程组装（Darcy + k(v) + 压缩项）
│   └── biot_coupling.h         # Biot 力计算 + AltMin 控制
├── source/
│   ├── main.cc                 # PhaseFieldFlow 入口
│   ├── phase_field_flow.cc     # 主类实现
│   ├── elastic_assembly.cc     # 组装 + Newton 求解
│   ├── damage_assembly.cc      # 组装 + 线性格式
│   ├── flow_assembly.cc        # 组装（从 DiffusionInterface 移植）
│   └── biot_coupling.cc        # ∇p 计算 + Biot 力
└── tests/
    └── terzaghi/               # Terzaghi 1D 固结验证
```

---

## 4. 里程碑（M1 → M4）

| 里程碑 | 交付 | 周期 | 依赖 |
|--------|------|------|------|
| **M1** | deal.II Phase Field 基础原型（仅 u + d 两个子问题） | 1 周 | deal.II 9.7 环境 |
| **M2** | 加第三个子问题 Darcy 压力场 + k(v) 渗透率 | 1 周 | M1 |
| **M3** | Biot 耦合完整 AltMin 循环 + Terzaghi 1D 验证 | 1 周 | M2 |
| **M4** | 自适应 AMR + 并行规模化 + KGD/圆形孔洞算例 | 2 周 | M3 |

**每个里程碑完成后都能跑通一个验证算例**——不会出现"代码写完但没编译过"的情况。

---

## 5. deal.II 入门资源（按优先级）

| 资源 | URL | 学什么 |
|------|-----|--------|
| Phase Field Gallery | https://dealii.org/developer/doxygen/code-gallery/Phase_field_fracture_model_in_3D/ | 我们直接的起点 |
| 官方 Tutorial 12 | https://dealii.org/developer/doxygen/deal.II/Tutorial.html | FEM 组装基础 |
| 官方 Tutorial 49 | https://dealii.org/developer/doxygen/deal.II/Tutorial.html | 多孔介质 Darcy（DG 方法） |
| IPACS 论文 | Wheeler+Wick+Lee, CMAME 2020 | Biot 固定应力迭代 |
| RACCOON | https://hugary1995.github.io/raccoon | MOOSE + deal.II 混合 AMR 参考 |

---

## 6. deal.II vs MEF90 关键差异（坑点预警）

| # | 差异 | MEF90 F90 写法 | deal.II C++ 写法 | 坑 |
|---|------|----------------|------------------|-----|
| 1 | **Gauss 积分** | `MEF90DiffusionOperatorSet` 预封装 | 手写 `fe_values.reinit(cell); for q: ...` | deal.II 更灵活但要自己写 Darcy |
| 2 | **非结构网格** | EXODUS II + DMplex stratum | deal.II 用 GridGenerator 或 GMSH 导入 | 我们的 Beam1x10-tri3 能直接转 GMSH |
| 3 | **历史变量** | `Vec` + VecAXPy | `CellDataStorage` 或 `QuadraturePointHistory` | deal.II 原生支持 AMR 后传递 ✅ |
| 4 | **Dirichlet BC** | `MEF90ISCreateCelltoVertex` + `VecSetValues` | `VectorTools::interpolate_boundary_values` | deal.II 自动处理 hanging nodes ✅ |
| 5 | **YAML** | PetscBag | ParameterHandler (.prm) | 参数格式不同，需要转换 |
| 6 | **Darcy 封装** | `MEF90DiffusionOperatorSet` | **无**，手写 | 这是最大工作量 — 参考 Tutorial 49 DG |
| 7 | **线性代数** | PETSc (强制) | Trilinos 或 PETSc 可选 | 我们 vDef 用 PETSc 后端最稳 |
