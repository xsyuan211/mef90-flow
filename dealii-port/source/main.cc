// main.cc — PhaseFieldFlow 入口
//
// 对应 vDef.F90 的主程序入口 + Init + AltMin 循环 + Cleanup
// deal.II 的 main 非常简洁，所有逻辑在 PhaseFieldFlow 类里

#include "phase_field_flow.h"

int main(int argc, char *argv[])
{
  try
  {
    // deal.II 自动处理 MPI 初始化（deal.II 9.4+）
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

    PhaseFieldFlow::PhaseFieldFlow<3> flow;
    flow.run();
  }
  catch (const std::exception &exc)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}
