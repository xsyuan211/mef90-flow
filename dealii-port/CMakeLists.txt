#
# deal.II CMakeLists — Phase Field Fracture + Darcy + Biot Coupling
#
# Requires: deal.II >= 9.7.1 with DEAL_II_WITH_TRILINOS or DEAL_II_WITH_PETSC
#           DEAL_II_WITH_P4EST = ON (for adaptive AMR)
#
# Build:
#   mkdir build && cd build
#   cmake -DDEAL_II_DIR=/path/to/dealii-9.7 ../
#   make -jN
#   mpirun -np 2 ./phase_field_flow parameters.prm
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

# ---- 依赖 deal.II ----
FIND_PACKAGE(deal.II 9.7.1
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n*** deal.II not found! ***\nSet DEAL_II_DIR env var or pass -DDEAL_II_DIR=/path/to/dealii-9.7")
ENDIF()

# ---- 校验必要组件 ----
IF(NOT ((DEAL_II_WITH_PETSC AND NOT DEAL_II_PETSC_WITH_COMPLEX) OR DEAL_II_WITH_TRILINOS))
  MESSAGE(FATAL_ERROR "deal.II needs DEAL_II_WITH_PETSC=ON or DEAL_II_WITH_TRILINOS=ON")
ENDIF()
IF(NOT DEAL_II_WITH_P4EST)
  MESSAGE(FATAL_ERROR "deal.II needs DEAL_II_WITH_P4EST=ON (for parallel AMR)")
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(phase_field_flow)

# ---- 源文件 ----
FILE(GLOB_RECURSE TARGET_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

SET(TARGET phase_field_flow)

ADD_EXECUTABLE(${TARGET} ${TARGET_SRC})

TARGET_INCLUDE_DIRECTORIES(${TARGET} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include)

DEAL_II_SETUP_TARGET(${TARGET})

# ---- 拷贝参数文件到 build ----
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/parameters.prm
  ${CMAKE_CURRENT_BINARY_DIR}/parameters.prm
  COPYONLY)
