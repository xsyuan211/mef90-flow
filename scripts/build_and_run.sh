#!/usr/bin/env bash
# ============================================================================
# vDef + m_Flow (Darcy + Biot) — One-Shot Build & Terzaghi Validation
# ============================================================================
# REQUIREMENTS:
#   - Bourdin patched PETSc-3.3 (gitlab.com/blaisebourdin/petsc.git branch mef90-3.3)
#   - gfortran, mpich/openmpi, exodusii, netcdf, sowing
#   - OR: Docker with bourdin/mef90centos7openmpio:latest
#
# USAGE (Docker, easiest):
#   docker run --rm -it -v $(pwd):/opt/MEF90 bourdin/mef90centos7openmpio:latest \
#     bash /opt/MEF90/scripts/build_and_run.sh
#
# USAGE (Linux with pre-installed PETSc):
#   export PETSC_DIR=/path/to/petsc-3.3-mef90
#   export PETSC_ARCH=Linux-gcc9-mpich32-O   # IMPORTANT: must contain a gcc version string!
#   cd mef90-master && bash scripts/build_and_run.sh
# ============================================================================

set -euo pipefail

MEF90_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$MEF90_DIR"

echo "============================================================"
echo "  vDef + Darcy Flow (Biot Coupling) — Build & Validate"
echo "  MEF90_DIR = $MEF90_DIR"
echo "  PETSC_DIR = ${PETSC_DIR:-<not set>}"
echo "  PETSC_ARCH = ${PETSC_ARCH:-<not set>}"
echo "============================================================"

# --- PETSc sanity check ---
if [ -z "${PETSC_DIR:-}" ] || [ ! -d "$PETSC_DIR" ]; then
  echo "ERROR: PETSC_DIR not set or invalid. Need Bourdin's patched PETSc-3.3."
  exit 1
fi
if [ -z "${PETSC_ARCH:-}" ]; then
  echo "ERROR: PETSC_ARCH not set. Example: Linux-gcc9-mpich32-O"
  exit 1
fi
if [ ! -f "$PETSC_DIR/lib/$PETSC_ARCH/libpetsc.a" ] && [ ! -f "$PETSC_DIR/lib/$PETSC_ARCH/libpetsc.so" ]; then
  echo "ERROR: PETSc library not found at $PETSC_DIR/lib/$PETSC_ARCH/"
  exit 1
fi
echo "PETSc OK ✓"

# --- Build ---
echo ""
echo ">>> Building MEF90 core..."
make MEF90 2>&1 | tail -5

echo ">>> Building m_DefMech..."
make m_DefMech 2>&1 | tail -5

echo ">>> Building m_HeatXfer..."
make m_HeatXfer 2>&1 | tail -5

echo ">>> Building m_Flow (NEW — Darcy + Biot)..."
make m_Flow 2>&1 | tail -20

echo ">>> Building vDef..."
make vDef 2>&1 | tail -10

echo ""
echo "============================================================"
echo "  BUILD SUCCESS ✓"
echo "  Binary: bin/$PETSC_ARCH/vDef"
echo "============================================================"

# --- Run Terzaghi 1D Consolidation Benchmark ---
echo ""
echo ">>> Running Terzaghi 1D Consolidation Benchmark..."
echo "    k_min=5e-11 m², μ=1e-3 Pa·s, S=4.6e-10 Pa⁻¹"
echo "    c_v ≈ 109 m²/s, T_max = 1.09 (at t=1s)"
echo ""

cd vDef
BINARY="../bin/$PETSC_ARCH/vDef"
PREFIX="data/Beam1x10-tri3"
OPTIONS="data/terzaghi1D.yaml"

if [ ! -f "$BINARY" ]; then
  echo "ERROR: Binary not found at $BINARY"
  exit 1
fi
if [ ! -f "$PREFIX.exo" ]; then
  echo "ERROR: Mesh not found at $PREFIX.exo"
  exit 1
fi
if [ ! -f "$OPTIONS" ]; then
  echo "ERROR: Options YAML not found at $OPTIONS"
  exit 1
fi

mpirun -np 2 "$BINARY" \
  -prefix "$PREFIX" \
  -options_file_yaml "$OPTIONS" \
  2>&1 | tee /tmp/terzaghi_output.log || {
    echo ""
    echo "WARNING: Solver exited with code $? — check /tmp/terzaghi_output.log"
  }

echo ""
echo "============================================================"
echo "  RUN COMPLETE"
echo "  Output: vDef/terzaghi1D-defmech.e-s.0000.*.exo"
echo "          vDef/terzaghi1D-flow.e-s.0000.*.exo"
echo "============================================================"
echo ""
echo "Terzaghi Analytical Validation:"
echo "  T(t) = c_v·t/H² = 109·t/100 = 1.09·t"
echo "  t=0.5s → T=0.54 → consolidation ≈ 42%"
echo "  t=1.0s → T=1.09 → consolidation ≈ 63%"
echo "  t=1.5s → T=1.63 → consolidation ≈ 80%"
echo ""
echo "To extract pressure time series from EXO:"
echo "  python3 -c \"import exodus; d=exodus.exodus('terzaghi1D-flow.e-s.0000.0.exo','r'); print(d.get_variable_names('n')); [print(f'Step {s}: p(z=0)={d.get_node_variable_values(\\\"Pressure\\\",s,ids=[1])[0]:.4e}') for s in range(d.num_steps()+1)]\""
echo "============================================================"
