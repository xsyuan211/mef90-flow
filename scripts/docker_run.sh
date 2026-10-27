#!/usr/bin/env bash
# ============================================================================
# DOCKER ONE-LINER: vDef + m_Flow Build & Terzaghi Validation
# ============================================================================
# This pulls Bourdin's pre-built Docker image (patched PETSc-3.3 + all deps)
# and mounts the local mef90-master directory inside to compile + run.
#
# PREREQUISITE: Docker Desktop (any platform)
#   Windows:  https://www.docker.com/products/docker-desktop/
#   macOS:    brew install --cask docker
#   Linux:    apt install docker.io
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MEF90_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Bourdin's official pre-built image (OpenMPI variant = fastest)
IMAGE="bourdin/mef90centos7openmpio:latest"
# Alternative: bourdin/mef90centos7mpichg:latest  (mpich debug variant)

echo "============================================================"
echo "  DOCKER BUILD + RUN"
echo "  Image:    $IMAGE"
echo "  Mount:    $MEF90_ROOT → /opt/MEF90"
echo "============================================================"

# --- Pull image if not cached ---
if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
  echo ">>> Pulling Docker image (first run — ~2GB)..."
  docker pull "$IMAGE"
fi

# --- Build & Run inside container ---
docker run --rm -it \
  -v "$MEF90_ROOT:/opt/MEF90" \
  -w /opt/MEF90 \
  "$IMAGE" \
  bash scripts/build_and_run.sh
