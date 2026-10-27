# ============================================================
#  vDef + Darcy Flow — Docker Build & Run (Windows PowerShell)
# ============================================================
#  PREREQUISITE: Install Docker Desktop from https://www.docker.com/
#
#  RIGHT-CLICK → Run with PowerShell
#  Or:  powershell -ExecutionPolicy Bypass -File run_docker.ps1
# ============================================================

$ErrorActionPreference = "Stop"

$Mef90Root = Split-Path -Parent $PSScriptRoot
$Image = "bourdin/mef90centos7openmpio:latest"

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "  vDef + m_Flow (Darcy + Biot) — Docker One-Liner"            -ForegroundColor Cyan
Write-Host "  Image:   $Image"                                             -ForegroundColor Yellow
Write-Host "  Mount:   $Mef90Root -> /opt/MEF90"                           -ForegroundColor Yellow
Write-Host "============================================================" -ForegroundColor Cyan

# --- Check Docker ---
try { docker --version | Out-Null }
catch {
    Write-Host "ERROR: Docker Desktop not found." -ForegroundColor Red
    Write-Host "Install from: https://www.docker.com/products/docker-desktop/" -ForegroundColor Red
    exit 1
}
Write-Host "Docker OK" -ForegroundColor Green

# --- Pull if needed ---
$exists = docker image inspect $Image 2>$null
if (-not $exists) {
    Write-Host ">>> Pulling Docker image (first run ~2GB)..." -ForegroundColor Yellow
    docker pull $Image
}

# --- Build & Run ---
Write-Host ">>> Compiling m_Flow + vDef + running Terzaghi benchmark..." -ForegroundColor Yellow
docker run --rm -it `
    -v "${Mef90Root}:/opt/MEF90" `
    -w /opt/MEF90 `
    $Image `
    bash scripts/build_and_run.sh

Write-Host ""
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "  DONE — check vDef/terzaghi1D-flow.e-s.0000.*.exo" -ForegroundColor Green
Write-Host "============================================================" -ForegroundColor Cyan
