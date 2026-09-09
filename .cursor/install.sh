#!/usr/bin/env bash
# Idempotent Cloud Agent setup for the bulkMAE R package.
#
# Layers:
#   1. System toolchain: R (>= 4.4, from the CRAN apt repo), pandoc, and the
#      C/Fortran dev libraries the Bioconductor stack compiles against.
#   2. R configuration: a user ~/.Rprofile pinning the Posit Public Package
#      Manager (P3M) binary repo for fast CRAN installs, plus a private user
#      library.
#   3. R packages: pak installs the package's hard dependencies and the core
#      analysis backends exercised by the default test suite and the executable
#      vignettes (mirrors .github/workflows/R-CMD-check.yaml extra-packages).
#
# Re-running is safe: apt, add-apt-repository and pak all converge without
# rewriting state. The heavier optional GitHub-only backends (MuSiC,
# immunedeconv, BayesPrism, WGCNA, ...) are intentionally left to
# full-backend-check; the default suite skips them cleanly when absent.
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive
R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"

log() { printf '\n=== %s ===\n' "$*"; }

apt_retry() {
  local i
  for i in 1 2 3 4; do
    if sudo -E apt-get "$@"; then
      return 0
    fi
    echo "apt-get $* failed (attempt $i); retrying..." >&2
    sleep $((i * 4))
  done
  return 1
}

# --- 1. System toolchain --------------------------------------------------
if ! command -v R >/dev/null 2>&1; then
  log "Installing base tooling and CRAN apt repository"
  apt_retry update -qq
  apt_retry install -y --no-install-recommends \
    software-properties-common dirmngr gnupg ca-certificates curl wget pandoc

  if [ ! -f /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc ]; then
    curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
      | sudo tee /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc >/dev/null
  fi
  . /etc/os-release
  sudo add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu ${UBUNTU_CODENAME}-cran40/"
  apt_retry update -qq

  log "Installing R"
  apt_retry install -y --no-install-recommends r-base r-base-dev
fi

log "Installing system libraries for the Bioconductor stack"
apt_retry update -qq
apt_retry install -y --no-install-recommends \
  build-essential gfortran pandoc \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  zlib1g-dev libbz2-dev liblzma-dev libpcre2-dev \
  libgit2-dev libssh2-1-dev libglpk-dev libgmp-dev \
  libblas-dev liblapack-dev cmake git

# --- 2. R configuration ---------------------------------------------------
log "Configuring R user profile and library"
mkdir -p "$R_LIBS_USER"
if ! grep -q "R_LIBS_USER=" "$HOME/.Renviron" 2>/dev/null; then
  echo "R_LIBS_USER=$R_LIBS_USER" >> "$HOME/.Renviron"
fi

cat > "$HOME/.Rprofile" <<'RPROFILE'
local({
  options(
    repos = c(P3M = "https://packagemanager.posit.co/cran/__linux__/noble/latest"),
    HTTPUserAgent = sprintf(
      "R/%s R (%s)",
      getRversion(),
      paste(getRversion(), R.version$platform, R.version$arch, R.version$os)
    ),
    Ncpus = max(1L, parallel::detectCores())
  )
})
RPROFILE

# --- 3. R packages --------------------------------------------------------
log "Installing pak and R package dependencies"
Rscript - <<'RSCRIPT'
lib <- Sys.getenv("R_LIBS_USER", unset = file.path(Sys.getenv("HOME"), "R", "library"))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(lib)

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

# Package hard dependencies (Depends / Imports / LinkingTo).
pak::local_install_deps(
  ".",
  dependencies = c("Depends", "Imports", "LinkingTo"),
  upgrade = FALSE
)

# Core analysis backends used by the default test suite and executable
# vignettes. Mirrors .github/workflows/R-CMD-check.yaml extra-packages.
core_backends <- c(
  "rcmdcheck", "testthat", "knitr", "rmarkdown", "devtools",
  "airway", "edgeR", "DESeq2", "limma",
  "clusterProfiler", "fgsea", "org.Hs.eg.db"
)
pak::pkg_install(core_backends, upgrade = FALSE)

# Install the package itself so `library(bulkMAE)` works out of the box.
pak::local_install(".", upgrade = FALSE)

cat("\nbulkMAE loadable:", requireNamespace("bulkMAE", quietly = TRUE), "\n")
RSCRIPT

log "Setup complete"
