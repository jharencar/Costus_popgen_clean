#!/bin/bash
# Build Dsuite from source for Farm.
# Uses conda env `dsuite` for a modern libstdc++.
# https://github.com/millanek/Dsuite
#
# Usage:
#   bash install_dsuite_farm.sh

set -euo pipefail

PREFIX="${PREFIX:-${HOME}/bin}"
SRC_DIR="${SRC_DIR:-${HOME}/src/Dsuite}"
CONDA_ENV="${CONDA_ENV:-dsuite}"
mkdir -p "${PREFIX}" "$(dirname "${SRC_DIR}")"

unset CONDA_EXE CONDA_PREFIX CONDA_PYTHON_EXE CONDA_SHLVL CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER || true
module load conda/latest
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV}"; then
  conda create -y -n "${CONDA_ENV}" -c conda-forge cxx-compiler zlib make
fi
conda activate "${CONDA_ENV}"

if [[ ! -d "${SRC_DIR}/.git" ]]; then
  git clone https://github.com/millanek/Dsuite.git "${SRC_DIR}"
else
  git -C "${SRC_DIR}" pull --ff-only || true
fi

cd "${SRC_DIR}"
make clean || true
make -j"${SLURM_CPUS_PER_TASK:-4}"

BIN_SRC=""
for cand in "${SRC_DIR}/Build/Dsuite" "${SRC_DIR}/Dsuite"; do
  if [[ -x "${cand}" ]]; then
    BIN_SRC="${cand}"
    break
  fi
done
if [[ -z "${BIN_SRC}" ]]; then
  echo "Dsuite binary not found after make" >&2
  exit 1
fi

install -m 0755 "${BIN_SRC}" "${PREFIX}/Dsuite.bin"
CONDA_LIB="${CONDA_PREFIX}/lib"
cat > "${PREFIX}/Dsuite" << WRAP
#!/bin/bash
export LD_LIBRARY_PATH="${CONDA_LIB}\${LD_LIBRARY_PATH:+:\${LD_LIBRARY_PATH}}"
exec "\${HOME}/bin/Dsuite.bin" "\$@"
WRAP
chmod +x "${PREFIX}/Dsuite"

echo "Installed ${PREFIX}/Dsuite"
"${PREFIX}/Dsuite" 2>&1 | head -20 || true
