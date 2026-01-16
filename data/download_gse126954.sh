#!/usr/bin/env bash
# data/download_gse126954.sh
#
# Minimal downloader for GEO supplementary files.
# - Uses curl if available, otherwise wget.
# - Does NOT require sha256sum (no checksum verification).

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data"
RAW_DIR="${DATA_DIR}/raw"

mkdir -p "${RAW_DIR}"

# GEO "download?acc=" endpoint redirects to the actual file
ANNOT_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126954&format=file&file=GSE126954%5Fcell%5Fannotation%2Ecsv%2Egz"
MTX_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126954&format=file&file=GSE126954%5Fgene%5Fby%5Fcell%5Fcount%5Fmatrix%2Etxt%2Egz"

ANNOT_OUT="${RAW_DIR}/GSE126954_cell_annotation.csv.gz"
MTX_OUT="${RAW_DIR}/GSE126954_gene_by_cell_count_matrix.txt.gz"

have_cmd() { command -v "$1" >/dev/null 2>&1; }

download() {
  local url="$1"
  local out="$2"

  if [[ -f "$out" ]]; then
    echo "[download] exists: ${out}"
    return 0
  fi

  echo "[download] ${url}"
  echo "           -> ${out}"

  if have_cmd curl; then
    curl -L --fail --retry 5 --retry-delay 2 -o "${out}.partial" "${url}"
  elif have_cmd wget; then
    wget -O "${out}.partial" "${url}"
  else
    echo "ERROR: need curl or wget." >&2
    exit 1
  fi

  mv "${out}.partial" "${out}"
}

echo "[setup] raw dir: ${RAW_DIR}"

download "${ANNOT_URL}" "${ANNOT_OUT}"
download "${MTX_URL}" "${MTX_OUT}"

cat <<EOF

[done] Files present:
  - ${ANNOT_OUT}
  - ${MTX_OUT}

[next] Run Phase-1 (example):
  python3 scripts/phase1_runner.py \\
    --annot data/raw/GSE126954_cell_annotation.csv.gz \\
    --mtx   data/raw/GSE126954_gene_by_cell_count_matrix.txt.gz \\
    --knn_k 20 --svd_components 50 --max_centers 5000 --global_perms 1000 --seed 0

EOF