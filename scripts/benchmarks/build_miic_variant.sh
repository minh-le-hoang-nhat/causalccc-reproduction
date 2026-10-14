#!/usr/bin/env bash
# =============================================================================
# Build a variant of miic (branch causalccc, the commit pinned in envs/install_R_packages.R) with zero or more
# patches, and install it in its own library, ~/R/miic_variants/<name>. The normal R library is not touched.
# Every variant, the unpatched "baseline" included, is compiled the same way, so build flags never differ.
#
# Usage (from the root of the repository):
#   bash scripts/benchmarks/build_miic_variant.sh baseline
#   bash scripts/benchmarks/build_miic_variant.sh no_rcout envs/patches/miic/01_no_rcout_in_omp.patch
# Then run R with this variant first in the library path:
#   R_LIBS=~/R/miic_variants/no_rcout Rscript scripts/run_miic.R ...
# =============================================================================

set -euo pipefail

MIIC_REPO="https://github.com/miicTeam/miic_R_package.git"
MIIC_COMMIT="5ca0b222cc36bee902aa87e5ecfc480428d98927"
SRC_CACHE="$HOME/miic_src/miic_R_package"

if [[ $# -lt 1 ]]; then
  echo "Usage: bash scripts/benchmarks/build_miic_variant.sh <name> [patch ...]" >&2
  exit 1
fi
name="$1"; shift
lib="$HOME/R/miic_variants/$name"

# Local clone of the miic repository, made once
if [[ ! -d "$SRC_CACHE" ]]; then
  mkdir -p "$(dirname "$SRC_CACHE")"
  git clone --quiet "$MIIC_REPO" "$SRC_CACHE"
fi
if ! git -C "$SRC_CACHE" cat-file -e "${MIIC_COMMIT}^{commit}" 2>/dev/null; then
  git -C "$SRC_CACHE" fetch --quiet origin causalccc
fi

# Clean copy of the pinned commit, then the patches in the order given
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
git -C "$SRC_CACHE" archive "$MIIC_COMMIT" | tar -x -C "$work"
for p in "$@"; do
  echo "Applying $p"
  patch --quiet -p1 -d "$work" < "$p"
done

rm -rf "$lib"
mkdir -p "$lib"
R CMD INSTALL --no-docs --no-html --library="$lib" "$work" > "$lib.build.log" 2>&1 || {
  tail -30 "$lib.build.log" >&2
  echo "ERROR: build failed, full log in $lib.build.log" >&2
  exit 1
}

# Record what this variant is made of
{
  echo "variant: $name"
  echo "miic commit: $MIIC_COMMIT"
  echo "built: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  for p in "$@"; do echo "patch: $p $(md5sum < "$p" | cut -d ' ' -f 1)"; done
} > "$lib/VARIANT.txt"

echo "Installed miic variant '$name' in $lib"
