#!/bin/bash
# =============================================================================
# Download NCHS Natality Detail Files from NBER (1968-1977)
# =============================================================================
#
# Downloads pre-built Stata .dta files (zipped) from NBER's data archive.
# These contain individual-level birth records for the entire US.
#
# Total download: ~327 MB (compressed), ~2 GB (uncompressed)
#
# Usage:
#   cd IV_dist/data/in/nber_natality && bash download_natality.sh
#
# =============================================================================

set -e

DEST_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_URL="https://data.nber.org/natality"

echo "Downloading NCHS natality data to: $DEST_DIR"
echo ""

for YEAR in $(seq 1968 1977); do
    FILE="natl${YEAR}.dta.zip"
    URL="${BASE_URL}/${YEAR}/${FILE}"
    OUT="${DEST_DIR}/${FILE}"
    DTA="${DEST_DIR}/natl${YEAR}.dta"

    if [ -f "$DTA" ]; then
        echo "  ${YEAR}: already exists (natl${YEAR}.dta), skipping"
        continue
    fi

    echo "  ${YEAR}: downloading ${FILE}..."
    curl -# -L -o "$OUT" "$URL"

    echo "  ${YEAR}: unzipping..."
    unzip -o -q "$OUT" -d "$DEST_DIR"
    rm -f "$OUT"

    echo "  ${YEAR}: done ($(du -h "$DTA" | cut -f1))"
done

echo ""
echo "All files downloaded."
ls -lh "$DEST_DIR"/natl*.dta
