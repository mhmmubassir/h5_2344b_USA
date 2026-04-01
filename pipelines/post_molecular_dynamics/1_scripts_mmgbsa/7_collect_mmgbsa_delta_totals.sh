#!/usr/bin/env bash
# collect_gbneck2_totals.sh
# Walk all case folders, find mmgbsa/mmgbsa_gbneck2*.dat (or no extension),
# extract "DELTA TOTAL" and its Std. Dev., and write a CSV.
#
# Usage:
#   ./collect_gbneck2_totals.sh             # writes mmgbsa_delta_totals.csv
#   ./collect_gbneck2_totals.sh -o out.csv  # custom output path

set -Eeuo pipefail

OUT="mmgbsa_delta_totals.csv"

# ---- parse flags ----
while getopts ":o:h" opt; do
  case "$opt" in
    o) OUT="$OPTARG";;
    h)
      echo "Usage: $0 [-o output.csv]"
      exit 0
      ;;
    \?) echo "Unknown option: -$OPTARG" >&2; exit 2;;
  esac
done

# CSV header
echo "top_folder,mmgbsa_dir,mmgbsa_file,delta_total_kcal_per_mol,stdev_kcal_per_mol" > "$OUT"

# Helper: extract "Delta TOTAL" and stdev from a file, print "VAL,SD" or nothing
extract_dt_sd() {
  local f="$1"
  # Try common MMPBSA/MMGBSA summary formats.
  # Accept lines like:
  #   DELTA TOTAL        -45.67   3.21
  #   DELTA G binding:   -45.67  +/- 3.21
  #   DELTA TOTAL (kcal/mol) = -45.67  3.21
  awk '
    BEGIN{ OFS="," }
    function isnum(x){ return (x ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) }
    # Preferred keys
    /DELTA[ _]TOTAL|DELTA G binding|DELTA[ _]G[ _]binding/ {
      n=0
      for (i=1;i<=NF;i++){
        if ($i=="+/-" && i<NF && isnum($(i+1))) { nums[++n]=$(i+1); i++; continue }
        if (isnum($i)) nums[++n]=$i
      }
      if (n>=2) { print nums[1], nums[2]; exit }
    }
    END { }' "$f"
}

# Find any mmgbsa_gbneck2 file within an mmgbsa/ subdir (handles with/without .dat, suffixes)
# Use -print0 to guard against odd chars
while IFS= read -r -d '' file; do
  rel="${file#./}"                    # strip leading ./ for prettier CSV
  mmgbsa_dir="$(dirname "$rel")"
  top_folder="${mmgbsa_dir%%/*}"      # first path component
  # If path starts with no slash after strip (unlikely), fallback:
  if [[ "$top_folder" == "$mmgbsa_dir" ]]; then
    top_folder="."
  fi

  parsed="$(extract_dt_sd "$file" || true)"

  if [[ -n "$parsed" ]]; then
    dt="${parsed%%,*}"
    sd="${parsed#*,}"
  else
    dt="NA"
    sd="NA"
    echo "WARN: Could not parse DELTA TOTAL from: $rel" >&2
  fi

  # Append CSV row
  printf "%s,%s,%s,%s,%s\n" \
    "$top_folder" "$mmgbsa_dir" "$rel" "$dt" "$sd" >> "$OUT"

done < <(find . -type f -path "*/mmgbsa/*" \( -name "mmgbsa_gbneck2" -o -name "mmgbsa_gbneck2.dat" -o -name "mmgbsa_gbneck2*.dat" \) -print0)

echo "Wrote: $OUT"
