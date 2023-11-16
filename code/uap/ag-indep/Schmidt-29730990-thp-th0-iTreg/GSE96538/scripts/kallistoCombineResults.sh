#!/bin/bash

nThreads=2
outDir="$*"

name=$(basename "$outDir")

# combine counts into one table
echo "[$name] Collecting results ..." >&2
tpmTemp=$(mktemp -d)
estcTemp=$(mktemp -d)
find "$outDir" -name abundance.tsv -print0 | xargs -0 -P $nThreads -n 1 -I{} \
    bash -c "sampName=\$(basename \$(dirname {}));cut -f 5 {} | sed \"1s/.*/\$sampName/\" > ${tpmTemp}/\$sampName"
find "$outDir" -name abundance.tsv -print0 | xargs -0 -P $nThreads -n 1 -I{} \
    bash -c "sampName=\$(basename \$(dirname {}));cut -f 4 {} | sed \"1s/.*/\$sampName/\" > ${estcTemp}/\$sampName"

echo "[$name] Combining results into one table ..." >&2
oneFile=$(find "$outDir" -name abundance.tsv | head -n1)
paste <(cut -f 1 "$oneFile") $tpmTemp/* > "${outDir}/tpm.tsv"
paste <(cut -f 1 "$oneFile") $estcTemp/* > "${outDir}/counts.tsv"
rm -r "$tpmTemp" "$estcTemp"

echo "[$name] Creating R objects ..." >&2
./makeRObjects.R "${outDir}/tpm.tsv" 2>&1 | xargs -n 1 -I{} echo "[$name] {}" >&2
./makeRObjects.R "${outDir}/counts.tsv" 2>&1 | xargs -n 1 -I{} echo "[$name] {}" >&2

echo "[$name] **DONE** finished successfully." >&2
