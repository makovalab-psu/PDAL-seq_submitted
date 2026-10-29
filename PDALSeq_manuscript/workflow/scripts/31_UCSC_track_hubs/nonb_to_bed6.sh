#Converts a non-B DNA annotation to BED6 for bedToBigBed
#$1 the non-B annotation (BED3, or BED6 with strand for G4 predictions)
#$2 the name written to column 4 (e.g. G4)
#Header lines are skipped; the strand is kept when column 6 is + or -

awk -v name="$2" 'BEGIN{OFS="\t"} $2 ~ /^[0-9]+$/ {
    strand = (NF >= 6 && ($6 == "+" || $6 == "-")) ? $6 : "."
    print $1, $2, $3, name, 0, strand
}' "$1"
