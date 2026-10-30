#Converts a non-B DNA annotation to BED6 for bedToBigBed
#$1 the non-B annotation (BED3, or BED6 with strand for G4 predictions)
#$2 the name written to column 4 (e.g. G4)
#$3 optional pair of chromosomes to swap, e.g. chr12,chr19 (matches chr12 and chr12_*)
#Header lines are skipped; the strand is kept when column 6 is + or -

awk -v name="$2" -v swap="${3:-}" 'BEGIN{OFS="\t"; split(swap, s, ",")} $2 ~ /^[0-9]+$/ {
    if (swap != "") {
        n = index($1, "_")
        head = n ? substr($1, 1, n - 1) : $1
        tail = n ? substr($1, n) : ""
        if (head == s[1]) $1 = s[2] tail
        else if (head == s[2]) $1 = s[1] tail
    }
    strand = (NF >= 6 && ($6 == "+" || $6 == "-")) ? $6 : "."
    print $1, $2, $3, name, 0, strand
}' "$1"
