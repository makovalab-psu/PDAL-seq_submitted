#Swaps two chromosome names in one column of a tab-separated stream
#-v pair=chr12,chr19  the two names to swap (chr12 also matches chr12_hap1 etc.); empty pair passes lines through
#-v col=1             the column holding the chromosome name (default 1)
#Used for siamang annotations made on mSymSyn1 v2.1, which only swapped the chr12 and chr19
#labels relative to v2.0 (GCA_028878055.2), the assembly the reads are mapped to

BEGIN { FS = OFS = "\t"; if (col == "") col = 1; split(pair, s, ",") }
pair != "" {
    n = index($col, "_")
    head = n ? substr($col, 1, n - 1) : $col
    tail = n ? substr($col, n) : ""
    if (head == s[1]) $col = s[2] tail
    else if (head == s[2]) $col = s[1] tail
}
{ print }
