#!/usr/bin/awk -f
BEGIN {
    OFS = "\t"  # Set the output field separator to a tab
}
{
    # Check if the line starts with 's' to ensure it's a sequence line
    if ($1 == "s") {
        # Assuming the fields are structured as follows:
        # s <chromosome> <start> <length> <strand> <score> <sequence>
        split($2, chrom_parts, ".")
        chrom = chrom_parts[2]
        strand = $5        # Strand information
        if (strand == "+"){
            start = $3 
            end = $3 + $4
        }
        if (strand == "-"){
            start = $6 - $3 - $4 
            end = $6 - $3 
        }
        # Print the fields in BED format
        print chrom, start, end, strand, $6, $7
    }
}

