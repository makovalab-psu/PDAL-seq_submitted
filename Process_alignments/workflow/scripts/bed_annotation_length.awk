#!/usr/bin/awk -f
BEGIN {
    OFS = "\t"  # Set the output field separator to a tab
}
{
    print $3 - $2 + 1

}

