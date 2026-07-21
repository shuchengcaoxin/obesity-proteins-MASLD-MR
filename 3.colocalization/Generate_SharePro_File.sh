#!/bin/bash

Trait1=(MASLD_MVP MASLD_MVP MASLD_MVP MASLD_MVP MASLD_MVP MASLD_MVP)
Trait2=(GUSB_DeCODE GUSB_UKBPPP HP_DeCODE NCAN_DeCODE NCAN_UKBPPP SEZ6L_DeCODE)
chrom=(7 7 16 19 19 22)
pos=(66771620 65767866 72067626 19219115 19215154 26302786)
start=(65771802 64768123 71068172 18219782 18219782 25302889)
end=(67771160 66767551 73067426 20218828 20211184 27302588)

# Get current values for this array job
t1=${Trait1[$i]}
t2=${Trait2[$i]}
chr=${chrom[$i]}
bp=${pos[$i]}

file1="${t1}_${chr}_${bp}_SharePro.txt"
file2="${t2}_${chr}_${bp}_SharePro.txt"
common_file="/scratch/butler/bangli.cao/BMI_driven_proteins_MASLD/3.colocalization/sharepro/predata/${t1}_${t2}_${chr}_${bp}_Common_Variants.txt"
output1_prefix="${t1}_${chr}_${bp}_SharePro_Final.txt"
output2_prefix="${t2}_${chr}_${bp}_SharePro_Final.txt"


# Create output with headline
head -n 1 "$file1" > "$output1"
head -n 1 "$file2" > "$output2"

# Keep other columns with matching SNPs and ensure uniqueness
awk 'NR==FNR {common[$1]; next} $1 in common && !seen[$1]++' "$common_variants" "$file1" >> "$output1"
awk 'NR==FNR {common[$1]; next} $1 in common && !seen[$1]++' "$common_variants" "$file2" >> "$output2"

echo "SharePro files created: $output1 and $output2"