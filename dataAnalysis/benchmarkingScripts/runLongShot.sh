fastqFileStr="$1";
refFileStr="$2";                # name of reference file to use

bamFileStr="tmp";               # name of output bam file
outVcfStr="tmp";                # name of output vcf file
outConsensusNameStr="variant";  # name of consensus to output
intVariant=1;                   # counter for number variants
strErr=0;                       # captures error message from bcftools
numZeroVar=0; # number times bcftools applied 0 varaints (prevent infinit loop)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Call variants
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

minimap2 \
    -ax map-ont \
    "$refFileStr" \
    "$fastqFileStr" |
  samtools sort |
  samtools \
    view \
    -F 0x04 \
    -b \
  > "$bamFileStr.bam";

samtools \
    index \
    "$bamFileStr.bam";

samtools \
    faidx \
    "$refFileStr";

longshot \
    --bam "$bamFileStr.bam" \
    --ref "$refFileStr" \
    --out "$outVcfStr.vcf" \
    --strand_bias_pvalue_cutoff 0.01; # cutoff for ONT data

rm \
    "$bamFileStr.bam.bai" \
    "$refFileStr.fai" \
    "$bamFileStr.bam"; # no longer need

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-5: Prepare vcf for building a consensus
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Get last line of file (check if header)
strErr="$(
    sed \
        -n \
        '$s/\t.*//p;' \
        "$outVcfStr.vcf" \
)"; # Get the 

bgzip "$outVcfStr.vcf"; # makes into vcf.gz

tabix \
    -p vcf \
    "$outVcfStr.vcf.gz"; # index vcf file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-6: Build a consensus for each variant
#    sec-6 sub-1: building a consensus when only one variant & is exact match
#    sec-6 sub-2: building a consensus when one variant not exact match
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-6 Sub-1: building a consensus when only one variant is exact match
#*******************************************************************************

if [[ "$strErr" == *"CHROM"* ]]; then
# if bulding a consensus for a single variant & vcf file has no variants
    bcftools \
        consensus \
        -f "$refFileStr" \
        -H "1" \
        "$outVcfStr.vcf.gz" \
      > "variant-$intVariant.fasta"; \

    printf "Applied 0 variants for only consensus\n";
fi # if bulding a consensus for a single variant & vcf file has no variants

#*******************************************************************************
# Sec-6 Sub-2: building a consensus when one variant not exact match
#*******************************************************************************

while [[ "$strErr" != *"CHROM"* ]]; do # CHROM means header line
# while there are variatns to apply
    strErr="$( \
        { \
            bcftools \
                consensus \
                -f "$refFileStr" \
                -H "$intVariant" \
                "$outVcfStr.vcf.gz" \
              > "variant-$intVariant.fasta"; \
        } \
        2>&1\
    )"; # build consensuses and get error when no more variants in vcf
    # Need the "{ command > file.str; } 2>&1" to capture the error message
    
    # check if bcftools errored out (no more variants)
    if [[ "$strErr" == *"Can"*  || "$strErr" == *"Failed"* ]]; then
        rm "variant-$intVariant.fasta";
        break; # no longer need to continue loop
    fi # if build a consensus for all variants

    # long shot doesn not seem to have this problem, still, just in case
    if [[ "$strErr" == *"Applied 0 variants"* ]]; then
        if [[ "$numZeroVar" -lt 1 ]]; then
            numZeroVar="$((numZeroVar + 1))";
        else
            rm "variant-$intVariant.fasta"; # remove duplicate consensus
            break; # no longer need to continue loop
        fi # check if just repeating the same variant
    fi # if there were no variations, make sure not infinite loop

    printf \
        "%s in consensus for variant %i\n" \
        "$strErr" \
	"$intVariant";

    intVariant=$((intVariant + 1));
done # while there are variatns to apply

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-7: Clean up
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# remove extra files
rm \
    -r \
    "$outVcfStr.vcf.gz" \
    "$outVcfStr.vcf.gz.tbi";

exit;
