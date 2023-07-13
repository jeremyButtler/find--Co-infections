#!/usr/bin/bash

################################################################################
# Use:
#    Does varint calling on a fastq file using Clair3 & reference
# Input
#    $1:
#        - Reads to do variant calling on    [Required: fastq]
#    $2:
#        - Reference to do variant calling with [Required: fasta]
# Output:
#    Files: variant-number.fasta holding cosenus for each variant
#           Clair3 detected
# Required:
#    samtools
#    minimap2
#    bcftools
#    Clair3 [hkubal/clair3:latest] (by docker), note mispelling
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    sec-1: Variable declerations
#    sec-3: Prepare data for variant calling
#    sec-4: Do variant calling with Clair3
#    sec-5: Build a consensus for each variant
#    sec-6: Final cleanup
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

readsFastqStr="$1";
refForVarCall="$2";

outConsensusNameStr="variant";
intVariant=1;      # will hold number variants clair found
numZeroVar=0; # number times bcftools applied 0 varaints (prevent infinit loop)

inDir="$(pwd)"; # location for docker to work out

strErr=""; # hold error message from bcftools consensus (so know when to quite)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Prepare data for variant calling
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# make sorted bam file with reads mapped to input reference
minimap2 \
    -ax map-ont \
    "$refForVarCall" \
    "$readsFastqStr" |
  samtools sort |
  samtools view \
    -b \
  > "tmp.bam";

# index the bam file
samtools index "tmp.bam";

# index the reference used to make the bam file
samtools faidx "$refForVarCall";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-4: Do variant calling with Clair3
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

docker \
    run \
    -it \
    --user $(id -u) \
    -v "$inDir:$inDir" \
    -v "$inDir:$inDir" \
    "hkubal/clair3:latest" \
    "/opt/bin/run_clair3.sh" \
    --bam_fn="$inDir/tmp.bam" \
    --ref_fn="$inDir/$refForVarCall" \
    --threads=3 \
    --platform="ont" \
    --output="$inDir" \
    --no_phasing_for_fa \
    --model_path="/opt/models/r941_prom_hac_g238" \
    --include_all_ctgs;
    #--haploid_precise; # only detects one haplotype
    #-u $(whoami) # runs docker as current user, instead of root, however user needs to be docker
    #-u $(id -u):$(id -g) # runs user id, but not system user, so does not have write permision
rm \
    "tmp.bam" \
    "tmp.bam.bai"; # no longer need the bam files

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-5: Build a consensus for each variant
#    sec-5 sub-1: building a consensus when only one variant & is exact match
#    sec-5 sub-2: building a consensus when one variant not exact match
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-5 Sub-1: building a consensus when only one variant is exact match
#*******************************************************************************

# unzip vcf to check if can keep
bgzip \
    -d \
    -k \
    "merge_output.vcf.gz";

# Get last line of file (check if header)
strErr="$(
    sed \
        -n \
        '$s/\t.*//p;' \
        "merge_output.vcf" \
)"; # Get the last line of the file

if [[ "$strErr" == *"CHROM"* ]]; then  # check if last lin was header
# if bulding a consensus for a single variant & vcf file has no variants
    bcftools \
        consensus \
        -f "$refForVarCall" \
        -H "1" \
        "$outVcfStr.vcf.gz" \
      > "variant-$intVariant.fasta"; \

    printf "Applied 0 variants for only consensus\n";
fi # if bulding a consensus for a single variant & vcf file has no variants

#*******************************************************************************
# Sec-6 Sub-2: building a consensus when one variant not exact match
#*******************************************************************************

while [[ "$strErr" != *"CHROM"* ]]; do
# while there are variatns to apply
    strErr="$( \
        { \
            bcftools \
                consensus \
                -f "$refForVarCall" \
                -H "$intVariant" \
                "merge_output.vcf.gz" \
              > "$outConsensusNameStr-$intVariant.fasta"; \
        } \
        2>&1\
    )"; # build consensuses and get error when no more variants in vcf
    # Need the "{ command > file.str; } 2>&1" to capture the error message
    
    # check if bcftools errored out (no more variants)
    if [[ "$strErr" == *"Can"*  || "$strErr" == *"Failed"* ]]; then
        rm "variant-$intVariant.fasta";
        break; # no longer need to continue loop
    fi # if build a consensus for all variants

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
# Sec-6: Final cleanup
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

rm \
    -r \
    "full_alignment.vcf.gz" \
    "full_alignment.vcf.gz.tbi" \
    "merge_output.vcf" \
    "merge_output.vcf.gz" \
    "merge_output.vcf.gz.tbi" \
    "pileup.vcf.gz" \
    "pileup.vcf.gz.tbi" \
    "run_clair3.log" \
    "log" \
    "tmp" \
    "$refForVarCall.fai"; # files clair made
