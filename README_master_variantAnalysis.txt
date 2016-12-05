Master script to run the variant calling pipeline, ahcg_pipeline.py, and perform a variant coverage analysis

Usage: ./master_variantAnalysis.sh [-h] -r hg19_reference -g gene_list -b BAM_file -v goldStandard_variants -c csv_file -o output_directory

Input Arguments:
    -h            Display this help and exit
    -r            hg19 reference file - MANDATORY
    -g            List of clinically relevant genes and NM numbers - MANDATORY
    -b            Merged BAM file - MANDATORY
    -v            Gold standard variants vcf file
    -o            output directory - MANDATORY
    -c            Use flag to remove all intermediate files

Outputs:
- the variants (final_variants) and coverage data (final_cov.bed) can be found in the specified output directory
