#!/bin/bash

show_help()
{
cat <<EOF

Usage: ./${0##*/} [-h] -r hg19_reference -g gene_list -b BAM_file -v goldStandard_variants -o output_directory -c

Input Arguments:
    -h            Display this help and exit
    -r 		  hg19 reference file - MANDATORY
    -g            List of clinically relevant genes and NM numbers - MANDATORY
    -b		  Merged BAM file - MANDATORY
    -v	          Gold standard variants vcf file
    -o		  output directory - MANDATORY
    -c		  Use flag to remove all intermediate files

EOF
}

#Variables
no_Gold_variants=0
clean=0

#Getopts
while getopts "h:r:g:b:v:o:c" opt
do 
	case $opt in
	h)
		show_help;
		exit
		;;
	r)
		hg19_ref=$OPTARG
		;;
	g)
		geneList=$OPTARG
		;;

	b)	
		bam=$OPTARG
		;;
	
	v)
		Gvariants=$OPTARG
		no_Gold_variants=1
		;;
	o)
		output_dir=$OPTARG
		;;
	c)	
		clean=1
		;;
	'?')
		show_help;
		exit
		;;
	:)
		show_help
		exit
		;;
	esac
done

#Mandatory files check
if [[ -z "$hg19_ref" ]]
then
	echo "No hg19 reference file found."
	exit 1
fi
if [[ -z "$geneList" ]]
then
        echo "No gene list found."
        exit 1
fi
if [[ -z "$bam" ]]
then
        echo "No BAM file found."
        exit 1
fi
if [[ -z "$output_dir" ]]
then
        echo "No output directory specified."
        exit 1
fi

echo "Creating FastQ files from BAM file"
bedtools bamtofastq -i $bam -fq read1.fastq -fq2 read2.fastq

echo "Running ahcg_pipeline on generated FastQ files"
python ahcg_pipeline.py -t ./lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b ./lib/bowtie2-2.2.9/bowtie2 -p ./lib/picard.jar -g ./lib/GenomeAnalysisTK.jar -i read1.fastq read2.fastq -w hg19 -d ./resources/dbsnp/dbsnp_138.hg19.vcf -r ./resources/genome/hg19.fa -a ./lib/Trimmomatic-0.36/adapters/TruSeq2-PE.fa -o ./output_temp/

echo "Recalibrating variants"
java -Xmx4g -jar /home/vagrant/ahcg_pipeline/lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R /home/vagrant/ahcg_pipeline/resources/genome/hg19.fa -input /home/vagrant/ahcg_pipeline/output_temp/variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./hapmap_3.3.hg19.sites.vcf.gz -resource:omni,known=false,training=true,truth=false,prior=12.0 ./1000G_omni2.5.hg19.sites.vcf.gz  -resource:1000G,known=false,training=true,truth=false,prior=10.0 ./1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /home/vagrant/ahcg_pipeline/resources/dbsnp/dbsnp_138.hg19.vcf.gz -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -recalFile output.recal -tranchesFile output.tranches

echo "Applying recalibration"
java -jar /home/vagrant/ahcg_pipeline/lib/GenomeAnalysisTK.jar -T ApplyRecalibration -R /home/vagrant/ahcg_pipeline/resources/genome/hg19.fa -input /home/vagrant/ahcg_pipeline/output_temp/variants.vcf --ts_filter_level 99.0 -tranchesFile output.tranches -recalFile output.recal -mode SNP  -o ./recal_variants.vcf

echo "Extracting genes of interest from hg19 reference file"
awk '{print "\\<" $2 "\\>" }' $geneList > nmNumbers.txt
grep -f nmNumbers.txt $hg19_ref > hg19_extracts.txt

echo "Creating BED file for genes of interest"
./BEDmaker.py -i hg19_extracts.txt -o genes.bed

echo "Extracting selected variants from all variants using BED file"
bedtools intersect -header -wa -a ./recal_variants.vcf -b genes.bed > foundVariants.vcf

if [ "$no_Gold_variants" == 1 ]
then
	echo "Comparing found variants to gold standard variants"
	bedtools intersect -header -wa -a $Gvariants -b genes.bed > goldVariants.vcf
	bedtools intersect -header -a foundVariants.vcf -b goldVariants.vcf > final_variants.vcf
	rm foundVariants.vcf
else
	mv foundVariants.vcf final_variants.vcf
fi

echo "Getting variant coverage"
samtools view -L genes.bed $bam -b > new.bam
bedtools genomecov -ibam new.bam -bga > coverage_output.bed
bedtools intersect -loj -split -a genes.bed -b coverage_output.bed >  cov.bed
awk '{printf("%s\t%s\t\%s\t%s\t%s\n", $1,$2,$3,$4,$10,$6)}' cov.bed > final_cov.bed

mv final_variants.vcf $output_dir
mv final_cov.bed $output_dir

# move intermediate files into intermediate directory
mv cov.bed ./output_temp/
mv coverage_output.bed ./output_temp/
mv new.bam ./output_temp/
mv read1.* ./output_temp/
mv read2.* ./output_temp/
mv genes.bed ./output_temp/
mv nmNumbers.txt ./output_temp/
mv hg19_extracts.txt ./output_temp/
mv recal_variants.vcf ./output_temp/
mv output.tranches ./output_temp/
mv output.recal* ./output_temp/

if [ "$clean" == 1 ]
then
	echo "Removing intermediate files"
	rm -r ./output_temp/
fi

echo "Program complete."
