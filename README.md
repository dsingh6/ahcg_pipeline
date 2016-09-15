# Variant calling pipeline

## Install Virtualbox and Basespace server
- https://www.virtualbox.org/wiki/Downloads
- https://da1s119xsxmu0.cloudfront.net/sites/developer/native/nativeappsvm/BaseSpace%20Native%20App%20VM%20(phix%20only)%20v9.ova
- https://developer.basespace.illumina.com/docs/content/documentation/native-apps/setup-dev-environment

## Virtualbox setup
- Login information:
	```{sh}
	username: vagrant
	password: vagrant
	```	

## SSH onto virtual machine (putty)
- Login information:
	```{sh}
	Hostname: vagrant@localhost port: 2222
	password: vagrant
	```	

## Program Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)

## Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

## Test data

Use the following protocol to download and prepare test dataset from NIST sample NA12878

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

## Help

To access help use the following command:

```{sh}
python3 ahcg_pipeline.py -h
```


## Set up pipeline dependencies
- download reference genome and snp db
	```{sh}
	wget www.prism.gatech.edu/~sravishankar9/resources.tar.gz
	tar -zxvf ./resources.tar.gz
	gunzip ./resources/dbsnp/dbsnp_138.hg19.vcf.gz
	```
	
- Build bowtie index in working directory
	```{sh}
	./lib/bowtie2-2.2.9/bowtie2-build ./resources/genome/hg19.fa hg19
	```
	
- Build fasta index in ./resources/genomes/ directory
	```{sh}
	sudo apt-get update
	sudo apt-get install samtools
	samtools faidx ./resources/genome/hg19.fa
	```
	
- Build genome dict file ./resources/genomes/ directory
	```{sh}
	java -jar ./lib/picard.jar CreateSequenceDictionary R=./hg19.fa O=./hg19.dict
	```
	
## Specify input parameter paths
	TrimPath = ./lib/Trimmomatic-0.36/trimmomatic-0.36.jar
	BowPath = ./lib/bowtie2-2.2.9/bowtie2
	PicPath = ./lib/picard.jar
	GATKPath = ./lib/GenomeAnalysisTK.jar 
	inPath = NIST7035_TAAGGCGA_L001_R1_001.fastq NIST7035_TAAGGCGA_L001_R2_001.fastq
	indexPath = hg19
	DBSNPpath = ./resources/dbsnp/dbsnp_138.hg19.vcf
	refPath = ./resources/genome/hg19.fa
	adapterPath  = ./lib/Trimmomatic-0.36/adapters/TruSeq2-PE.fa
	output = ./output/

## Run pipeline
```{sh}
python ahcg_pipeline.py -t ./lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b ./lib/bowtie2-2.2.9/bowtie2 -p ./lib/picard.jar -g ./lib/GenomeAnalysisTK.jar -i NIST7035_TAAGGCGA_L001_R1_001.fastq NIST7035_TAAGGCGA_L001_R2_001.fastq -w hg19 -d ./resources/dbsnp/dbsnp_138.hg19.vcf -r ./resources/genome/hg19.fa -a ./lib/Trimmomatic-0.36/adapters/TruSeq2-PE.fa -o ./output/
```

## Set up github
- Change remote path
	```{sh}
	vim .git/config
	url = https://github.com/dsingh6/ahcg_pipeline
	```

- Add files to ignore to gitignore
	```{sh}
	vim .gitignore
	add directories to be ignored
	```

- Add file to git repository
	```{sh}
	git add <filename>
	git config --global user.email dsingh6@gatech.edu
 	git commit -m <"commit message">
 	git push origin master
	```

## Extract sequences for the gene of interest: BRCA1
- download gene coordinates file for hg19
	```{sh}
	wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt
	```

- download bedtools
	```{sh}
	$ wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
	$ tar -zxvf bedtools-2.25.0.tar.gz
	$ cd bedtools2
	$ make
	```

- Find BRCA1 gene in gene coordinate file
	```{sh}
	grep BRCA1 hg19_refGene.txt
	```

- Select BRCA1 variant
	```{sh}
	NM_007294
	https://dnasu.org/DNASU/AdvancedSearchOptions.do
	```
	
- Write script to create bed file from extracted variant
	```{sh}
	python pullCoordinates.py hg19_refGene.txt allExomes.bed
	grep NM_007294 allExomes.bed > exomes007294.bed
	```
	
- extract sequences from original reference and bed file
	```{sh}
	bedtools getfasta -s -fi ./resources/genome/hg19.fa -bed exomes007294.bed -fo nm_007294.out
	```

## Set up to call variants on NIST's NA12878
- Find the .bam file for NIST NA12878 genome
	```{sh}
	wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam
	wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam
	wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_1_NA12878.bwa.markDuplicates.bam
	wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_2_NA12878.bwa.markDuplicates.bam
	```

- merge BAM files
	```{sh}
	samtools merge <merged.bam> <in.1.bam in.2.bam ...>
	```

- extract region of interest from whole bam file (whole bam -> extracted bam) using BED file
	```{sh}
	samtools view <BAM in> -L <input.bed> -b -o <output BAM>
	```
	
- Convert BAM file to FASTQ file for region of interest
	```{sh}
	bedtools bamtofastq [OPTIONS] -i <BAM> -fq <FASTQ> -fq1 <READ2>
	```
## Download gold standard vcf files for NA12878
	```{sh}
	Genome in a Bottle Lastest Release
	wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.gz
	illumina
	ftp://ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz
	Genome in a Bottle
	wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST.hc.snps.indels.vcf
	```
## Find the variants in the generated file 
	```{sh}
	tabix method (finds ALL variants in region):
	bgzip <file.vcf>
	tabix -p vcf <file.vcf>
	tabix <file.vcf.gz> chr17:41195311-41278500
	Script method (finds variants only in the exome region):
	extractVariants.py <infile> > <outfile>
	```
## Find the variants in the correct locations (exome) of the whole genome vcf file
	```{sh}
	extractVariants.py <infile> > <outfile>
	```
