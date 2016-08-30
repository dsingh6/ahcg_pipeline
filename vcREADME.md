# Variant calling pipeline

## Virtualbox setup
- Login information:
	username: vagrant
	password: vagrant

## SSH onto virtual machine (putty)
- Login information:
	Hostname: vagrant@localhost port: 2222
	password: vagrant

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
	.git/config
	```

- Add files to ignore to gitignore
	.gitignore

- Add file to git repository 
	git add vcREADME.txt
	git commit -m "ReadMe file describing the pipeline"
	git config --global user.email dsingh6@gatech.edu
 	git commit -m "ReadMe file describing the pipeline"
 	git push origin master
