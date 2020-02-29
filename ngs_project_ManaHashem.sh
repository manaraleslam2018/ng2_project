


###########identify the most common germline ERBB-family SNPs in HER-2 positive BC patients by via high depth NGS ########





##******************Download Sampl Fastq file********

#SRX4211716: DNA-seq of Adult Female HER2 Breast Cancer
#(Library:
#Name: SOL6717_R1
#Instrument: Illumina MiSeq
#Strategy: WGS
#Source: GENOMIC
#Selection: other
#Layout: PAIRED)


mkdir -p ~/BC_sample/fqData && cd ~/BC_sample/fqData


#Download sample (1 ILLUMINA (Illumina MiSeq) run: 1.2M spots, 183.6M bases, 72.4Mb downloads)  
Wget http://127.0.0.1:5001/files/BC_sample/fqData/SRR7309325?_xsrf=2%7C83f45831%7C0ff84e8a76f41efa570bfad14bae4d5c%7C1582896007


#split the two reads from the paired end SRA sample 



#This command converts the interlaced fastq file into 8-column tsv file, cuts columns 1-4 (read 1 lines), changes from tsv to fastq format (by replacing tabs with newlines) and redirects the output to read1.fq. In the same STDOUT stream (for speed), using tee, it cuts columns 5-8 (read 2 lines), etc, and redirects the output to read2.fq. 


paste - - - - - - - - < sra_data.fixed.fastq | tee >(cut -f 1-4 | tr "\t" "\n" > SRR7309325_1.fastq)  | cut -f 5-8 | tr "\t" "\n" > SRR7309325_2.fastq  



##Rename the sample file name

mv  SRR7309325_1.fastq  SRR7309325_CAGATC_L002_R1.fastq
mv  SRR7309325_2.fastq  SRR7309325_CAGATC_L002_R2.fastq

****************QUALITY CONTROL*********************** 

mkdir  -p ~/BC_sample/FASTQC && cd ~/BC_sample/FASTQC
conda install -c bioconda fastqc 
conda install -c bioconda multiqc 

for f in ~/BC_sample/fqData/SRR7309325_CAGATC_L002_*.fastq;do fastqc -t 1 -f fastq -noextract $f;done





*********************Alignment*************


mkdir -p ~/BC_sample/fqData && cd ~/BC_sample/fqData

#Download reference file 


mkdir -p ~/workdir/sample_data && cd ~/workdir/sample_data 
(Download  ch2 ,ch7 ,ch12,ch17 from NCBI) 

#combine the sequence of 4 chromosomes 

awk '/^>/ {if (c++ == 0) {print;}next}/^$/ {next}{printf "%s", $0}END {print ""}' 'NC_000002.12[1..242193529].fa' 'NC_000007.14[54999756..55230889].fa' 'NC_000012.12[4783517..47911381].fa' 'NC_000017.11[7666487..7689465].fa' > chrom2_7_12_17.fa 


conda install -c bioconda -y bwa

#index the genome

mkdir -p ~/BC_sample/BWA_align/bwaIndex 
cd ~/BC_sample/BWA_align/bwaIndex
ln -s ~/BC_sample/fqData/ chrom2_7_12_17.fa.
bwa index -a bwtsw chrom2_7_12_17.fa 





#sequence alignment with read group information


cd ~/BC_sample/BWA_align

 

for R1 in ~/BC_sample/fqData/*R1.fastq.gz;do
    SM=$(basename $R1 | cut -d"_" -f1)                                          
    LB=$(basename $R1 | cut -d"_" -f1,2)                                        
    PL="Illumina"                                      
    RGID=$(zcat $R1 | head -n1 | sed 's/ /_/g' |cut -d "_" -f1)                 
    PU=$RGID.$LB                                                                
    echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"

    R2=$(echo $R1 | sed 's/_R1./_R2./')
    echo $R1 $R2
    index=”/home/BC_sample/bwa_align/bwaIndex/chrom2_7_12_17.fa”
    bwa mem -t 4 -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" $index $R1 $R2 > $(basename $R1 _R1_001.pe.fq.gz).sam
done



# statistics of the alligment  
samtools flagstat SRR7309325_CAGATC_L002_R1.fastq.gz.sam> bwa_alignment_sample_stats.out 



#Generate bam file 

samtools view -hbo SRR7309325_CAGATC_L002_R1.fq.bam SRR7309325_CAGATC_L002_R1.fastq.gz.sam  

samtools sort SRR7309325_CAGATC_L002_R1.fq.bam -o SRR7309325_CAGATC_L002_R1.fq.sorted.bam 



####################mapping QC############## 
  samtools depth SRR7309325_CAGATC_L002_R1.fq.sorted.bam  | awk '{{sum+=$3}} END {{print "Average = ",sum/NR, "No of covered Nuc = ", NR}}' > SRR7309325_CAGATC_L002_R1.fq.sorted.cov 

  samtools flagstat SRR7309325_CAGATC_L002_R1.fq.sorted.bam > SRR7309325_CAGATC_L002_R1.fq.sorted.stat






################### Mark duplicate ############## 
conda install -c bioconda picard 
picard_path=$CONDA_PREFIX/share/picard-* ## 2.21.7-0


java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=SRR7309325_CAGATC_L002_R1.fq.sorted.bam OUTPUT=SRR7309325_CAGATC_L002_R1.fq.sorted.dedup.bam METRICS_FILE=SRR7309325_CAGATC_L002_R1.fq.sorted.bam.metrics.txt

  

################ indexing to gtak variant calling######## 
#sample 
java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT= SRR7309325_CAGATC_L002_R1.fq.sorted.dedup.bam


#reference
ln -s /home/ngs/BC_sample/BWA_align/bwaIndex/chrom2_7_12_17.fa .
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R= chrom2_7_12_17.fa O=chrom2_7_12_17.fa.dict
samtools faidx chrom2_7_12_17.fa 



################Download known varinats ################## 

wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr2.vcf.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr7.vcf.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr12.vcf.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr17.vcf.gz  

##concatenat to one vcf file and index this file 

grep "^#" homo_sapiens-chr2.vcf > merged_variants.vcf 
grep "^2" homo_sapiens-chr2.vcf >> merged_variants.vcf
grep "^7" homo_sapiens-chr7.vcf >> merged_variants.vcf 
grep "^12" homo_sapiens-chr12.vcf >> merged_variants.vcf
grep "^17" homo_sapiens-chr17.vcf >> merged_variants.vcf  

gatk IndexFeatureFile -I merged_variants.vcf  









