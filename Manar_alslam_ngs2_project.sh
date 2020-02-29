1- download data 
'''
mkdir /home/manar/ngs2_project/Breast_cancer_samples && cd /home/manar/ngs2_project/Breast_cancer_samples

wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR730/SRR7309332/SRR7309332.sra
wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR730/SRR7309338/SRR7309338.sra
fastq-dump --outdir . --gzip --split-3 SRR7309332.sra
fastq-dump --outdir . --gzip --split-3 SRR7309338.sra
'''

2- download reference:
 
download chr2, 7, 12, 17 after that concatenate them in one fasta file (https://www.biostars.org/p/256796/) (https://unix.stackexchange.com/questions/158941/how-to-combine-gunzipped-fastq-files)

'''
wget -c ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.7.fa.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.12.fa.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.17.fa.gz

gunzip -k *.fa.gz 

cat Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa Homo_sapiens.GRCh38.dna_sm.chromosome.7.fa Homo_sapiens.GRCh38.dna_sm.chromosome.12.fa Homo_sapiens.GRCh38.dna_sm.chromosome.17.fa > Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa
'''


3- check the quality of data

'''
mkdir ~/ngs2_project/Breast_cancer_samples/FASTQC_sample1_2 && cd ~/ngs2_project/Breast_cancer_samples/FASTQC_sample1_2
cp  /home/manar/ngs2_project/Breast_cancer_samples/SRR7309332_1.fastq.gz .
cp  /home/manar/ngs2_project/Breast_cancer_samples/SRR7309332_2.fastq.gz .

cp /home/manar/ngs2_project/Breast_cancer_samples/SRR7309338_1.fastq.gz .
cp /home/manar/ngs2_project/Breast_cancer_samples/SRR7309338_2.fastq.gz .


for f in ~/ngs2_project/Breast_cancer_samples/FASTQC_sample1_2/*.fastq.gz;do fastqc -t 1 -f fastq -noextract $f;done
'''

3-BWA aligment 
 
a) indixing the the genome refernce 

'''
mkdir -p ~/ngs2_project/bwa_align/bwaIndex && cd ~/ngs2_project/bwa_align/bwaIndex
ln -s /home/manar/ngs2_project/Breast_cancer_samples/Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa .
bwa index -a bwtsw Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa
'''

#aligment trial to check if the chromosoms are correct as a reference 
'''
#cd ~/ngs2_project/bwa_align/

#R1="/home/manar/ngs2_project/Breast_cancer_samples/SRR7309332_1.fastq.gz"
#R2="/home/manar/ngs2_project/Breast_cancer_samples/SRR7309332_2.fastq.gz"

#/usr/bin/time -v bwa mem bwaIndex/Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa $R1 $R2 > SRR7309332.sam

#samtools flagstat /home/manar/ngs2_project/bwa_align/SRR7309332.sam > bwa_alignment_sample${r}_stats.out

rename the sanples to make alligment 
'''
B) renam the sample name to be seutable to gtak 
'''
mv SRR7309332_1.fastq.gz SRR7309332_TGACCA_L001_R1.fastq.gz
mv SRR7309332_2.fastq.gz SRR7309332_TGACCA_L001_R2.fastq.gz
mv SRR7309338_1.fastq.gz SRR7309338_TAGCTT_L002_R1.fastq.gz   
mv SRR7309338_2.fastq.gz SRR7309338_TAGCTT_L002_R2.fastq.gz
'''
C) Alligment 
'''
mkdir -p ~/ngs2_project/GATK_varient_calling && cd ~/ngs2_project/GATK_varient_calling
for R1 in ~/ngs2_project/Breast_cancer_samples/*_R1.fastq.gz;do
    SM=$(basename $R1 | cut -d"_" -f1)                                          ##sample ID
    LB=$(basename $R1 | cut -d"_" -f1,2)                                        ##library ID
    PL="Illumina"                                                               ##platform (e.g. illumina, solid)
    RGID=$(zcat $R1 | head -n1 | sed 's/ /_/g' |cut -d "_" -f1)                 ##read group identifier 
    PU=$RGID.$LB                                                                ##Platform Unit
    echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"

    R2=$(echo $R1 | sed 's/_R1./_R2./')
    echo $R1 $R2
    index="/home/manar/ngs2_project/bwa_align/bwaIndex/Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa"
    bwa mem -t 4 -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" $index $R1 $R2 > $(basename $R1 _R1_001.pe.fq.gz).sam
done
'''
#########from this step I worked each sample separately without for loop as the lab didnt work for a long time after aligment####
D) statistics to sam file after aligment 
'''
samtools flagstat /home/manar/ngs2_project/GATK_varient_calling/SRR7309332_TGACCA_L001_R1.fastq.gz.sam > SRR7309332_TGACCA_L001_sam_stats.out
samtools flagstat /home/manar/ngs2_project/GATK_varient_calling/SRR7309338_TAGCTT_L002_R1.fastq.gz.sam > SRR7309338_TAGCTT_L002_sam_stats.out
'''
E) rename file name 
'''
mv SRR7309332_TGACCA_L001_R1.fastq.gz.sam SRR7309332_TGACCA_L001.sam
mv SRR7309338_TAGCTT_L002_R1.fastq.gz.sam SRR7309338_TAGCTT_L002.sam
'''


4) generate & sort BAM file
'''
samtools view -hbo SRR7309332_TGACCA_L001.bam SRR7309332_TGACCA_L001.sam
samtools sort SRR7309332_TGACCA_L001.bam -o SRR7309332_TGACCA_L001.sorted.bam
samtools view -hbo SRR7309338_TAGCTT_L002.bam SRR7309338_TAGCTT_L002.sam
samtools sort SRR7309338_TAGCTT_L002.bam -o SRR7309338_TAGCTT_L002.sorted.bam
'''

5) mark deduplicates
A) sample 1
'''
picard_path=$CONDA_PREFIX/share/picard-* ## 2.21.7-0
java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=SRR7309332_TGACCA_L001.sorted.bam OUTPUT= SRR7309332_TGACCA_L001.sorted.dedup.bam METRICS_FILE= SRR7309332_TGACCA_L001.sorted.bam.metrics.txt
samtools flagstat SRR7309332_TGACCA_L001.sorted.dedup.bam > SRR7309332_TGACCA_L001.sorted.dedup.stat
'''
A) sample 2
'''
java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=SRR7309338_TAGCTT_L002.sorted.bam OUTPUT=SRR7309338_TAGCTT_L002.sorted.dedup.bam METRICS_FILE= SRR7309332_TGACCA_L001.sorted.bam.metrics.txt
samtools flagstat SRR7309338_TAGCTT_L002.sorted.dedup.bam > SRR7309338_TAGCTT_L002.sorted.dedup.stat
'''

6) indexing to gtak variant calling 
A) sample 1
'''
  java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR7309332_TGACCA_L001.sorted.dedup.bam
'''
B) sample 2
'''
  java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR7309338_TAGCTT_L002.sorted.dedup.bam
'''
b) reference 
'''
ln -s ~/ngs2_project/bwa_align/bwaIndex/Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa O=Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.dict
samtools faidx Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa
'''
7) Download known varinats & concatenat them in  a onr filr and index this file: 
'''
wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr2.vcf.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr7.vcf.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr12.vcf.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr17.vcf.gz
gunzip -k *.vcf.gz

###bcftools concat --output homo_sapiens-chr2_7_12_17.vcf homo_sapiens-chr2.vcf homo_sapiens-chr7.vcf homo_sapiens-chr12.vcf homo_sapiens-chr17.vcf 

grep "^#" homo_sapiens-chr2.vcf > homo_sapiens-chr2_7_12_17.vcf
grep "^2" homo_sapiens-chr2.vcf >> homo_sapiens-chr2_7_12_17.vcf
grep "^7" homo_sapiens-chr7.vcf >> homo_sapiens-chr2_7_12_17.vcf
grep "^12" homo_sapiens-chr12.vcf >> homo_sapiens-chr2_7_12_17.vcf
grep "^17" homo_sapiens-chr17.vcf >> homo_sapiens-chr2_7_12_17.vcf

gatk IndexFeatureFile -I homo_sapiens-chr2_7_12_17.vcf
'''
8) Recalibration 
A) Sample 1
'''
 gatk --java-options "-Xmx2G" BaseRecalibrator \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa -I SRR7309332_TGACCA_L001.sorted.dedup.bam --known-sites homo_sapiens-chr2_7_12_17.vcf \
-O SRR7309332_TGACCA_L001.sorted.report

 gatk --java-options "-Xmx2G" ApplyBQSR \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa -I SRR7309332_TGACCA_L001.sorted.dedup.bam -bqsr SRR7309332_TGACCA_L001.sorted.report \
-O SRR7309332_TGACCA_L001.sorted.bqsr.bam --add-output-sam-program-record --emit-original-quals
'''
B) Sample 2
'''
gatk --java-options "-Xmx2G" BaseRecalibrator \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa -I SRR7309338_TAGCTT_L002.sorted.dedup.bam --known-sites homo_sapiens-chr2_7_12_17.vcf \
-O SRR7309338_TAGCTT_L002.sorted.report

 gatk --java-options "-Xmx2G" ApplyBQSR \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa -I SRR7309338_TAGCTT_L002.sorted.dedup.bam -bqsr SRR7309338_TAGCTT_L002.sorted.report \
-O SRR7309338_TAGCTT_L002.sorted.dedup.bqsr.bam --add-output-sam-program-record --emit-original-quals
'''

9) Joint variant calling using HaplotypeCaller

Call germline SNPs and indels via local re-assembly of haplotypes
A) assess genotype likelihood per-sample
A.1) sample 1
'''  
  gatk --java-options "-Xmx2G" HaplotypeCaller \
  -R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa -I SRR7309332_TGACCA_L001.sorted.bqsr.bam \
  --emit-ref-confidence GVCF \
  --pcr-indel-model NONE \
  -O SRR7309332_TGACCA_L001.sorted.gvcf
'''
A.2) sample 2
'''
  gatk --java-options "-Xmx2G" HaplotypeCaller \
  -R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa -I SRR7309338_TAGCTT_L002.sorted.dedup.bqsr.bam \
  --emit-ref-confidence GVCF \
  --pcr-indel-model NONE \
  -O SRR7309338_TAGCTT_L002.sorted.gvcf
'''
B) combine samples
'''
gatk --java-options "-Xmx2G" CombineGVCFs \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa \
-V SRR7309332_TGACCA_L001.sorted.gvcf \
-V SRR7309338_TAGCTT_L002.sorted.gvcf \
-O raw_variants.gvcf
'''

C) Joint Genotyping
'''
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa \
-V raw_variants.gvcf \
--max-alternate-alleles 3 \
-O raw_variants.vcf
'''
D) annotated output
'''
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa \
-V raw_variants.gvcf \
--max-alternate-alleles 3 \
--dbsnp homo_sapiens-chr2_7_12_17.vcf \
-O raw_variants_ann.vcf
'''

10) Split SNPs and indels
'''
gatk --java-options "-Xmx2G" SelectVariants \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa \
-V raw_variants_ann.vcf \
--select-type-to-include SNP \
-O raw_variants_ann_SNP.vcf

gatk --java-options "-Xmx2G" SelectVariants \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa \
-V raw_variants_ann.vcf \
--select-type-to-include INDEL \
-O raw_variants_ann_INDEL.vcf
'''
11) VCF statitics
'''
bgzip -c raw_variants_ann.vcf > raw_variants_ann.vcf.gz
tabix -p vcf raw_variants_ann.vcf.gz
rtg vcfstats raw_variants_ann.vcf.gz > stats.txt
'''
11) Assess the different filters in both known and novel

'''
for var in "SNP" "INDEL";do
 input="raw_variants_ann_"$var".vcf"
 for filter in "QD" "MQ" "MQRankSum" "FS" "SOR" "ReadPosRankSum" "AN" "DP" "InbreedingCoeff";do
  filterValues=$var.$filter
  awk -v k="$filter=" '!/#/{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~"^"k) {sub(k,$3" ",a[i]); print a[i]}}' $input > $filterValues
  grep -v "^\." $filterValues > known.$var.$filter
  grep "^\." $filterValues > novel.$var.$filter
done;

mkdir filters && cd filters
mv ../{*.SNP.*,SNP.*,*.INDEL.*,INDEL.*} .

wget https://raw.githubusercontent.com/dib-lab/dogSeq/master/scripts/densityCurves.R
Rscript -e "install.packages('ggplot2', contriburl=contrib.url('http://cran.r-project.org/'))"
for f in SNP.* INDEL.*;do
  Rscript densityCurves.R "$f"
done
'''

12) SNP Variant filteration

'''
cd ~/ngs2_project/GATK_varient_calling
gatk --java-options "-Xmx2G" VariantFiltration \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa \
-V raw_variants_ann_SNP.vcf \
--filter-name "snpQD" \
--filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
--filter-name "snpMQ" \
--filter-expression "vc.hasAttribute('MQ') && MQ < 40.0" \
--filter-name "snpMQRankSum" \
--filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
--filter-name "snpFS" \
--filter-expression "vc.hasAttribute('FS') && FS > 60.0" \
--filter-name "snpSOR" \
--filter-expression "vc.hasAttribute('SOR') && SOR > 4.0" \
--filter-name "snpReadPosRankSum" \
--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
--filter-name "snpDP" \
--filter-expression "vc.hasAttribute('DP') && DP > 3105" \
-O raw_variants_ann_SNP_clean.vcf
'''
13) extract all passed record in a new vcf

'''
vcftools --vcf raw_variants_ann_SNP_clean.vcf --remove-filtered-all --recode --recode-INFO-all --stdout > raw_variants_ann_SNP_clean_all_passed.vcf
'''
14) VCF statitics to the passed file
'''
bgzip -c raw_variants_ann_SNP_clean_all_passed.vcf > raw_variants_ann_SNP_clean_all_passed.vcf.gz
tabix -p vcf raw_variants_ann_SNP_clean_all_passed.vcf.gz

rtg vcfstats raw_variants_ann_SNP_clean_all_passed.vcf.gz > raw_variants_ann_SNP_clean_all_passed.stats.txt
'''
15) extract known snps for further analysis in a text file:
'''
grep -v "^#" raw_variants_ann_SNP_clean_all_passed.vcf | awk '{print $3}' | grep "^rs" > raw_variants_ann_SNP_clean_all_passed.stats.txt
'''

############## I tried to download metadata and make functional annotation using FilterFuncotations in Gatk howewer;
there is no metadata to SRA data & no enogh memory in lap & no enough time######
download metadata
'''
sudo apt install ncbi-entrez-direct
esearch -db sra -query SRR7309332 | efetch -format runinfo
'''
functional annotation : For germline data sources: (https://gatkforums.broadinstitute.org/gatk/discussion/11193/funcotator-information-and-tutorial)
'''
./gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download
'''
##grep "^#"  raw_variants_ann_SNP_clean.vcf > raw_variants_ann_SNP_clean_passed_filter.vcf
##grep -v "^#" raw_variants_ann_SNP_clean.vcf | awk '{if($7="PASS")print $0}' >> raw_variants_ann_SNP_clean_passed_filter.vcf

##grep -v "^#" raw_variants_ann_SNP_clean_all_passed.vcf | awk'{if($3="rs")print $0}'
##grep -v "^#" raw_variants_ann_SNP_clean_all_passed.vcf | awk '{print $3}' | grep "^rs" | wc -l
