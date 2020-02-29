# ng2_project

ABSTRACT

Background
Breast cancer (BC) is the most common malignancy in women with over 25% of all cancers being diagnosed as BC in 2018. HER2-positive BC accounts for around 20% of all human BCs and HER2 overexpression is associated with poor prognosis and an aggressive phenotype. Trastuzumab is a monoclonal antibody targeted to HER2, has well established efficacy in the treatment of HER2-positive BC. However, a significant proportion of patients with the disease have tumors that initially do not respond or that acquire resistance to trastuzumab after an initial period of response . ERBB-family genes which encode the HER family of proteins EGFR, HER2, HER3 and HER4 are commonly studied in HER2-positive BC, and some studies have identified the role of HER2 SNPs in response to trastuzumab .

Aim :
Here we want to identify the most common germline ERBB-family SNPs in HER-2 positive BC patients by via high depth NGS.

Methods (Sequencing analysis) :
1. Download two samples whole exome sequencing of Irish HER2+ breast cancer patients from SRA (paired end)
   Library: Name: SOL6717_R1 (The custom library designed comprised of 132 regions)
2. Check QC and Trimming in case of bad quality
3. Reads will be aligned with BWA.
4. Duplicate reads will be marked by Picard tools
5. local realignment and base recalibration will be conducted with GATK
6. Joint variant calling using HaplotypeCaller
7- Split SNPs and indels
8- Assess the different filters in both known and novel to decide the threshould for the filteration in each filter 
9-SNP Variant filteration
10- separate the passed SNPs in a vcf file
