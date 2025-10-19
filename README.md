# 1000 Genomes in silico allo-HCT simulation

Modified_mHag is a WDL workflow designed to identify minor histocompatibility antigens (mHAgs) from paired donor-recipient sequencing data. 

The WGS bam files for the 1000 Genomes Project (1000G), mapped to the GRCH37 reference genome, were obtained from the Google Brain Genomics repository. Suitable D-R pairs were selected based on the HLA typing information available for each individual in the 1000G dataset. The criterion for selection was that greater or equal to 5 alleles in the HLA class I genes (HLA-A, HLA-B and HLA-C) had to be matched between donor and recipient. 

WGS bam files from the selected individuals were reduced to include only the coding region of the genes included in the AML and hematopoietic filters. Variant calling was performed on the reduced bam files using DeepVariant (v1.1.0), generating germline variant call files (VCFs). Allele frequency of each GvL SNP was calculated for the 1000G cohort and compared to the gnomAD allele frequency which was incorporated into the Funcotator task. 
