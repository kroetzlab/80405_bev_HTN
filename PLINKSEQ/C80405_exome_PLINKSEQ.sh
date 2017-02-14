# Run single-variant association tests on common (MAF >0.10) variants

# Dependencies ***WILL UPDATE WITH SCRIPTS USED TO CREATE FILES***
# PLINK/SEQ library (https://atgu.mgh.harvard.edu/plinkseq/)
# VCF (.vcf) file containing genetic data
# Phenotype (.phe) file containing phenotypes and covariates

# Create new project
pseq C80405_candidate new-project

# Load .vcf data
pseq C80405_candidate load-vcf --vcf C80405_exome_candidate.vcf

# Load phenotype data
pseq C80405_candidate load-pheno --file C80405_pheno.phe

# Perform association test with no covariates (this will nicely summarize the allele counts and homo/hetero-zygosity between cases and controls for each tested variant)
pseq C80405_candidate v-assoc --phenotype casecontrol --mask maf=0.1-1 > output/C80405_PLINKSEQ_candidate_results_vassoc.txt

# Run logistic regression with covariates
pseq C80405_candidate glm --phenotype casecontrol --covar sex age BMI_Cat25 PRHPT Diabetes --mask maf=0.1-1 | sort -k9 -g > output/C80405_PLINKSEQ_candidate_results_glm.txt

 
