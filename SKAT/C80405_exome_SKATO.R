# Run gene-level association test on rare (MAF <0.03) variants

# Dependencies ***WILL UPDATE WITH SCRIPTS USED TO CREATE FILES***
# SKAT package (https://CRAN.R-project.org/package=SKAT)
# Binary PLINK files ('.bed', '.bim') containing genetic data converted from .vcf file using PLINK2
# SetID file ('.setID') which indicates which variants belong to which gene
# Phenotype file ('.fam')
# Covariate file ('.cov')

# Load library
library(SKAT)

# Set file names
File.Bed <- "C80405_exome_candidate.bed"
File.Bim <- "C80405_exome_candidate.bim"
File.SetID <- "C80405_exome_candidate_rare03.setID"
File.Fam <- "C80405_Caucasian.fam"
File.Cov <- "C80405_Caucasian.cov"

# Read in phenotype and covariate files
FAM_Cov <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary = TRUE, flag1 = 0, cov_header = TRUE)

# Set variable names
y = FAM_Cov$Phenotype
sex = FAM_Cov$Sex
age = FAM_Cov$AGE
bmi = FAM_Cov$BMI_CAT25
prhpt = FAM_Cov$PRHPT
diabetes = FAM_Cov$DIABETES

# Create a SNP set ID ('.SSD') file and info ('.info') file
File.SSD <- "C80405_exome_candidate_rare03.SSD"
File.Info <- "C80405_exome_candidate_rare03.info"
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)

# Open .SSD file
SSD.INFO <- Open_SSD(File.SSD, File.Info)

# Create null model based on phenotypes and covariates
null <- SKAT_Null_Model(y ~ sex + age + bmi + prhpt + diabetes, out_type = "D", Adjustment = TRUE)

# Perform SKAT-O
SKAT_data <- SKAT.SSD.All(SSD.INFO, null, method = "optimal.adj")
results <- SKAT_data$results

# Sort results by P-value
results_sorted <- (results[order(results$P.value),])

# Save results
write.table(results_sorted, "output/SKATO_results.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Close .SSD file
Close_SSD()
