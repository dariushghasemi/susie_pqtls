
#--------------------------------#
#------      Read GWAS     ------#
#--------------------------------#

#seqid <- "seq.3054.3"
#locus <- "16_69449233_73516227"

# GWAS file
#file_pwas <- paste0("/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/glm_model/out/seq.3054.3/16_69449233_73516227.regenie")
file_pwas <- paste0("/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/glm_model/out/seq.12593.33/19_54219624_54427869.regenie")

pwas <- data.table::fread(file_pwas)


#--------------------------------#
#------      Read PGEN     ------#
#--------------------------------#

# File paths
pgen_path <- "/exchange/healthds/pQTL/INTERVAL/Genetic_QC_files/pgen/qc_recoded_harmonised/impute_recoded_selected_sample_filter_hq_var_new_id_alleles_19.pgen"
pvar_path <- "/exchange/healthds/pQTL/INTERVAL/Genetic_QC_files/pgen/qc_recoded_harmonised/impute_recoded_selected_sample_filter_hq_var_new_id_alleles_19.pvar"
psam_path <- "/exchange/healthds/pQTL/INTERVAL/Genetic_QC_files/pgen/qc_recoded_harmonised/impute_recoded_selected_sample_filter_hq_var_new_id_alleles_19.psam"

# Read the PVAR and PSAM files
pvar <- data.table::fread(pvar_path, header = TRUE)
psam <- data.table::fread(psam_path, header = TRUE)

# Initialize PGEN object
pvar_obj <- pgenlibr::NewPvar(pvar_path)
pgen <- pgenlibr::NewPgen(pgen_path, pvar = pvar_obj)

# Extract genotype data from the PGEN file
num_samples <- nrow(psam)
num_variants <- nrow(pvar)
variant_subset <- c(1:num_variants)

# Read genotype matrix
geno_mat <- pgenlibr::ReadList(pgen, variant_subset, meanimpute = FALSE)
if (length(which(apply(geno_mat, 1, function(row) any(is.na(row))))) != 0) {
  stop("NAs are present!")
}

# Convert geno_mat to numeric matrix
geno_mat <- as.matrix(geno_mat)
geno_mat <- apply(geno_mat, 2, as.numeric)

rownames(geno_mat) <- as.character(psam$IID)
colnames(geno_mat) <- as.character(pvar$ID)


#--------------------------------#
#------    In-sample LD    ------#
#--------------------------------#

# compute in-sample LD from the genotype/dosages
# cor_matrix <- dosage %>% 
#   select(- c(IID, FID, PAT, MAT, SEX, PHENOTYPE)) %>%
#   #summarise(across(where(is.numeric), ~ sum(is.na(.x))))
#   cor(use = "pairwise.complete.obs") #"complete.obs"


# convert2matrix <- function(df){
# 
# # Extract unique variant IDs
# ids <- unique(c(df[[1]], df[[2]]))
# 
# # Create an empty matrix
# ld_matrix <- matrix(NA, nrow = length(ids), ncol = length(ids))
# rownames(ld_matrix) <- ids
# colnames(ld_matrix) <- ids
# 
# # Fill the matrix with LD values
# for (i in 1:nrow(df)) {
#   id_a <- df[[i,1]]
#   id_b <- df[[i,2]]
#   ld_matrix[id_a, id_b] <- df[[i,3]]
#   ld_matrix[id_b, id_a] <- df[[i,3]]  # Since the matrix is symmetric
# }
# 
# # Replace NA with 1 for diagonal elements (self correlation)
# diag(ld_matrix) <- 1
# 
# # Print the matrix
# return(ld_matrix)
# }



#--------------------------------#
#------   Pre-computed LD  ------#
#--------------------------------#

# script: /scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/glm_model/02_subset_gwas.sh
# read ld created by plink2 r-unphased
#file_ld <- "/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/lz_plt/16_69449233_73516227.ld"
file_ld <- "/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/glm_model/out/seq.12593.33/19_54219624_54427869.unphased.vcor1"

file_ld_snps <- paste0(file_ld, ".vars")

ld_raw  <- data.table::fread(file_ld)
ld_snps <- data.table::fread(file_ld_snps, header = F)

colnames(ld_raw) <- ld_snps %>% unlist()
rownames(ld_raw) <- ld_snps %>% unlist() %>% unique()

# subset LD to only snps in PWAS
col_snps <- which(colnames(ld_raw) %in% pwas$SNPID)
row_snps <- which(rownames(ld_raw) %in% pwas$SNPID)

# to be used by susie
ld_matrix <- ld_raw[row_snps] %>% select(pwas$SNPID) %>% as.matrix()



#--------------------------------#
#------     LocusZoom      ------#
#--------------------------------#

pwas %>% 
  ggplot(aes(POS, MLOG10P)) + #-log10(p.value)
  geom_point(size = 3, color = "grey30", show.legend = F) + 
  theme_classic()



#--------------------------------#
#------       SUSIE        ------#
#--------------------------------#


#install.packages("susieR")  v0.12.35

#------------#

#n_sample <- nrow(dosage)
#var_prot <- somagen %>% summarise_at(seqid, ~ sd(.x, na.rm = T))

set.seed(777)

fit_model <- susieR::susie_rss(
  bhat = pwas$BETA,
  shat = pwas$SE,
  #z = (pwas$BETA/pwas$SE)/10,
  n = 9251,
  R = ld_matrix,
  var_y = 1,
  L = 10, 
  estimate_residual_variance = F,
  max_iter = 100
)

# results attributes
fit_model$sets
summary(fit_model)$vars
summary(fit_model)$cs
susieR::susie_get_niter(fit_model)

# plot of posterior priorities
susieR::susie_plot(fit_model, y="PIP", b=pwas$BETA, xlab="Variants")


positive_semi_definite <- TRUE
tryCatch({
  chol(ld_matrix)
}, error = function(e) {
  positive_semi_definite <- FALSE
})


if (!positive_semi_definite) {
  stop("The LD matrix is not positive semi-definite.")
}



