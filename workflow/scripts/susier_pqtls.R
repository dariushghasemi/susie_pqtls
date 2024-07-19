
#--------------------------------#
#------      Read GWAS     ------#
#--------------------------------#

seqid <- "seq.3054.3"
locus <- "16_69549233_73416227" #"16_69449233_73516227"

# GWAS file
path_base <- "/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/glm_model/out/"
file_pwas <- paste0(path_base, seqid, "/", locus, ".regenie")
file_dose <- paste0(path_base, seqid, "/", locus, "_dosage.raw")

#file_pwas <- paste0("/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/glm_model/out/seq.12593.33/19_54219624_54427869.regenie")

pwas <- data.table::fread(file_pwas)
dose <- data.table::fread(file_dose)




#--------------------------------#
#------    In-sample LD    ------#
#--------------------------------#

# compute in-sample LD from the genotype/dosages
cor_matrix <- dose[,1:100] %>%
  select(- c(IID, FID, PAT, MAT, SEX, PHENOTYPE)) %>%
  #summarise(across(where(is.numeric), ~ sum(is.na(.x))))
  cor(use = "pairwise.complete.obs") #"complete.obs"


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
#file_ld <- "/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/glm_model/out/seq.12593.33/19_54219624_54427869.unphased.vcor1"
file_ld <- "/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/glm_model/out/seq.3054.3/16_69449233_73516227.ld"

# file_ld_snps <- paste0(file_ld, ".vars")
# 
# ld_raw  <- data.table::fread(file_ld)
# ld_snps <- data.table::fread(file_ld_snps, header = F)
# 
# colnames(ld_raw) <- ld_snps %>% unlist()
# rownames(ld_raw) <- ld_snps %>% unlist() %>% unique()
# 
# # subset LD to only snps in PWAS
# col_snps <- which(colnames(ld_raw) %in% pwas$SNPID)
# row_snps <- which(rownames(ld_raw) %in% pwas$SNPID)
# 
# # to be used by susie
# ld_matrix <- ld_raw[row_snps] %>% select(pwas$SNPID) %>% as.matrix()



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

n_sample <- nrow(dose)
#var_prot <- somagen %>% summarise_at(seqid, ~ sd(.x, na.rm = T))

set.seed(777)

fit_model <- susieR::susie_rss(
  bhat = pwas$BETA,
  shat = pwas$SE,
  #z = (pwas$BETA/pwas$SE)/10,
  n = n_sample,
  R = cor_matrix,
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


#--------------------------------#
#------  Positive controls ------#
#--------------------------------#

library(tidyverse)
library(gt) #The package "webshot2" is required to save gt tables as images.

# seq.6556.5
tibble(
  "rsid1"  = c("rs144538541","rs35942390", "rs1047153", NA),
  "snpid1" = c("6:46090938", "6:46110854", "6:46128745", NA),
  "rsid2"  = c("rs12201262", "rs1047153",  "rs34109856", "rs77746021"),
  "snpid2" = c("6:46070642", "6:46128745", "6:46135884", "6:46138158"),
  "Meta_Analysis" = c("6:46101453", "6:46142439", NA, NA)
) %>%
  gt() %>%
  tab_header(title = "seq.6556.5", subtitle = "ENPP5") %>%
  tab_footnote(
    footnote = "Missing in MA.",
    locations = cells_body(
      columns = snpid1, rows = 3
    )) %>%
  tab_spanner(
    label = "INTERVAL",
    columns = c(rsid1, snpid1)
  ) %>%
  tab_spanner(
    label = "deCode",
    columns = c(rsid2, snpid2)
  )

  
# seq.3054.3
tibble(
  "rsid1"  = c("rs217184", "rs9302635", "rs34042070", "rs34309753", NA, NA),
  "snpid1" = c("16:72105965", "16:72144174", "16:72101525", "16:72163932", NA, NA),
  "deCode" = c("rs217184", "rs9302635","SNP #12", "SNP #11", "SNP #10", "SNP #9"),
  "rsid2"  = c("rs61733122", "rs58609735", "rs75203664", "rs763665", "rs34042070", "rs217181"),
  "snpid2" = c("16:71955314", "16:71993165", "16:72067332", "16:72078043", "16:72101525", "16:72114002")
) %>%
  gt() %>%
  tab_header(title = "seq.3054.3", subtitle = "Haptoglobin HP") %>%
  tab_footnote(
    footnote = "Available in MA.",
    locations = cells_body(
      columns = snpid1, rows = 1:3
    )) %>%
  tab_footnote(
    footnote = "Missing in MA.",
    locations = cells_body(
      columns = snpid1, rows = 4
    )) %>%
  tab_footnote(
    footnote = md("**Replicated**"),
    locations = cells_body(
      columns = snpid2, rows = 5
    )) %>%
  tab_spanner(
    label = "INTERVAL",
    columns = c(rsid1, snpid1)
  ) %>%
  tab_spanner(
    label = "Meta_Analysis",
    columns = c(rsid2, snpid2)
  ) #%>% gtsave("18-Jul-24_seq.3054.3.jpeg") #, expand = 10










