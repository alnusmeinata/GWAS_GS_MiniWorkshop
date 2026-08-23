if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
R.utils::gunzip("geno_imputed.vcf.gz", destname = "geno_imputed.vcf", remove = FALSE, overwrite = TRUE)
