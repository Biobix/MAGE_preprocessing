#library("devtools")
#library(roxygen2)
#devtools::load_all()

#document()
#check()

#setwd("..")
#build("testimpr")
#install("testimpr")

library(testimpr)

chromosomes <- c(15, 11, 14, 6, 7, 8, 9, 10, 12, 13, 16, 17, 18, 19, 20, 21, 22, 1, 2, 3, 4, 5)
wd_samples <- "/data/jeroeng/TCGAkidney/bamfiles/"
wd_seq <-  "/data/jeroeng/TCGAkidney/bamfiles/"
wd_seqem <- "/data/tineg/Package/kidney/imprinted/SeqEM_data/" 
wd_res <- "/data/tineg/Package/kidney/imprinted/"
file_samples <- paste(wd_samples, "kidney_samples.txt", sep = "")

#######
#START#
#######
sample_info <- read.table(file_samples, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
names(sample_info) <- c("sample_nr", "sample")
samples_all <- sample_info$sample
nr_samples_all <- length(samples_all)

for (chr in chromosomes) {
  #READ COUNT FILES
  sample_data <- read.table(paste(wd_seq, "counts_chr", chr, ".txt", sep=""), header = FALSE, stringsAsFactors = F, sep = "\t")
  #sample_data <- read.table(paste(wd_seq, "counts_SNP_chr", chr, ".txt", sep=""), header = FALSE, stringsAsFactors = F, sep = "\t")
  names(sample_data) <- c("chromosome", "position", "ref_alleles", "variationtype", "dbSNP_ref", "gene", "A", "T", "C", "G", "sample")
  sample_data$sample_nr <- unlist(lapply(sample_data$sample, function(x) sample_info[which(sample_info$sample == x), "sample_nr"]))
  sample_data$position <- as.character(sample_data$position)

  #DETERMINE AMOUNT OF LOCI
  positions <- unique(sample_data$position)
  loci <- length(positions)

  #CREATE HASH WITH DATA FRAME FOR EVERY POSITION
  data <- hash::hash()
  for (i in positions) {
    data[[ i ]] <- sample_data[which(sample_data$position==i), ]
  }
  
  #DELETE INDELS AND SUBSTITUTIONS
  for (z in positions) {
    if (data[[z]]$variationtype[1] != "SNV"){
      data[[z]] <- NULL
    }
  }
  positions <- hash::keys(data); loci <- length(positions)

  #DETERMINE STANDARD ALLELES AND REFERENCE/VARIANT COUNTS
  data_c <- hash::hash()
  for (z in positions) {
    data[[z]] <- testimpr::standard_alleles(data[[z]])
    #SAVE CONTROL DATA BEFORE FILTERING
    data_c[[z]] <- data[[z]]
  }

  #PRIOR FILTER OF LOCI
  for (z in positions) {
    data[[z]] <- testimpr::prior_filter(data[[z]], samples_filter = 30)
  }
  positions <- hash::keys(data); loci <- length(positions)

  #ESIMATE PARAMETERS WITH SeqEM (SEQUENCING ERROR RATE, ALLELE FREQUENCY AND INBREEDING), CALCUTE SYMMETRY GOF PER SNP AND FILTER SNPS
  SE_all <- NULL; F_all <- NULL
  for (z in positions) {
    seqem_results <- testimpr::estimate_parameters(data[[z]]$ref_count, data[[z]]$var_count, wd_seqem, ref_allele = data[[z]][1, "ref"], var_allele = data[[z]][1, "var"], position = z, dbSNP = data[[z]][1, "dbSNP_ref"], chr = chr)
    data[[z]]$allelefreq <- seqem_results$allelefreq
    data[[z]]$est_SE <- seqem_results$SE
    SE_all <- c(SE_all, seqem_results$SE)
    F_all <- c(F_all, seqem_results$inbr)

    data[[z]]$sym <- testimpr::symmetry_gof(data[[z]]$ref_count, data[[z]]$var_count, seqem_results$allelefreq)

    if (seqem_results$allelefreq <= 0.15 || seqem_results$allelefreq >= 0.85 || seqem_results$SE > 0.035 ||data[[z]]$sym[1] <= 0.05) {
      data[[z]] <- NULL
    }
  }
  f_inb_chr <- median(F_all)
  SE <- median(SE_all)
  positions <- hash::keys(data); loci <- length(positions)

  #PERFORM IMPRINTING ANALYSIS PER SNP AND CREATE RESULTS DATA FRAME
  results <- data.frame()
  for (z in positions) {
    lrt_results <- testimpr::lrt_i(data[[z]]$ref_count, data[[z]]$var_count, allelefreq = data[[z]][1, "allelefreq"], SE = SE, inbr = f_inb_chr)
    med_imp <- testimpr::median_imprinting(data[[z]]$ref_count, data[[z]]$var_count, allelefreq = data[[z]][1, "allelefreq"], inbr = f_inb_chr)
    results_z <- data.frame("position" = z, "gene" = data[[z]]$gene[1], "LRT" = lrt_results$LRT, "p" = lrt_results$p_value, "estimated.i" = lrt_results$est_i, "allele.frequency" = data[[z]]$allelefreq[1], "dbSNP" = data[[z]]$dbSNP_ref[1], "reference" = data[[z]]$ref[1], "variant" = data[[z]]$var[1], "est_SE" = data[[z]]$est_SE[1], "coverage" = data[[z]]$coverage[1], "nr_samples" = nrow(data[[z]]), "GOF" = lrt_results$GOF_likelihood, "symmetry" = data[[z]]$sym[1], "med_impr" = med_imp, "inbreeding" = f_inb_chr, stringsAsFactors = FALSE)
    results <- rbind(results, results_z)
  }

  #FINAL FILTERING OF INTERESTING AND IMPRINTED SNPS AND WRITE RESULTS FILES
  #results <- final_filter(chr, data, results, wd_res, gof_filt = 0.1, med_impr_filt = 0.1, i_filt = 0.2, file_all = TRUE, file_impr = TRUE, file_all_counts = FALSE, file_impr_counts = TRUE)
  results <- testimpr::final_filter(chr, data, results, wd_res, gof_filt = 1.2, med_impr_filt = 0.8, i_filt = 0.6, file_all = TRUE, file_impr = TRUE, file_all_counts = FALSE, file_impr_counts = TRUE)

  positions_impr <- as.character(results$position)

  #PLOT SIGNIFICANTLY IMPRINTED GENES#
  if (length(positions_impr) > 0) {
    for (z in positions_impr) {
      testimpr::plot_imprinting(data[[z]]$ref_count, data[[z]]$var_count, allelefreq = results[which(results$position==z), "allele.frequency"], impr = results[which(results$position==z), "estimated.i"], SE = SE, wd_res = wd_res, chr = chr, position = z, gene  = results[which(results$position==z), "gene"], inbr = f_inb_chr, coverage = 40, plot_hwe = FALSE)
    }
  }
}


