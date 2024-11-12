#### 1 - Main analysis ####
 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)

vte <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_VTE.gz")
    vte <- vte[vte$pval < 5e-4,]
    vte <- vte[vte$af_alt > 0.01 & vte$af_alt < 0.99,]
    vte$"#chrom" <- ifelse(vte$"#chrom"==23, "X", vte$"#chrom")
  
longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)
longcovid_strict <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_strictcontrols.gz") 
    longcovid_strict$"#chrom" <- ifelse(longcovid_strict$"#chrom"==23, "X", longcovid_strict$"#chrom")
    colnames(longcovid_strict)[which(colnames(longcovid_strict)=="rsid")] <- "rsids"
    colnames(longcovid_strict)[which(colnames(longcovid_strict)=="stderr_beta")] <- "sebeta"
    colnames(longcovid_strict)[which(colnames(longcovid_strict)=="alt_allele_freq")] <- "af_alt"
    longcovid_strict$pval <- 10^(-longcovid_strict$neg_log_pvalue)

 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
 for (i in c("vte")){
   exp_u <- get(i)  

   for (pvalue in c(5e-4, 5e-6, 5e-8, 5e-10)) {
   exp <- exp_u[exp_u$pval<pvalue,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   exp$"#chrom" <- ifelse(exp$"#chrom"==23, "X", exp$"#chrom")                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   if (is.null(exp) || nrow(exp) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     
  for (outc in c("longcovid", "longcovid_strict")) {
    outcome <- get(outc)
     outcome_overlap <- outcome[outcome$pos %in% exp$pos,]                                                 # Only selecting the variants that are overlapping
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste(outc)
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     exp_overlap <- exp[exp$pos %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     exp_overlap_2 <- exp_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
     exp_overlap_2$beta <- exp_overlap_2$beta * -1
     exp_overlap_2$af_alt  <- 1 - exp_overlap_2$af_alt
     colnames(exp_overlap_2)[colnames(exp_overlap_2) %in% c("ref", "alt")] <- c("alt", "ref")
     exp_overlap <- rbind(exp_overlap, exp_overlap_2)
     exp_overlap$ID <- paste(exp_overlap$"#chrom", exp_overlap$pos, exp_overlap$alt, exp_overlap$ref, sep=":")
     exp_overlap$phen <- i
     exp_overlap <- as.data.frame(exp_overlap)
     exp_overlap <- format_data(exp_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     rm(exp_overlap_2)
     
     dat_u <- harmonise_data(exposure_dat=exp_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
     dat_u <- dat_u[dat_u$mr_keep==T,]

     dat_u$chrom_pos <- paste(dat_u$chr.exposure, dat_u$pos.exposure, sep="_")
     outcome_rsid$chrom_pos <- paste(outcome_rsid$"#chrom", outcome_rsid$pos, sep="_")
     dat_u <- merge(dat_u, outcome_rsid, by="chrom_pos", all.x=T)
     
     for (rsq in c(0.1, 0.01, 0.001, 0.0001)) {
     dat <- dat_u
    
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     dat$pval_thresh <- pvalue
     dat$rsq_thresh <- rsq
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     rm(exp_overlap, outcome_overlap, outcome_rsid, clump)
     
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outc), nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval)
      } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_egger <- MendelianRandomization::mr_egger(dat1)
       output_mr_median <- MendelianRandomization::mr_median(dat1, weighting = "weighted")
       output_mr_mbe <- MendelianRandomization::mr_mbe(dat1, weighting = "weighted")
       results0 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue)
       results1 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Estimate, 
                             se=output_mr_egger@StdError.Est, pval=output_mr_egger@Pvalue.Est)
       results2 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger intercept", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Intercept, 
                             se=output_mr_egger@StdError.Int, pval=output_mr_egger@Pvalue.Int)
       results3 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_median@SNPs, 
                             method="Median", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_median@Estimate, 
                             se=output_mr_median@StdError, pval=output_mr_median@Pvalue)
       results4 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Mode-based", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_mbe@Estimate, 
                             se=output_mr_mbe@StdError, pval=output_mr_mbe@Pvalue)
       results <- rbind(results0, results1, results2, results3, results4)
       rm(results0, results1, results2, results3, results4)
       
      }
      }
     
  if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)
   }
   }
   }
   }
   }
   }
 
   write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1a_main_thromb_longcovid.csv") 
   write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1b_main_thromb_longcovid_instr.csv") 
}

 

   # 1.a - Main analysis - other phenotypes ####
 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)

  vte <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_VTE.gz")
    vte <- vte[vte$pval < 5e-8,]
    vte <- vte[vte$af_alt > 0.01 & vte$af_alt < 0.99,]
    vte$"#chrom" <- ifelse(vte$"#chrom"==23, "X", vte$"#chrom")
  dvt <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_PHLETHROMBDVTLOW.gz")
    dvt <- dvt[dvt$pval < 5e-8,]
    dvt <- dvt[dvt$af_alt > 0.01 & dvt$af_alt < 0.99,]
    dvt$"#chrom" <- ifelse(dvt$"#chrom"==23, "X", dvt$"#chrom")
  pulmem <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_PULMEMB.gz")
    pulmem <- pulmem[pulmem$pval < 5e-8,]
    pulmem <- pulmem[pulmem$af_alt > 0.01 & pulmem$af_alt < 0.99,]
    pulmem$"#chrom" <- ifelse(pulmem$"#chrom"==23, "X", pulmem$"#chrom") 
  afib <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_AF.gz")
    afib <- afib[afib$pval < 5e-8,]
    afib <- afib[afib$af_alt > 0.01 & afib$af_alt < 0.99,]
    afib$"#chrom" <- ifelse(afib$"#chrom"==23, "X", afib$"#chrom")
  asth <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_J10_ASTHMA_EXMORE.gz")
    asth <- asth[asth$pval < 5e-8,]
    asth <- asth[asth$af_alt > 0.01 & asth$af_alt < 0.99,]
    asth$"#chrom" <- ifelse(asth$"#chrom"==23, "X", asth$"#chrom")
  ckd <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_N14_CHRONKIDNEYDIS.gz")
    ckd <- ckd[ckd$pval < 5e-8,]
    ckd <- ckd[ckd$af_alt > 0.01 & ckd$af_alt < 0.99,]
    ckd$"#chrom" <- ifelse(ckd$"#chrom"==23, "X", ckd$"#chrom")
  chron_low <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_J10_LOWCHRON.gz")
    chron_low <- chron_low[chron_low$pval < 5e-8,]
    chron_low <- chron_low[chron_low$af_alt > 0.01 & chron_low$af_alt < 0.99,]
    chron_low$"#chrom" <- ifelse(chron_low$"#chrom"==23, "X", chron_low$"#chrom")
  cad <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_REVASC.gz")
    cad <- cad[cad$pval < 5e-8,]
    cad <- cad[cad$af_alt > 0.01 & cad$af_alt < 0.99,]
    cad$"#chrom" <- ifelse(cad$"#chrom"==23, "X", cad$"#chrom")
  dement <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_F5_DEMENTIA.gz")
    dement <- dement[dement$pval < 5e-8,]
    dement <- dement[dement$af_alt > 0.01 & dement$af_alt < 0.99,]
    dement$"#chrom" <- ifelse(dement$"#chrom"==23, "X", dement$"#chrom")
  diab <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_T2D.gz")
    diab <- diab[diab$pval < 5e-8,]
    diab <- diab[diab$af_alt > 0.01 & diab$af_alt < 0.99,]
    diab$"#chrom" <- ifelse(diab$"#chrom"==23, "X", diab$"#chrom")
  htn <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_HYPTENS.gz")
    htn <- htn[htn$pval < 5e-8,]
    htn <- htn[htn$af_alt > 0.01 & htn$af_alt < 0.99,]
    htn$"#chrom" <- ifelse(htn$"#chrom"==23, "X", htn$"#chrom")
  migr <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_MIGRAINE_TRIPTAN.gz")
    migr <- migr[migr$pval < 5e-8,]
    migr <- migr[migr$af_alt > 0.01 & migr$af_alt < 0.99,]
    migr$"#chrom" <- ifelse(migr$"#chrom"==23, "X", migr$"#chrom")
  ms <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_G6_MS.gz")
    ms <- ms[ms$pval < 5e-8,]
    ms <- ms[ms$af_alt > 0.01 & ms$af_alt < 0.99,]
    ms$"#chrom" <- ifelse(ms$"#chrom"==23, "X", ms$"#chrom")
  osteoarthr <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_M13_ARTHROSIS.gz")
    osteoarthr <- osteoarthr[osteoarthr$pval < 5e-8,]
    osteoarthr <- osteoarthr[osteoarthr$af_alt > 0.01 & osteoarthr$af_alt < 0.99,]
    osteoarthr$"#chrom" <- ifelse(osteoarthr$"#chrom"==23, "X", osteoarthr$"#chrom")
  rheuma <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_M13_RHEUMA.gz")
    rheuma <- rheuma[rheuma$pval < 5e-8,]
    rheuma <- rheuma[rheuma$af_alt > 0.01 & rheuma$af_alt < 0.99,]
    rheuma$"#chrom" <- ifelse(rheuma$"#chrom"==23, "X", rheuma$"#chrom")

longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)

 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
 for (i in c("vte", "dvt", "pulmem", "afib", "asth", "ckd", "chron_low", "cad", "dement", "diab", "htn", "migr", "ms", "osteoarthr", "rheuma")){
   exp_u <- get(i)  

   for (pvalue in c(5e-8)) {
   exp <- exp_u[exp_u$pval<pvalue,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   exp$"#chrom" <- ifelse(exp$"#chrom"==23, "X", exp$"#chrom")                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   if (is.null(exp) || nrow(exp) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     
  for (outc in c("longcovid")) {
    outcome <- get(outc)
     outcome_overlap <- outcome[outcome$pos %in% exp$pos,]                                                 # Only selecting the variants that are overlapping
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste(outc)
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     exp_overlap <- exp[exp$pos %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     exp_overlap_2 <- exp_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
     exp_overlap_2$beta <- exp_overlap_2$beta * -1
     exp_overlap_2$af_alt  <- 1 - exp_overlap_2$af_alt
     colnames(exp_overlap_2)[colnames(exp_overlap_2) %in% c("ref", "alt")] <- c("alt", "ref")
     exp_overlap <- rbind(exp_overlap, exp_overlap_2)
     exp_overlap$ID <- paste(exp_overlap$"#chrom", exp_overlap$pos, exp_overlap$alt, exp_overlap$ref, sep=":")
     exp_overlap$phen <- i
     exp_overlap <- as.data.frame(exp_overlap)
     exp_overlap <- format_data(exp_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     rm(exp_overlap_2)
     
     dat_u <- harmonise_data(exposure_dat=exp_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
     dat_u <- dat_u[dat_u$mr_keep==T,]

     dat_u$chrom_pos <- paste(dat_u$chr.exposure, dat_u$pos.exposure, sep="_")
     outcome_rsid$chrom_pos <- paste(outcome_rsid$"#chrom", outcome_rsid$pos, sep="_")
     dat_u <- merge(dat_u, outcome_rsid, by="chrom_pos", all.x=T)
     
     for (rsq in c(0.001)) {
     dat <- dat_u
    
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     dat$pval_thresh <- pvalue
     dat$rsq_thresh <- rsq
     
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     rm(exp_overlap, outcome_overlap, outcome_rsid, clump)
     
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outc), nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
      } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
      } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_egger <- MendelianRandomization::mr_egger(dat1)
       output_mr_median <- MendelianRandomization::mr_median(dat1, weighting = "weighted")
       output_mr_mbe <- MendelianRandomization::mr_mbe(dat1, weighting = "weighted")
       results0 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results1 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Estimate, 
                             se=output_mr_egger@StdError.Est, pval=output_mr_egger@Pvalue.Est, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results2 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger intercept", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Intercept, 
                             se=output_mr_egger@StdError.Int, pval=output_mr_egger@Pvalue.Int, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results3 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_median@SNPs, 
                             method="Median", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_median@Estimate, 
                             se=output_mr_median@StdError, pval=output_mr_median@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results4 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Mode-based", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_mbe@Estimate, 
                             se=output_mr_mbe@StdError, pval=output_mr_mbe@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results <- rbind(results0, results1, results2, results3, results4)
       rm(results0, results1, results2, results3, results4)
       
      }
      }
     
  if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw, output_mr_egger)
   }
   }
   }
   }
   }
   }
 
   write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1a_main_otherphenos_longcovid.csv") 
   write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1b_main_otherphenos_longcovid_instr.csv") 
}

 

     # 1.a.1 - Main analysis - Steiger ####
 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)

  vte <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_VTE.gz")
    vte <- vte[vte$pval < 5e-8,]
    vte <- vte[vte$af_alt > 0.01 & vte$af_alt < 0.99,]
    vte$"#chrom" <- ifelse(vte$"#chrom"==23, "X", vte$"#chrom")

longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)

 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
 for (i in c("vte")){
   exp_u <- get(i)  

   for (pvalue in c(5e-8)) {
   exp <- exp_u[exp_u$pval<pvalue,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   exp$"#chrom" <- ifelse(exp$"#chrom"==23, "X", exp$"#chrom")                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   if (is.null(exp) || nrow(exp) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     
  for (outc in c("longcovid")) {
    outcome <- get(outc)
     outcome_overlap <- outcome[outcome$pos %in% exp$pos,]                                                 # Only selecting the variants that are overlapping
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste(outc)
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     exp_overlap <- exp[exp$pos %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     exp_overlap_2 <- exp_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
     exp_overlap_2$beta <- exp_overlap_2$beta * -1
     exp_overlap_2$af_alt  <- 1 - exp_overlap_2$af_alt
     colnames(exp_overlap_2)[colnames(exp_overlap_2) %in% c("ref", "alt")] <- c("alt", "ref")
     exp_overlap <- rbind(exp_overlap, exp_overlap_2)
     exp_overlap$ID <- paste(exp_overlap$"#chrom", exp_overlap$pos, exp_overlap$alt, exp_overlap$ref, sep=":")
     exp_overlap$phen <- i
     exp_overlap <- as.data.frame(exp_overlap)
     exp_overlap <- format_data(exp_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     rm(exp_overlap_2)
     
     dat_u <- harmonise_data(exposure_dat=exp_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
     dat_u <- dat_u[dat_u$mr_keep==T,]

     dat_u$chrom_pos <- paste(dat_u$chr.exposure, dat_u$pos.exposure, sep="_")
     outcome_rsid$chrom_pos <- paste(outcome_rsid$"#chrom", outcome_rsid$pos, sep="_")
     dat_u <- merge(dat_u, outcome_rsid, by="chrom_pos", all.x=T)
     
     for (rsq in c(0.001)) {
     dat <- dat_u
    
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     dat$pval_thresh <- pvalue
     dat$rsq_thresh <- rsq
     
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     
     dat$units.exposure <- "log odds"
     dat$ncase.exposure <- 21021
     dat$ncontrol.exposure <- 391160
     dat$prevalence.exposure <- 21021 / (21021+391160)
     dat$units.outcome <- "log odds"
     dat$ncase.outcome <- 3018
     dat$ncontrol.outcome <- 994582
     dat$prevalence.outcome <- 3018 / (3018+994582)
     
     dat <- steiger_filtering(dat)
     
     rm(exp_overlap, outcome_overlap, outcome_rsid, clump)
     
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outc), nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
      } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
      } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_egger <- MendelianRandomization::mr_egger(dat1)
       output_mr_median <- MendelianRandomization::mr_median(dat1, weighting = "weighted")
       output_mr_mbe <- MendelianRandomization::mr_mbe(dat1, weighting = "weighted")
       results0 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results1 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Estimate, 
                             se=output_mr_egger@StdError.Est, pval=output_mr_egger@Pvalue.Est, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results2 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger intercept", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Intercept, 
                             se=output_mr_egger@StdError.Int, pval=output_mr_egger@Pvalue.Int, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results3 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_median@SNPs, 
                             method="Median", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_median@Estimate, 
                             se=output_mr_median@StdError, pval=output_mr_median@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results4 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Mode-based", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_mbe@Estimate, 
                             se=output_mr_mbe@StdError, pval=output_mr_mbe@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results <- rbind(results0, results1, results2, results3, results4)
       rm(results0, results1, results2, results3, results4)
       
      }
      }
     
  if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw, output_mr_egger)
   }
   }
   }
   }
   }
   }
 
   write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1a_main_steiger_longcovid.csv") 
   write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1b_main_steiger_longcovid_instr.csv") 
}

 

   # 1.b - Main analysis - other phenotypes (COVID) ####
 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)

  vte <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_VTE.gz")
    vte <- vte[vte$pval < 5e-8,]
    vte <- vte[vte$af_alt > 0.01 & vte$af_alt < 0.99,]
    vte$"#chrom" <- ifelse(vte$"#chrom"==23, "X", vte$"#chrom")

 severe_covid <- fread("https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz")
 hospitalized_covid <- fread("https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz")
 sars_cov <- fread("https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz")

 for (i in c("severe_covid", "hospitalized_covid", "sars_cov")) {
    outc <- get(i)
    outc <- outc[outc$all_meta_AF > 0.01 & outc$all_meta_AF < 0.99,]
    outc <- outc[outc$all_inv_var_het_p > 0.05,]
    colnames(outc)[which(colnames(outc)=="#CHR")] <- "#chrom"
    colnames(outc)[which(colnames(outc)=="POS")] <- "pos"
    colnames(outc)[which(colnames(outc)=="rsid")] <- "rsid"
    colnames(outc)[which(colnames(outc)=="REF")] <- "ref"
    colnames(outc)[which(colnames(outc)=="ALT")] <- "alt"
    colnames(outc)[which(colnames(outc)=="all_inv_var_meta_p")] <- "pval"
    colnames(outc)[which(colnames(outc)=="all_inv_var_meta_beta")] <- "beta"
    colnames(outc)[which(colnames(outc)=="all_inv_var_meta_sebeta")] <- "stderr_beta"
    colnames(outc)[which(colnames(outc)=="all_meta_AF")] <- "alt_allele_freq"
    outc <- outc[,c("#chrom", "pos", "rsid", "ref", "alt", "pval", "beta", "stderr_beta", "alt_allele_freq")]
    outc$"#chrom" <- ifelse(outc$"#chrom"==23, "X", outc$"#chrom")
    colnames(outc)[which(colnames(outc)=="rsid")] <- "rsids"
    colnames(outc)[which(colnames(outc)=="stderr_beta")] <- "sebeta"
    colnames(outc)[which(colnames(outc)=="alt_allele_freq")] <- "af_alt"
    assign(i, outc)
  }

 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
 for (i in c("vte")){
   exp_u <- get(i)  

   for (pvalue in c(5e-8)) {
   exp <- exp_u[exp_u$pval<pvalue,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   exp$"#chrom" <- ifelse(exp$"#chrom"==23, "X", exp$"#chrom")                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   if (is.null(exp) || nrow(exp) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     
  for (outc in c("severe_covid", "hospitalized_covid", "sars_cov")) {
    outcome <- get(outc)
     outcome_overlap <- outcome[outcome$pos %in% exp$pos,]                                                 # Only selecting the variants that are overlapping
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste(outc)
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     exp_overlap <- exp[exp$pos %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     exp_overlap_2 <- exp_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
     exp_overlap_2$beta <- exp_overlap_2$beta * -1
     exp_overlap_2$af_alt  <- 1 - exp_overlap_2$af_alt
     colnames(exp_overlap_2)[colnames(exp_overlap_2) %in% c("ref", "alt")] <- c("alt", "ref")
     exp_overlap <- rbind(exp_overlap, exp_overlap_2)
     exp_overlap$ID <- paste(exp_overlap$"#chrom", exp_overlap$pos, exp_overlap$alt, exp_overlap$ref, sep=":")
     exp_overlap$phen <- i
     exp_overlap <- as.data.frame(exp_overlap)
     exp_overlap <- format_data(exp_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     rm(exp_overlap_2)
     
     dat_u <- harmonise_data(exposure_dat=exp_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
     dat_u <- dat_u[dat_u$mr_keep==T,]

     dat_u$chrom_pos <- paste(dat_u$chr.exposure, dat_u$pos.exposure, sep="_")
     outcome_rsid$chrom_pos <- paste(outcome_rsid$"#chrom", outcome_rsid$pos, sep="_")
     dat_u <- merge(dat_u, outcome_rsid, by="chrom_pos", all.x=T)
     
     for (rsq in c(0.001)) {
     dat <- dat_u
    
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     dat$pval_thresh <- pvalue
     dat$rsq_thresh <- rsq
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     rm(exp_overlap, outcome_overlap, outcome_rsid, clump)
     
     #dat <- dat[dat$SNP != ]
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outc), nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
      } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
      } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_egger <- MendelianRandomization::mr_egger(dat1)
       output_mr_median <- MendelianRandomization::mr_median(dat1, weighting = "weighted")
       output_mr_mbe <- MendelianRandomization::mr_mbe(dat1, weighting = "weighted")
       results0 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results1 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Estimate, 
                             se=output_mr_egger@StdError.Est, pval=output_mr_egger@Pvalue.Est, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results2 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger intercept", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Intercept, 
                             se=output_mr_egger@StdError.Int, pval=output_mr_egger@Pvalue.Int, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results3 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_median@SNPs, 
                             method="Median", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_median@Estimate, 
                             se=output_mr_median@StdError, pval=output_mr_median@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results4 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Mode-based", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_mbe@Estimate, 
                             se=output_mr_mbe@StdError, pval=output_mr_mbe@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results <- rbind(results0, results1, results2, results3, results4)
       rm(results0, results1, results2, results3, results4)
       
      }
      }
     
  if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw, output_mr_egger)
   }
   }
   }
   }
   }
   }
 
   write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1a_main_otherphenoscovid_longcovid.csv") 
   write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1b_main_otherphenoscovid_longcovid_instr.csv") 
}

 

     # 1.b.1 - Main analysis - adjusted for other phenotypes (COVID) ####
 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)
 library(MVMR)

  vte <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_VTE.gz")
    list_vte <- vte$rsids[vte$pval < 5e-8]
    vte <- vte[vte$af_alt > 0.01 & vte$af_alt < 0.99,]
    vte$"#chrom" <- ifelse(vte$"#chrom"==23, "X", vte$"#chrom")

 severe_covid <- fread("https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz")
 hospitalized_covid <- fread("https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz")
 sars_cov <- fread("https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz")

 for (i in c("severe_covid", "hospitalized_covid", "sars_cov")) {
    outc <- get(i)
    outc <- outc[outc$all_meta_AF > 0.01 & outc$all_meta_AF < 0.99,]
    outc <- outc[outc$all_inv_var_het_p > 0.05,]
    colnames(outc)[which(colnames(outc)=="#CHR")] <- "#chrom"
    colnames(outc)[which(colnames(outc)=="POS")] <- "pos"
    colnames(outc)[which(colnames(outc)=="rsid")] <- "rsid"
    colnames(outc)[which(colnames(outc)=="REF")] <- "ref"
    colnames(outc)[which(colnames(outc)=="ALT")] <- "alt"
    colnames(outc)[which(colnames(outc)=="all_inv_var_meta_p")] <- "pval"
    colnames(outc)[which(colnames(outc)=="all_inv_var_meta_beta")] <- "beta"
    colnames(outc)[which(colnames(outc)=="all_inv_var_meta_sebeta")] <- "stderr_beta"
    colnames(outc)[which(colnames(outc)=="all_meta_AF")] <- "alt_allele_freq"
    outc <- outc[,c("#chrom", "pos", "rsid", "ref", "alt", "pval", "beta", "stderr_beta", "alt_allele_freq")]
    outc$"#chrom" <- ifelse(outc$"#chrom"==23, "X", outc$"#chrom")
    colnames(outc)[which(colnames(outc)=="rsid")] <- "rsids"
    colnames(outc)[which(colnames(outc)=="stderr_beta")] <- "sebeta"
    colnames(outc)[which(colnames(outc)=="alt_allele_freq")] <- "af_alt"
    assign(i, outc)
    assign(paste0("list_", i), outc$rsids[outc$pval<5e-8])
 }
 
 longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)

 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
 for (i in c("vte")){
   exp_u <- get(i)  

  for (exp_2 in c("severe_covid", "hospitalized_covid", "sars_cov")) {
   list_12 <- c(list_vte, get(paste0("list_", exp_2))) 
   list_12 <- unlist(strsplit(as.character(list_12), ","))
   
   exp2 <- get(exp_2)
   exp2 <- exp2 %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
   exp2 <- exp2[!is.na(exp2$rsids) & exp2$rsids %in% list_12,]
   exp2$"#chrom" <- ifelse(exp2$"#chrom"==23, "X", exp2$"#chrom")                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   for (pvalue in c(5e-8)) {
   exp <- exp_u %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
   exp <- exp[!is.na(exp$rsids) & exp$rsids %in% list_12,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   exp$"#chrom" <- ifelse(exp$"#chrom"==23, "X", exp$"#chrom")                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   if (is.null(exp) || nrow(exp) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome <- get("longcovid")
     outcome_overlap <- outcome[outcome$pos %in% exp$pos & outcome$pos %in% exp2$pos ,]                                                 # Only selecting the variants that are overlapping
    
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste("longcovid")
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     exp_overlap <- exp[exp$pos %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     exp_overlap_2 <- exp_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
     exp_overlap_2$beta <- exp_overlap_2$beta * -1
     exp_overlap_2$af_alt  <- 1 - exp_overlap_2$af_alt
     colnames(exp_overlap_2)[colnames(exp_overlap_2) %in% c("ref", "alt")] <- c("alt", "ref")
     exp_overlap <- rbind(exp_overlap, exp_overlap_2)
     exp_overlap$ID <- paste(exp_overlap$"#chrom", exp_overlap$pos, exp_overlap$alt, exp_overlap$ref, sep=":")
     exp_overlap$phen <- i
     exp_overlap <- as.data.frame(exp_overlap)
     exp_overlap <- format_data(exp_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")

     exp2_overlap <- exp2[exp2$pos %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     exp2_overlap_2 <- exp2_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
     exp2_overlap_2$beta <- exp2_overlap_2$beta * -1
     exp2_overlap_2$af_alt  <- 1 - exp2_overlap_2$af_alt
     colnames(exp2_overlap_2)[colnames(exp2_overlap_2) %in% c("ref", "alt")] <- c("alt", "ref")
     exp2_overlap <- rbind(exp2_overlap, exp2_overlap_2)
     exp2_overlap$ID <- paste(exp2_overlap$"#chrom", exp2_overlap$pos, exp2_overlap$alt, exp2_overlap$ref, sep=":")
     exp2_overlap$phen <- i
     exp2_overlap <- as.data.frame(exp2_overlap)
     exp2_overlap <- format_data(exp2_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     rm(exp_overlap_2, exp2_overlap_2)
     
     exp_overlap$chrom_pos <- paste(exp_overlap$chr.exposure, exp_overlap$pos.exposure, sep="_")
     outcome_rsid$chrom_pos <- paste(outcome_rsid$"#chrom", outcome_rsid$pos, sep="_")
     exp_overlap <- merge(exp_overlap, outcome_rsid, by="chrom_pos", all.x=T)

     exp2_overlap$chrom_pos <- paste(exp2_overlap$chr.exposure, exp2_overlap$pos.exposure, sep="_")
     exp2_overlap <- merge(exp2_overlap, outcome_rsid, by="chrom_pos", all.x=T)
     outcome_overlap$chrom_pos <- paste(outcome_overlap$chr.outcome, outcome_overlap$pos.outcome, sep="_")
     outcome_overlap <- merge(outcome_overlap, outcome_rsid, by="chrom_pos", all.x=T)
     
     exp_overlap <- exp_overlap[exp_overlap$SNP %in% outcome_overlap$SNP & exp_overlap$SNP %in% exp2_overlap$SNP ,]
       exp_overlap <- exp_overlap[order(exp_overlap$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
       exp_overlap <- exp_overlap[!duplicated(exp_overlap$SNP),]
       exp_overlap <- exp_overlap[!duplicated(exp_overlap$rsids),]
       exp_overlap <- exp_overlap[order(exp_overlap$SNP),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     exp2_overlap <- exp2_overlap[exp2_overlap$SNP %in% exp_overlap$SNP & exp2_overlap$SNP %in% outcome_overlap$SNP ,]
       exp2_overlap <- exp2_overlap[order(exp2_overlap$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
       exp2_overlap <- exp2_overlap[!duplicated(exp2_overlap$SNP),]
       exp2_overlap <- exp2_overlap[order(exp2_overlap$SNP),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     outcome_overlap <- outcome_overlap[outcome_overlap$SNP %in% exp_overlap$SNP & outcome_overlap$SNP %in% exp2_overlap$SNP,]
       outcome_overlap <- outcome_overlap[order(outcome_overlap$pval.outcome),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
       outcome_overlap <- outcome_overlap[!duplicated(outcome_overlap$SNP),]
       outcome_overlap <- outcome_overlap[order(outcome_overlap$SNP),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])

     for (rsq in c(0.001)) {
     colnames(exp_overlap)[colnames(exp_overlap) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     colnames(exp2_overlap)[colnames(exp2_overlap) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     colnames(outcome_overlap)[colnames(outcome_overlap) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column

     clump <- ld_clump(dplyr::tibble(rsid=exp_overlap$SNP, pval=exp_overlap$pval.exposure, id=exp_overlap$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     clump <- clump[clump$pval<pvalue,]
     clump2 <- ld_clump(dplyr::tibble(rsid=exp2_overlap$SNP, pval=exp2_overlap$pval.exposure, id=exp2_overlap$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     clump2 <- clump2[clump2$pval<pvalue,]
     clump_list <- c(clump$rsid, clump2$rsid)
     
     colnames(exp_overlap)[colnames(exp_overlap) %in% c("pos_id", "SNP")] <- c("SNP", "rsids")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     colnames(exp2_overlap)[colnames(exp2_overlap) %in% c("pos_id", "SNP")] <- c("SNP", "rsids")
     colnames(outcome_overlap)[colnames(outcome_overlap) %in% c("pos_id", "SNP")] <- c("SNP", "rsids")
     exp_overlap <- exp_overlap[exp_overlap$rsids %in% clump_list,]
     exp2_overlap <- exp2_overlap[exp2_overlap$rsids %in% clump_list,]
     outcome_overlap <- outcome_overlap[outcome_overlap$rsids %in% clump_list,]

     dat <- format_mvmr(BXGs = data.frame(exp_overlap$beta.exposure, exp2_overlap$beta.exposure),
                           BYG = outcome_overlap$beta.outcome,
                           seBXGs = data.frame(exp_overlap$se.exposure, exp2_overlap$se.exposure),
                           seBYG = outcome_overlap$se.outcome,
                           RSID = outcome_overlap$rsids)
     fstat <- strength_mvmr(dat)

     dat_mrmv <- mr_mvinput(bx=as.matrix(data.frame(exp_overlap$beta.exposure, exp2_overlap$beta.exposure)), 
                            bxse=as.matrix(data.frame(exp_overlap$se.exposure, exp2_overlap$se.exposure)), 
                            by=outcome_overlap$beta.outcome, 
                            byse=outcome_overlap$se.outcome, 
                            exposure=c(i, exp_2), 
                            outcome="long covid", 
                            snps=dat$SNP)
     output_mvmr <- mr_mvivw(dat_mrmv)
     results <- data.frame(outc=paste("long covid"), exp=output_mvmr@Exposure, nsnp=output_mvmr@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mvmr@Estimate, 
                             se=output_mvmr@StdError, pval=output_mvmr@Pvalue, fstats=c(fstat$exposure1, fstat$exposure2))
     
   if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping"))
   } else {
      df_sum <- rbind(df_sum, results)
      rm(dat, results, dat_mrmv, fstat, output_mvmr)
   }

   }
   }
   }
   }
   }
    
   write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1a_add_multivariable_covid.csv") 
   write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1b_add_multivariable_covid_instr.csv") 
}

 

   # 1.c - Main analysis - reverse MR ####
 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)

vte <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_VTE.gz")
    vte <- vte[vte$af_alt > 0.01 & vte$af_alt < 0.99,]
    vte$"#chrom" <- ifelse(vte$"#chrom"==23, "X", vte$"#chrom")
  
longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)
    longcovid <- longcovid[longcovid$pval < 5e-4,]
    
 severe_covid <- fread("https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz")
 hospitalized_covid <- fread("https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz")
 sars_cov <- fread("https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz")
 for (i in c("severe_covid", "hospitalized_covid", "sars_cov")) {
    outc <- get(i)
    outc <- outc[outc$all_meta_AF > 0.01 & outc$all_meta_AF < 0.99,]
    outc <- outc[outc$all_inv_var_het_p > 0.05,]
    colnames(outc)[which(colnames(outc)=="#CHR")] <- "#chrom"
    colnames(outc)[which(colnames(outc)=="POS")] <- "pos"
    colnames(outc)[which(colnames(outc)=="rsid")] <- "rsid"
    colnames(outc)[which(colnames(outc)=="REF")] <- "ref"
    colnames(outc)[which(colnames(outc)=="ALT")] <- "alt"
    colnames(outc)[which(colnames(outc)=="all_inv_var_meta_p")] <- "pval"
    colnames(outc)[which(colnames(outc)=="all_inv_var_meta_beta")] <- "beta"
    colnames(outc)[which(colnames(outc)=="all_inv_var_meta_sebeta")] <- "stderr_beta"
    colnames(outc)[which(colnames(outc)=="all_meta_AF")] <- "alt_allele_freq"
    outc <- outc[,c("#chrom", "pos", "rsid", "ref", "alt", "pval", "beta", "stderr_beta", "alt_allele_freq")]
    outc$"#chrom" <- ifelse(outc$"#chrom"==23, "X", outc$"#chrom")
    colnames(outc)[which(colnames(outc)=="rsid")] <- "rsids"
    colnames(outc)[which(colnames(outc)=="stderr_beta")] <- "sebeta"
    colnames(outc)[which(colnames(outc)=="alt_allele_freq")] <- "af_alt"
    assign(i, outc)
  }

 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
 for (i in c("longcovid", "severe_covid", "hospitalized_covid", "sars_cov")){
   exp_u <- get(i)  

   for (pvalue in c(5e-8)) {
   exp <- exp_u[exp_u$pval<pvalue,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   exp$"#chrom" <- ifelse(exp$"#chrom"==23, "X", exp$"#chrom")                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   if (is.null(exp) || nrow(exp) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     
  for (outc in c("vte")) {
    outcome <- get(outc)
     outcome_overlap <- outcome[outcome$pos %in% exp$pos,]                                                 # Only selecting the variants that are overlapping
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste(outc)
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     exp_overlap <- exp[exp$pos %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     exp_overlap_2 <- exp_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
     exp_overlap_2$beta <- exp_overlap_2$beta * -1
     exp_overlap_2$af_alt  <- 1 - exp_overlap_2$af_alt
     colnames(exp_overlap_2)[colnames(exp_overlap_2) %in% c("ref", "alt")] <- c("alt", "ref")
     exp_overlap <- rbind(exp_overlap, exp_overlap_2)
     exp_overlap$ID <- paste(exp_overlap$"#chrom", exp_overlap$pos, exp_overlap$alt, exp_overlap$ref, sep=":")
     exp_overlap$phen <- i
     exp_overlap <- as.data.frame(exp_overlap)
     exp_overlap <- format_data(exp_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     rm(exp_overlap_2)
     
     dat_u <- harmonise_data(exposure_dat=exp_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
     dat_u <- dat_u[dat_u$mr_keep==T,]

     dat_u$chrom_pos <- paste(dat_u$chr.exposure, dat_u$pos.exposure, sep="_")
     outcome_rsid$chrom_pos <- paste(outcome_rsid$"#chrom", outcome_rsid$pos, sep="_")
     dat_u <- merge(dat_u, outcome_rsid, by="chrom_pos", all.x=T)
     
     for (rsq in c(0.001)) {
     dat <- dat_u
    
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     dat$pval_thresh <- pvalue
     dat$rsq_thresh <- rsq
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     rm(exp_overlap, outcome_overlap, outcome_rsid, clump)
     
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(outc=paste(outc), exp=paste(i), nsnp=results_mr$nsnp, pval_thresh=pvalue, rsq_thresh=rsq, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
      } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_egger <- MendelianRandomization::mr_egger(dat1)
       output_mr_median <- MendelianRandomization::mr_median(dat1, weighting = "weighted")
       output_mr_mbe <- MendelianRandomization::mr_mbe(dat1, weighting = "weighted")
       results0 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results1 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Estimate, 
                             se=output_mr_egger@StdError.Est, pval=output_mr_egger@Pvalue.Est, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results2 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger intercept", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Intercept, 
                             se=output_mr_egger@StdError.Int, pval=output_mr_egger@Pvalue.Int, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results3 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_median@SNPs, 
                             method="Median", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_median@Estimate, 
                             se=output_mr_median@StdError, pval=output_mr_median@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results4 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Mode-based", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_mbe@Estimate, 
                             se=output_mr_mbe@StdError, pval=output_mr_mbe@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       results <- rbind(results0, results1, results2, results3, results4)
       rm(results0, results1, results2, results3, results4)
       
      }
      }
     
  if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)
   }
   }
   }
   }
   }
   }
 
   write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1a_main_thromb_longcovid_reverse.csv") 
   write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1b_main_thromb_longcovid_reverse_instr.csv") 
}

 


   # 1.d - Main analysis - VTE (MVP) and COVID-resembling conditions ####
 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)

  vte <- fread("/medpop/esp2/projects/summary_statistics/MVP/release/submission/sub20200416/anno.CLEANED.MVP.EUR.VTE.results.txt.gz")
    vte <- rename(vte, rsids = SNPID, "#chrom" = CHROM, pos_hg19 = POS, alt = EFFECT_ALLELE, ref = OTHER_ALLELE, 
                       af_alt = EAF, beta = BETA, sebeta = SE, pval = PVAL)
    vte <- vte[vte$pval < 5e-4,]
    vte <- vte[vte$af_alt > 0.01 & vte$af_alt < 0.99,]
    vte$"#chrom" <- ifelse(vte$"#chrom"==23, "X", vte$"#chrom")

  vte_2 <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_VTE.gz")
    vte_2 <- vte_2[vte_2$af_alt > 0.01 & vte_2$af_alt < 0.99,]
    vte_2$"#chrom" <- ifelse(vte_2$"#chrom"==23, "X", vte_2$"#chrom")
  dvt <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_PHLETHROMBDVTLOW.gz")
    dvt <- dvt[dvt$af_alt > 0.01 & dvt$af_alt < 0.99,]
    dvt$"#chrom" <- ifelse(dvt$"#chrom"==23, "X", dvt$"#chrom")
  pulmem <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_PULMEMB.gz")
    pulmem <- pulmem[pulmem$af_alt > 0.01 & pulmem$af_alt < 0.99,]
    pulmem$"#chrom" <- ifelse(pulmem$"#chrom"==23, "X", pulmem$"#chrom") 
  bleeding <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_BLEEDING.gz")
    bleeding <- bleeding[bleeding$af_alt > 0.01 & bleeding$af_alt < 0.99,]
    bleeding$"#chrom" <- ifelse(bleeding$"#chrom"==23, "X", bleeding$"#chrom")
  chronic_lower_resp <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_J10_LOWCHRON.gz")
    chronic_lower_resp <- chronic_lower_resp[chronic_lower_resp$af_alt > 0.01 & chronic_lower_resp$af_alt < 0.99,]
    chronic_lower_resp$"#chrom" <- ifelse(chronic_lower_resp$"#chrom"==23, "X", chronic_lower_resp$"#chrom")
  dyspepsia <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_K11_FUNCDYSP.gz")
    dyspepsia <- dyspepsia[dyspepsia$af_alt > 0.01 & dyspepsia$af_alt < 0.99,]
    dyspepsia$"#chrom" <- ifelse(dyspepsia$"#chrom"==23, "X", dyspepsia$"#chrom")
  ibs <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_K11_IBS.gz")
    ibs <- ibs[ibs$af_alt > 0.01 & ibs$af_alt < 0.99,]
    ibs$"#chrom" <- ifelse(ibs$"#chrom"==23, "X", ibs$"#chrom")
  endometriosis <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_N14_ENDOMETRIOSIS.gz")
    endometriosis <- endometriosis[endometriosis$af_alt > 0.01 & endometriosis$af_alt < 0.99,]
    endometriosis$"#chrom" <- ifelse(endometriosis$"#chrom"==23, "X", endometriosis$"#chrom")
  fibromyalgia <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_M13_FIBROMYALGIA.gz")
    fibromyalgia <- fibromyalgia[fibromyalgia$af_alt > 0.01 & fibromyalgia$af_alt < 0.99,]
    fibromyalgia$"#chrom" <- ifelse(fibromyalgia$"#chrom"==23, "X", fibromyalgia$"#chrom")
  anxiety <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_F5_ALLANXIOUS.gz")
    anxiety <- anxiety[anxiety$af_alt > 0.01 & anxiety$af_alt < 0.99,]
    anxiety$"#chrom" <- ifelse(anxiety$"#chrom"==23, "X", anxiety$"#chrom")
  memloss <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_MEMLOSS.gz")
    memloss <- memloss[memloss$af_alt > 0.01 & memloss$af_alt < 0.99,]
    memloss$"#chrom" <- ifelse(memloss$"#chrom"==23, "X", memloss$"#chrom")
  sleep <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_SLEEP.gz")
    sleep <- sleep[sleep$af_alt > 0.01 & sleep$af_alt < 0.99,]
    sleep$"#chrom" <- ifelse(sleep$"#chrom"==23, "X", sleep$"#chrom")
  mood <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_F5_MOOD.gz")
    mood <- mood[mood$af_alt > 0.01 & mood$af_alt < 0.99,]
    mood$"#chrom" <- ifelse(mood$"#chrom"==23, "X", mood$"#chrom")
  persist_mood <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_F5_PERSMOOD.gz")
    persist_mood <- persist_mood[persist_mood$af_alt > 0.01 & persist_mood$af_alt < 0.99,]
    persist_mood$"#chrom" <- ifelse(persist_mood$"#chrom"==23, "X", persist_mood$"#chrom")
  depress <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_F5_DEPRESSIO.gz")
    depress <- depress[depress$af_alt > 0.01 & depress$af_alt < 0.99,]
    depress$"#chrom" <- ifelse(depress$"#chrom"==23, "X", depress$"#chrom")
  postvir_fatigue <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_G6_POSTVIRFAT.gz")
    postvir_fatigue <- postvir_fatigue[postvir_fatigue$af_alt > 0.01 & postvir_fatigue$af_alt < 0.99,]
    postvir_fatigue$"#chrom" <- ifelse(postvir_fatigue$"#chrom"==23, "X", postvir_fatigue$"#chrom")

      
 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
 for (i in c("vte")){
   exp_u <- get(i)  

   for (pvalue in c(5e-8)) {
   exp <- exp_u[exp_u$pval<pvalue,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   exp$"#chrom" <- ifelse(exp$"#chrom"==23, "X", exp$"#chrom")                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   if (is.null(exp) || nrow(exp) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     
  for (outc in c("vte_2", "dvt", "pulmem", "bleeding", "chronic_lower_resp", "dyspepsia", 
                 "ibs", "endometriosis", "fibromyalgia", "anxiety", "memloss",
                "sleep", "mood", "persist_mood", "depress", "postvir_fatigue")) {
    outcome <- get(outc)
     outcome_overlap <- outcome[outcome$rsids %in% exp$rsids,]                                                 # Only selecting the variants that are overlapping
     exp_overlap <- exp[exp$rsids %in% outcome_overlap$rsids,]  
     exp_overlap <- merge(exp_overlap, outcome_overlap[,c("rsids", "pos")], by="rsids", all.x=T, all.y=F)
     outcome_overlap$alt_ad <- ifelse(nchar(outcome_overlap$alt) > nchar(outcome_overlap$ref), "I", 
                                      ifelse(nchar(outcome_overlap$alt) < nchar(outcome_overlap$ref), "D", outcome_overlap$alt))
     outcome_overlap$ref_ad <- ifelse(nchar(outcome_overlap$ref) > nchar(outcome_overlap$alt), "I", 
                                      ifelse(nchar(outcome_overlap$ref) < nchar(outcome_overlap$alt), "D", outcome_overlap$ref))
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste(outc)
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt_ad, outcome_overlap$ref_ad, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt_ad", other_allele_col="ref_ad", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     exp_overlap <- exp_overlap[exp_overlap$pos %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     exp_overlap_2 <- exp_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
     exp_overlap_2$beta <- exp_overlap_2$beta * -1
     exp_overlap_2$af_alt  <- 1 - exp_overlap_2$af_alt
     colnames(exp_overlap_2)[colnames(exp_overlap_2) %in% c("alt", "ref")] <- c("ref", "alt")
     exp_overlap <- rbind(exp_overlap, exp_overlap_2)
     exp_overlap$ID <- paste(exp_overlap$"#chrom", exp_overlap$pos, exp_overlap$alt, exp_overlap$ref, sep=":")
     exp_overlap$phen <- i
     exp_overlap <- as.data.frame(exp_overlap)
     exp_overlap <- format_data(exp_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     rm(exp_overlap_2)
     
     dat_u <- harmonise_data(exposure_dat=exp_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
     dat_u <- dat_u[dat_u$mr_keep==T,]

     dat_u$chrom_pos <- paste(dat_u$chr.exposure, dat_u$pos.exposure, sep="_")
     outcome_rsid$chrom_pos <- paste(outcome_rsid$"#chrom", outcome_rsid$pos, sep="_")
     dat_u <- merge(dat_u, outcome_rsid, by="chrom_pos", all.x=T)
     
     for (rsq in c(0.001)) {
     dat <- dat_u
    
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     dat$pval_thresh <- pvalue
     dat$rsq_thresh <- rsq
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     rm(exp_overlap, outcome_overlap, outcome_rsid, clump)
     
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outc), nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval)
      } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_egger <- MendelianRandomization::mr_egger(dat1)
       output_mr_median <- MendelianRandomization::mr_median(dat1, weighting = "weighted")
       output_mr_mbe <- MendelianRandomization::mr_mbe(dat1, weighting = "weighted")
       results0 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue)
       results1 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Estimate, 
                             se=output_mr_egger@StdError.Est, pval=output_mr_egger@Pvalue.Est)
       results2 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger intercept", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Intercept, 
                             se=output_mr_egger@StdError.Int, pval=output_mr_egger@Pvalue.Int)
       results3 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_median@SNPs, 
                             method="Median", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_median@Estimate, 
                             se=output_mr_median@StdError, pval=output_mr_median@Pvalue)
       results4 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Mode-based", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_mbe@Estimate, 
                             se=output_mr_mbe@StdError, pval=output_mr_mbe@Pvalue)
       results <- rbind(results0, results1, results2, results3, results4)
       rm(results0, results1, results2, results3, results4)
       
      }
      }
     
  if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)
   }
   }
   }
   }
   }
   }
 
   write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1a_main_thromb_longcovid_covidresemblingconditions.csv") 
   write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/1b_main_thromb_longcovid_covidresemblingconditions_instr.csv") 
}

 


#### 2 - Replication analysis - VTE in MVP ####
 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)

vte <- fread("/medpop/esp2/projects/summary_statistics/MVP/release/submission/sub20200416/anno.CLEANED.MVP.EUR.VTE.results.txt.gz")
    vte <- rename(vte, rsids = SNPID, "#chrom" = CHROM, pos_hg19 = POS, alt = EFFECT_ALLELE, ref = OTHER_ALLELE, 
                       af_alt = EAF, beta = BETA, sebeta = SE, pval = PVAL)
    vte <- vte[vte$pval < 5e-4,]
    vte <- vte[vte$af_alt > 0.01 & vte$af_alt < 0.99,]
    vte$"#chrom" <- ifelse(vte$"#chrom"==23, "X", vte$"#chrom")
  
longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)

 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
 for (i in c("vte")){
   exp_u <- get(i)  

   for (pvalue in c(5e-4, 5e-6, 5e-8, 5e-10)) {
   exp <- exp_u[exp_u$pval<pvalue,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   exp$"#chrom" <- ifelse(exp$"#chrom"==23, "X", exp$"#chrom")                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   if (is.null(exp) || nrow(exp) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     
  for (outc in c("longcovid")) {
    outcome <- get(outc)
     outcome_overlap <- outcome[outcome$rsids %in% exp$rsids,]                                                 # Only selecting the variants that are overlapping
     exp_overlap <- exp[exp$rsids %in% outcome_overlap$rsids,]  
     exp_overlap <- merge(exp_overlap, outcome_overlap[,c("rsids", "pos")], by="rsids", all.x=T, all.y=F)
     outcome_overlap$alt_ad <- ifelse(nchar(outcome_overlap$alt) > nchar(outcome_overlap$ref), "I", 
                                      ifelse(nchar(outcome_overlap$alt) < nchar(outcome_overlap$ref), "D", outcome_overlap$alt))
     outcome_overlap$ref_ad <- ifelse(nchar(outcome_overlap$ref) > nchar(outcome_overlap$alt), "I", 
                                      ifelse(nchar(outcome_overlap$ref) < nchar(outcome_overlap$alt), "D", outcome_overlap$ref))
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste(outc)
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt_ad, outcome_overlap$ref_ad, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt_ad", other_allele_col="ref_ad", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     exp_overlap <- exp_overlap[exp_overlap$pos %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     exp_overlap_2 <- exp_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
     exp_overlap_2$beta <- exp_overlap_2$beta * -1
     exp_overlap_2$af_alt  <- 1 - exp_overlap_2$af_alt
     colnames(exp_overlap_2)[colnames(exp_overlap_2) %in% c("alt", "ref")] <- c("ref", "alt")
     exp_overlap <- rbind(exp_overlap, exp_overlap_2)
     exp_overlap$ID <- paste(exp_overlap$"#chrom", exp_overlap$pos, exp_overlap$alt, exp_overlap$ref, sep=":")
     exp_overlap$phen <- i
     exp_overlap <- as.data.frame(exp_overlap)
     exp_overlap <- format_data(exp_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     rm(exp_overlap_2)
     
     dat_u <- harmonise_data(exposure_dat=exp_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
     dat_u <- dat_u[dat_u$mr_keep==T,]

     dat_u$chrom_pos <- paste(dat_u$chr.exposure, dat_u$pos.exposure, sep="_")
     outcome_rsid$chrom_pos <- paste(outcome_rsid$"#chrom", outcome_rsid$pos, sep="_")
     dat_u <- merge(dat_u, outcome_rsid, by="chrom_pos", all.x=T)
     
     for (rsq in c(0.1, 0.01, 0.001, 0.0001)) {
     dat <- dat_u
    
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     dat$pval_thresh <- pvalue
     dat$rsq_thresh <- rsq
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     rm(exp_overlap, outcome_overlap, outcome_rsid, clump)
     
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outc), nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval)
      } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_egger <- MendelianRandomization::mr_egger(dat1)
       output_mr_median <- MendelianRandomization::mr_median(dat1, weighting = "weighted")
       output_mr_mbe <- MendelianRandomization::mr_mbe(dat1, weighting = "weighted")
       results0 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue)
       results1 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Estimate, 
                             se=output_mr_egger@StdError.Est, pval=output_mr_egger@Pvalue.Est)
       results2 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger intercept", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Intercept, 
                             se=output_mr_egger@StdError.Int, pval=output_mr_egger@Pvalue.Int)
       results3 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_median@SNPs, 
                             method="Median", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_median@Estimate, 
                             se=output_mr_median@StdError, pval=output_mr_median@Pvalue)
       results4 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Mode-based", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_mbe@Estimate, 
                             se=output_mr_mbe@StdError, pval=output_mr_mbe@Pvalue)
       results <- rbind(results0, results1, results2, results3, results4)
       rm(results0, results1, results2, results3, results4)
       
      }
      }
     
  if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)
   }
   }
   }
   }
   }
   }
 
   write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/2a_repl_thromb_longcovid_mvp.csv") 
   write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/2b_repl_thromb_longcovid_mvp_instr.csv") 
}

 


#### 3 - Replication analysis - long COVID in PHOSP ####
 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)

 expand_alleles <- function(row) {
    alleles <- strsplit(row$V4, ":")[[1]]
    prefix <- paste(alleles[1:3], collapse = ":")
    new_alleles <- strsplit(alleles[4], ",")[[1]]
    data.frame(V1 = row$V1,
               V2 = row$V2,
               V3 = row$V3,
               V4 = paste0(prefix, ":", new_alleles),
               stringsAsFactors = FALSE)
  }
 
vte <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_VTE.gz")
    vte <- vte[vte$pval < 5e-4,]
    vte <- vte[vte$af_alt > 0.01 & vte$af_alt < 0.99,]
    vte$"#chrom" <- ifelse(vte$"#chrom"==23, "X", vte$"#chrom")
  
 longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/PHOSP-COVID_WAIN_NQ1.3_F6_ALL_ALL_EUR_697_400_REGENIE_20231010.txt.gz")  
    longcovid <- longcovid[longcovid$MAF>0.01,]
    colnames(longcovid)[which(colnames(longcovid)=="CHROM")] <- "#chrom"
    colnames(longcovid)[which(colnames(longcovid)=="GENPOS")] <- "pos"
    colnames(longcovid)[which(colnames(longcovid)=="ALLELE0")] <- "ref"
    colnames(longcovid)[which(colnames(longcovid)=="ALLELE1")] <- "alt"
    colnames(longcovid)[which(colnames(longcovid)=="LOG10P")] <- "neg_log_pvalue"
    colnames(longcovid)[which(colnames(longcovid)=="BETA")] <- "beta"
    colnames(longcovid)[which(colnames(longcovid)=="SE")] <- "stderr_beta"
    colnames(longcovid)[which(colnames(longcovid)=="A1FREQ")] <- "alt_allele_freq"
    longcovid <- longcovid[,c("#chrom", "pos", "ref", "alt", "neg_log_pvalue", "beta", "stderr_beta", "alt_allele_freq")]
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)

 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   
 rsid_u <- fread("/medpop/esp2/btruong/Tools/hg38_common_splitMULTI.db151.vcf.gz")
  rsid_u <- rsid_u[,-c("QUAL", "FILTER", "INFO")]
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
 for (i in c("vte")){
   exp_u <- get(i)  

   for (pvalue in c(5e-4, 5e-6, 5e-8, 5e-10)) {
   exp <- exp_u[exp_u$pval<pvalue,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   exp$"#chrom" <- ifelse(exp$"#chrom"==23, "X", exp$"#chrom")                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   rsid <- rsid_u[rsid_u$POS %in% exp$pos,]
   rsid_2 <- rsid
   colnames(rsid_2)[which(colnames(rsid_2) %in% c("REF", "ALT"))] <- c("ALT", "REF")
   rsid <- rbind(rsid, rsid_2)
   rm(rsid_2)
   
   if (is.null(exp) || nrow(exp) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     
  for (outc in c("longcovid")) {
   outcome <- get(outc)
     outcome_overlap <- outcome[outcome$pos %in% exp$pos,]                                                 # Only selecting the variants that are overlapping
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_overlap$phen <- paste(outc)
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     exp_overlap <- exp[exp$pos %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     exp_overlap_2 <- exp_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
     exp_overlap_2$beta <- exp_overlap_2$beta * -1
     exp_overlap_2$af_alt  <- 1 - exp_overlap_2$af_alt
     colnames(exp_overlap_2)[colnames(exp_overlap_2) %in% c("ref", "alt")] <- c("alt", "ref")
     exp_overlap <- rbind(exp_overlap, exp_overlap_2)
     exp_overlap$ID <- paste(exp_overlap$"#chrom", exp_overlap$pos, exp_overlap$alt, exp_overlap$ref, sep=":")
     exp_overlap$phen <- i
     exp_overlap <- as.data.frame(exp_overlap)
     exp_overlap <- format_data(exp_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     rm(exp_overlap_2)
     
     dat_u <- harmonise_data(exposure_dat=exp_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
     dat_u <- dat_u[dat_u$mr_keep==T,]

     dat_u$chrom_pos <- paste(dat_u$chr.exposure, dat_u$pos.exposure, sep="_")
     rsid$chrom_pos <- paste(rsid$"#CHROM", rsid$POS, sep="_")
     dat_u <- merge(dat_u, rsid[,c("chrom_pos", "ID")], by="chrom_pos", all.x=T) 
     colnames(dat_u)[colnames(dat_u) %in% c("ID")] <- c("rsids")
     
     for (rsq in c(0.1, 0.01, 0.001, 0.0001)) {
     dat <- dat_u
    
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     dat$pval_thresh <- pvalue
     dat$rsq_thresh <- rsq
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     rm(exp_overlap, outcome_overlap, outcome_rsid, clump)
     
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outc), nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval)
      } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_egger <- MendelianRandomization::mr_egger(dat1)
       output_mr_median <- MendelianRandomization::mr_median(dat1, weighting = "weighted")
       output_mr_mbe <- MendelianRandomization::mr_mbe(dat1, weighting = "weighted")
       results0 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue)
       results1 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Estimate, 
                             se=output_mr_egger@StdError.Est, pval=output_mr_egger@Pvalue.Est)
       results2 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Egger intercept", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger@Intercept, 
                             se=output_mr_egger@StdError.Int, pval=output_mr_egger@Pvalue.Int)
       results3 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_median@SNPs, 
                             method="Median", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_median@Estimate, 
                             se=output_mr_median@StdError, pval=output_mr_median@Pvalue)
       results4 <- data.frame(outc=paste(outc), exp=paste(i), nsnp=output_mr_egger@SNPs, 
                             method="Mode-based", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_mbe@Estimate, 
                             se=output_mr_mbe@StdError, pval=output_mr_mbe@Pvalue)
       results <- rbind(results0, results1, results2, results3, results4)
       rm(results0, results1, results2, results3, results4)
       
      }
      }
     
  if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)
   }
   }
   }
   }
   }
   }
 
   write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/2a_repl_thromb_longcovid_phosp.csv") 
   write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/2b_repl_thromb_longcovid_phosp_instr.csv") 
}

 


#### 4 - Downstream analysis - per coagulation and complement gene (only variants within 200 kilobases) ####
 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)
 
 sumstats_info <- fread("/medpop/esp2/aschuerm/coagulation_longcovid/input/map_gtex.csv")
    targets <- c('A2M', 'BDKRB1', 'BDKRB2', 'C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S', 'C2', 'C3', 'C3AR1', 'C4A', 'C4B', 
                 'C4BPA', 'C4BPB', 'C5', 'C5AR1', 'C6', 'C7', 'C8A', 'C8B', 'C8G', 'C9', 'CD46', 'CD55', 'CD59', 'CFB', 'CFD', 
                 'CFH', 'CFI', 'CPB2', 'CR1', 'CR2', 'F10', 'F11', 'F12', 'F13A1', 'F13B', 'F2', 'F2R', 'F3', 'F5', 'F7', 'F8', 
                 'F9', 'FGA', 'FGB', 'FGG', 'KLKB1', 'KNG1', 'MASP1', 'MASP2', 'MBL2', 'PLAT', 'PLAU', 'PLAUR', 'PLG', 'PROC', 
                 'PROS1', 'SERPINA1', 'SERPINA5', 'SERPINC1', 'SERPIND1', 'SERPINE1', 'SERPINF2', 'SERPING1', 'TFPI', 'THBD', 'VWF')
    sumstats_info <- sumstats_info[sumstats_info$external_gene_name %in% targets,]
    rm(targets)

  vte <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_VTE.gz")
      vte <- vte[vte$pval < 5e-4,]
      vte <- vte[vte$af_alt > 0.01 & vte$af_alt < 0.99,]
      vte$"#chrom" <- ifelse(vte$"#chrom"==23, "X", vte$"#chrom")
    
  longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
      longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
      colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
      colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
      colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
      longcovid$pval <- 10^(-longcovid$neg_log_pvalue)

 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
 for (i in c("vte")){

   exp_u <- get(i)
   exp_genes <- exp_u[exp_u$"#chrom" == sumstats_info$chromosome_name[sumstats_info$external_gene_name==sumstats_info$external_gene_name[1]]]
   exp_genes <- exp_genes[exp_genes$pos > (sumstats_info$start_position[sumstats_info$external_gene_name==sumstats_info$external_gene_name[1]] - 200000) &                              # Selecting the cis region only (here defined as 1Mb before or after the protein-encoding region)
                          exp_genes$pos < (sumstats_info$end_position[sumstats_info$external_gene_name==sumstats_info$external_gene_name[1]] + 200000), ]
   exp_genes$gene <- sumstats_info$external_gene_name[1]
   
   for (nr in 2:nrow(sumstats_info)){
     exp_genes_add <- exp_u[exp_u$"#chrom" == sumstats_info$chromosome_name[sumstats_info$external_gene_name==sumstats_info$external_gene_name[nr]]]
     exp_genes_add <- exp_genes_add[exp_genes_add$pos > (sumstats_info$start_position[sumstats_info$external_gene_name==sumstats_info$external_gene_name[nr]] - 200000) &                              # Selecting the cis region only (here defined as 1Mb before or after the protein-encoding region)
                                    exp_genes_add$pos < (sumstats_info$end_position[sumstats_info$external_gene_name==sumstats_info$external_gene_name[nr]] + 200000), ]
     exp_genes_add$gene <- sumstats_info$external_gene_name[nr]
     exp_genes <- rbind(exp_genes, exp_genes_add)
   } 
   
   for (gene_unique in unique(exp_genes$gene)) {
   for (pvalue in c(5e-4, 5e-6, 5e-8)) {
     
   exp <- exp_genes[exp_genes$gene==gene_unique,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   exp <- exp[exp$pval<pvalue,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   exp$"#chrom" <- ifelse(exp$"#chrom"==23, "X", exp$"#chrom")                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   if (is.null(exp) || nrow(exp) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     
  for (outc in c("longcovid")) {
     
    outcome <- get(outc)
     outcome_overlap <- outcome[outcome$pos %in% exp$pos,]                                                 # Only selecting the variants that are overlapping
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste(outc)
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     exp_overlap <- exp[exp$pos %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     exp_overlap_2 <- exp_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
     exp_overlap_2$beta <- exp_overlap_2$beta * -1
     exp_overlap_2$af_alt  <- 1 - exp_overlap_2$af_alt
     colnames(exp_overlap_2)[colnames(exp_overlap_2) %in% c("ref", "alt")] <- c("alt", "ref")
     exp_overlap <- rbind(exp_overlap, exp_overlap_2)
     exp_overlap$ID <- paste(exp_overlap$"#chrom", exp_overlap$pos, exp_overlap$alt, exp_overlap$ref, sep=":")
     exp_overlap$phen <- i
     exp_overlap <- as.data.frame(exp_overlap)
     exp_overlap <- format_data(exp_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     rm(exp_overlap_2)
     dat <- harmonise_data(exposure_dat=exp_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
     dat <- dat[dat$mr_keep==T,]
     
     dat$chrom_pos <- paste(dat$chr.exposure, dat$pos.exposure, sep="_")
     outcome_rsid$chrom_pos <- paste(outcome_rsid$"#chrom", outcome_rsid$pos, sep="_")
     dat <- merge(dat, outcome_rsid, by="chrom_pos", all.x=T)
     
     rsq <- 0.001
     
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     rm(exp_overlap, outcome_overlap, outcome_rsid, clump)
     
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=i, outc=paste(outc), gene=gene_unique, nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval, pval_thresh=pvalue, rsq_thresh=rsq, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
     } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), gene=gene_unique, nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
      } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), gene=gene_unique, nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       }
      }
     
  if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)
   }
   }
   }
   }
   }
   }
   
   write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/3a_downstr_thromb_longcovid_indivgenes.csv") 
   write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_coag_compl_mr/3b_downstr_thromb_longcovid_indivgenes_instr.csv") 
}

 


#### 5 - Downstream analysis - evaluation of proteins on long COVID risk (UKB-PPP) ####

 library(synapser) 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)
 
 setwd("/broad/hptmp/aschuerm/synapser")
  synLogin(email='artschuermans', authToken="eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIl0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTcyMTExOTM2MSwiaWF0IjoxNzIxMTE5MzYxLCJqdGkiOiIxMDA3MyIsInN1YiI6IjM0NzQzMjEifQ.vDm7V43ilOyV0UM9DYSmDepfqU7S0M_yrdBbDy38MVzT-fJY2jEzlNXfcV6W5o3y2-mhOLd6jGV8H6sURi9xMBGejd5E_Rk4u-YJdwIiVnu-FH24J-pr0zDQSC7M6zLbngvNNcmSLHalJ4Ra4AlIAq0mxrSEKt0f91QL3ZxBYwKssyJpq9fbPNcBNGDSsUoUzJKwF__MYwZYkzv69hvFKttq0yekaq1-fHgnkDadH6x4eyyIXqywGW8gFD6M9xmHNuoByuUapDQnUC66e5AMzVbV1HOpc_0zyuFLgAdT587A7T0Q0xB5BXh3L2GVf68kif9S0RY4PCFL57iQy6GWsA")
 
 list <- fread("/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/3a_downstr_thromb_longcovid_indivgenes.csv")
    list <- list[list$method %in% c("Wald ratio", "Inverse variance weighted (correlation inc)"),]
    list <- list[list$pval < 0.1 | list$b > 0.5,]
 sumstats_info <- fread("/medpop/esp2/aschuerm/ukb_proteomics_cvd/input_files/olink_protein_map_3k_v1.tsv")  
    sumstats_info <- sumstats_info[sumstats_info$Assay %in% list$gene,]
    rm(list)
    
 longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)
    
 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, pval_thresh=NA, rsq_thresh=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA, rsq=NA, pval=NA)[-1,]

 
 for (i in sumstats_info$Code){

   syn_code <- synGet(entity=i, downloadLocation = paste(getwd(), "sumstat_prot", sep="/")) 
   untar(paste(syn_code$path), list=F, exdir=paste(syn_code$cacheDir))
   chrom_u <- fread(paste0(syn_code$cacheDir, "/", gsub(".tar", "", sumstats_info[sumstats_info$Code==i,]$Docname[1]), "/", 
                  "discovery_chr", sumstats_info[sumstats_info$Code==i,]$chr[1], "_", sumstats_info[sumstats_info$Code==i,]$UKBPPP_ProteinID[1], 
                  ":", sumstats_info[sumstats_info$Code==i,]$Panel[1], ".gz"))
   
   chrom_u <- chrom_u[chrom_u$GENPOS > (sumstats_info[sumstats_info$Code==i,]$gene_start[1] - 200000) &
                    chrom_u$GENPOS < (sumstats_info[sumstats_info$Code==i,]$gene_end[1] + 200000), ]
   chrom_u$P <- 10^-chrom_u$LOG10P
   
   for (pvalue in c(5e-4, 5e-6, 5e-8, 5e-10)) {
     
   chrom <- chrom_u[chrom_u$P<pvalue,]
   chrom$CHROM <- ifelse(chrom$CHROM==23, "X", chrom$CHROM)
   
   if (is.null(chrom) || nrow(chrom) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
    
   for (outc in "longcovid") {
     outcome <- get(outc)
     outcome_overlap <- outcome[outcome$`#chrom` == sumstats_info[sumstats_info$Code==i,]$chr[1],]                                      
     outcome_overlap <- outcome_overlap[outcome_overlap$pos %in% chrom$GENPOS,]
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste("longcovid")
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     chrom_overlap <- chrom[chrom$GENPOS %in% outcome_overlap$pos.outcome,]
     chrom_overlap_2 <- chrom_overlap
     chrom_overlap_2$BETA <- chrom_overlap_2$BETA * -1
     chrom_overlap_2$A1FREQ  <- 1 - chrom_overlap_2$A1FREQ
     colnames(chrom_overlap_2)[colnames(chrom_overlap_2) %in% c("ALLELE0", "ALLELE1")] <- c("ALLELE1", "ALLELE0")
     chrom_overlap <- rbind(chrom_overlap, chrom_overlap_2)
     chrom_overlap$ID <- paste(chrom_overlap$CHROM, chrom_overlap$GENPOS, chrom_overlap$ALLELE1, chrom_overlap$ALLELE0, sep=":")
     chrom_overlap$phen <- sumstats_info[sumstats_info$Code==i,]$Assay[1]
     chrom_overlap <- as.data.frame(chrom_overlap)
     chrom_overlap <- format_data(chrom_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                                effect_allele_col="ALLELE1", other_allele_col="ALLELE0", pval_col="LOG10P", chr_col="CHROM", samplesize_col="N", pos_col="GENPOS", log_pval=T)
     rm(chrom_overlap_2, chrom)
     
    for (rsq in c(0.1, 0.01, 0.001, 0.0001)) {

     dat <- harmonise_data(exposure_dat=chrom_overlap, outcome_dat=outcome_overlap)
     dat <- dat[dat$mr_keep==T,]
     
     if (sumstats_info[sumstats_info$Code==i,]$chr[1]=="X"){
       dat <- merge(dat, ref_rsid[,c("V1", "V2", "V3")], by.x="pos.exposure", by.y="V2", all.x=T) 
       colnames(dat)[colnames(dat) %in% c("V1", "V3")] <- c("#chrom", "rsids")
     } else {
      dat <- merge(dat, outcome_rsid, by.x="pos.exposure", by.y="pos", all.x=T)
     }
     
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")
     dat <- dat[order(dat$pval.exposure),]
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     dat$rsq <- rsq
     dat$pval <- pvalue
     rm(clump)
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste("longcovid"), nsnp=results_mr$nsnp, 
                             method=results_mr$method, pval_thresh=pvalue, rsq_thresh=rsq, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval)
     } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       ld <- ld_matrix(dat$SNP, bfile="/medpop/esp2/aschuerm/tools/g1000_eur", plink_bin=genetics.binaRies::get_plink_binary())
       rownames(ld) <- gsub("\\_.*", "", rownames(ld))
       dat_dupl <- dat[!(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        colnames(dat_dupl)[which(colnames(dat_dupl) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
        dat_dupl$beta.exposure <- dat_dupl$beta.exposure*-1
        dat_dupl$beta.outcome <- dat_dupl$beta.outcome*-1
        dat_dupl$eaf.exposure <- 1-dat_dupl$eaf.exposure
        dat_dupl$eaf.outcome <- 1-dat_dupl$eaf.outcome
        dat <- dat[(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        dat <- rbind(dat, dat_dupl)
        dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       dat2 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome,
                                                 correlation = ld)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_ivw_corr <- MendelianRandomization::mr_ivw(dat2, correl = TRUE)
       results0 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue)
       results1 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw_corr@SNPs, 
                             method="Inverse variance weighted (correlation inc)", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw_corr@Estimate, 
                             se=output_mr_ivw_corr@StdError, pval=output_mr_ivw_corr@Pvalue)
       results <- rbind(results0, results1)
       rm(results0, results1)
      } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       ld <- ld_matrix(dat$SNP, bfile="/medpop/esp2/aschuerm/tools/g1000_eur", plink_bin=genetics.binaRies::get_plink_binary())
       rownames(ld) <- gsub("\\_.*", "", rownames(ld))
       dat_dupl <- dat[!(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        colnames(dat_dupl)[which(colnames(dat_dupl) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
        dat_dupl$beta.exposure <- dat_dupl$beta.exposure*-1
        dat_dupl$beta.outcome <- dat_dupl$beta.outcome*-1
        dat_dupl$eaf.exposure <- 1-dat_dupl$eaf.exposure
        dat_dupl$eaf.outcome <- 1-dat_dupl$eaf.outcome
        dat <- dat[(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        dat <- rbind(dat, dat_dupl)
        dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       dat2 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome,
                                                 correlation = ld)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_ivw_corr <- MendelianRandomization::mr_ivw(dat2, correl = TRUE)
       output_mr_egger_corr <- MendelianRandomization::mr_egger(dat2, correl = TRUE)
       results0 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue)
       results1 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw_corr@SNPs, 
                             method="Inverse variance weighted (correlation inc)", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw_corr@Estimate, 
                             se=output_mr_ivw_corr@StdError, pval=output_mr_ivw_corr@Pvalue)
       results2 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_egger_corr@SNPs, 
                             method="Egger (correlation inc)", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger_corr@Estimate, 
                             se=output_mr_egger_corr@StdError.Est, pval=output_mr_egger_corr@Pvalue.Est)
       results3 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_egger_corr@SNPs, 
                             method="Egger intercept (correlation inc)", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger_corr@Intercept, 
                             se=output_mr_egger_corr@StdError.Int, pval=output_mr_egger_corr@Pvalue.Int)
       results <- rbind(results0, results1, results2, results3)
       rm(results0, results1, results2, results3)
       
        }
        }
       
    if (is.null(results) || nrow(results) == 0) {
          print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
     } else {
        df_sum <- rbind(df_sum, results)
        df_instr <- rbind(df_instr, dat[,-c(34)])
        rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)
     }
     }
     }
     }
     }
   }  
   
     write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/4a_downstr_proteins_longcovid_ukbppp.csv") 
     write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/4b_downstr_proteins_longcovid_ukbppp.csv") 
     print(paste(sumstats_info[sumstats_info$Code==i,]$Assay[1], "done"))
     unlink(paste0(gsub("\\\\[^\\\\]*$", "", syn_code$cacheDir)), recursive=T)
 }
 
  
#### 6 - Downstream analysis - evaluation of proteins on long COVID risk (UKB-PPP; only intragenic variants) ####

 library(synapser) 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)
 
 setwd("/broad/hptmp/aschuerm/synapser")
  synLogin(email='artschuermans', authToken="eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIl0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTcyMTExOTM2MSwiaWF0IjoxNzIxMTE5MzYxLCJqdGkiOiIxMDA3MyIsInN1YiI6IjM0NzQzMjEifQ.vDm7V43ilOyV0UM9DYSmDepfqU7S0M_yrdBbDy38MVzT-fJY2jEzlNXfcV6W5o3y2-mhOLd6jGV8H6sURi9xMBGejd5E_Rk4u-YJdwIiVnu-FH24J-pr0zDQSC7M6zLbngvNNcmSLHalJ4Ra4AlIAq0mxrSEKt0f91QL3ZxBYwKssyJpq9fbPNcBNGDSsUoUzJKwF__MYwZYkzv69hvFKttq0yekaq1-fHgnkDadH6x4eyyIXqywGW8gFD6M9xmHNuoByuUapDQnUC66e5AMzVbV1HOpc_0zyuFLgAdT587A7T0Q0xB5BXh3L2GVf68kif9S0RY4PCFL57iQy6GWsA")
 
 list <- fread("/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/3a_downstr_thromb_longcovid_indivgenes.csv")
    list <- list[list$method %in% c("Wald ratio", "Inverse variance weighted (correlation inc)"),]
    list <- list[list$pval < 0.1 | list$b > 0.5,]
 sumstats_info <- fread("/medpop/esp2/aschuerm/ukb_proteomics_cvd/input_files/olink_protein_map_3k_v1.tsv")  
    sumstats_info <- sumstats_info[sumstats_info$Assay %in% list$gene,]
    rm(list)
    
 longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)
    
 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, pval_thresh=NA, rsq_thresh=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA, rsq=NA, pval=NA)[-1,]

 
 for (i in sumstats_info$Code){

   syn_code <- synGet(entity=i, downloadLocation = paste(getwd(), "sumstat_prot", sep="/")) 
   untar(paste(syn_code$path), list=F, exdir=paste(syn_code$cacheDir))
   chrom_u <- fread(paste0(syn_code$cacheDir, "/", gsub(".tar", "", sumstats_info[sumstats_info$Code==i,]$Docname[1]), "/", 
                  "discovery_chr", sumstats_info[sumstats_info$Code==i,]$chr[1], "_", sumstats_info[sumstats_info$Code==i,]$UKBPPP_ProteinID[1], 
                  ":", sumstats_info[sumstats_info$Code==i,]$Panel[1], ".gz"))
   
   chrom_u <- chrom_u[chrom_u$GENPOS > (sumstats_info[sumstats_info$Code==i,]$gene_start[1] - 0) &
                    chrom_u$GENPOS < (sumstats_info[sumstats_info$Code==i,]$gene_end[1] + 0), ]
   chrom_u$P <- 10^-chrom_u$LOG10P
   
   for (pvalue in c(5e-4, 5e-6, 5e-8, 5e-10)) {
     
   chrom <- chrom_u[chrom_u$P<pvalue,]
   chrom$CHROM <- ifelse(chrom$CHROM==23, "X", chrom$CHROM)
   
   if (is.null(chrom) || nrow(chrom) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
    
   for (outc in "longcovid") {
     outcome <- get(outc)
     outcome_overlap <- outcome[outcome$`#chrom` == sumstats_info[sumstats_info$Code==i,]$chr[1],]                                      
     outcome_overlap <- outcome_overlap[outcome_overlap$pos %in% chrom$GENPOS,]
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste("longcovid")
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     chrom_overlap <- chrom[chrom$GENPOS %in% outcome_overlap$pos.outcome,]
     chrom_overlap_2 <- chrom_overlap
     chrom_overlap_2$BETA <- chrom_overlap_2$BETA * -1
     chrom_overlap_2$A1FREQ  <- 1 - chrom_overlap_2$A1FREQ
     colnames(chrom_overlap_2)[colnames(chrom_overlap_2) %in% c("ALLELE0", "ALLELE1")] <- c("ALLELE1", "ALLELE0")
     chrom_overlap <- rbind(chrom_overlap, chrom_overlap_2)
     chrom_overlap$ID <- paste(chrom_overlap$CHROM, chrom_overlap$GENPOS, chrom_overlap$ALLELE1, chrom_overlap$ALLELE0, sep=":")
     chrom_overlap$phen <- sumstats_info[sumstats_info$Code==i,]$Assay[1]
     chrom_overlap <- as.data.frame(chrom_overlap)
     chrom_overlap <- format_data(chrom_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                                effect_allele_col="ALLELE1", other_allele_col="ALLELE0", pval_col="LOG10P", chr_col="CHROM", samplesize_col="N", pos_col="GENPOS", log_pval=T)
     rm(chrom_overlap_2, chrom)
     
    for (rsq in c(0.1, 0.01, 0.001, 0.0001)) {

     dat <- harmonise_data(exposure_dat=chrom_overlap, outcome_dat=outcome_overlap)
     dat <- dat[dat$mr_keep==T,]
     
     if (sumstats_info[sumstats_info$Code==i,]$chr[1]=="X"){
       dat <- merge(dat, ref_rsid[,c("V1", "V2", "V3")], by.x="pos.exposure", by.y="V2", all.x=T) 
       colnames(dat)[colnames(dat) %in% c("V1", "V3")] <- c("#chrom", "rsids")
     } else {
      dat <- merge(dat, outcome_rsid, by.x="pos.exposure", by.y="pos", all.x=T)
     }
     
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")
     dat <- dat[order(dat$pval.exposure),]
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     dat$rsq <- rsq
     dat$pval <- pvalue
     rm(clump)
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste("longcovid"), nsnp=results_mr$nsnp, 
                             method=results_mr$method, pval_thresh=pvalue, rsq_thresh=rsq, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval)
     } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       ld <- ld_matrix(dat$SNP, bfile="/medpop/esp2/aschuerm/tools/g1000_eur", plink_bin=genetics.binaRies::get_plink_binary())
       rownames(ld) <- gsub("\\_.*", "", rownames(ld))
       dat_dupl <- dat[!(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        colnames(dat_dupl)[which(colnames(dat_dupl) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
        dat_dupl$beta.exposure <- dat_dupl$beta.exposure*-1
        dat_dupl$beta.outcome <- dat_dupl$beta.outcome*-1
        dat_dupl$eaf.exposure <- 1-dat_dupl$eaf.exposure
        dat_dupl$eaf.outcome <- 1-dat_dupl$eaf.outcome
        dat <- dat[(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        dat <- rbind(dat, dat_dupl)
        dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       dat2 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome,
                                                 correlation = ld)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_ivw_corr <- MendelianRandomization::mr_ivw(dat2, correl = TRUE)
       results0 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue)
       results1 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw_corr@SNPs, 
                             method="Inverse variance weighted (correlation inc)", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw_corr@Estimate, 
                             se=output_mr_ivw_corr@StdError, pval=output_mr_ivw_corr@Pvalue)
       results <- rbind(results0, results1)
       rm(results0, results1)
      } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       ld <- ld_matrix(dat$SNP, bfile="/medpop/esp2/aschuerm/tools/g1000_eur", plink_bin=genetics.binaRies::get_plink_binary())
       rownames(ld) <- gsub("\\_.*", "", rownames(ld))
       dat_dupl <- dat[!(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        colnames(dat_dupl)[which(colnames(dat_dupl) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
        dat_dupl$beta.exposure <- dat_dupl$beta.exposure*-1
        dat_dupl$beta.outcome <- dat_dupl$beta.outcome*-1
        dat_dupl$eaf.exposure <- 1-dat_dupl$eaf.exposure
        dat_dupl$eaf.outcome <- 1-dat_dupl$eaf.outcome
        dat <- dat[(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        dat <- rbind(dat, dat_dupl)
        dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       dat2 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome,
                                                 correlation = ld)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       output_mr_ivw_corr <- MendelianRandomization::mr_ivw(dat2, correl = TRUE)
       output_mr_egger_corr <- MendelianRandomization::mr_egger(dat2, correl = TRUE)
       results0 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue)
       results1 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw_corr@SNPs, 
                             method="Inverse variance weighted (correlation inc)", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_ivw_corr@Estimate, 
                             se=output_mr_ivw_corr@StdError, pval=output_mr_ivw_corr@Pvalue)
       results2 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_egger_corr@SNPs, 
                             method="Egger (correlation inc)", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger_corr@Estimate, 
                             se=output_mr_egger_corr@StdError.Est, pval=output_mr_egger_corr@Pvalue.Est)
       results3 <- data.frame(outc=paste(outc), exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_egger_corr@SNPs, 
                             method="Egger intercept (correlation inc)", pval_thresh=pvalue, rsq_thresh=rsq, b=output_mr_egger_corr@Intercept, 
                             se=output_mr_egger_corr@StdError.Int, pval=output_mr_egger_corr@Pvalue.Int)
       results <- rbind(results0, results1, results2, results3)
       rm(results0, results1, results2, results3)
       
        }
        }
       
    if (is.null(results) || nrow(results) == 0) {
          print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
     } else {
        df_sum <- rbind(df_sum, results)
        df_instr <- rbind(df_instr, dat[,-c(34)])
        rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)
     }
     }
     }
     }
     }
   }  
   
     write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/4a_downstr_proteins_longcovid_ukbppp_intragenic.csv") 
     write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/4b_downstr_proteins_longcovid_ukbppp_intragenic_instr.csv") 
     print(paste(sumstats_info[sumstats_info$Code==i,]$Assay[1], "done"))
     unlink(paste0(gsub("\\\\[^\\\\]*$", "", syn_code$cacheDir)), recursive=T)
 }
 
  
#### 7 - Downstream analysis - evaluation of proteins on long COVID risk (VTE; coloc) ####
 
 library(synapser) 
 library(coloc) 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)
 
 sumstats_info <- fread("/medpop/esp2/aschuerm/ukb_proteomics_cvd/input_files/olink_protein_map_3k_v1.tsv")                   # Information on the sum stats / protein codes etc.; functions as a linker file
    targets <- c('F2R') 
    sumstats_info <- sumstats_info[sumstats_info$Assay %in% targets,]
    rm(targets)
    
 longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)
    
 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   # This is a linker file to match RSIDs with chromosome positions for the X chromosome
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

  vte <- fread("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_VTE.gz")
    vte <- vte[vte$af_alt > 0.01 & vte$af_alt < 0.99,]
    vte$"#chrom" <- ifelse(vte$"#chrom"==23, "X", vte$"#chrom")
    colnames(vte)[which(colnames(vte)=="pos")] <- 'GENPOS'
    colnames(vte)[which(colnames(vte)=="af_alt")] <- 'A1FREQ'
    colnames(vte)[which(colnames(vte)=="alt")] <- 'ALLELE1'
    colnames(vte)[which(colnames(vte)=="ref")] <- 'ALLELE0'
    colnames(vte)[which(colnames(vte)=="pval")] <- 'P'
    colnames(vte)[which(colnames(vte)=="#chrom")] <- '#CHROM'
    colnames(vte)[which(colnames(vte)=="beta")] <- 'BETA'
    colnames(vte)[which(colnames(vte)=="sebeta")] <- 'SE'
    
 chrom_get <- "vte"
 outcome_get <- "longcovid"
 
    for (window in c(200000)) {
 
    outcome <- get(outcome_get)
    chrom_u <- get(chrom_get)
    
    chrom <- chrom_u[chrom_u$GENPOS > (sumstats_info[sumstats_info$Assay=="F2R",]$gene_start[1] - window) &                              # Selecting the cis region only (here defined as 1Mb before or after the protein-encoding region)
                    chrom_u$GENPOS < (sumstats_info[sumstats_info$Assay=="F2R",]$gene_end[1] + window) &
                    chrom_u$"#CHROM" == sumstats_info[sumstats_info$Assay=="F2R",]$chr[1]  , ]
    
    outcome_overlap <- outcome[outcome$`#chrom` == sumstats_info[sumstats_info$Assay=="F2R",]$chr[1],]                                      # Only selecting the chromosome of interest to speed stuff up downstream from here
    outcome_overlap <- outcome_overlap[outcome_overlap$pos %in% chrom$GENPOS,]                                                 # Only selecting the variants that are overlapping
    outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
    outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
    outcome_overlap$phen <- paste(outcome_get)
    outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
    outcome_overlap <- as.data.frame(outcome_overlap)
    outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                              effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
    
    chrom_overlap <- chrom[chrom$GENPOS %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
    chrom_overlap_2 <- chrom_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
    chrom_overlap_2$BETA <- chrom_overlap_2$BETA * -1
    chrom_overlap_2$A1FREQ  <- 1 - chrom_overlap_2$A1FREQ
    colnames(chrom_overlap_2)[colnames(chrom_overlap_2) %in% c("ALLELE0", "ALLELE1")] <- c("ALLELE1", "ALLELE0")
    chrom_overlap <- rbind(chrom_overlap, chrom_overlap_2)
    chrom_overlap$ID <- paste(chrom_overlap$"#CHROM", chrom_overlap$GENPOS, chrom_overlap$ALLELE1, chrom_overlap$ALLELE0, sep=":")
    chrom_overlap$phen <- sumstats_info[sumstats_info$Assay=="F2R",]$Assay[1]
    chrom_overlap <- as.data.frame(chrom_overlap)
    chrom_overlap <- format_data(chrom_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                              effect_allele_col="ALLELE1", other_allele_col="ALLELE0", pval_col="P", chr_col="#CHROM", pos_col="GENPOS", log_pval=F)
    rm(chrom_overlap_2)
    dat <- harmonise_data(exposure_dat=chrom_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
    
     if (sumstats_info[sumstats_info$Assay=="F2R",]$chr[1]=="X"){                                                                    # This little if-else-statement just makes sure that you get the appropriate RSIDs for each variant; because the X-chromosome requires an additional file, this one is in a separate loop
       dat <- merge(dat, ref_rsid[,c("V1", "V2", "V3")], by.x="pos.exposure", by.y="V2", all.x=T) 
       colnames(dat)[colnames(dat) %in% c("V1", "V3")] <- c("#chrom", "rsids")
     } else {
      dat <- merge(dat, outcome_rsid, by.x="pos.exposure", by.y="pos", all.x=T)
     }
     
    colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
    dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
    dat <- dat[!duplicated(dat$SNP),]
    dat <- dat[!duplicated(dat$pos_id),]
       ld <- ld_matrix(dat$SNP, bfile="/medpop/esp2/aschuerm/tools/g1000_eur", plink_bin=genetics.binaRies::get_plink_binary())
       rownames(ld) <- gsub("\\_.*", "", rownames(ld))
       dat_dupl <- dat[!(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        colnames(dat_dupl)[which(colnames(dat_dupl) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
        dat_dupl$beta.exposure <- dat_dupl$beta.exposure*-1
        dat_dupl$beta.outcome <- dat_dupl$beta.outcome*-1
        dat_dupl$eaf.exposure <- 1-dat_dupl$eaf.exposure
        dat_dupl$eaf.outcome <- 1-dat_dupl$eaf.outcome
        dat <- dat[(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        dat <- rbind(dat, dat_dupl)
        dat <- dat[dat$SNP %in% rownames(ld),]
        dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
       colnames(ld) <- gsub("\\_.*", "", colnames(ld))
   
    dat_exp <- list(beta=dat$beta.exposure, varbeta=dat$se.exposure^2, snp=dat$SNP, position=dat$pos.exposure, type="quant", sdY=1, LD=ld, N=unique(dat$samplesize.exposure))
    dat_outc <- list(beta=dat$beta.outcome, varbeta=dat$se.outcome^2, snp=dat$SNP, position=dat$pos.outcome, type="cc", LD=ld, N=3018+994582, s=3018/(3018+994582))
    
    coloc_outc <- coloc.abf(dataset1=dat_exp, dataset2=dat_outc)
    
    results <- data.frame(exp=sumstats_info[sumstats_info$Assay=="F2R",]$Assay[1], outc=paste(outcome_get), nsnps=coloc_outc$summary["nsnps"],
                      pp_h0=coloc_outc$summary["PP.H0.abf"], pp_h1=coloc_outc$summary["PP.H1.abf"], pp_h2=coloc_outc$summary["PP.H2.abf"],
                      pp_h3=coloc_outc$summary["PP.H3.abf"], pp_h4=coloc_outc$summary["PP.H4.abf"],
                      pp_h3_cond=coloc_outc$summary["PP.H3.abf"]/(coloc_outc$summary["PP.H3.abf"]+coloc_outc$summary["PP.H4.abf"]), 
                      pp_h4_cond=coloc_outc$summary["PP.H4.abf"]/(coloc_outc$summary["PP.H3.abf"]+coloc_outc$summary["PP.H4.abf"]), 
                      window=window)
    df_sum <- rbind(df_sum, results)
    rm(results, outcome, outcome_overlap, outcome_rsid, chrom_overlap, dat, dat_exp, dat_outc, coloc_outc)

    }
  
    write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/4a_downstr_proteins_longcovid_ukbppp_coloc.csv") 
    print(paste(sumstats_info[sumstats_info$Assay=="F2R",]$Assay[1], "done"))
    unlink(paste0(gsub("\\\\[^\\\\]*$", "", syn_code$cacheDir)), recursive=T)
  
 
     
     

#### 8 - Downstream analysis - evaluation of proteins on long COVID risk (UKB-PPP; coloc) ####
 
 library(synapser) 
 library(coloc) 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)
 
 setwd("/broad/hptmp/aschuerm/synapser")
  synLogin(email='artschuermans', authToken="eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIl0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTcyMTExOTM2MSwiaWF0IjoxNzIxMTE5MzYxLCJqdGkiOiIxMDA3MyIsInN1YiI6IjM0NzQzMjEifQ.vDm7V43ilOyV0UM9DYSmDepfqU7S0M_yrdBbDy38MVzT-fJY2jEzlNXfcV6W5o3y2-mhOLd6jGV8H6sURi9xMBGejd5E_Rk4u-YJdwIiVnu-FH24J-pr0zDQSC7M6zLbngvNNcmSLHalJ4Ra4AlIAq0mxrSEKt0f91QL3ZxBYwKssyJpq9fbPNcBNGDSsUoUzJKwF__MYwZYkzv69hvFKttq0yekaq1-fHgnkDadH6x4eyyIXqywGW8gFD6M9xmHNuoByuUapDQnUC66e5AMzVbV1HOpc_0zyuFLgAdT587A7T0Q0xB5BXh3L2GVf68kif9S0RY4PCFL57iQy6GWsA") # This is to log in to the Synapse platform, needed to gain access to the protein sum stats

 sumstats_info <- fread("/medpop/esp2/aschuerm/ukb_proteomics_cvd/input_files/olink_protein_map_3k_v1.tsv")                   # Information on the sum stats / protein codes etc.; functions as a linker file
    targets <- c('F2R') 
    sumstats_info <- sumstats_info[sumstats_info$Assay %in% targets,]
    rm(targets)
    
 longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)
    
 ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   # This is a linker file to match RSIDs with chromosome positions for the X chromosome
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 i <- sumstats_info$Code
   syn_code <- synGet(entity=i, downloadLocation = paste(getwd(), "sumstat_prot", sep="/")) # Downloading the summary statistics for the protein of interest
   untar(paste(syn_code$path), list=F, exdir=paste(syn_code$cacheDir))
   chrom_olink <- fread(paste0(syn_code$cacheDir, "/", gsub(".tar", "", sumstats_info[sumstats_info$Code==i,]$Docname[1]), "/", 
                  "discovery_chr", sumstats_info[sumstats_info$Code==i,]$chr[1], "_", sumstats_info[sumstats_info$Code==i,]$UKBPPP_ProteinID[1], 
                  ":", sumstats_info[sumstats_info$Code==i,]$Panel[1], ".gz"))
   chrom_olink <- chrom_olink[chrom_olink$GENPOS > (sumstats_info[sumstats_info$Code==i,]$gene_start[1] - 200000) &                              # Selecting the cis region only (here defined as 1Mb before or after the protein-encoding region)
                    chrom_olink$GENPOS < (sumstats_info[sumstats_info$Code==i,]$gene_end[1] + 200000), ]
   chrom_olink$P <- 10^-chrom_olink$LOG10P
   chrom_olink$CHROM <- ifelse(chrom_olink$CHROM==23, "X", chrom_olink$CHROM)                                                                     # Transferring 23rd chromosome to "X" for consistency
   chrom_olink <- chrom_olink[chrom_olink$A1FREQ > 0.01 & chrom_olink$A1FREQ < 0.99,]
   
 chrom_get <- "chrom_olink"
 outcome_get <- "longcovid"
 
    for (window in c(200000)) {
 
    outcome <- get(outcome_get)
    chrom_u <- get(chrom_get)
    
    chrom <- chrom_u[chrom_u$GENPOS > (sumstats_info[sumstats_info$Code==i,]$gene_start[1] - window) &                              # Selecting the cis region only (here defined as 1Mb before or after the protein-encoding region)
                    chrom_u$GENPOS < (sumstats_info[sumstats_info$Code==i,]$gene_end[1] + window), ]
    
    outcome_overlap <- outcome[outcome$`#chrom` == sumstats_info[sumstats_info$Code==i,]$chr[1],]                                      # Only selecting the chromosome of interest to speed stuff up downstream from here
    outcome_overlap <- outcome_overlap[outcome_overlap$pos %in% chrom$GENPOS,]                                                 # Only selecting the variants that are overlapping
    outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
    outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
    outcome_overlap$phen <- paste(outcome_get)
    outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
    outcome_overlap <- as.data.frame(outcome_overlap)
    outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                              effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
    
    chrom_overlap <- chrom[chrom$GENPOS %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
    chrom_overlap_2 <- chrom_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
    chrom_overlap_2$BETA <- chrom_overlap_2$BETA * -1
    chrom_overlap_2$A1FREQ  <- 1 - chrom_overlap_2$A1FREQ
    colnames(chrom_overlap_2)[colnames(chrom_overlap_2) %in% c("ALLELE0", "ALLELE1")] <- c("ALLELE1", "ALLELE0")
    chrom_overlap <- rbind(chrom_overlap, chrom_overlap_2)
    chrom_overlap$ID <- paste(chrom_overlap$CHROM, chrom_overlap$GENPOS, chrom_overlap$ALLELE1, chrom_overlap$ALLELE0, sep=":")
    chrom_overlap$phen <- sumstats_info[sumstats_info$Code==i,]$Assay[1]
    chrom_overlap <- as.data.frame(chrom_overlap)
    chrom_overlap <- format_data(chrom_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                              effect_allele_col="ALLELE1", other_allele_col="ALLELE0", pval_col="P", chr_col="CHROM", samplesize_col="N", pos_col="GENPOS", log_pval=F)
    rm(chrom_overlap_2)
    dat <- harmonise_data(exposure_dat=chrom_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
    
     if (sumstats_info[sumstats_info$Code==i,]$chr[1]=="X"){                                                                    # This little if-else-statement just makes sure that you get the appropriate RSIDs for each variant; because the X-chromosome requires an additional file, this one is in a separate loop
       dat <- merge(dat, ref_rsid[,c("V1", "V2", "V3")], by.x="pos.exposure", by.y="V2", all.x=T) 
       colnames(dat)[colnames(dat) %in% c("V1", "V3")] <- c("#chrom", "rsids")
     } else {
      dat <- merge(dat, outcome_rsid, by.x="pos.exposure", by.y="pos", all.x=T)
     }
     
    colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
    dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
    dat <- dat[!duplicated(dat$SNP),]
    dat <- dat[!duplicated(dat$pos_id),]
       ld <- ld_matrix(dat$SNP, bfile="/medpop/esp2/aschuerm/tools/g1000_eur", plink_bin=genetics.binaRies::get_plink_binary())
       rownames(ld) <- gsub("\\_.*", "", rownames(ld))
       dat_dupl <- dat[!(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        colnames(dat_dupl)[which(colnames(dat_dupl) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
        dat_dupl$beta.exposure <- dat_dupl$beta.exposure*-1
        dat_dupl$beta.outcome <- dat_dupl$beta.outcome*-1
        dat_dupl$eaf.exposure <- 1-dat_dupl$eaf.exposure
        dat_dupl$eaf.outcome <- 1-dat_dupl$eaf.outcome
        dat <- dat[(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        dat <- rbind(dat, dat_dupl)
        dat <- dat[dat$SNP %in% rownames(ld),]
        dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
       colnames(ld) <- gsub("\\_.*", "", colnames(ld))
   
    dat_exp <- list(beta=dat$beta.exposure, varbeta=dat$se.exposure^2, snp=dat$SNP, position=dat$pos.exposure, type="quant", sdY=1, LD=ld, N=unique(dat$samplesize.exposure))
    dat_outc <- list(beta=dat$beta.outcome, varbeta=dat$se.outcome^2, snp=dat$SNP, position=dat$pos.outcome, type="cc", LD=ld, N=3018+994582, s=3018/(3018+994582))
    
    coloc_outc <- coloc.abf(dataset1=dat_exp, dataset2=dat_outc)
    
    results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outcome_get), nsnps=coloc_outc$summary["nsnps"],
                      pp_h0=coloc_outc$summary["PP.H0.abf"], pp_h1=coloc_outc$summary["PP.H1.abf"], pp_h2=coloc_outc$summary["PP.H2.abf"],
                      pp_h3=coloc_outc$summary["PP.H3.abf"], pp_h4=coloc_outc$summary["PP.H4.abf"],
                      pp_h3_cond=coloc_outc$summary["PP.H3.abf"]/(coloc_outc$summary["PP.H3.abf"]+coloc_outc$summary["PP.H4.abf"]), 
                      pp_h4_cond=coloc_outc$summary["PP.H4.abf"]/(coloc_outc$summary["PP.H3.abf"]+coloc_outc$summary["PP.H4.abf"]), 
                      window=window)
    df_sum <- rbind(df_sum, results)
    rm(results, outcome, outcome_overlap, outcome_rsid, chrom_overlap, dat, dat_exp, dat_outc, coloc_outc)

    }
  
    write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/4a_downstr_proteins_longcovid_ukbppp_coloc.csv") 
    print(paste(sumstats_info[sumstats_info$Code==i,]$Assay[1], "done"))
    unlink(paste0(gsub("\\\\[^\\\\]*$", "", syn_code$cacheDir)), recursive=T)
  
 
     
     

#### 9 - Downstream analysis - evaluation of transcripts on long COVID ####
   # 9.a - GTEx preparation (EUR ancestry) ####
  
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(coloc)
  library(TwoSampleMR)
  
  sumstats_info <- fread("/medpop/esp2/aschuerm/coagulation_longcovid/input/map_gtex.csv")

  list_gtex_tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", 
                         "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", 
                         "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", 
                         "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", 
                         "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", 
                         "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", 
                         "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", 
                         "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", 
                         "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", 
                         "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", 
                         "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", 
                         "Whole_Blood")
  
  for (gtex in list_gtex_tissues) {
       object_name <- paste0("gtex_", gtex)
       gtex_sumstats <- fread(paste0("/medpop/esp2/aschuerm/resources/gtex/GTEx_Analysis_v8_eQTL_EUR/eqtls/", gtex, ".v8.EUR.signif_pairs.txt.gz"))
       gtex_sumstats <- gtex_sumstats[,-c("pval_nominal_threshold", "min_pval_nominal", "pval_beta")]
       gtex_sumstats$gene_id <- sub("\\..*", "", gtex_sumstats$phenotype_id)
       gtex_sumstats <- gtex_sumstats[gtex_sumstats$gene_id %in% sumstats_info$ensembl_gene_id,]
       gtex_sumstats <- separate(gtex_sumstats, variant_id, into = c("chromosome", "position", "ref", "alt", "build"), sep = "_")
       gtex_sumstats$chromosome <- sub("chr", "", gtex_sumstats$chromosome)
       write.csv(gtex_sumstats, paste0("/medpop/esp2/aschuerm/coagulation_longcovid/input/gtex/eur/", gtex, "_coagulation_complement_gtex_eur.csv"), row.names=F)
  }
    
  
   # 9.b - GTEx MR ####

  library(data.table)
  library(tidyr)
  library(dplyr)
  library(coloc)
  library(TwoSampleMR)
  library(ieugwasr)

  sumstats_info <- fread("/medpop/esp2/aschuerm/coagulation_longcovid/input/map_gtex.csv")

  list_gtex_tissues <- c("Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Lung")
   
longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)
  
  ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   # This is a linker file to match RSIDs with chromosome positions for the X chromosome
  df_sum <- data.frame(exp=NA)[-1,]
  df_instr <- data.frame(pos.exposure=NA)[-1,]
   
  
  for (gtex in list_gtex_tissues) {
    object_name <- paste0("gtex_", gtex)
    gtex_sumstats <- fread(paste0("/medpop/esp2/aschuerm/coagulation_longcovid/input/gtex/eur/", gtex, "_coagulation_complement_gtex_eur.csv"))
    assign(object_name, gtex_sumstats)
  }


 for (i in "F2R"){

  for (gtex in list_gtex_tissues) {
      gtex_cisregion <- get(paste0("gtex_", gtex))
      gtex_cisregion <- gtex_cisregion[gtex_cisregion$gene_id==sumstats_info$ensembl_gene_id[sumstats_info$external_gene_name==i],]
      if (nrow(gtex_cisregion)==0) {
        print("Skipping a GTEX tissue for a protein")
      } else {
      
      for (window in c(0, 50000, 100000, 200000)) {
      gtex_cisregion_2 <- gtex_cisregion[gtex_cisregion$position > sumstats_info$start_position[sumstats_info$external_gene_name==i]-window & gtex_cisregion$position < sumstats_info$end_position[sumstats_info$external_gene_name==i]+window, ]
 
    if (nrow(gtex_cisregion_2)==0) {
        print("Skipping")
      } else {
        
    for (pvalue in c(5e-4, 5e-6, 5e-8)){
      gtex_cisregion_3 <- data.frame(beta=gtex_cisregion_2$slope, se=gtex_cisregion_2$slope_se, pval=gtex_cisregion_2$pval_nominal, maf=gtex_cisregion_2$maf, 
                                   trait=gtex, pos=gtex_cisregion_2$position, id=paste(gtex_cisregion_2$chromosome, gtex_cisregion_2$position, 
                                   gtex_cisregion_2$alt, gtex_cisregion_2$ref, sep=":"), chrom=gtex_cisregion_2$chromosome, ref=gtex_cisregion_2$ref, alt=gtex_cisregion_2$alt)
      gtex_cisregion_3 <- gtex_cisregion_3[!is.na(gtex_cisregion_3$se),]
      gtex_cisregion_3 <- gtex_cisregion_3[gtex_cisregion_3$pval < pvalue,]                                                                                                # Selecting "region-wide" significant cis-pQTs

     if (nrow(gtex_cisregion_3)==0) {
        print("Skipping")
      } else {
      
  for (outc in c("longcovid")){
    outcome <- get(outc)
   
    outcome_overlap <- outcome[outcome$`#chrom` == sumstats_info$chromosome_name[sumstats_info$external_gene_name==i],]                                      # Only selecting the chromosome of interest to speed stuff up downstream from here
    outcome_overlap <- outcome_overlap[outcome_overlap$pos %in% gtex_cisregion_2$pos,]                                                 # Only selecting the variants that are overlapping
      outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
      outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
      outcome_rsid <- outcome_rsid[!duplicated(outcome_rsid$pos),]
    outcome_overlap$phen <- paste(outc)
    outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")

    gtex_cisregion_overlap <- gtex_cisregion_3[gtex_cisregion_3$pos %in% outcome_overlap$pos,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     
    if (nrow(gtex_cisregion_overlap)==0) {
        print("Skipping")
      
      } else {
    gtex_cisregion_overlap_2 <- gtex_cisregion_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
    gtex_cisregion_overlap_2$beta <- gtex_cisregion_overlap_2$beta * -1
    colnames(gtex_cisregion_overlap_2)[colnames(gtex_cisregion_overlap_2) %in% c("ref", "alt")] <- c("alt", "ref")
    gtex_cisregion_overlap <- rbind(gtex_cisregion_overlap, gtex_cisregion_overlap_2)
    gtex_cisregion_overlap$id <- paste(gtex_cisregion_overlap$chrom, gtex_cisregion_overlap$pos, gtex_cisregion_overlap$alt, gtex_cisregion_overlap$ref, sep=":")
    gtex_cisregion_overlap$phen <- paste0(sumstats_info$external_gene_name[sumstats_info$external_gene_name==i], "_", gtex)
    gtex_cisregion_overlap <- gtex_cisregion_overlap[gtex_cisregion_overlap$id %in% outcome_overlap$id,]
    outcome_overlap <- outcome_overlap[outcome_overlap$id %in% gtex_cisregion_overlap$id,]
    
    gtex_cisregion_overlap <- merge(gtex_cisregion_overlap, outcome_overlap[,c("id", "af_alt")], by="id")
    
    rm(gtex_cisregion_overlap_2)
    
     if (sumstats_info$chromosome_name[sumstats_info$external_gene_name==i]=="X"){                                                                    # This little if-else-statement just makes sure that you get the appropriate RSIDs for each variant; because the X-chromosome requires an additional file, this one is in a separate loop
       outcome_overlap <- outcome_overlap[,-(c("#chrom", "rsids"))]
       outcome_overlap <- merge(outcome_overlap, ref_rsid[,c("V1", "V2", "V3")], by.x="pos", by.y="V2", all.x=T) 
       colnames(outcome_overlap)[colnames(outcome_overlap) %in% c("V1", "V3")] <- c("#chrom", "rsids")
       gtex_cisregion_overlap <- merge(gtex_cisregion_overlap, ref_rsid[,c("V1", "V2", "V3")], by.x="pos", by.y="V2", all.x=T) 
       colnames(gtex_cisregion_overlap)[colnames(gtex_cisregion_overlap) %in% c("V1", "V3")] <- c("#chrom", "rsids")
     } else {
       outcome_overlap <- outcome_overlap[,-(c("#chrom", "rsids"))]
       outcome_overlap <- merge(outcome_overlap, outcome_rsid, by.x="pos", by.y="pos", all.x=T)
       gtex_cisregion_overlap <- merge(gtex_cisregion_overlap, outcome_rsid, by.x="pos", by.y="pos", all.x=T)
     }

    outcome_overlap <- outcome_overlap[!duplicated(outcome_overlap$rsids),]
    outcome_overlap <- outcome_overlap[!is.na(outcome_overlap$rsids),]
    gtex_cisregion_overlap <- gtex_cisregion_overlap[!duplicated(gtex_cisregion_overlap$rsids),]
    gtex_cisregion_overlap <- gtex_cisregion_overlap[!is.na(gtex_cisregion_overlap$rsids),]
    
    outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
    chrom_overlap <- format_data(gtex_cisregion_overlap, type="exposure", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="se", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval_nominal", chr_col="chrom", pos_col="pos")
    
    for (rsq in c(0.001)){
    dat <- harmonise_data(exposure_dat=chrom_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
    dat <- dat[dat$mr_keep==T,]
    
     if (sumstats_info$chromosome_name[sumstats_info$external_gene_name==i]=="X"){                                                                    # This little if-else-statement just makes sure that you get the appropriate RSIDs for each variant; because the X-chromosome requires an additional file, this one is in a separate loop
       dat <- merge(dat, ref_rsid[,c("V1", "V2", "V3")], by.x="pos.exposure", by.y="V2", all.x=T) 
       colnames(dat)[colnames(dat) %in% c("V1", "V3")] <- c("#chrom", "rsids")
     } else {
      dat <- merge(dat, outcome_rsid, by.x="pos.exposure", by.y="pos", all.x=T)
     }
     
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     rm(clump)
     
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info$external_gene_name[sumstats_info$external_gene_name==i]))
       results <- NULL
      } else {
    
     if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=i, outc=paste(outc), tissue=gtex, nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval, pval_thresh=pvalue, rsq_thresh=rsq, window=window, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
     } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), tissue=gtex, nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, window=window, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
      } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), tissue=gtex, nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, window=window, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       }
      }
     
  if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ",sumstats_info$external_gene_name[sumstats_info$external_gene_name==i]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)

  }
  }
  }
  }
  }
  }
  }
  }
  }
  }
   
    write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/5a_downstr_transcr_longcovid_gtex.csv") 
    write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/5b_downstr_transcr_longcovid_gtex_instr.csv") 
    print(paste(sumstats_info$external_gene_name[sumstats_info$external_gene_name==i], "done"))
 }
  
  
  
   # 9.c - Platelet eQTLs ####

  library(data.table)
  library(tidyr)
  library(dplyr)
  library(coloc)
  library(TwoSampleMR)
  library(ieugwasr)

  sumstats_info <- fread("/medpop/esp2/aschuerm/coagulation_longcovid/input/map_gtex.csv")

  longcovid <- fread("/medpop/esp2/aschuerm/resources/other_sumstats/long_covid_lammi_2023_strictcases_broadcontrols.gz") 
    longcovid$"#chrom" <- ifelse(longcovid$"#chrom"==23, "X", longcovid$"#chrom")
    colnames(longcovid)[which(colnames(longcovid)=="rsid")] <- "rsids"
    colnames(longcovid)[which(colnames(longcovid)=="stderr_beta")] <- "sebeta"
    colnames(longcovid)[which(colnames(longcovid)=="alt_allele_freq")] <- "af_alt"
    longcovid$pval <- 10^(-longcovid$neg_log_pvalue)
  
  ref_rsid <- fread("/medpop/esp2/aschuerm/tools/hg38_common_chrpos_X.txt")                                                   # This is a linker file to match RSIDs with chromosome positions for the X chromosome
  df_sum <- data.frame(exp=NA)[-1,]
  df_instr <- data.frame(pos.exposure=NA)[-1,]
   
  gtex_cisregion_u <- fread("/medpop/esp2/aschuerm/coagulation_longcovid/input/platelet_eqtls_coagcompl_keramati.csv")
  

 for (i in "F2R"){
  gtex_cisregion <- gtex_cisregion_u
  gtex_cisregion <- merge(gtex_cisregion, sumstats_info[,c("ensembl_gene_id", "external_gene_name")], by.x="gene", by.y="external_gene_name", all.x=T, all.y=F)
  colnames(gtex_cisregion)[colnames(gtex_cisregion)=="ensembl_gene_id"] <- "gene_id"
  
      for (window in c(0, 50000, 100000, 200000)) {
      gtex_cisregion_2 <- gtex_cisregion[gtex_cisregion$position > sumstats_info$start_position[sumstats_info$external_gene_name==i]-window & gtex_cisregion$position < sumstats_info$end_position[sumstats_info$external_gene_name==i]+window, ]
 
    if (nrow(gtex_cisregion_2)==0) {
        print("Skipping")
      } else {
        
   for (pvalue in c(5e-4, 5e-6, 5e-8)) {
      gtex_cisregion_3 <- data.frame(beta=gtex_cisregion_2$slope, se=gtex_cisregion_2$slope_se, pval=gtex_cisregion_2$pval_nominal, 
                                   trait="platelets", pos=gtex_cisregion_2$position, id=paste(gtex_cisregion_2$chromosome, gtex_cisregion_2$position, 
                                   gtex_cisregion_2$alt, gtex_cisregion_2$ref, sep=":"), chrom=gtex_cisregion_2$chromosome, ref=gtex_cisregion_2$ref, alt=gtex_cisregion_2$alt)
      gtex_cisregion_3 <- gtex_cisregion_3[!is.na(gtex_cisregion_3$se),]
      gtex_cisregion_3 <- gtex_cisregion_3[gtex_cisregion_3$pval < pvalue,]                                                                                                # Selecting "region-wide" significant cis-pQTs

     if (nrow(gtex_cisregion_3)==0) {
        print("Skipping")
      } else {
      
  for (outc in c("longcovid")){
    outcome <- get(outc)
   
    outcome_overlap <- outcome[outcome$`#chrom` == sumstats_info$chromosome_name[sumstats_info$external_gene_name==i],]                                      # Only selecting the chromosome of interest to speed stuff up downstream from here
    outcome_overlap <- outcome_overlap[outcome_overlap$pos %in% gtex_cisregion_3$pos,]                                                 # Only selecting the variants that are overlapping
      outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
      outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
      outcome_rsid <- outcome_rsid[!duplicated(outcome_rsid$pos),]
    outcome_overlap$phen <- paste(outc)
    outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")

    gtex_cisregion_overlap <- gtex_cisregion_3[gtex_cisregion_3$pos %in% outcome_overlap$pos,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     
    if (nrow(gtex_cisregion_overlap)==0) {
        print("Skipping")
      
      } else {
    gtex_cisregion_overlap_2 <- gtex_cisregion_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
    gtex_cisregion_overlap_2$beta <- gtex_cisregion_overlap_2$beta * -1
    colnames(gtex_cisregion_overlap_2)[colnames(gtex_cisregion_overlap_2) %in% c("ref", "alt")] <- c("alt", "ref")
    gtex_cisregion_overlap <- rbind(gtex_cisregion_overlap, gtex_cisregion_overlap_2)
    gtex_cisregion_overlap$id <- paste(gtex_cisregion_overlap$chrom, gtex_cisregion_overlap$pos, gtex_cisregion_overlap$alt, gtex_cisregion_overlap$ref, sep=":")
    gtex_cisregion_overlap$phen <- paste0(sumstats_info$external_gene_name[sumstats_info$external_gene_name==i], "_", "platelets")
    gtex_cisregion_overlap <- gtex_cisregion_overlap[gtex_cisregion_overlap$id %in% outcome_overlap$id,]
    outcome_overlap <- outcome_overlap[outcome_overlap$id %in% gtex_cisregion_overlap$id,]
    
    gtex_cisregion_overlap <- merge(gtex_cisregion_overlap, outcome_overlap[,c("id", "af_alt")], by="id")
    
    rm(gtex_cisregion_overlap_2)
    
     if (sumstats_info$chromosome_name[sumstats_info$external_gene_name==i]=="X"){                                                                    # This little if-else-statement just makes sure that you get the appropriate RSIDs for each variant; because the X-chromosome requires an additional file, this one is in a separate loop
       outcome_overlap <- outcome_overlap[,-(c("#chrom", "rsids"))]
       outcome_overlap <- merge(outcome_overlap, ref_rsid[,c("V1", "V2", "V3")], by.x="pos", by.y="V2", all.x=T) 
       colnames(outcome_overlap)[colnames(outcome_overlap) %in% c("V1", "V3")] <- c("#chrom", "rsids")
       gtex_cisregion_overlap <- merge(gtex_cisregion_overlap, ref_rsid[,c("V1", "V2", "V3")], by.x="pos", by.y="V2", all.x=T) 
       colnames(gtex_cisregion_overlap)[colnames(gtex_cisregion_overlap) %in% c("V1", "V3")] <- c("#chrom", "rsids")
     } else {
       outcome_overlap <- outcome_overlap[,-(c("#chrom", "rsids"))]
       outcome_overlap <- merge(outcome_overlap, outcome_rsid, by.x="pos", by.y="pos", all.x=T)
       gtex_cisregion_overlap <- merge(gtex_cisregion_overlap, outcome_rsid, by.x="pos", by.y="pos", all.x=T)
     }

    outcome_overlap <- outcome_overlap[!duplicated(outcome_overlap$rsids),]
    outcome_overlap <- outcome_overlap[!is.na(outcome_overlap$rsids),]
    gtex_cisregion_overlap <- gtex_cisregion_overlap[!duplicated(gtex_cisregion_overlap$rsids),]
    gtex_cisregion_overlap <- gtex_cisregion_overlap[!is.na(gtex_cisregion_overlap$rsids),]
    
    outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
    chrom_overlap <- format_data(gtex_cisregion_overlap, type="exposure", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="se", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval_nominal", chr_col="chrom", pos_col="pos")
    
    for (rsq in c(0.001)){
    dat <- harmonise_data(exposure_dat=chrom_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
    dat <- dat[dat$mr_keep==T,]
    
     if (sumstats_info$chromosome_name[sumstats_info$external_gene_name==i]=="X"){                                                                    # This little if-else-statement just makes sure that you get the appropriate RSIDs for each variant; because the X-chromosome requires an additional file, this one is in a separate loop
       dat <- merge(dat, ref_rsid[,c("V1", "V2", "V3")], by.x="pos.exposure", by.y="V2", all.x=T) 
       colnames(dat)[colnames(dat) %in% c("V1", "V3")] <- c("#chrom", "rsids")
     } else {
      dat <- merge(dat, outcome_rsid, by.x="pos.exposure", by.y="pos", all.x=T)
     }
     
     colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[!duplicated(dat$pos_id),]
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     rm(clump)
     
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info$external_gene_name[sumstats_info$external_gene_name==i]))
       results <- NULL
      } else {
    
      if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=paste(i), outc=paste(outc), tissue="platelets", nsnp=results_mr$nsnp, method=results_mr$method, 
                             b=results_mr$b, pval_thresh=pvalue, rsq_thresh=rsq, window=window, 
                           se=results_mr$se, pval=results_mr$pval, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
      } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), tissue="platelets", nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, window=window, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
       dat1 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome)
       output_mr_ivw <- MendelianRandomization::mr_ivw(dat1)
       results <- data.frame(outc=paste(outc), exp=paste(i), tissue="platelets", nsnp=output_mr_ivw@SNPs, 
                             method="Inverse variance weighted", pval_thresh=pvalue, rsq_thresh=rsq, window=window, b=output_mr_ivw@Estimate, 
                             se=output_mr_ivw@StdError, pval=output_mr_ivw@Pvalue, fstat_mean=mean(dat$beta.exposure^2/dat$se.exposure^2))
       
      }
      }
     
  if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ",sumstats_info$external_gene_name[sumstats_info$external_gene_name==i]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)

  }
  }
  }
  }
  }
  }
  }
  }

    write.csv(df_sum, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/5a_downstr_transcr_longcovid_platelet_qtls.csv") 
    write.csv(df_instr, "/broad/hptmp/aschuerm/synapser/longcovid_thrombo_mr/5b_downstr_transcr_longcovid_platelet_qtls_instr.csv") 
    print(paste(sumstats_info$external_gene_name[sumstats_info$external_gene_name==i], "done"))
 }
  
 
 
