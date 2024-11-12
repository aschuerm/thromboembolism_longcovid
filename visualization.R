# 0 - Libraries ####

  library(data.table)
  library(ggplot2)
  library(ggtext)
  library(ggpubr)
  library(R.utils)


# 1 - Main analyses ####
  # 1.a - Main analysis, other genetic instruments, and long COVID ####

  other <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1a_main_otherphenos_longcovid.csv")
    other <- other[other$method %in% c("Inverse variance weighted", "Wald ratio"),]
    other <- other[,-c("fstat_mean")]

    other$or <- exp(other$b)
    other$or_lci <- exp(other$b + qnorm(0.025)*other$se)
    other$or_uci <- exp(other$b + qnorm(0.975)*other$se)
    other$exposure <- c("Venous\nthromboembolism", "Deep vein thrombosis", "Pulmonary embolism", "Atrial fibrillation", "Asthma", 
                        "Chronic kidney\ndisease", "Chronic lower\nrespiratory disease", "Coronary artery\ndisease", "Dementia", "Diabetes mellitus", "Hypertension", "Migraine", 
                        "Multiple sclerosis", "Osteoarthritis","Rheumatoid arthritis")
  other$exposure <- factor(other$exposure, levels=c("Venous\nthromboembolism", "Deep vein thrombosis", "Pulmonary embolism", "Atrial fibrillation", "Asthma", 
                        "Chronic kidney\ndisease", "Chronic lower\nrespiratory disease", "Coronary artery\ndisease", "Dementia", "Diabetes mellitus", "Hypertension", "Migraine", 
                        "Multiple sclerosis", "Osteoarthritis","Rheumatoid arthritis"))
  other <- other[order(other$exposure),]
  other$rank <- 1:15
  other$group <- c("main", rep("coag", 2), rep("other", 12))
  
  f <- ggplot(other, aes(x=rank, y=or, color=group)) +
          geom_hline(yintercept=1, color="grey90") +
          geom_vline(xintercept=1.5, color="grey90", linetype="dashed") +
          geom_vline(xintercept=3.5, color="grey90", linetype="dashed") +
          geom_hline(yintercept=other$or[other$exposure=="Venous\nthromboembolism"], color="#B3B3FF", linetype="solid", alpha=0.5) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = other$or_lci[other$exposure=="Venous\nthromboembolism"], ymax = other$or_uci[other$exposure=="Venous\nthromboembolism"],  
                   fill = "#B3B3FF", alpha=.1) +
          geom_hline(yintercept=other$or_lci[other$exposure=="Venous\nthromboembolism"], color="#B3B3FF", linetype="solid", alpha=0.2) +
          geom_hline(yintercept=other$or_uci[other$exposure=="Venous\nthromboembolism"], color="#B3B3FF", linetype="solid", alpha=0.2) +
          geom_errorbar(aes(ymax=or_uci, ymin=or_lci), width=0.1, size=0.8) +
          geom_point(size=2) +
          scale_y_continuous(name="OR (95% CI)", limits=c(0.69, 1.37)) +
          scale_x_continuous(name="", breaks=1:15, labels=as.character(other$exposure)) +
          ggtitle("Associations of different genetically predicted conditions with long COVID:\n")+
          scale_color_manual(values=c("#C6C6FF", "#6666FF", "grey70")) +
          theme_classic() +
          theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.position="",
                axis.text.x = element_text(angle = 45, hjust = 1), plot.title=element_text(size=10), axis.title=element_text(size=10))
  
  rm(other)
  

  # 1.b - Other SARS-CoV-2-related phenotypes ####
  
  main <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1a_main_thromb_longcovid.csv")
    main <- main[main$method %in% c("Inverse variance weighted", "Wald ratio"),]
    main <- main[main$rsq_thresh==0.001 & main$pval_thresh==5e-8,]
    main <- main[main$outc=="longcovid",]
  covid <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1a_main_otherphenoscovid_longcovid.csv")
    covid <- covid[covid$method=="Inverse variance weighted",]
    covid <- covid[,-c("fstat_mean")]
  covid <- rbind(main, covid)
  main$or <- exp(main$b)
    main$or_lci <- exp(main$b + qnorm(0.025)*main$se)
    main$or_uci <- exp(main$b + qnorm(0.975)*main$se)
  
  covid$outcome <- c("Long COVID", "Critical COVID-19", "           Hospitalized\nCOVID-19", "SARS-CoV-2\ninfection")
  covid$rank <- c(1, 4, 3, 2)
    covid$or <- exp(covid$b)
    covid$or_lci <- exp(covid$b + qnorm(0.025)*covid$se)
    covid$or_uci <- exp(covid$b + qnorm(0.975)*covid$se)
    covid$sign <- ifelse(covid$pval<0.05, "sign", "nonsign")
    covid$in_ci <- ifelse(covid$outc=="longcovid", "longcovid", ifelse(covid$or>main$or_lci, "in", "not_in"))
    covid <- covid[order(covid$rank),]
  
  g <- ggplot(covid, aes(x=rank, y=or, color=in_ci)) +
          geom_vline(xintercept=1.5, color="grey90", linetype="dashed") +
          geom_hline(yintercept=1, color="grey90") +
          geom_hline(yintercept=main$or, color="#B3B3FF", linetype="solid", alpha=0.5) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = main$or_lci, ymax = main$or_uci,  
                   fill = "#B3B3FF", alpha=.1) +
          geom_hline(yintercept=main$or_lci, color="#B3B3FF", linetype="solid", alpha=0.2) +
          geom_hline(yintercept=main$or_uci, color="#B3B3FF", linetype="solid", alpha=0.2) +
          geom_errorbar(aes(ymax=or_uci, ymin=or_lci), width=0.1, size=0.8) +
          geom_point(size=2) +
          scale_color_manual(values=c("#C6C6FF", "#6666FF", "grey70")) +
          scale_y_continuous(name="OR (95% CI)", limits=c(0.75, 1.36)) +
          scale_x_continuous(name="", breaks=1:4, labels=as.character(covid$outcome), limits=c(0.5, 4.5)) +
          ggtitle("Associations of genetically predicted venous\nthromboembolism with acute COVID-19:")+
          theme_classic() +
          theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.position="",
                axis.text.x = element_text(angle = 45, hjust = 1), plot.title=element_text(size=10), axis.title=element_text(size=10))
    rm(covid, main)
  
  
  # 1.c - Other long COVID-resembling phenotypes ####
  
  main <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1a_main_thromb_longcovid.csv")
    main <- main[main$method %in% c("Inverse variance weighted", "Wald ratio"),]
    main <- main[main$rsq_thresh==0.001 & main$pval_thresh==5e-8,]
    main <- main[main$outc=="longcovid",]
    main$or <- exp(main$b)
    main$or_lci <- exp(main$b + qnorm(0.025)*main$se)
    main$or_uci <- exp(main$b + qnorm(0.975)*main$se)
  longcovid_resembling <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1a_main_thromb_longcovid_covidresemblingconditions.csv")
    longcovid_resembling <- longcovid_resembling[longcovid_resembling$method=="Inverse variance weighted",]
    longcovid_resembling <- longcovid_resembling[!(longcovid_resembling$outc %in% c("vte_2", "dvt", "pulmem", "bleeding", "persist_mood")),]

  longcovid_resembling$outcome <- c("Chronic lower\nrespiratory disease","Functional dyspepsia", "Irritable bowel\nsyndrome", "Endometriosis", "Fibromyalgia", "Anxiety disorders", 
                                    "Memory loss", "Sleep disorders", "Mood disturbances", "Depression", "Postviral fatigue")
  longcovid_resembling$rank <- c(1:nrow(longcovid_resembling))
    longcovid_resembling$or <- exp(longcovid_resembling$b)
    longcovid_resembling$or_lci <- exp(longcovid_resembling$b + qnorm(0.025)*longcovid_resembling$se)
    longcovid_resembling$or_uci <- exp(longcovid_resembling$b + qnorm(0.975)*longcovid_resembling$se)
    longcovid_resembling$in_ci <- ifelse(longcovid_resembling$or>main$or_lci, "in", "not_in")
    longcovid_resembling <- longcovid_resembling[order(longcovid_resembling$rank),]

  g_bis <- ggplot(longcovid_resembling, aes(x=rank, y=or, color=in_ci)) +
          geom_hline(yintercept=1, color="grey90") +
          geom_hline(yintercept=main$or, color="#B3B3FF", linetype="solid", alpha=0.5) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = main$or_lci, ymax = main$or_uci,  
                   fill = "#B3B3FF", alpha=.1) +
          geom_hline(yintercept=main$or_lci, color="#B3B3FF", linetype="solid", alpha=0.2) +
          geom_hline(yintercept=main$or_uci, color="#B3B3FF", linetype="solid", alpha=0.2) +
          geom_errorbar(aes(ymax=or_uci, ymin=or_lci), width=0.1, size=0.8) +
          geom_point(size=2) +
          scale_color_manual(values=c("grey70")) +
          scale_y_continuous(name="OR (95% CI)", limits=c(0.75, 1.36)) +
          scale_x_continuous(name="", labels=as.character(longcovid_resembling$outcome), breaks=1:nrow(longcovid_resembling)) +
          ggtitle("Associations of genetically predicted venous thromboembolism with long COVID-resembling conditions:\n")+
          theme_classic() +
          theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.position="",
                axis.text.x = element_text(angle = 45, hjust = 1), plot.title=element_text(size=10), axis.title=element_text(size=10))
    rm(covid, main)
  
  
# 2 - Sensitivity and replication analyses ####
  # 1.a - Different replication panels ####
  
  sens_mvp <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/2a_repl_thromb_longcovid_mvp.csv")
    sens_mvp$group <- "3_mvp"
  sens_phosp <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/2a_repl_thromb_longcovid_phosp.csv")
    sens_phosp$group <- "1_phosp"
  sens_thresh <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1a_main_thromb_longcovid.csv")
    sens_thresh$group <- "4_thresh"
    sens_thresh <- sens_thresh[sens_thresh$outc=="longcovid",]
  sens_strict <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1a_main_thromb_longcovid.csv")
    sens_strict$group <- "2_strict"
    sens_strict <- sens_strict[sens_strict$outc=="longcovid_strict",]
  
  sens_all <- rbindlist(list(sens_mvp, sens_phosp, sens_thresh, sens_strict))
    sens_all$analysis_group <- paste(sens_all$group, sens_all$pval_thresh, sens_all$rsq_thresh, sep=" / ")
    sens <- sens_all[sens_all$method %in% c("Inverse variance weighted", "Wald ratio"),]
    sens$or <- exp(sens$b)
    sens$or_lci <- exp(sens$b + qnorm(0.025)*sens$se)
    sens$or_uci <- exp(sens$b + qnorm(0.975)*sens$se)
    sens$analysis <- ifelse(sens$pval_thresh==5e-08 & sens$rsq_thresh==0.001, "main", "sens")
    sens$rsq_thresh <- formatC(sens$rsq_thresh, format = "e", digits = 0)

  sens_exp_1 <- sens[sens$group %in% c("3_mvp"),]
  sens_exp_2 <- sens[sens$group %in% c("4_thresh"),]
  sens_exp_3 <- sens[sens$group %in% c("2_strict"),]
  sens_exp_4 <- sens[sens$group %in% c("1_phosp"),]
  rm(sens_mvp, sens_strict, sens_thresh, sens_phosp, sens)

  e <- ggplot(sens_exp_1, aes(x=as.factor(pval_thresh), y=or, color=as.factor(rsq_thresh), alpha=analysis)) +
          geom_hline(yintercept=1, color="grey90") +
          geom_errorbar(aes(ymax=or_uci, ymin=or_lci), na.rm=TRUE, width=0.3, size=0.8, position=position_dodge(width = 0.7)) +
          geom_point(na.rm=TRUE, position=position_dodge(width = 0.7)) +
          scale_y_continuous(name="OR (95% CI) for long COVID", breaks=c(0.9, 1, 1.1, 1.2), limits=c(0.88, 1.25)) +
          scale_x_discrete(name=" ") +
          scale_color_manual(values=c( "#CC80CC", "pink2", "#9999FF","#ADD8E6"), name="*R*<sup>2</sup> threshold:") +
          scale_shape_manual(values=c(15,16,17,18), name=" ") +
          scale_alpha_manual(values=c(1, 0.3), name="Analysis:", labels=c("Main", "Sensitivity")) +
          theme_classic() +
          coord_flip() +
          ggtitle("Replication using genetic instruments from the MVP:") +
          theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.title=element_markdown(), plot.title=element_markdown(size = rel(0.8)))
  e_bis <- ggplot(sens_exp_2, aes(x=as.factor(pval_thresh), y=or, color=as.factor(rsq_thresh), alpha=analysis)) +
          geom_hline(yintercept=1, color="grey90") +
          geom_errorbar(aes(ymax=or_uci, ymin=or_lci), na.rm=TRUE, width=0.3, size=0.8, position=position_dodge(width = 0.7)) +
          geom_point(na.rm=TRUE, position=position_dodge(width = 0.7)) +
          scale_y_continuous(name="OR (95% CI) for long COVID", limits=c(0.85, 1.38), breaks=c(0.9, 1, 1.1, 1.2, 1.3)) +
          scale_x_discrete(name="*P*-value threshold") +
          scale_color_manual(values=c( "#CC80CC", "pink2", "#9999FF","#ADD8E6"), name="*R*<sup>2</sup> threshold:") +
          scale_shape_manual(values=c(15,16,17,18), name="*P*-value threshold:") +
          scale_alpha_manual(values=c(1, 0.3), name="Analysis:", labels=c("Main", "Sensitivity")) +
          theme_classic() +
          coord_flip() +
          ggtitle("Primary analyses: ") +
          theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.title=element_markdown(), plot.title=element_markdown(size = rel(0.8)))
  e_tris <- ggplot(sens_exp_3, aes(x=as.factor(pval_thresh), y=or, color=as.factor(rsq_thresh), alpha=analysis)) +
          geom_hline(yintercept=1, color="grey90") +
          geom_errorbar(aes(ymax=or_uci, ymin=or_lci), na.rm=TRUE, width=0.3, size=0.8, position=position_dodge(width = 0.7)) +
          geom_point(na.rm=TRUE, position=position_dodge(width = 0.7)) +
          scale_y_continuous(name="OR (95% CI) for long COVID", limits=c(0.9,1.29)) +
          scale_x_discrete(name=" ") +
          scale_color_manual(values=c( "#CC80CC", "pink2", "#9999FF","#ADD8E6"), name="*R*<sup>2</sup> threshold:") +
          scale_shape_manual(values=c(15,16,17,18), name=" ") +
          scale_alpha_manual(values=c(1, 0.3), name="Analysis:", labels=c("Main", "Sensitivity")) +
          theme_classic() +
          coord_flip() +
          ggtitle("Replication using stricter long COVID criteria:") +
          theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.title=element_markdown(), plot.title=element_markdown(size = rel(0.8)))
  e_quadr <- ggplot(sens_exp_4, aes(x=as.factor(pval_thresh), y=or, color=as.factor(rsq_thresh), alpha=analysis)) +
          geom_hline(yintercept=1, color="grey90") +
          geom_errorbar(aes(ymax=or_uci, ymin=or_lci), na.rm=TRUE, width=0.3, size=0.8, position=position_dodge(width = 0.7)) +
          geom_point(na.rm=TRUE, position=position_dodge(width = 0.7)) +
          scale_y_continuous(name="OR (95% CI) for long COVID", limits=c(0.58, 2.25), breaks=c(0.6, 1, 1.4, 1.8, 2.2)) +
          scale_x_discrete(name=" ") +
          scale_color_manual(values=c( "#CC80CC", "pink2", "#9999FF","#ADD8E6"), name="*R*<sup>2</sup> threshold:") +
          scale_shape_manual(values=c(15,16,17,18), name=" ") +
          scale_alpha_manual(values=c(1, 0.3), name="Analysis:", labels=c("Main", "Sensitivity")) +
          theme_classic() +
          coord_flip() +
          ggtitle("Replication in PHOSP-COVID:") +
          theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.title=element_markdown(), plot.title=element_markdown(size = rel(0.8)))
  rm(sens_all, sens_exp_1, sens_exp_2, sens_exp_3, sens_exp_4)
  
  
  # 1.b - Scatter plot ####
  
  regr <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1a_main_otherphenos_longcovid.csv")
    regr <- regr[regr$exp=="vte", ]
  
  instruments <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1b_main_otherphenos_longcovid_instr.csv")
    instruments <- instruments[instruments$exp=="vte", ]
    instruments$b_exp_use <- ifelse(instruments$beta.exposure<0, -instruments$beta.exposure, instruments$beta.exposure)
    instruments$b_outc_use <- ifelse(instruments$beta.exposure<0, -instruments$beta.outcome, instruments$beta.outcome)
    instruments$b_exp_min <- instruments$b_exp_use + qnorm(0.025)*instruments$se.exposure
    instruments$b_exp_max <- instruments$b_exp_use + qnorm(0.975)*instruments$se.exposure
    instruments$b_outc_min <- instruments$b_outc_use + qnorm(0.025)*instruments$se.outcome
    instruments$b_outc_max <- instruments$b_outc_use + qnorm(0.975)*instruments$se.outcome

  d <- ggplot(data=instruments, aes(x=b_exp_use, y=b_outc_use)) +
           geom_hline(yintercept=0, color="grey90") +
           geom_vline(xintercept=0, color="grey90") +
           geom_errorbar(aes(xmin=b_exp_min, xmax=b_exp_max), color="#C4CCD8", alpha=0.2) +
           geom_errorbar(aes(ymin=b_outc_min, ymax=b_outc_max), color="#C4CCD8", alpha=0.2) +
          geom_point(color="#6666FF", alpha=0.3) +
          geom_abline(aes(intercept=0, slope=regr$b[regr$method=="Median"], color="Median"), size=1) +
          geom_abline(aes(intercept=0, slope=regr$b[regr$method=="Mode-based"], color="Mode-based"), size=1) +
          geom_abline(aes(intercept=0, slope=regr$b[regr$method=="Inverse variance weighted"], color="IVW"), size=1) +
          geom_abline(aes(intercept=regr$b[regr$method=="Egger intercept"], slope=regr$b[regr$method=="Egger"], color="Egger"), size=1) +
          labs(x = "*\u03B2* for venous thromboembolism", y = "*\u03B2* for long COVID") + 
          scale_color_manual(name="Method:", breaks=c('IVW', 'Egger', 'Median', 'Mode-based'), values=c('IVW'='#6666FF', 'Egger'='#B3B3FF', 'Median'="grey60", 'Mode-based'="grey40"))+
          theme_classic()+
          theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.position="top", 
                legend.key = element_rect(colour = NA, fill = NA)) + 
          coord_cartesian(ylim=c(-0.34,0.34), xlim=c(-0.09,0.9))
  rm(regr, instruments)
  
  
  
# 3 - Identification of molecular drivers ####
  # 3.a - Partitioned genetic risk score for VTE ####
  
  regions <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/3a_downstr_thromb_longcovid_indivgenes.csv")
    regions <- regions[regions$method %in% c("Inverse variance weighted", "Wald ratio"),]
    regions_1 <- regions[regions$pval_thresh==5e-8,]
    regions_2 <- regions[regions$pval_thresh==5e-6 & (!(regions$gene %in% regions_1$gene)),]
    regions_3 <- regions[regions$pval_thresh==5e-4 & (!(regions$gene %in% c(regions_1$gene, regions_2$gene))),]
    regions <- rbind(regions_1, regions_2, regions_3)
    rm(regions_1, regions_2, regions_3)
    
  main <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1a_main_otherphenos_longcovid.csv")
    main <- main[main$exp=="vte" & main$method=="Inverse variance weighted",]
    main$gene <- NA
  regions <- rbind(regions, main)
  regions <- regions[regions$method %in% c("Wald ratio", "Inverse variance weighted (correlation inc)", "Inverse variance weighted"),]
  regions$exp <- ifelse(regions$exp=="C1QA", "C1q", ifelse(regions$exp=="C1R", "C1r", ifelse(regions$exp=="C1S", "C1s", 
                    ifelse(regions$exp=="C4BPB", "C4BP-\u03B2", ifelse(regions$exp=="C8B", "C8-\u03B2", ifelse(regions$exp=="CD46", "MCP / CD46", 
                    ifelse(regions$exp=="CD55", "DAF / CD55", ifelse(regions$exp=="CFB", "FB", ifelse(regions$exp=="CFD", "FD", 
                    ifelse(regions$exp=="CFH", "FH", ifelse(regions$exp=="CFI", "FI", ifelse(regions$exp=="CR1", "CR1 / CD35", 
                    ifelse(regions$exp=="CR2", "CR2 / CD21", ifelse(regions$exp=="F10", "FX", ifelse(regions$exp=="F11", "FXI", 
                    ifelse(regions$exp=="F12", "FXII", ifelse(regions$exp=="F13B", "FXIII-\u03B2", ifelse(regions$exp=="F2", "FII", 
                    ifelse(regions$exp=="F3", "TF / FIII", ifelse(regions$exp=="F7", "FVII", ifelse(regions$exp=="F2R", "PAR-1", 
                    ifelse(regions$exp=="F9", "FIX", ifelse(regions$exp=="FGA", "Fg-\u03B1", ifelse(regions$exp=="KLKB1", "PKa", 
                    ifelse(regions$exp=="MASP1", "MASP-1", ifelse(regions$exp=="MBL2", "MBL", ifelse(regions$exp=="PLAT", "tPA", 
                    ifelse(regions$exp=="PLAU", "uPA", ifelse(regions$exp=="PLAUR", "uPAR", ifelse(regions$exp=="PLG", "Plg", 
                    ifelse(regions$exp=="PROC", "PC", ifelse(regions$exp=="PROS1", "PS", ifelse(regions$exp=="SERPINA1", "\u03B11AT", 
                    ifelse(regions$exp=="SERPINA5", "PCI", ifelse(regions$exp=="SERPINC1", "AT", ifelse(regions$exp=="SERPIND1", "HCII", 
                    ifelse(regions$exp=="SERPINE1", "PAI-1", ifelse(regions$exp=="SERPINF2", "\u03B12AP", ifelse(regions$exp=="SERPING1", "C1-INH", 
                    ifelse(regions$exp=="THBD", "TM", regions$exp))))))))))))))))))))))))))))))))))))))))
  regions$b_lci <- (regions$b + qnorm(0.025)*regions$se)
  regions$b_uci <- (regions$b + qnorm(0.975)*regions$se)
  regions$sign <- ifelse(is.na(regions$gene), "prim", ifelse(regions$pval < 0.05, "sign", "nonsign"))
  regions$gene <- ifelse(is.na(regions$gene), "Genome-wide", regions$gene)
  regions$gene_form <- ifelse(regions$gene=="Genome-wide", "Genome-wide", paste0("*", regions$gene, "*"))
  regions <- regions[order(regions$b),]
  regions <- rbind(regions[regions$gene!="Genome-wide", ], regions[regions$gene=="Genome-wide", ])
  regions$rank <- 1:nrow(regions)
  write.csv(regions, "C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/3a_downstr_thromb_longcovid_indivgenes_modified.csv", row.names=F)
  
  h_tris <- ggplot(data=regions, aes(y=b, ymax=b_lci, ymin=b_uci, x=rank, color=sign, alpha=sign)) + 
    geom_hline(yintercept=0, color="grey90") +
    geom_hline(yintercept=regions$b[regions$gene=="Genome-wide"], color="#B3B3FF", linetype="solid", alpha=0.5) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = regions$b_lci[regions$gene=="Genome-wide"], ymax = regions$b_uci[regions$gene=="Genome-wide"],  
             fill = "#B3B3FF", alpha=.1) +
    geom_hline(yintercept=regions$b_lci[regions$gene=="Genome-wide"], color="#B3B3FF", linetype="solid", alpha=0.2) +
    geom_hline(yintercept=regions$b_uci[regions$gene=="Genome-wide"], color="#B3B3FF", linetype="solid", alpha=0.2) +
    geom_vline(xintercept=nrow(regions)-0.5, color="grey90", linetype="dashed") +
    geom_point() +
    geom_errorbar(width=0.3) +
    scale_y_continuous(name="Association with long<br>COVID (*\u03B2* [95% CI])") +
    scale_x_continuous(breaks=1:nrow(regions), labels=regions$gene_form, name="") +
    scale_color_manual(values=c("#FFA9D4", "#6666FF", "#FF69B4")) +
    scale_alpha_manual(values=c(0.7, 1, 1)) +
    theme_classic() + 
    theme(axis.text.x = element_markdown(angle = 45, hjust = 1), axis.title.y = element_markdown(), 
          legend.position = "none", legend = NULL)
  
  h <- ggplot(data=regions[regions$gene_form!="Genome-wide",], aes(y=(b), x=rank, color=b, size=-log10(pval))) + 
    geom_hline(yintercept=0, color="grey90") +
    geom_point() +
    scale_y_continuous(name="Association with<br>long COVID (*\u03B2*)", breaks=c(-1, -0.5, 0, 0.5, 1, 1.5), limits=c(-1.2, 1.7)) +
    scale_x_continuous(breaks=1:nrow(regions[regions$gene_form!="Genome-wide",]), labels=regions[regions$gene_form!="Genome-wide",]$gene_form, name="") +
    scale_colour_gradientn(colors=c("#3A3A98", "#7C7CBA", "#BDBDDD", "#FFD4EA", "#FFA9D4", "#D67D99", scales::muted("red")),values=c(0,0.12,0.25,0.37,0.5,0.67,0.8, 1), name="*\u03B2*:", breaks=c(-1, 0, 1)) +
    scale_size_continuous(breaks = c(0.5, 1, 1.5), limits=c(0.004, 2), name="-Log<sub>10</sub>(*P*):") + 
    ggtitle("Associations of gene-specific genetic instruments for venous thromboembolism with long COVID:") +
    theme_classic() + 
    theme(axis.text.x = element_markdown(angle = 45, hjust = 1), axis.title.y = element_markdown(),
          legend.title=element_markdown(), legend.position="right", legend.box="horizontal", plot.title=element_markdown(size = 10), axis.title=element_markdown(size = 10))
  
  rm(main, regions)
  
  
  # 3.b - Protein associations ####
  
  prot_ukb <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/4a_downstr_proteins_longcovid_ukbppp.csv")
  prot_ukb$window <- 200000
  prot_ukb_intr <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/4a_downstr_proteins_longcovid_ukbppp_intragenic.csv")
  prot_ukb_intr$window <- 0
  prot_ukb <- rbind(prot_ukb, prot_ukb_intr)
  rm(prot_ukb_intr)
  
  sens_all <- prot_ukb
    sens_all <- sens_all[sens_all$exp %in% c("F2R"),]
    sens_all$analysis_group <- paste(sens_all$group, sens_all$pval_thresh, sens_all$rsq_thresh, sep=" / ")
    sens <- sens_all[sens_all$method %in% c("Inverse variance weighted", "Wald ratio"),]
    sens$or <- exp(sens$b)
    sens$or_lci <- exp(sens$b + qnorm(0.025)*sens$se)
    sens$or_uci <- exp(sens$b + qnorm(0.975)*sens$se)
    sens$b_lci <- (sens$b + qnorm(0.025)*sens$se)
    sens$b_uci <- (sens$b + qnorm(0.975)*sens$se)
    sens$analysis <- ifelse(sens$pval_thresh==5e-08 & sens$rsq_thresh==0.001 & sens$window==200000, "main", "sens")
    sens$rsq_thresh <- formatC(sens$rsq_thresh, format = "e", digits = 0)
    sens$window_words <- ifelse(sens$window==0, "1_Intragenic variants", "0_Variants within\n200 kilobases")
    
  rm(prot_ukb, sens_all)

  i <- ggplot(sens, aes(x=as.factor(pval_thresh), y=b, color=as.factor(rsq_thresh), shape=window_words, alpha=analysis)) +
          geom_hline(yintercept=0, color="grey90") +
          geom_errorbar(aes(ymax=b_uci, ymin=b_lci), na.rm=TRUE, width=0.3, size=0.8, position=position_dodge(width = 0.7)) +
          geom_point(na.rm=TRUE, position=position_dodge(width = 0.7), size=1.5, fill="white") +
          scale_y_continuous(name="*\u03B2* (95% CI)", limits=c(-1.15, 1.15)) +
          scale_x_discrete(name="*P*-value threshold") +
          scale_color_manual(values=c( "#CC80CC", "pink2", "#832424","#ADD8E6"), name="*R*<sup>2</sup> threshold:") +
          scale_shape_manual(values=c(21,16), name="Variants:", labels=c("Variants within\n200 kilobases", "Intragenic\nvariants")) +
          scale_alpha_manual(values=c(1, 0.3), name="Analysis:", labels=c("Main", "Sensitivity"), guide="none") +
          theme_classic() +
          coord_flip() +
          ggtitle("Associations of circulating PAR-1 with long COVID:") +
          theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.title=element_markdown(), plot.title=element_markdown(size = 10), axis.title=element_markdown(size = 10))
  rm(sens)
  
  
  # 3.c - Identification of molecular drivers ####

  tissues <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/5a_downstr_transcr_longcovid_gtex.csv")
    platelets <- fread("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/5a_downstr_transcr_longcovid_platelet_qtls.csv")
    tissues <- rbind(tissues, platelets)
    tissues <- tissues[tissues$window %in% c(0, 200000),]
    tissues$tissue_window <- paste0(tissues$tissue, "_", tissues$window)
    tissues_1 <- tissues[tissues$pval_thresh == 5e-8 & tissues$rsq_thresh == 1e-3,]
    tissues_2 <- tissues[tissues$pval_thresh == 5e-6 & tissues$rsq_thresh == 1e-3 & (!(tissues$tissue_window %in% tissues_1$tissue_window)),]
    tissues_3 <- tissues[tissues$pval_thresh == 5e-4 & tissues$rsq_thresh == 1e-3 & (!(tissues$tissue_window %in% tissues_1$tissue_window)) & (!(tissues$tissue_window %in% tissues_2$tissue_window)),]
    tissues <- rbind(tissues_1, tissues_2, tissues_3)
    rm(tissues_1, tissues_2, tissues_3, platelets)
    
    tissues$b_lci <- (tissues$b + qnorm(0.025)*tissues$se)
    tissues$b_uci <- (tissues$b + qnorm(0.975)*tissues$se)
    tissues$or <- exp(tissues$b)
    tissues$or_lci <- exp(tissues$b + qnorm(0.025)*tissues$se)
    tissues$or_uci <- exp(tissues$b + qnorm(0.975)*tissues$se)
    tissues$sign <- ifelse(tissues$pval < 0.05, "sign", "nonsign")
    tissues <- tissues[order(tissues$b),]
    tissues$rank <- c(1, 1, 2, 5, 2, 3, 3, 4, 5)
    tissues$window_words <- ifelse(tissues$window==0, "1_Intragenic variants", "0_Variants within\n200 kilobases")
    tissues$tissue_form <- DescTools::StrCap(tolower(gsub("_", " ", tissues$tissue)))
    tissues$tissue_form <- ifelse(tissues$tissue_form=="Artery tibial", "Artery\n(tibial)",
                                  ifelse(tissues$tissue_form=="Muscle skeletal", "Muscle\n(skeletal)",
                                  ifelse(tissues$tissue_form=="Esophagus muscularis", "Esophagus\n(muscularis)",
                                  ifelse(tissues$tissue_form=="Small intestine terminal ileum", "Small intestine\n(terminal ileum)",
                                  ifelse(tissues$tissue_form=="Adipose subcutaneous", "Adipose\n(subcutaneous)",
                                  ifelse(tissues$tissue_form=="Colon sigmoid", "Colon\n(sigmoid)",
                                  ifelse(tissues$tissue_form=="Adipose visceral omentum", "Adipose\n(visceral omentum)",
                                  ifelse(tissues$tissue_form=="Nerve tibial", "Nerve\n(tibial)",
                                  ifelse(tissues$tissue_form=="Esophagus mucosa", "Esophagus\n(mucosa)",
                                  ifelse(tissues$tissue_form=="Heart left ventricle", "Heart\n(left ventricle)",
                                  ifelse(tissues$tissue_form=="Skin not sun exposed suprapubic", "Skin\n(suprapubic)",
                                  ifelse(tissues$tissue_form=="Colon transverse", "Colon\n(transverse)",
                                  ifelse(tissues$tissue_form=="Skin sun exposed lower leg", "Skin\n(lower leg)",
                                  ifelse(tissues$tissue_form=="Breast mammary tissue", "Breast\n(mammary tissue)",
                                  ifelse(tissues$tissue_form=="Esophagus gastroesophageal junction", "Esophagus\n(GE junction)",
                                  ifelse(tissues$tissue_form=="Breast mammary tissue", "Breast\n(mammary tissue)",
                                  ifelse(tissues$tissue_form=="Artery aorta", "Artery\n(aorta)",
                                  ifelse(tissues$tissue_form=="Heart atrial appendage", "Heart\n(atrial appendage)", tissues$tissue_form))))))))))))))))))
    write.csv(tissues, "C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/5a_downstr_transcr_longcovid_gtex.csv", row.names=F)
    
    j_bis <- ggplot(data=tissues, aes(y=b, x=rank, color=window_words, shape=window_words, size=-log10(pval))) + 
      geom_hline(yintercept=0, color="grey90") +
      geom_point(position=position_dodge(0.5), stroke=1.15) +
      scale_y_continuous(name="Association with<br>long COVID (*\u03B2*)", limits=c(-0.8,0.8)) +
      scale_x_continuous(name="Tissues", breaks=1:5, labels=unique(tissues$tissue_form[order(tissues$rank)])) +
      scale_colour_manual(values=c("grey20", "grey20"), name="Variants:", labels=c("Variants within\n200 kilobases", "Intragenic\nvariants")) +
      scale_shape_manual(values=c(1,16), name="Variants:", labels=c("Variants within\n200 kilobases", "Intragenic\nvariants")) +
      scale_size_continuous(breaks = c(0.5, 1, 1.5), limits=c(0.01, 1.5), name="-Log<sub>10</sub>(*P*):") + 
      ggtitle("Associations of *F2R* expression with long COVID:<br>") +
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title=element_markdown(), axis.title.y = element_markdown(), 
            legend.position="right", plot.title=element_markdown(size = 10), axis.title=element_markdown(size = 10))

    
  
# 4 - PDFs ####
  
  pdf("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1_primary_analyses.pdf", width=10, height=3)
    f
    dev.off()
    ggsave("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/1_primary_analyses.svg", plot=f, width=10, height=3)
  
  
  pdf("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/2_different_outcomes.pdf", width=12, height=3)
    ggarrange(NULL, g, NULL, NULL, g_bis, widths=c(0.4, 7, 0.2, 0.4, 16), labels=c("a", "", "", "b", ""), nrow=1, ncol=5)
    dev.off()
    ggsave("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/2_different_outcomes.svg", 
           plot=ggarrange(NULL, g, NULL, NULL, g_bis, widths=c(0.4, 7, 0.2, 0.4, 16), labels=c("a", "", "", "b", ""), nrow=1, ncol=5), width=12, height=3)
  
  pdf("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/3_sensitivity_analyses.pdf", width=15, height=4)
    ggarrange(NULL, e_bis, NULL, e, NULL, e_tris, NULL, e_quadr, nrow=1, ncol=8, widths=rep(c(0.2,4), 4), labels=c("a", "", "b", "", "c", "", "d", ""), common.legend = T)
    dev.off()
    ggsave("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/3_sensitivity_analyses.svg", ggarrange(NULL, e_bis, NULL, e, NULL, e_tris, NULL, e_quadr, nrow=1, ncol=8, widths=rep(c(0.2,4), 4), labels=c("a", "", "b", "", "c", "", "d", ""), common.legend = T), width=15, height=4)

  pdf("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/4_molecular_analyses.pdf", width=10, height=7)
    ggarrange(NULL, ggarrange(NULL, h, widths=c(0.4, 8.9), nrow=1, ncol=2), NULL, ggarrange(NULL, NULL, NULL, NULL, NULL, i, NULL, j_bis, ncol=4, nrow=2, widths=c(0.4, 4.25 ,0.4, 4.25), heights=c(0.4, 5), common.legend = F, labels=c("b", "", "c", "", "", "", "", "")), 
            nrow=4, ncol=1, heights=c(0.4, 5, 0.2, 6), labels=c("a", "", "", ""))
    dev.off()
    ggsave("C:/Users/artsc/Documents/- Onderzoek/AS/14 - Coagulation in long COVID/output/4_molecular_analyses.svg", 
           plot=ggarrange(NULL, ggarrange(NULL, h, widths=c(0.4, 8.9), nrow=1, ncol=2), NULL, ggarrange(NULL, NULL, NULL, NULL, NULL, i, NULL, j_bis, ncol=4, nrow=2, widths=c(0.4, 4.25 ,0.4, 4.25), heights=c(0.4, 5), common.legend = F, labels=c("b", "", "c", "", "", "", "", "")), 
            nrow=4, ncol=1, heights=c(0.4, 5, 0.2, 6), labels=c("a", "", "", "")), width=10, height=7)
  
  
