
suppressMessages(library("foreach"))
suppressMessages(library(doParallel))
suppressMessages(library(optparse))
source("../bcr/RIMA/RIMA_pipeline/src/immune_repertoire/trust4_metric_functions.R")

suppressMessages(library(tidyverse))
option_list = list(
make_option(c("-s", "--sample"), type="character", default=NULL,
            help="samplex info", metavar="character"))
samples <- read_csv("../samples.csv")
n.cores <- parallel::detectCores() - 1
registerDoParallel(cores=n.cores)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

outdir_main <- "./R/"
createResDir <- function(dir_name  = ""){
  if (!dir.exists(paste0( dir_name))){
    dir.create(paste0( dir_name))
  }else{
    print("dir exists") 
  }}
createResDir(outdir_main)
  s <- samples[samples$library == opt$sample,]
  outdir <- paste0(outdir_main, s$library)
  
  cdr3 <- read.table(file = paste0("../bcr_take2_trust4/output_reads/", s$library, "_report.tsv"), sep = "\t", header = T, comment.char = "$", stringsAsFactors = FALSE)
    cdr3 <- cdr3 %>% mutate(count = X.count)
    cdr3 <- cdr3 %>% select( count, everything()) %>% select(-X.count)
  phenotype <- paste0(s$condition, "_", s$cell_type)        
  cdr3 <- subset(cdr3, count > 0) %>% 
    mutate(V = as.character(V), J = as.character(J), C = as.character(C), CDR3aa = as.character(CDR3aa)) %>% 
    mutate(clinic = as.character(phenotype))
  cdr3$is_complete <- sapply(cdr3$CDR3aa, function(x) ifelse(x != "partial" && x != "out_of_frame" && !grepl("^_",x) && !grepl("^\\?", x),"Y","N")) 
  
  #exact the TCR and BCR
  cdr3.bcr <- subset(cdr3, grepl("^IG",V) | grepl("^IG",J) | grepl("^IG",C))

  #add lib size and clinic traits
  cdr3.bcr <- cdr3.bcr %>% mutate(lib.size = sum(count)) 

  #split BCR into heavy chain and light chain
  cdr3.bcr.heavy <- subset(cdr3.bcr, grepl("^IGH",V) | grepl("^IGH",J) | grepl("^IGH",C))
  cdr3.bcr.light <- subset(cdr3.bcr, grepl("^IG[K|L]",V) | grepl("^IG[K|L]",J) | grepl("^IG[K|L]",C))
  
  #save BCR and TCR info for downsteam use
  outdir <- paste0(outdir_main, s$library)
  sample <- s$library 
  cdr3.bcr.light <- cdr3.bcr.light %>% 
    mutate(lib.size = sum(count))      
  cdr3.bcr.heavy$sample <- sample
  cdr3.bcr.light$sample <- sample
  cdr3.bcr.heavy <- cdr3.bcr.heavy %>% 
    mutate(lib.size = sum(count))  
  #save(cdr3.bcr.light, file = paste(outdir, "_TRUST4_BCR_light.Rdata",sep = ""))
  #save(cdr3.bcr.heavy, file = paste(outdir, "_TRUST4_BCR_heavy.Rdata",sep = ""))
  write_csv(cdr3.bcr.light, file = paste(outdir, "_TRUST4_BCR_light.csv",sep = ""))
  write_csv(cdr3.bcr.heavy, file = paste(outdir, "_TRUST4_BCR_heavy.csv",sep = ""))
  #BCR clustering 
  #Note that not every sample have BCR cluster




