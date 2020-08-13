####################################################################################
# environmental variables, strato or local (niekverw) 
f.mapping.proteincoding = "data/mapping.proteincoding.tsv"

if (Sys.getenv("USER")=="niek"){
  # print("loading strato variables")
  # f_static="/Users/niek/data/ecg/snpindexes.Rdata"
  # f_ecg_unadjusted="/Users/niek/data/ecg/all.SignalAverage.unadjusted.mean.QRS.tsv.gz.mean.tsv"
  # f_ecg_stretch="/Users/niek/data/ecg/all.SignalAverage.stretch.mean.QRS.tsv.gz.mean.tsv"
  # dir_data = "/Users/niek/data/ecg/" # contains "all_stats.unadjusted" and "all_stats.stretch"
  # dir_recombination = "/data_work/databases/genetic_recombination_rates/hapmap/"
  datadir="/Users/niek/data/ecg/data" #(300001-1) %/% 300000
  dir_recombination = "/Users/niek/data/ecg/data/recombination_rates/"
  
} else if (Sys.getenv("USER") == "benzhi") {
  print("loading strato variables")
  f_static="/mnt/linux_data/data/ecgenetics/snpindexes.Rdata"
  f_ecg_unadjusted="/mnt/linux_data/data/ecgenetics/all.SignalAverage.unadjusted.mean.QRS.tsv.gz.mean.tsv"
  f_ecg_stretch="/mnt/linux_data/data/ecgenetics/all.SignalAverage.stretch.mean.QRS.tsv.gz.mean.tsv"
  dir_data = "/mnt/linux_data/data/ecgenetics/"
  dir_recombination = "/mnt/linux_data/data/ecgenetics/recombination_rates/"
} else {
  # f_static="/data/datasets/ecg_morph_stats/snpindexes.Rdata"
  # f_ecg_unadjusted="/data/datasets/ecg_morph_stats/all.SignalAverage.unadjusted.mean.QRS.tsv.gz.mean.tsv"
  # f_ecg_stretch="/data/datasets/ecg_morph_stats/all.SignalAverage.stretch.mean.QRS.tsv.gz.mean.tsv"
  # dir_data = "/data/datasets/ecg_morph_stats/" # contains "all_stats.unadjusted" and "all_stats.stretch"
  # dir_recombination = "/data/datasets/recombination_rates/"
  datadir="/data/datasets/ecg_morph_stats/data"
}
################################################
# #https://github.com/Automattic/node-canvas/wiki/Installation:-Mac-OS-X
# #https://plot.ly/r/shinyapp-linked-click/
# install.packages("Cairo")
library(ggplot2)
library(Cairo)   # For nicer ggplot2 output when deployed on Linux
library(ggthemes)
library(plotly)
library(shiny)

library(data.table)
library(ggplot2)
library(matrixStats)
library(heatmaply)

library("bedr")
library(seqminer)
library(shinyjs)
library(MendelianRandomization)
### 
# https://www.r-bloggers.com/long-running-tasks-with-shiny-challenges-and-solutions/
# https://cran.r-project.org/web/packages/promises/vignettes/shiny.html
library(promises)
library(future)
plan(future::multiprocess)
#plan("multisession")
# [ONE-TIME WARNING] Forked processing ('multicore') is disabled in future (>= 1.13.0) 
# when running R from RStudio, because it is considered unstable. Because of this, plan("multicore") 
# will fall back to plan("sequential"), and plan("multiprocess") will fall back to plan("multisession") - 
# not plan("multicore") as in the past. For more details, how to control forked processing or not, 
# and how to silence this warning in future R sessions, see ?future::supportsMulticore 
# > # genes do not always show on regional plot when th
#   
#library(future.callr)
#plan(callr)
##################
#### DATAFRAMES
##################

### MAKE REVISED DATAFRAMES
##### RSID:
#df.static <- fread(paste0(dir_data,"/stretch.logP.outfile.index.tsv.gz.tophits.gz"))
# names(df.static) <- c("i","uniqid","SNP","CHR","BP","EFAL","NEFAL","EAF","INFO")
# df.static.rsid <- df.static[startsWith(df.static$SNP,"rs"),c("CHR", "SNP","i")]
# df.static.rsid$SNP <- as.numeric(sub("rs","",df.static.rsid$SNP))
# setkey(df.static.rsid,"SNP")
# ##### POSITIONAL
# df.static.pos <- df.static[,c("CHR","BP","i")]
# setkey(df.static.pos,CHR,BP)

# save(df.static.rsid,df.static.pos,file=f_static)
# df.static<-NULL

if (!exists("df.static.rsid")){load(paste0(datadir,"/various_data/snpindexes.Rdata"))}

# check.size(df.static)
# check.size(df.static.rsid)
# check.size(df.static.pos)
#############
if (!exists("mapping.proteincoding")){
  try(mapping.proteincoding <- data.frame(fread(f.mapping.proteincoding)))
}

# 
#df.static$chromosome <- paste0("chr",df.static$CHR)

if (!exists("df_ecg_unadjusted")){df_ecg_unadjusted<-fread(paste0(datadir,"/various_data/all.SignalAverage.unadjusted.mean.QRS.tsv.gz.mean.tsv"))}
if (!exists("df_ecg_stretch")){df_ecg_stretch<-fread(paste0(datadir,"/various_data/all.SignalAverage.stretch.mean.QRS.tsv.gz.mean.tsv"))}


#### MULTIPLE SNPS
rsids=c("rs540538026","rs367896724","rs201106462")
rsids=c("rs201106462")

