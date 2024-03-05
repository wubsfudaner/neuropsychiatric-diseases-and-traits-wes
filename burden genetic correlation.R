library(bhr)
library(tidyr)
library(data.table)
###
single_result_path="/bhr/input/rare/"
###
setwd(single_result_path)
filelist <- list.files(single_result_path)
filelist<-grep("_plof.txt",filelist,value = T)
pheno_list<-c(0)
for (i in 1:length(filelist)) {
  pheno_path<-filelist[i]
  pheno_name<-gsub("_plof.txt","",pheno_path)
  
  BHR_input_path<-paste("/bhr/input/rare/",pheno_name,"_plof.txt",sep = "")
  BHR_input_rare<-as.data.frame(fread(BHR_input_path))
  baseline <- read.table("ms_baseline_oe5.txt")
  assign(pheno_name,BHR_input_rare)
  pheno_list<-c(pheno_list,pheno_name)
}
pheno_list<-pheno_list[-1]

result_rg<-data.frame(pheno1_name=as.character(0),pheno2_name=as.character(0),rg=as.numeric(0),rg_se=as.numeric(0))
for (i in 1:ncol(combn(pheno_list, 2))){
  bivariate <- BHR(mode = "bivariate",
                   trait1_sumstats = get(combn(pheno_list, 2)[,i][1]),
                   trait2_sumstats = get(combn(pheno_list, 2)[,i][2]),
                   annotations = list(baseline))
  pheno1_name<-combn(pheno_list, 2)[,i][1]
  pheno2_name<-combn(pheno_list, 2)[,i][2]
  rg<-bivariate$rg$rg_mixed
  rg_se<-bivariate$rg$rg_mixed_se
  result<-cbind(pheno1_name,pheno2_name,rg,rg_se)
  colnames(result)<-c("pheno1_name","pheno2_name","rg","rg_se")
  result_rg<-rbind(result_rg,result)
}

write.table(result_rg,"/bhr/output/bhr_rg_plof.txt",sep = "\t",row.names = F,quote = F)