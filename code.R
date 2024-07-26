# immunce cell trait exposure download.
immid <- paste0('ebi-a-GCST9000',c(1391:2121))
# obtain ID to trait file
library(MRInstruments)
library(plyr)
library(dplyr)
library(TwoSampleMR)
idTotrait <- available_outcomes()
rownames(idTotrait) <- idTotrait$id
ID <- idTotrait[immid,]
ID_Trait <- list(immID <-immid,ID<- ID)
save(ID_Trait,file = "TraitID.RData")

# dowload 731 immune cell trait exposure data.
dir.create("/Users/llls2012163.com/GWAS/Immune cell trait")
setwd("/Users/llls2012163.com/GWAS/Immune cell trait")
immucell_ID <- as.vector(immid)
for(i in immucell_ID){
  expo_data <- TwoSampleMR::extract_instruments(outcome=i,
                                   p1= 1e-5,
                                   clump = T,
                                   p2= 1e-5,
                                   r2 = 0.1,
                                   kb = 500)
  filename = i
  write.table(expo_data,file=paste0(filename,'.txt'),row.names = F,sep = '\t',quote = F)
}
  
  
  

##download outcome data. immune cell trait to Lung cancer

# ukb-a-54 :Lung cancer
# ieu-a-984 :Lung adenocarcinoma
# ieu-a-988 :small cell lung carcinoma
# ieu-a-989 : squamous cell lung cancer
# finn-b-C3_LUNG_NONSMALL : no small lung cancer
# ukb-b-14521 : lung cancer father
# ukb-b-20176 : lung cancer mother
# ebi-a-GCST90013972 : lung cancer any parental
# ukb-b-15826: lung cancer sibilings
# ebi-a-GCST004749 : lung cancer smoker
# ebi-a-GCST004747: lung cancer no smoker


File_name <- list.files("/Users/llls2012163.com/GWAS/Immune cell trait")
library(MRInstruments)
library(plyr)
library(dplyr)
library(TwoSampleMR)
#immune cell trait to Lung Cancer: outcome download function.--------
Download.outcome <- function(x,y,z){# x: File_name y: GWAS_ID z: cancer type
for(i in x){
path_promt <- paste0("/Users/llls2012163.com/GWAS/Immune cell trait/",i)
data_promt <- read.delim(path_promt)
data_outcome <- extract_outcome_data(snps = data_promt$SNP,outcomes = y)
outcome_path <- paste0("/Users/llls2012163.com/Lung cancer/CancerTypeOutcome/",z,"/",i)
write.table(data_outcome,file=outcome_path,row.names = F,sep = '\t',quote = F)
}
}
## immune cell trait to Lung Cancer: data check function.
Datacheck <- function(x,y){# x: File_name, y: cancer type
path_file <- paste0('/Users/llls2012163.com/Lung cancer/CancerTypeOutcome/',y)
Files_downloaded <- list.files(path_file)
File_udownload <- setdiff(File_name,Files_downloaded)
File_size <- c()
for(i in Files_downloaded){
path_promt <- paste0(path_file,'/',i)
size_promt <- file.size(path_promt)
if(length(File_size) == 0){
File_size <- size_promt
}else{
File_size <- c(File_size,size_promt)
}
}
File_size <- cbind(Files_downloaded,File_size)
## The empty file for server error: file size = 1 byte which suggested empty file and was redownloaded
File_name_empty <- File_size[File_size[,2]==1,][,1]
File_name_empty <- c(File_udownload,File_name_empty)
return(File_name_empty)
}
### immune cell trait to Lung Cancer: data harmonise function.
Harmose_Fun <- function(x,y,z){# x: File_name # y: exposure directory # z: output directory
for(i in x){
path_exposure <- paste0("/Users/llls2012163.com/GWAS/Immune cell trait/",i)
data_exposure <- read.delim(path_exposure)
path_outcome  <- paste0('/Users/llls2012163.com/Lung cancer/',y,'/',i)
data_outcome <- read.delim(path_outcome)
data_promt <- harmonise_data(data_exposure,data_outcome)
path_output <- paste0("/Users/llls2012163.com/Lung cancer/Harmnoise_immToCancer",z,'/',i)
write.table(data_promt,file = path_output,row.names = F,sep = '\t',quote = F)
}
}
#### detect the file name which do no map the outcome information
disctinct_otherwise_empty_function <- function(x,y){# x:File_name y:aim directionary
  path_file <- paste0('/Users/llls2012163.com/Lung cancer/Harmnoise_immToCancer',y)
  File_size <- c()
  
  for(i in x){
    path_promt <- paste0(path_file,'/',i)
    size_promt <- file.size(path_promt)
    if(length(File_size) == 0){
      File_size <- size_promt
    }else{
      File_size <- c(File_size,size_promt)
    }
  }
  File_size <- cbind(File_name,File_size)
  empty <- File_size[File_size[,2] == 1,][1]
  return(empty)
}
##### remove the weak IVs
# This function calculate the F value and reamin SNP with F > 10 in the harmosed file.
calculate_F_lager_tha10_exist <- function(x){
  file_path <- paste0('/Users/llls2012163.com/Lung cancer/Harmnoise_immToCancer',x)
  file_name <- list.files(file_path)
  for(i in file_name){
    path_promt <- paste0(file_path,'/',i)
    data_promt <- read.delim(path_promt)
    R_squre <- TwoSampleMR::get_r_from_bsen(data_promt$beta.exposure,data_promt$se.exposure,data_promt$samplesize.exposure)
    F_value <- as.numeric(R_squre)*(as.numeric(data_promt$samplesize.exposure) - 2)/(1 - as.numeric(R_squre))
    data_promt <- cbind(data_promt, R2 = as.numeric(R_squre), F_value = F_value)
    data_finall <- data_promt[as.numeric(data_promt$F_value) > 10,]
    write.table(data_finall,file = path_promt,row.names = F,sep = '\t',quote = F)
  }
} 
###### MR analysis function
MR_analysis <- function(x){
  result_mr <- data.frame()
  file_path <- paste0('/Users/llls2012163.com/Lung cancer/Harmnoise_immToCancer',x)
  file_name <- list.files(file_path)
  for(j in file_name){
    path_promt <- paste0(file_path,'/',j)
    data_promt <- read.delim(path_promt)
    res <- mr(data_promt,method_list = c('mr_ivw_fe','mr_ivw_mre','mr_egger_regression','mr_weighted_median','mr_weighted_mode','mr_simple_mode'))
    if(dim(result_mr)[1]==0){
      result_mr=res
    }else{
      result_mr=rbind(result_mr,res)
    }
    outcome_path <- paste0('/Users/llls2012163.com/Lung cancer/result_file/MR_result',x,'.txt')
    write.table(result_mr,file = outcome_path,row.names = F,sep = '\t',quote = F)
  }
}
####### MR heterogeneity analysis
MR_heterogeneity_analysis <- function(x){
  mr_het <- data.frame()
  file_path <- paste0('/Users/llls2012163.com/Lung cancer/Harmnoise_immToCancer',x)
  file_name <- list.files(file_path)
  for(j in file_name){
    path_promt <- paste0(file_path,'/',j)
    data_promt <- read.delim(path_promt)
    het <- mr_heterogeneity(data_promt)
    if(dim(mr_het)[1]==0){
      mr_het <- het
    }else{
      mr_het <- rbind(mr_het,het)
    }
  }
  outcome_path <- paste0('/Users/llls2012163.com/Lung cancer/result_file/MR_heterogeneity/',x,'.txt')
  write.table(mr_het,file = outcome_path,row.names = F,sep = '\t',quote = F)
}
######## MR pleiotropy
MR_pleiotropy_analysis <- function(x){
  mr_ple <- data.frame()
  file_path <- paste0('/Users/llls2012163.com/Lung cancer/Harmnoise_immToCancer',x)
  file_name <- list.files(file_path)
  for(j in file_name){
    path_promt <- paste0(file_path,'/',j)
    data_promt <- read.delim(path_promt)
    ple <- mr_pleiotropy_test(data_promt)
    if(dim(mr_ple)[1]==0){
      mr_ple <- ple
    }else{
      mr_ple <- rbind(mr_ple,ple)
    }
    outcome_path <- paste0('/Users/llls2012163.com/Lung cancer/result_file/MR_pleiotropy/',x,'.txt')
    write.table(mr_ple,file = outcome_path,row.names = F,sep = '\t',quote = F)
  }
}
######### MR presso
doMR_presso <- function(x){# dir_name
  File_name_condi <- vector()
  print(x)
  file_path <- paste0('/Users/llls2012163.com/Lung cancer/Harmnoise_immToCancer',x)
  file_name <- list.files(file_path)
  for(j in file_name){
    File_path <- paste0(file_path,'/',j)
    data_promt <- read.delim(File_path)
    if(nrow(data_promt) >= 4){File_name_promt <- j} else{
      File_name_promt <- 'NA'
    }
    if(length(File_name_condi)==0){
      File_name_condi = File_name_promt
    }else{
      File_name_condi = c(File_name_condi,File_name_promt)
    }
  }
  File_name <- File_name_condi[File_name_condi !='NA']
  presso_mr <- data.frame()
  for(m in File_name){
    path <- paste0('/Users/llls2012163.com/Lung cancer/Harmnoise_immToCancer',x,'/',m)
    print(paste0('processing:',m))
    data_promt <- read.delim(path)
    het_adb_ada <- MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome',BetaExposure = "beta.exposure",SdOutcome = 'se.outcome',SdExposure = "se.exposure",
                                       OUTLIERtest = TRUE,DISTORTIONtest = TRUE,data = data_promt,NbDistribution = 1000,SignifThreshold = 0.05)
    result_promt <- cbind(id.exposure = data_promt$id.exposure[2],het_adb_ada$`Main MR results`)
    result_promt$global_pval <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$Pvalue,'NA')
    result_promt$global_RSSobs <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$RSSobs,'NA')
    result_promt$distortion <- c(het_adb_ada$`MR-PRESSO results`$`Distortion Test`$Pvalue,'NA')
    if(dim(presso_mr)[1]==0){
      presso_mr=result_promt
    }else{
      presso_mr=rbind(presso_mr,result_promt)
    }
    outcome_path <- paste0('/Users/llls2012163.com/Lung cancer/result_file/MR_presso/',x,'.txt')
    write.table(presso_mr,file = outcome_path,row.names = F,sep = '\t',quote = F)
  }
}
#This function resolve the problem aimed to the doMR_presso function.
doMR_presso_continue <- function(x,y){# dir_name
  File_name_condi <- vector()
  print(x)
  outcome_path <- paste0('/Users/llls2012163.com/Lung cancer/result_file/MR_presso/',x,'.txt')
  if(!file.exists(outcome_path)){
    file_promt <- c("id.exposure","Exposure","MR.Analysis", "Causal.Estimate","Sd","T.stat","P.value","global_pval",
                    "global_RSSobs","distortion")
    writeLines(file_promt,con =  outcome_path,sep = '\t',useBytes = FALSE)
  }
  file_path <- paste0('/Users/llls2012163.com/Lung cancer/Harmnoise_immToCancer',x)
  file_name <- list.files(file_path)
  for(j in file_name){
    File_path <- paste0(file_path,'/',j)
    data_promt <- read.delim(File_path)
    if(nrow(data_promt) >= 4){File_name_promt <- j} else{
      File_name_promt <- 'NA'
    }
    if(length(File_name_condi)==0){
      File_name_condi = File_name_promt
    }else{
      File_name_condi = c(File_name_condi,File_name_promt)
    }
  }
  File_name_condi <- File_name_condi[File_name_condi !='NA']
  path_target <- paste0('/Users/llls2012163.com/Lung cancer/result_file/MR_presso/',x,'.txt')
  result_exist <- read.delim(path_target)
  File_name_exist <- result_exist[,1][!duplicated(result_exist[,1])]
  File_name_exist <- paste0(File_name_exist,'.txt')
  File_name <- setdiff(File_name_condi,File_name_exist)
  presso_mr <- data.frame()
  for(m in File_name){
    path <- paste0('/Users/llls2012163.com/Lung cancer/Harmnoise_immToCancer',x,'/',m)
    print(paste0('processing:',m))
    data_promt <- read.delim(path)
    if(dim(data_promt)[1] > 300){
      data_promt <- data_promt[data_promt$pval.exposure < 5*10^-8 & data_promt$F_value > y,]
    }else{
      data_promt <- data_promt
    }
    het_adb_ada <- MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome',BetaExposure = "beta.exposure",SdOutcome = 'se.outcome',SdExposure = "se.exposure",
                                       OUTLIERtest = TRUE,DISTORTIONtest = TRUE,data = data_promt,NbDistribution = 1000,SignifThreshold = 0.05)
    result_promt <- cbind(id.exposure = data_promt$id.exposure[2],het_adb_ada$`Main MR results`)
    result_promt$global_pval <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$Pvalue,'NA')
    result_promt$global_RSSobs <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$RSSobs,'NA')
    if(dim(data_promt)[1]>300){
      result_promt$distortion <- c('NA','NA')
    }else{
      result_promt$distortion <- c(het_adb_ada$`MR-PRESSO results`$`Distortion Test`$Pvalue,'NA')
    }
    presso_mr <- result_promt
    write.table(presso_mr,file = outcome_path,row.names = F,sep = '\t',quote = F,col.names = FALSE,append = TRUE)
  }
}

#--------------------------------------------------From cancertype to immunecell trait
######## Get IVs with cancer types
Get_IVsCancertype <- function(x){ # x dirNames # y gwasID
  for(i in x){
    path_promt_pre <- paste0("/Users/llls2012163.com/Lung cancer/CancerTypeExposure/",i,'/')
    print(i)
    dirNames <- c("Lung adenocarcinoma","Lung cancer","lung cancer_any_parental","lung cancer_ever_smoker",
                  "lung cancer_father","lung cancer_mother","lung cancer_no_smoker","lung cancer_sibilings",
                  "Small-cell lung carcinoma","Squamous cell lung cancer",'No small cell lung carcinoma')
    place <- which(dirNames == i)
    gwasID <- c('ieu-a-984','ukb-a-54','ebi-a-GCST90013972','ebi-a-GCST004749','ukb-b-14521',
                'ukb-b-20176','ebi-a-GCST004747','ukb-b-15826','ieu-a-988','ieu-a-989','finn-b-C3_LUNG_NOSMALL')
    id <- gwasID[place]
    path_promt_finall <- paste0(path_promt_pre,id,'.txt')
    data_promt <- TwoSampleMR::extract_instruments(outcome=id,
                                                   p1= 5e-8,
                                                   clump = T,
                                                   p2= 5e-8,
                                                   r2 = 0.001,
                                                   kb = 10000)
    write.table(data_promt,file = path_promt_finall,row.names = F,sep = '\t',quote = F)
  }
}
######### Download outcome data for cancer type on immune traits
getOutcome_cancerType <- function(x){# x dirNames
  for(j in x){
    # detect target directory
    print(paste0('detecting: ',j))
    outcome_path <- paste0('/Users/llls2012163.com/Lung cancer/LungToimmunetrait','/',j)
    File_downloaded <- list.files(outcome_path)
    File_name <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
    File_undownloaded <- setdiff(File_name,File_downloaded)
    
    File_empty <- c()
    for(i in File_downloaded){
      path_promt <- paste0(outcome_path,'/',i)
      size_promt <- file.size(path_promt)
      if(length(File_empty) == 0){
        File_empty <- size_promt
      }else{
        File_empty <- c(File_empty,size_promt)
      }
    }
    
    File_empty <- cbind(trait = File_downloaded,File_empty)
    File_empty <- as.data.frame(File_empty)
    File_empty <- File_empty[File_empty[,2]==1,][,1]
    File_empty <- c(File_undownloaded,File_empty)
    if(!length(File_empty) ==0){
      File_empty <- limma::strsplit2(File_empty,'\\.')[,1]
      path_exposure <- paste0('/Users/llls2012163.com/Lung cancer/CancerTypeExposure','/',j)
      Gwas_ID_exposure <- list.files(path_exposure)
      path_exposure_promt <- paste0(path_exposure,'/',Gwas_ID_exposure)
      exposure_IV <- read.delim(path_exposure_promt)
      print(Gwas_ID_exposure)
      for(k in File_empty){
        data_promt <- TwoSampleMR::extract_outcome_data(snps = exposure_IV$SNP,outcomes = k)
        print(k)
        outcome_path_promt <- paste0(outcome_path,'/',k,'.txt')
        write.table(data_promt,file = outcome_path_promt,row.names = F,sep = '\t',quote = F)
      }
    }
  }
}
#### harmonise data
Harmose_cancerType_to_immune_fun <- function(x){# x dirName
  for(i in x){
    print(i)
    path_exposure <- paste0("/Users/llls2012163.com/Lung cancer/CancerTypeExposure/",i)
    path_outcome <- paste0("/Users/llls2012163.com/Lung cancer/LungToimmunetrait/",i,'/')
    Target_path <- paste0("/Users/llls2012163.com/Lung cancer/Harmonise_cancerType_To_immune/",i,'/')
    path_file <- list.files(path_exposure)
    path_promt <- paste0(path_exposure,'/',path_file)
    data_exposure <- read.delim(path_promt)
    File_name <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
    for(j in File_name){
      print(j)
      path_outcome_promt <- paste0(path_outcome,j)
      if(!file.size(path_outcome_promt)==0){
        data_outcome <- read.delim(path_outcome_promt)
        data_promt <- TwoSampleMR::harmonise_data(data_exposure,data_outcome)
        Target_path_promt <- paste0(Target_path,j)
        write.table(data_promt,file = Target_path_promt,row.names = F,sep = '\t',quote = F)
      }else{
        Target_path_promt <- paste0(Target_path,j)
        file.create(Target_path_promt)
      }
      
    }
  }
}
#### calculate F value and remove the weak IVs
calculate_F_lager_than10_exist <- function(x){# dirName
  for(j in x){
    print(j)
    file_path <- paste0('/Users/llls2012163.com/Lung cancer/Harmonise_cancerType_To_immune/',j)
    file_name <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
    for(i in file_name){
      print(i)
      path_promt <- paste0(file_path,'/',i)
      if(!file.size(path_promt) ==0){
        data_promt <- read.delim(path_promt)
        R_squre <- TwoSampleMR::get_r_from_bsen(data_promt$beta.exposure,data_promt$se.exposure,data_promt$samplesize.exposure)
        F_value <- as.numeric(R_squre)*(as.numeric(data_promt$samplesize.exposure) - 2)/(1 - as.numeric(R_squre))
        data_promt <- cbind(data_promt, R2 = as.numeric(R_squre), F_value = F_value)
        data_finall <- data_promt[as.numeric(data_promt$F_value) > 10,]
        write.table(data_finall,file = path_promt,row.names = F,sep = '\t',quote = F)
      }else{
        file.create(path_promt)
      }
      
    }
  }
}




#------------------------------------------------------------------------------
# ukb-a-54 :Lung cancer.     1----------------
dir.create("/Users/llls2012163.com/Lung cancer/Lung cancer")
File_name_empty <- Datacheck(File_name,'Lung cancer')
Download.outcome(File_name_empty,'ukb-a-54','Lung cancer')
dir.create("/Users/llls2012163.com/Lung cancer/ImmToLung")
Harmose_Fun(File_name,'Lung cancer','ImmToLung')

# ieu-a-984 :Lung adenocarcinoma.    2-------------------------
dir.create("/Users/llls2012163.com/Lung cancer/Lung adenocarcinoma")
File_name_empty <- Datacheck(File_name,'Lung adenocarcinoma')
Download.outcome(File_name_empty,'ieu-a-984','Lung adenocarcinoma')
dir.create("/Users/llls2012163.com/Lung cancer/ImmToLungAdenocarinoma")
Harmose_Fun(File_name,'Lung adenocarcinoma','ImmToLungAdenocarinoma')

# ieu-a-988 :small cell lung carcinoma.  3----------------------
dir.create("/Users/llls2012163.com/Lung cancer/Small-cell lung carcinoma")
File_name_empty <- Datacheck(File_name,'Small-cell lung carcinoma')
Download.outcome(File_name_empty,'ieu-a-988','Small-cell lung carcinoma')
dir.create("/Users/llls2012163.com/Lung cancer/immToSmallcell")
Harmose_Fun(File_name,'Small-cell lung carcinoma','immToSmallcell')

# ieu-a-989 : squamous cell lung cancer.  4--------------------
dir.create("/Users/llls2012163.com/Lung cancer/Squamous cell lung cancer")
File_name_empty <- Datacheck(File_name,'Squamous cell lung cancer')
Download.outcome(File_name_empty,'ieu-a-989','Squamous cell lung cancer')
dir.create("/Users/llls2012163.com/Lung cancer/immToSquamous")
Harmose_Fun(File_name,'Squamous cell lung cancer','immToSquamous')

# finn-b-C3_LUNG_NONSMALL : no small lung cancer.     5---------------
dir.create("/Users/llls2012163.com/Lung cancer/NosmallCell")
File_name_empty <- Datacheck(File_name,'NosmallCell')
Download.outcome(File_name_empty,'ieu-a-989','NosmallCell')
dir.create("/Users/llls2012163.com/Lung cancer/immToNosmallCell")
Harmose_Fun(File_name,'NosmallCell','immToNosmallCell')

# ukb-b-14521 : lung cancer father.    6------------------
dir.create("/Users/llls2012163.com/Lung cancer/lungCancer_father")
File_name_empty <- Datacheck(File_name,'lungCancer_father')
Download.outcome(File_name_empty,'ukb-b-14521','lungCancer_father')
dir.create("/Users/llls2012163.com/Lung cancer/immTolungCancer_father")
Harmose_Fun(File_name,'lungCancer_father','immTolungCancer_father')

# ukb-b-20176 : lung cancer mother.    7-------------------
dir.create("/Users/llls2012163.com/Lung cancer/lungCancer_mother")
File_name_empty <- Datacheck(File_name,'lungCancer_mother')
Download.outcome(File_name_empty,'ukb-b-20176','lungCancer_mother')
dir.create("/Users/llls2012163.com/Lung cancer/immTolungCancer_mother")
Harmose_Fun(File_name,'lungCancer_mother','immTolungCancer_mother')

# ebi-a-GCST90013972 : lung cancer any parental.  8-----------------
dir.create("/Users/llls2012163.com/Lung cancer/lungCancer_AnyParental")
File_name_empty <- Datacheck(File_name,'lungCancer_AnyParental')
Download.outcome(File_name_empty,'ebi-a-GCST90013972','lungCancer_AnyParental')
dir.create("/Users/llls2012163.com/Lung cancer/immTolungCancer_AnyParental")
Harmose_Fun(File_name,'lungCancer_AnyParental','immTolungCancer_AnyParental')

# ukb-b-15826: lung cancer sibilings.  9-------------------------
dir.create("/Users/llls2012163.com/Lung cancer/lungCancer_sibilings")
File_name_empty <- Datacheck(File_name,'lungCancer_sibilings')
Download.outcome(File_name_empty,'ukb-b-15826','lungCancer_sibilings') # ebi-a-GCST90002092 do not map the outcome information and skiped
empty <- disctinct_otherwise_empty_function(File_name,'lungCancer_sibilings')
dir.create("/Users/llls2012163.com/Lung cancer/immTolungCancer_sibilings")
Harmose_Fun(File_name[-which(File_name == empty)],'lungCancer_sibilings','immTolungCancer_sibilings')

# ebi-a-GCST004747: lung cancer no smoker.    10--------------------------
dir.create("/Users/llls2012163.com/Lung cancer/Lung cancerNosmok")
File_name_empty <- Datacheck(File_name,'Lung cancerNosmok')
Download.outcome(File_name_empty,'ebi-a-GCST004747','Lung cancerNosmok')
dir.create("/Users/llls2012163.com/Lung cancer/immTolung cancerNosmok")
Harmose_Fun(File_name,'Lung cancerNosmok','immTolung cancerNosmok')

# ebi-a-GCST004749 : lung cancer smoker.   11----------------------
dir.create("/Users/llls2012163.com/Lung cancer/Lung cancerSmok")
File_name_empty <- Datacheck(File_name,'Lung cancerSmok')
Download.outcome(File_name_empty,'ebi-a-GCST004749','Lung cancerSmok')
dir.create("/Users/llls2012163.com/Lung cancer/immTolung cancerSmok")
Harmose_Fun(File_name,'Lung cancerSmok','immTolung cancerSmok')


#-------------------------------------------------------------------------------
## remove the weak IVs
dir_name <- list.files("/Users/llls2012163.com/Lung cancer")[3:13]
for(i in dir_name){
  print(i)
  calculate_F_lager_tha10_exist(i)
}

### MR analysis
dir.create('/Users/llls2012163.com/Lung cancer/result_file/MR_result')
for(i in dir_name){
  print(i)
  MR_analysis(i)
}

#### MR heterogeneity analysis
dir.create('/Users/llls2012163.com/Lung cancer/result_file/MR_heterogeneity')
for(i in dir_name){
  print(i)
  MR_heterogeneity_analysis(i)
}

##### MR pleiotropy analysis
dir.create('/Users/llls2012163.com/Lung cancer/result_file/MR_pleiotropy')
for(i in dir_name){
  print(i)
  MR_pleiotropy_analysis(i)
}

###### MR presso
dir.create('/Users/llls2012163.com/Lung cancer/result_file/MR_presso')
dir_name <- list.files("/Users/llls2012163.com/Lung cancer")[3:13][length(file):13]
for(i in dir_name){
  doMR_presso_continue(i,300)
}


#####----------------------------------------------------from cancer type to immune cell traits.

###### Get IVs from cancer type to immune cell traits
dir.create("/Users/llls2012163.com/Lung cancer/CancerTypeExposure")
dirNames <- c("Lung adenocarcinoma","Lung cancer","lung cancer_any_parental","lung cancer_ever_smoker",
              "lung cancer_father","lung cancer_mother","lung cancer_no_smoker","lung cancer_sibilings",
              "Small-cell lung carcinoma","Squamous cell lung cancer",'No small cell lung carcinoma')
gwasID <- c('ieu-a-984','ukb-a-54','ebi-a-GCST90013972','ebi-a-GCST004749','ukb-b-14521',
            'ukb-b-20176','ebi-a-GCST004747','ukb-b-15826','ieu-a-988','ieu-a-989','finn-b-C3_LUNG_NONSMALL')

for(i in dirNames){
  path_promt <- paste0("/Users/llls2012163.com/Lung cancer/CancerTypeExposure/",i)
  dir.create(path_promt)
}
Get_IVsCancertype(dirNames)

## The finngen_R10_C3_LUNG_NONSMALL_EXALLC can not download from the IEU GWAS project and downloaded from the FinnG database.
## https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_C3_LUNG_NONSMALL_EXALLC.gz
library(data.table)
data_orignial <- fread('~/Lung cancer/finngen_R10_C3_LUNG_NONSMALL_EXALLC.gz',header = TRUE)
data_orignial <- data_orignial[data_orignial$pval < 5e-6,]
library(TwoSampleMR)
exposure <- format_data(data_orignial,
                        type = 'exposure',
                        snp_col = 'rsids',
                        phenotype_col = 'phenotypes',
                        beta_col = 'beta',
                        se_col = 'sebeta',
                        eaf_col = 'af_alt',
                        effect_allele_col = 'alt',
                        other_allele_col = 'ref',
                        pval_col = 'pval')
exposure$exposure <-  c(rep('C3_LUNG_NONSMALL_EXALLC',dim(exposure)[1]))
## clump data
exposure_data <- TwoSampleMR::clump_data(exposure,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 5e-6,clump_p2 = 5e-8)
write.table(exposure_data,file = '~/Lung cancer/CancerTypeExposure/No small cell lung carcinoma/finn_C3_LUNG_NONSMALL_EXALLC.txt',row.names = F,sep = '\t',quote = F)

#######Extract outcome with cancer type of immune cell traits
dir.create('/Users/llls2012163.com/Lung cancer/LungToimmunetrait')
dircancertypeToimm <- paste0('/Users/llls2012163.com/Lung cancer/LungToimmunetrait','/',dirNames)
for(i in dirNames){
  path_promt <- paste0('/Users/llls2012163.com/Lung cancer/LungToimmunetrait','/',i)
  dir.create(path_promt)
}
getOutcome_cancerType(dirNames)

## Harminose GWAS cancerType (exposure) on immune cell traits (outcome).

dir.create("/Users/llls2012163.com/Lung cancer/Harmonise_cancerType_To_immune")
for(i in dirNames){
  path_promt <- paste0("/Users/llls2012163.com/Lung cancer/Harmonise_cancerType_To_immune/",i)
  dir.create(path_promt)
}
Harmose_cancerType_to_immune_fun <- function(x = dirNames){# x dirName
  for(i in x){
    print(i)
    path_exposure <- paste0("/Users/llls2012163.com/Lung cancer/CancerTypeExposure/",i)
    path_outcome <- paste0("/Users/llls2012163.com/Lung cancer/LungToimmunetrait/",i,'/')
    Target_path <- paste0("/Users/llls2012163.com/Lung cancer/Harmonise_cancerType_To_immune/",i,'/')
    path_file <- list.files(path_exposure)
    path_promt <- paste0(path_exposure,'/',path_file)
    data_exposure <- read.delim(path_promt)
    File_name <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
    for(j in File_name){
      print(j)
      path_outcome_promt <- paste0(path_outcome,j)
      if(!file.size(path_outcome_promt)==0){
        data_outcome <- read.delim(path_outcome_promt)
        data_promt <- TwoSampleMR::harmonise_data(data_exposure,data_outcome)
        Target_path_promt <- paste0(Target_path,j)
        write.table(data_promt,file = Target_path_promt,row.names = F,sep = '\t',quote = F)
      }else{
        Target_path_promt <- paste0(Target_path,j)
        file.create(Target_path_promt)
      }
      
    }
 }
}

Harmose_cancerType_to_immune_fun(dirNames)

# add the samplesize.exposure to the dirName[11]. This information were obtained from finnG. 218792 samples.
file_name <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
for(i in file_name){
  file_path <- paste0('/Users/llls2012163.com/Lung cancer/Harmonise_cancerType_To_immune/',dirNames[11],'/',i)
  data_promt <- read.delim(file_path)
  data_promt$samplesize.exposure <- rep(218792,dim(data_promt)[1])
  write.table(data_promt,file = file_path,row.names = F,sep = '\t',quote = F)
}


## Remove IVs with F value < 10

calculate_F_lager_than10_exist <- function(x = dirNames){# dirName
  for(j in x){
    print(j)
    file_path <- paste0('/Users/llls2012163.com/Lung cancer/Harmonise_cancerType_To_immune/',j)
    file_name <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
    for(i in file_name){
      print(i)
      path_promt <- paste0(file_path,'/',i)
      if(!file.size(path_promt) ==0){
        data_promt <- read.delim(path_promt)
        R_squre <- TwoSampleMR::get_r_from_bsen(data_promt$beta.exposure,data_promt$se.exposure,data_promt$samplesize.exposure)
        F_value <- as.numeric(R_squre)*(as.numeric(data_promt$samplesize.exposure) - 2)/(1 - as.numeric(R_squre))
        data_promt <- cbind(data_promt, R2 = as.numeric(R_squre), F_value = F_value)
        data_finall <- data_promt[as.numeric(data_promt$F_value) > 10,]
        write.table(data_finall,file = path_promt,row.names = F,sep = '\t',quote = F)
      }else{
        file.create(path_promt)
      }
      
      }
    }
  }
calculate_F_lager_than10_exist(dirNames)

## MR analysis from cancerType to immune cell traits
dir.create("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_result")
MR_analysis_CancerType_to_immune <- function(x = dirNames){
  for(i in x){
    print(i)
    path_promt <- paste0('/Users/llls2012163.com/Lung cancer/Harmonise_cancerType_To_immune/',i)
    immid <- paste0('ebi-a-GCST9000',c(1391:2121))
    output_name <- paste0(i,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_result/",output_name)
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.outcome[!duplicated( output_promt$id.outcome)]
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j,'.txt')
      print(j)
      
      if(!file.size(path_promt_pre)==0){
        data_promt <- read.delim(path_promt_pre)
        res <- TwoSampleMR::mr(data_promt,method_list = c('mr_ivw_fe','mr_ivw_mre','mr_egger_regression','mr_weighted_median','mr_weighted_mode','mr_simple_mode'))
      }  
      
      if(j == immid[1]){
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
      }else{
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
      }
    }
    
  }
}

MR_analysis_CancerType_to_immune(dirNames)

## MR heterogeneity from cancerType to immune cell traits
dir.create("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_heterogeneity")
MR_heterogeneity_CancerType_to_immune <- function(x = dirNames){
  for(i in x){
    print(i)
    path_promt <- paste0('/Users/llls2012163.com/Lung cancer/Harmonise_cancerType_To_immune/',i)
    immid <- paste0('ebi-a-GCST9000',c(1391:2121))
    output_name <- paste0(i,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_heterogeneity/",output_name)
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.outcome[!duplicated( output_promt$id.outcome)]
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j,'.txt')
      print(j)
      
      if(!file.size(path_promt_pre)==0){
        data_promt <- read.delim(path_promt_pre)
        res <- TwoSampleMR::mr_heterogeneity(data_promt)
      }  
      
      if(j == immid[1]){
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
      }else{
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
      }
    }
    
  }
}

MR_heterogeneity_CancerType_to_immune(dirNames)

## MR pleiotropy from cancerType to immune cell traits
dir.create("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_pleiotropy")
MR_pleiotropy_CancerType_to_immune <- function(x = dirNames){
  for(i in x){
    print(i)
    path_promt <- paste0('/Users/llls2012163.com/Lung cancer/Harmonise_cancerType_To_immune/',i)
    immid <- paste0('ebi-a-GCST9000',c(1391:2121))
    output_name <- paste0(i,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_pleiotropy/",output_name)
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.outcome[!duplicated( output_promt$id.outcome)]
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j,'.txt')
      print(j)
      
      if(!file.size(path_promt_pre)==0){
        data_promt <- read.delim(path_promt_pre)
        res <- TwoSampleMR::mr_pleiotropy_test(data_promt)
      }  
      
      if(j == immid[1]){
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
      }else{
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
      }
    }
    
  }
}

MR_pleiotropy_CancerType_to_immune(dirNames)

## MR presso from cancerType to immune cell traits
dir.create('/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_presso')



MR_presso_CancerType_to_immune <- function(x = dirNames){
  for(i in x){
    print(i)
    path_promt <- paste0('/Users/llls2012163.com/Lung cancer/Harmonise_cancerType_To_immune/',i)
    immid <- paste0('ebi-a-GCST9000',c(1391:2121))
    output_name <- paste0(i,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_presso/",output_name)
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.outcome[!duplicated(output_promt$id.outcome)]
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j,'.txt')
      print(j)
      if(file.size(path_promt_pre)==0){message('Null file: the outcome do not mapped SNPs')}
      
      if(!file.size(path_promt_pre)==0){
        data_promt <- read.delim(path_promt_pre)
        if(dim(data_promt)[1] >=4){
          het_adb_ada <- MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome',BetaExposure = "beta.exposure",SdOutcome = 'se.outcome',SdExposure = "se.exposure",
                                             OUTLIERtest = TRUE,DISTORTIONtest = TRUE,data = data_promt,NbDistribution = 1000,SignifThreshold = 0.05)
          result_promt <- cbind(id.exposure = data_promt$id.exposure[2],het_adb_ada$`Main MR results`)
          result_promt$global_pval <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$Pvalue,'NA')
          result_promt$global_RSSobs <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$RSSobs,'NA')
          result_promt$id.outcome <- rep(j,dim(result_promt)[1])
          if(j == immid[1]){
            write.table(result_promt,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
          }else{
            write.table(result_promt,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
          }
        }else{
          message('Number of SNP < 4')
        }
        
      }  
      
    }
    
  }
}
MR_presso_CancerType_to_immune(dirNames)

## multithreaded perform MR Presso
library(snow)
c1 <- makeCluster(8,type = 'SOCK')
clusterEvalQ(c1,library(snowfall))
clusterApply(c1,dirNames,MR_presso_CancerType_to_immune)

#### 41 Circulating inflammatory cytokines GWAS were downloaded from University of BRISTOL health sciences database
#### https://data.bris.ac.uk/data/dataset/3g3i5smgghp0s2uvm1doflkx9x
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure")
File_name <- list.files("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inflammation cytosines")#--------------------------------------------------####
for(i in File_name){
  print(i)
  input_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inflammation cytosines/",i)
  data_promt <- read.delim(input_path)
  data_promt <- data_promt[data_promt$P.value < 5e-6,]
  output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure/",i)
  write.table(data_promt,file = output_path,row.names = F,sep = '\t',quote = F,append = FALSE,col.names = TRUE)
}

## clumping data
for(i in File_name){
  print(i)
  path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure/",i)
  data_promt <- read.delim(path_promt)
  phenotype <- limma::strsplit2(i,'\\.')[,1]
  data_promt$phenotype <- rep(phenotype,dim(data_promt)[1])
  data_promt <- TwoSampleMR::format_data(data_promt,
                          type = 'exposure',
                          snp_col = 'MarkerName',
                          phenotype_col = 'phenotype',
                          beta_col = 'Effect',
                          se_col = 'StdErr',
                          eaf_col = 'Freq1',
                          effect_allele_col = 'Allele1',
                          other_allele_col = 'Allele2',
                          pval_col = 'P.value')
  clump_data <- TwoSampleMR::clump_data(data_promt,clump_p1 = 5e-6,clump_kb = 10000,clump_p2 = 5e-8,clump_r2 = 0.001)
  write.table(clump_data,file = path_promt,row.names = F,sep = '\t',quote = F,append = FALSE,col.names = TRUE)
}
## map IVs with inflammatory cytokines on the immune cell traits-------------------------------------
# we first create a IV poor for inflammatory cytokines on 731 immune cell traits
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/infToimm_IVpoor")
Read_SNP_inf <- function(x = 'File_name'){
  SNP <- c()
  for(i in x){
    path_exposure <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure/",i)
    data_promt <- read.delim(path_exposure)
    SNP_promt <- data_promt$SNP
    if(length(SNP)==0){
      SNP <- SNP_promt
    }else{
      SNP <- c(SNP,SNP_promt)
    }
    SNP <- SNP[!duplicated(SNP)]
  }
  return(SNP)
}
snp_poor <- Read_SNP_inf(File_name)
Map_IVs <- function(x = 'infToimm_IVpoor'){
  # detecting target directory
  imm_ID <- paste0('ebi-a-GCST9000',c(1391:2121))
  path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/",x)
  file <- list.files(path)
  if(!length(file)==0){
    File_downloaded <- limma::strsplit2(file,'\\.')[,1]
  }else{
    File_downloaded <- 'NA'
  }
  
  file_ud <- setdiff(imm_ID,File_downloaded)
  
  File_empty <- c()
  for(i in File_downloaded){
    path_promt <- paste0(path,'/',i,'.txt')
    size_promt <- file.size(path_promt)
    if(length(File_empty) == 0){
      File_empty <- size_promt
    }else{
      File_empty <- c(File_empty,size_promt)
    }
  }
  
  File_empty <- cbind(trait = File_downloaded,File_empty)
  File_empty <- as.data.frame(File_empty)
  File_empty <- File_empty[File_empty[,2]==1,][,1]
  File_empty <- c(file_ud,File_empty)
  
  for(j in File_empty){
    data_promt <- TwoSampleMR::extract_outcome_data(snps = snp_poor,outcomes = j)
    print(paste0('Extracting outcome:',j))
    output_path <- paste0(path,'/',j,'.txt')
    write.table(data_promt,file = output_path,row.names = F,sep = '\t',quote = F)
  }
}
Map_IVs("infToimm_IVpoor")  
## Then map the IVs for each inflammatory cytokines
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome")
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/infToimmune")
for(i in File_name){
  File <- limma::strsplit2(i,'\\.')[,1]
  path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/infToimmune/",File)
  dir.create(path_promt)
}
Select_IV_form_poor <- function(x = 'File_name'){
  imm_ID <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
  for(i in x){
    exposure_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure/",i)
    exposure_data <- read.delim(exposure_path)
    SNP_promt <- exposure_data$SNP
    print(i)
    cytokine <- limma::strsplit2(i,'\\.')[,1]
    output_path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/infToimmune/",cytokine,'/')
    for(j in imm_ID){
      IVpoor_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/infToimm_IVpoor/",j)
      outcome_promt <- read.delim(IVpoor_path)
      outcome_promt <- outcome_promt[outcome_promt$SNP %in% SNP_promt,]
      output_path <- paste0(output_path_promt,j)
      print(j)
      write.table(outcome_promt,output_path,row.names = F,sep = '\t',quote = F)
    }
  }
}
Select_IV_form_poor(File_name)

#### Method 2 dowload data one by one
getOutcome_inf_immFun <- function(x = File_name){# File_name
  for(j in x){
    # detect target directory
    print(paste0('detecting: ',j))
    id <- limma::strsplit2(j,'\\.')[,1]
    outcome_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/infToimmune/",id)
    File_downloaded <- list.files(outcome_path)
    File_outcome <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
    File_undownloaded <- setdiff(File_outcome,File_downloaded)
    
    File_empty <- c()
    for(i in File_downloaded){
      path_promt <- paste0(outcome_path,'/',i)
      size_promt <- file.size(path_promt)
      if(length(File_empty) == 0){
        File_empty <- size_promt
      }else{
        File_empty <- c(File_empty,size_promt)
      }
    }
    
    File_empty <- cbind(trait = File_downloaded,File_empty)
    File_empty <- as.data.frame(File_empty)
    File_empty <- File_empty[File_empty[,2]==1,][,1]
    File_empty <- c(File_undownloaded,File_empty)
    if(!length(File_empty) ==0){
      File_empty <- limma::strsplit2(File_empty,'\\.')[,1]
      path_exposure <- paste0('/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure/',j)
      exposure_IV <- read.delim(path_exposure)
      print(j)
      for(k in File_empty){
        data_promt <- TwoSampleMR::extract_outcome_data(snps = exposure_IV$SNP,outcomes = k)
        print(k)
        outcome_path_promt <- paste0(outcome_path,'/',k,'.txt')
        write.table(data_promt,file = outcome_path_promt,row.names = F,sep = '\t',quote = F)
      }
    }
  } 
}
getOutcome_inf_immFun(File_name[1])

library(snow)
c1 <- makeCluster(8,type = 'SOCK')
clusterEvalQ(c1,library(snowfall))
clusterApply(c1,File_name[4:12],getOutcome_inf_immFun)

# Harmonise data (inflammatory cytokines on the immune cell traits)
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/infToimmune_harmnoise")
for(i in File_name){
  File <- limma::strsplit2(i,'\\.')[,1]
  path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/infToimmune_harmnoise/",File)
  dir.create(path_promt)
}
Harmose_infToimme <- function(x = 'File_name'){
  imm_ID <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
  for(i in x){
    File <- limma::strsplit2(i,'\\.')[,1]
    exposure_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure/",i)
    outcome_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/infToimmune/",File)
    output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/infToimmune_harmnoise/",File)
    exposure_data <- read.delim(exposure_path)
    message(paste0('Processing:',i))
    for(j in imm_ID){
      outcome_path_promt <- paste0(outcome_path,'/',j)
      output_path_promt <- paste0(output_path,'/',j)
      if(!file.size(outcome_path_promt)== 1){
        outcome_data <- read.delim(outcome_path_promt)
        res <- TwoSampleMR::harmonise_data(exposure_data,outcome_data)
        write.table(res,output_path_promt,row.names = F,sep = '\t',quote = F)
        print(j)
      }else{
        message(paste0(i,'did not mapped any IV'))
        print(j)
      }
      
    }
  }
}
Harmose_infToimme(File_name)

# Add the exposure.sample.size to harmonise data

Add_samplesizeFun <- function(x = 'File_name'){
  for(i in x){
    path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inflammation cytosines/",i)
    print(i)
    Sample_size <- readLines(path_promt,n = 3)[2]
    Sample_size <- limma::strsplit2(Sample_size,'\t')[,16]
    Sample_size <- as.numeric(Sample_size)
    File <- limma::strsplit2(i,'\\.')[,1]
    path_target <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/infToimmune_harmnoise/",File)
    imm_ID <- list.files(path_target)
    for(j in imm_ID){
      path_target_promt <- paste0(path_target,'/',j)
      print(j)
      if(file.size(path_target_promt)== 1){
        message('No IVs were mapped')
      }else{
        data_promt <- read.delim(path_target_promt)
        data_promt <- cbind(data_promt,samplesize.exposure = rep(Sample_size,dim(data_promt)[1]))
        write.table(data_promt,file = path_target_promt,row.names = F,sep = '\t',quote = F)
      }
    }
  }
}
Add_samplesizeFun(File_name)

# caculate F value
calculate_F_lager_than10_exist_infToimmFun <- function(x = File_name){# dirName
  for(i in x){
    print(i)
    File <- limma::strsplit2(i,'\\.')[,1]
    file_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/infToimmune_harmnoise/",File)
    imm_id <- list.files(file_path)
    for(j in imm_id){
      path_promt <- paste0(file_path,'/',j)
      if(!file.size(path_promt) ==0){
        data_promt <- read.delim(path_promt)
        R_squre <- TwoSampleMR::get_r_from_bsen(data_promt$beta.exposure,data_promt$se.exposure,data_promt$samplesize.exposure)
        F_value <- as.numeric(R_squre)*(as.numeric(data_promt$samplesize.exposure) - 2)/(1 - as.numeric(R_squre))
        data_promt <- cbind(data_promt, R2 = as.numeric(R_squre), F_value = F_value)
        data_finall <- data_promt[as.numeric(data_promt$F_value) > 10,]
        write.table(data_finall,file = path_promt,row.names = F,sep = '\t',quote = F)
        print(j)
      }else{
        message('No IVs were mapped')
      }
      
    }
  }
}

calculate_F_lager_than10_exist_infToimmFun(File_name)

# MR analysis
# create result file
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results")
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_resultinfToimm")

MR_analysis_infToimmFun <- function(x = File_name){
  for(i in x){
    print(i)
    File <- limma::strsplit2(i,'\\.')[,1]
    path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/infToimmune_harmnoise/",File)
    immid <- list.files(path_promt)
    output_name <- paste0(File,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_resultinfToimm/",output_name)
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.outcome[!duplicated( output_promt$id.outcome)]
      MR_done <- paste0(MR_done,'.txt')
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j)
      print(j)
     
      if(!file.size(path_promt_pre)==0){
          data_promt <- read.delim(path_promt_pre)
          res <- TwoSampleMR::mr(data_promt,method_list = c('mr_ivw_fe','mr_ivw_mre','mr_egger_regression','mr_weighted_median','mr_weighted_mode','mr_simple_mode'))
      }else{
          message('No ivs were mapped')
        }  
      

      if(j == immid[1]){
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
      }else{
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
      }
    }
    
  }
}
library(snow)
c1 <- makeCluster(8,type = 'SOCK')
clusterEvalQ(c1,library(snowfall))
clusterApply(c1,File_name,MR_analysis_infToimmFun)
MR_analysis_infToimmFun(File_name)

# heterogeneity
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_HeterogeneityinfToimm")
MR_heterogeneity_inf_to_immune <- function(x = File_name){
  for(i in x){
    print(i)
    File <- limma::strsplit2(i,'\\.')[,1]
    path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/infToimmune_harmnoise/",File)
    immid <- paste0('ebi-a-GCST9000',c(1391:2121))
    output_name <- paste0(File,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_HeterogeneityinfToimm/",output_name)
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.outcome[!duplicated( output_promt$id.outcome)]
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j,'.txt')
      print(j)
      
      if(file.exists(path_promt_pre)){
        data_promt <- read.delim(path_promt_pre)
        if(dim(data_promt)[1] >1){
          res <- TwoSampleMR::mr_heterogeneity(data_promt)
        }else{
          message(paste0(i,': Not enough SNPs'))
        }
        
      }  
      
      if(j == immid[1]){
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
      }else{
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
      }
    }
    
  }
}
library(snow)
c1 <- makeCluster(8,type = 'SOCK')
clusterEvalQ(c1,library(snowfall))
clusterApply(c1,File_name,MR_heterogeneity_inf_to_immune)
MR_heterogeneity_inf_to_immune()

# MR pleiotropy
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_pleiotrpyinfToimm")
MR_pleiotropy_inf_to_immune <- function(x = File_name){
  for(i in x){
    print(i)
    File <- limma::strsplit2(i,'\\.')[,1]
    path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/infToimmune_harmnoise/",File)
    immid <- paste0('ebi-a-GCST9000',c(1391:2121))
    output_name <- paste0(i,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_pleiotrpyinfToimm/",output_name)
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.outcome[!duplicated( output_promt$id.outcome)]
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j,'.txt')
      print(j)
      
      if(file.exists(path_promt_pre)){
        data_promt <- read.delim(path_promt_pre)
        if(dim(data_promt)[1] >1){
          res <- TwoSampleMR::mr_pleiotropy_test(data_promt)
        }else{
          message(j,': not enoguh snps')
        }
        
      }else{
        message(j,': do not exist')
      }  
      
      if(j == immid[1]){
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
      }else{
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
      }
    }
    
  }
}

library(snow)
c1 <- makeCluster(8,type = 'SOCK')
clusterEvalQ(c1,library(snowfall))
clusterApply(c1,File_name,MR_pleiotropy_inf_to_immune)

MR_pleiotropy_inf_to_immune(File_name)

# MR presso 
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_pressoinfToimm")

MR_presso_inf_to_immune <- function(x = File_name){
  for(i in x){
    print(i)
    File <- limma::strsplit2(i,'\\.')[,1]
    path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/infToimmune_harmnoise/",File)
    immid <- paste0('ebi-a-GCST9000',c(1391:2121))
    output_name <- paste0(i,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_pressoinfToimm/",output_name)
    
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.outcome[!duplicated(output_promt$id.outcome)]
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j,'.txt')
      print(j)
      if(!file.exists(path_promt_pre)){message(paste0(j,': file do not exist'))}
      
      if(file.exists(path_promt_pre)){
        data_promt <- read.delim(path_promt_pre)
        if(dim(data_promt)[1] >=4){
          het_adb_ada <- MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome',BetaExposure = "beta.exposure",SdOutcome = 'se.outcome',SdExposure = "se.exposure",
                                             OUTLIERtest = TRUE,DISTORTIONtest = TRUE,data = data_promt,NbDistribution = 1000,SignifThreshold = 0.05)
          result_promt <- cbind(id.exposure = data_promt$id.exposure[2],het_adb_ada$`Main MR results`)
          result_promt$global_pval <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$Pvalue,'NA')
          result_promt$global_RSSobs <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$RSSobs,'NA')
          result_promt$id.outcome <- rep(j,dim(result_promt)[1])
          if(j == immid[1]){
            write.table(result_promt,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
          }else{
            write.table(result_promt,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
          }
        }else{
          message('Number of SNP < 4')
        }
        
      }  
      
    }
    
  }
}

library(snow)
c1 <- makeCluster(8,type = 'SOCK')
clusterEvalQ(c1,library(snowfall))
clusterApply(c1,File_name[2:41],MR_presso_inf_to_immune)#-----------------------


MR_presso_inf_to_immune()



##---------------------------------------------


## map IVs with immune cell traits on the inflammatory cytokines
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/immuneToinf")
for(i in File_name){
  File <- limma::strsplit2(i,'\\.')[,1]
  path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/immuneToinf/",File)
  dir.create(path_promt)
}

getOutcome_imm_infFun <- function(x = File_name){
  File_imm <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
  for(i in x){
    File <- limma::strsplit2(i,'\\.')[,1]
    output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/immuneToinf/",File)
    Map_down <- list.files(output_path)
    un_do <- setdiff(File_imm,Map_down)
    print(i)
    if(!length(un_do) == 0){
      path_outcome <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inflammation cytosines/",i)
      data_outcome <- read.delim(path_outcome)
      for(j in un_do){
        path_exposure <- paste0("/Users/llls2012163.com/GWAS/Immune cell trait/",j)
        data_exposure <- read.delim(path_exposure)
        data_promt <- data_outcome[data_outcome[,1] %in% data_exposure$SNP,]
        path_promt <- paste0(output_path,'/',j)
        print(j)
        write.table(data_promt,file = path_promt,row.names = F,sep = '\t',quote = F)
    }
    }
  }
}

getOutcome_imm_infFun()

## Rename the column names for files with immune cell traits on the inflammatory cytokines

Re_name_Fun <- function(x = File_name){
  File_imm <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
  for(i in x){
    File <- limma::strsplit2(i,'\\.')[,1]
    path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/immuneToinf/",File)
    print(i)
    for(j in File_imm){
      path_promt_finall <- paste0(path_promt,'/',j)
      data_promt <- read.delim(path_promt_finall)
      cols <- c('SNP','effect_allele.outcome',"other_allele.outcome","eaf.outcome",'seEAF.outcome','minaf.outcome',
                'maxaf.outcome','beta.outcome','se.outcome','pval.outcome','direction.outcome',"HetISq.outcome",
                "HetChiSq.outcome","HetDf.outcome","HetPVal.outcome","samplesize.outcome")
      colnames(data_promt) <- cols
      data_promt$id.outcome <- rep(File,dim(data_promt)[1])
      data_promt$outcome <- rep(File,dim(data_promt)[1])
      print(j)
      write.table(data_promt,file = path_promt_finall,row.names = F,sep = '\t',quote = F)
    }
  }
  
}
Re_name_Fun()

## Harmonise data exposure immune cell traits on outcome inflammatory cytokines
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie")
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/immuneToinf_harmnoise")
for(i in File_name){
  File <- limma::strsplit2(i,'\\.')[,1]
  path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/immuneToinf_harmnoise/",File)
  dir.create(path_promt)
}
  
Harmnose_imm_infFun <- function(x= File_name){
  File_imm <- paste0('ebi-a-GCST9000',c(1391:2121),'.txt')
  for(i in x){
    File <- limma::strsplit2(i,'\\.')[,1]
    path_outcome_pre <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_outcome/immuneToinf/",File)
    print(i)
    for(j in File_imm){
      print(j)
      path_outcome <- paste0(path_outcome_pre,'/',j)
      data_promt_outcome <- read.delim(path_outcome)
      path_exposure <- paste0("/Users/llls2012163.com/GWAS/Immune cell trait/",j)
      data_promt_exposure <- read.delim(path_exposure)
      data_promt <- TwoSampleMR::harmonise_data(data_promt_exposure,data_promt_outcome,action = 2)
      output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/immuneToinf_harmnoise/",File,'/',j)
      write.table(data_promt,file = output_path,row.names = F,sep = '\t',quote = F)
    }
  }
}
Harmnose_imm_infFun(File_name)

calculate_F_lager_than10_exist_immToinfFun <- function(x = 'File_name'){# dirName
  for(i in x){
    print(i)
    File <- limma::strsplit2(i,'\\.')[,1]
    file_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/immuneToinf_harmnoise/",File)
    imm_id <- list.files(file_path)
    for(j in imm_id){
      path_promt <- paste0(file_path,'/',j)
      if(!file.size(path_promt) ==0){
        data_promt <- read.delim(path_promt)
        R_squre <- TwoSampleMR::get_r_from_bsen(data_promt$beta.exposure,data_promt$se.exposure,data_promt$samplesize.exposure)
        F_value <- as.numeric(R_squre)*(as.numeric(data_promt$samplesize.exposure) - 2)/(1 - as.numeric(R_squre))
        data_promt <- cbind(data_promt, R2 = as.numeric(R_squre), F_value = F_value)
        data_finall <- data_promt[as.numeric(data_promt$F_value) > 10,]
        write.table(data_finall,file = path_promt,row.names = F,sep = '\t',quote = F)
        print(j)
      }else{
        message('No IVs were mapped')
      }
      
    }
  }
}
calculate_F_lager_than10_exist_immToinfFun(File_name)

### MR analysis
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_resultimmToinf")
MR_analysis_immToinfFun <- function(x = File_name){
  for(i in x){
    print(i)
    File <- limma::strsplit2(i,'\\.')[,1]
    path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/immuneToinf_harmnoise/",File)
    immid <- list.files(path_promt)
    output_name <- paste0(File,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_resultimmToinf/",output_name)
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.exposure[!duplicated(output_promt$id.exposure)]
      MR_done <- paste0(MR_done,'.txt')
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j)
      print(j)
      
      if(!file.size(path_promt_pre)==0){
        data_promt <- read.delim(path_promt_pre)
        res <- TwoSampleMR::mr(data_promt,method_list = c('mr_ivw_fe','mr_ivw_mre','mr_egger_regression','mr_weighted_median','mr_weighted_mode','mr_simple_mode'))
      }else{
        message('No ivs were mapped')
      }  
      
      
      if(!file.exists(output_path)){
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
      }else{
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
      }
    }
    
  }
}
MR_analysis_immToinfFun(File_name)

library(snow)
c1 <- makeCluster(8,type = 'SOCK')
clusterEvalQ(c1,library(snowfall))
clusterApply(c1,File_name,MR_analysis_immToinfFun)

## MR heterogeneity___________________________________________________________________________________________________________
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_HeterogeneityimmToinf")

MR_heterogeneity_imm_to_inf <- function(x = File_name){
  for(i in x){
    print(i)
    File <- limma::strsplit2(i,'\\.')[,1]
    path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/immuneToinf_harmnoise/",File)
    immid <- paste0('ebi-a-GCST9000',c(1391:2121))
    output_name <- paste0(File,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_HeterogeneityimmToinf/",output_name)
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.exposure[!duplicated( output_promt$id.exposure)]
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j,'.txt')
      print(j)
      
      if(file.exists(path_promt_pre)){
        data_promt <- read.delim(path_promt_pre)
        if(dim(data_promt)[1] >1){
          res <- TwoSampleMR::mr_heterogeneity(data_promt)
        }else{
          message(paste0(i,': Not enough SNPs'))
        }
        
      }  
      
      if(j == immid[1]){
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
      }else{
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
      }
    }
    
  }
}
MR_heterogeneity_imm_to_inf()

library(snow)
c1 <- makeCluster(8,type = 'SOCK')
clusterEvalQ(c1,library(snowfall))
clusterApply(c1,File_name,MR_heterogeneity_imm_to_inf)

## MR pleiotropy analysis

dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_pleiotrpyimmToinf")
MR_pleiotropy_immune_to_inf <- function(x = File_name){
  for(i in x){
    print(i)
    File <- limma::strsplit2(i,'\\.')[,1]
    path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/immuneToinf_harmnoise/",File)
    immid <- paste0('ebi-a-GCST9000',c(1391:2121))
    output_name <- paste0(i,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_pleiotrpyimmToinf/",output_name)
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.exposure[!duplicated( output_promt$id.exposure)]
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j,'.txt')
      print(j)
      
      if(file.exists(path_promt_pre)){
        data_promt <- read.delim(path_promt_pre)
        if(dim(data_promt)[1] >1){
          res <- TwoSampleMR::mr_pleiotropy_test(data_promt)
        }else{
          message(j,': not enoguh snps')
        }
        
      }else{
        message(j,': do not exist')
      }  
      
      if(j == immid[1]){
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
      }else{
        write.table(res,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
      }
    }
    
  }
}
MR_pleiotropy_immune_to_inf()
clusterApply(c1,File_name,MR_pleiotropy_immune_to_inf)

## MR presso analysis
dir.create("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_pressoimmToinf")

MR_presso_immune_to_inf <- function(x = File_name){
  for(i in x){
    print(i)
    File <- limma::strsplit2(i,'\\.')[,1]
    path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_harmnosie/immuneToinf_harmnoise/",File)
    immid <- paste0('ebi-a-GCST9000',c(1391:2121))
    output_name <- paste0(i,'_','To','imm.txt')
    output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/MR_pressoimmToinf/",output_name)
    
    
    if(file.exists(output_path)){
      output_promt <- read.delim(output_path)
      MR_done <- output_promt$id.exposure[!duplicated(output_promt$id.exposure)]
      MR_undo <- setdiff(immid,MR_done)
    }else{
      MR_undo <- immid
    }
    
    for(j in MR_undo){
      path_promt_pre <- paste0(path_promt,'/',j,'.txt')
      print(j)
      if(!file.exists(path_promt_pre)){message(paste0(j,': file do not exist'))}
      
      if(file.exists(path_promt_pre)){
        data_promt <- read.delim(path_promt_pre)
        if(dim(data_promt)[1] >= 1000){
          message('Too large to compute')
        }
        
        if(dim(data_promt)[1] >=4 & dim(data_promt)[1] < 1000){
          het_adb_ada <- MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome',BetaExposure = "beta.exposure",SdOutcome = 'se.outcome',SdExposure = "se.exposure",
                                             OUTLIERtest = TRUE,DISTORTIONtest = TRUE,data = data_promt,NbDistribution = 1000,SignifThreshold = 0.05)
          result_promt <- cbind(id.exposure = data_promt$id.exposure[2],het_adb_ada$`Main MR results`)
          result_promt$global_pval <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$Pvalue,'NA')
          result_promt$global_RSSobs <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$RSSobs,'NA')
          result_promt$id.outcome <- rep(j,dim(result_promt)[1])
          if(j == immid[1]){
            write.table(result_promt,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = TRUE)
          }else{
            write.table(result_promt,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = FALSE)
          }
        }else{
          message('Number of SNP < 4')
        }
        
      }  
      
    }
    
  }
}

MR_presso_immune_to_inf(File_name)

clusterApply(c1,File_name[29:33],MR_presso_immune_to_inf)



####Circulating inflammtory cytokines to Lung cancer
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung")
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung")
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/IVs_poor")

Lung_ID <- c('Lung cancer','Lung adenocarcinoma','small cell lung carcinoma','squamous cell lung cancer',
             'lung cancer father','lung cancer mother','lung cancer any parental',
             'lung cancer sibilings','lung cancer smoker','lung cancer no smoker','no small lung cancer')
access_ID <- c('ukb-a-54','ieu-a-984','ieu-a-988','ieu-a-989','ukb-b-14521','ukb-b-20176','ebi-a-GCST90013972',
               'ukb-b-15826','ebi-a-GCST004749','ebi-a-GCST004747','finn-b-C3_LUNG_NONSMALL')

outcome_inf_LungIVpoor <- function(x=c(1:11)){
  Lung_ID <- c('Lung cancer','Lung adenocarcinoma','small cell lung carcinoma','squamous cell lung cancer',
               'lung cancer father','lung cancer mother','lung cancer any parental',
               'lung cancer sibilings','lung cancer smoker','lung cancer no smoker','no small lung cancer')
  access_ID <- c('ukb-a-54','ieu-a-984','ieu-a-988','ieu-a-989','ukb-b-14521','ukb-b-20176','ebi-a-GCST90013972',
                 'ukb-b-15826','ebi-a-GCST004749','ebi-a-GCST004747','finn-b-C3_LUNG_NONSMALL')
  
  exposure_dir <- "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure"
  inf_name <- list.files(exposure_dir)
  SNPs <- c()
  for(i in inf_name){
    path_promt <- paste0(exposure_dir,'/',i)
    data_promt <- read.delim(path_promt)
    if(length(SNPs)==0){
      SNPs <- data_promt$SNP
    }else{
      SNPs <- c(SNPs,data_promt$SNP)
    }
  }
  
  SNPs <- SNPs[!duplicated(SNPs)]
  if(x= c(1:10)){
    for(j in x){
      data_outcome <- TwoSampleMR::extract_outcome_data(snps = SNPs,outcomes = access_ID[j])
      outcome_path_promt <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/IVs_poor/",Lung_ID[j],'.txt')
      write.table(data_outcome,file=outcome_path_promt,row.names = F,sep = '\t',quote = F)
    }
  }
  if(x= 11){
    data_orignial <- data.table::fread('~/Lung cancer/finngen_R10_C3_LUNG_NONSMALL_EXALLC.gz',header = TRUE)
    data_orignial <- TwoSampleMR::format_data(data_orignial,
                                         type = 'outcome',
                                         snp_col = 'rsids',
                                         phenotype_col = 'phenotypes',
                                         beta_col = 'beta',
                                         se_col = 'sebeta',
                                         eaf_col = 'af_alt',
                                         effect_allele_col = 'alt',
                                         other_allele_col = 'ref',
                                         pval_col = 'pval')
    data_outcome <- data_orignial[data_orignial$SNP %in% SNPs,]
    outcome_path_promt <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/IVs_poor/",Lung_ID[11],'.txt')
    write.table(data_outcome,file=outcome_path_promt,row.names = F,sep = '\t',quote = F)
  }
  
    
  }
outcome_inf_LungIVpoor()

dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/outcome")
MapIVs_inf <- function(x = Lung_ID){
  exposure_dir <- "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure"
  inf_name <- list.files(exposure_dir)
  for(i in x){
    dir_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/outcome/",i)
    dir.create(dir_path)
    snp_poor_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/IVs_poor/",i,'.txt')
    snp_poor <- read.delim(snp_poor_path)
    print(i)
    for(j in inf_name){
      path_promt <- paste0(exposure_dir,'/',j)
      data_promt <- read.delim(path_promt)
      data_outcome <- snp_poor[snp_poor$SNP %in% data_promt$SNP,]
      outcome_path <- paste0(dir_path,'/',j)
      write.table(data_outcome,file=outcome_path,row.names = F,sep = '\t',quote = F)
      print(j)
    }
  }
}
MapIVs_inf()

dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/Harmonise_data")
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/Harmonise_data/Harmonise_inf_Lung")

Harmonise_inf_Lung <- function(x = Lung_ID){
  exposure_dir <- "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure"
  inf_name <- list.files(exposure_dir)
  for(i in x){
    dir_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/Harmonise_data/Harmonise_inf_Lung/",i)
    dir.create(dir_path)
    
    outcome_path_pre <- paste0('/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/outcome/',i,'/')
    print(i)
    
    for(j in inf_name){
      exposure_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure/",j)
      exposure_data <- read.delim(exposure_path)
      
      outcome_path <- paste0(outcome_path_pre,j)
      outcome_data <- read.delim(outcome_path)
      
      output_path <- paste0(dir_path,'/',j)
      data_promt <- TwoSampleMR::harmonise_data(exposure_dat = exposure_data,outcome_dat = outcome_data)
      write.table(data_promt,file=output_path,row.names = F,sep = '\t',quote = F)
      print(j)
    }
  }
}
Harmonise_inf_Lung()


Add_samplesizeFun <- function(x = Lung_ID){
  for(i in x){
    path_target <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/Harmonise_data/Harmonise_inf_Lung/",i)
    inf_name <- list.files(path_target)
    for(j in inf_name){
      path_promt <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inflammation cytosines/",j)
      Sample_size <- readLines(path_promt,n = 3)[2]
      Sample_size <- limma::strsplit2(Sample_size,'\t')[,16]
      Sample_size <- as.numeric(Sample_size)
      print(j)
      path_target_promt <- paste0(path_target,'/',j)
      data_promt <- read.delim(path_target_promt)
      data_promt <- cbind(data_promt,samplesize.exposure = rep(Sample_size,dim(data_promt)[1]))
      write.table(data_promt,file = path_target_promt,row.names = F,sep = '\t',quote = F)
    }
  }
  
}
Add_samplesizeFun(Lung_ID)

calculate_F_lager_than10 <- function(x = Lung_ID){# dirName
  for(i in x){
    print(i)
    file_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/Harmonise_data/Harmonise_inf_Lung/",i)
    File <- list.files(file_path)
    for(j in File){
      path_promt <- paste0(file_path,'/',j)
      if(!file.size(path_promt) ==0){
        data_promt <- read.delim(path_promt)
        R_squre <- TwoSampleMR::get_r_from_bsen(data_promt$beta.exposure,data_promt$se.exposure,data_promt$samplesize.exposure)
        F_value <- as.numeric(R_squre)*(as.numeric(data_promt$samplesize.exposure) - 2)/(1 - as.numeric(R_squre))
        data_promt <- cbind(data_promt, R2 = as.numeric(R_squre), F_value = F_value)
        data_finall <- data_promt[as.numeric(data_promt$F_value) > 10,]
        write.table(data_finall,file = path_promt,row.names = F,sep = '\t',quote = F)
        print(j)
      }else{
        message('No IVs were mapped')
      }
      
    }
  }
}
calculate_F_lager_than10()


# MR analysis
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/result")
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_resultinfToLung")

MR_anlaysis_inf_Lung <- function(x = Lung_ID){
  exposure_dir <- "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure"
  inf_name <- list.files(exposure_dir)
  for(i in x){
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_resultinfToLung/",i,'.txt')
    data_path_pre <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/Harmonise_data/Harmonise_inf_Lung/",i)
    print(i)
    
    result_promt <- data.frame()
    for(j in inf_name){
      
      data_path <- paste0(data_path_pre,'/',j)
      data_promt <- read.delim(data_path)
      res <- TwoSampleMR::mr(data_promt,method_list = c('mr_ivw_fe','mr_ivw_mre','mr_egger_regression','mr_weighted_median','mr_weighted_mode','mr_simple_mode'))
      inf <- limma::strsplit2(j,'\\.')[,1]
      result <- cbind(exposure_inf = rep(inf,dim(res)[1]),res)
      if(dim(result_promt)[1]==0){
        result_promt <- result
      }else{
        result_promt <- rbind(result_promt,result)
      }
      
      
      write.table(result_promt,file=output_path,row.names = F,sep = '\t',quote = F)
      print(j)
    }
  }
}
MR_anlaysis_inf_Lung(Lung_ID)

library(snow)
c1 <- makeCluster(8,type = 'SOCK')
clusterEvalQ(c1,library(snowfall))
clusterApply(c1,Lung_ID[2:11],MR_anlaysis_inf_Lung)

# heterogenity
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_HeteroinfToLung")
MR_heterogenity_inf_Lung <- function(x = Lung_ID){
  exposure_dir <- "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure"
  inf_name <- list.files(exposure_dir)
  for(i in x){
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_HeteroinfToLung/",i,'.txt')
    data_path_pre <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/Harmonise_data/Harmonise_inf_Lung/",i)
    print(i)
    
    result_promt <- data.frame()
    for(j in inf_name){
      
      data_path <- paste0(data_path_pre,'/',j)
      data_promt <- read.delim(data_path)
      res <- TwoSampleMR::mr_heterogeneity(data_promt)
      inf <- limma::strsplit2(j,'\\.')[,1]
      result <- cbind(exposure_inf = rep(inf,dim(res)[1]),res)
      if(dim(result_promt)[1]==0){
        result_promt <- result
      }else{
        result_promt <- rbind(result_promt,result)
      }
      
      
      write.table(result_promt,file=output_path,row.names = F,sep = '\t',quote = F)
      print(j)
    }
  }
}
MR_heterogenity_inf_Lung()

# MR pleiotropy analysis
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_pleiotropyinfToLung")
MR_pleiotropy_inf_Lung <- function(x = Lung_ID){
  exposure_dir <- "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure"
  inf_name <- list.files(exposure_dir)
  for(i in x){
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_pleiotropyinfToLung/",i,'.txt')
    data_path_pre <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/Harmonise_data/Harmonise_inf_Lung/",i)
    print(i)
    
    result_promt <- data.frame()
    for(j in inf_name){
      
      data_path <- paste0(data_path_pre,'/',j)
      data_promt <- read.delim(data_path)
      res <- TwoSampleMR::mr_pleiotropy_test(data_promt)
      inf <- limma::strsplit2(j,'\\.')[,1]
      result <- cbind(exposure_inf = rep(inf,dim(res)[1]),res)
      if(dim(result_promt)[1]==0){
        result_promt <- result
      }else{
        result_promt <- rbind(result_promt,result)
      }
      
      
      write.table(result_promt,file=output_path,row.names = F,sep = '\t',quote = F)
      print(j)
    }
  }
}
MR_pleiotropy_inf_Lung()

# MR presso
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_pressoinfToLung")
MR_presso_inf_Lung <- function(x = Lung_ID){
  exposure_dir <- "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inf_exposure"
  inf_name <- list.files(exposure_dir)
  for(i in x){
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_pressoinfToLung/",i,'.txt')
    data_path_pre <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_inf_to_Lung/Harmonise_data/Harmonise_inf_Lung/",i)
    print(i)
    
    if(file.exists(outcome_path)){
      result_down <- read.delim(output_path)
      inf_down <- result_down[,1][!duplicated(result_down[,1])]
      inf_down <- paste0(inf_down,'.txt')
      inf_undo <- setdiff(inf_name,inf_down)
    }else{
      inf_undo <- inf_name
    }
    result_promt <- data.frame()
    for(j in inf_undo){
      data_path <- paste0(data_path_pre,'/',j)
      data_promt <- read.delim(data_path)
      
      if(dim(data_promt)[1] >=4 & dim(data_promt)[1] < 1000){
        het_adb_ada <- MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome',BetaExposure = "beta.exposure",SdOutcome = 'se.outcome',SdExposure = "se.exposure",
                                           OUTLIERtest = TRUE,DISTORTIONtest = TRUE,data = data_promt,NbDistribution = 1000,SignifThreshold = 0.05)
        result <- cbind(id.exposure = data_promt$id.exposure[2],het_adb_ada$`Main MR results`)
        result$global_pval <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$Pvalue,'NA')
        result$global_RSSobs <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$RSSobs,'NA')
        inf <- limma::strsplit2(j,'\\.')[,1]
        result <- cbind(exposure_inf = rep(inf,dim(result)[1]),result)
        if(dim(result_promt)[1]==0){
        result_promt <- result
      }else{
        result_promt <- rbind(result_promt,result)
      }
      
      
      write.table(result_promt,file=output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
      print(j)
      }else{
      message('nSNPs out of the range')
    }
}
  }
}
MR_presso_inf_Lung(Lung_ID[1])
clusterApply(c1,Lung_ID[3:11],MR_presso_inf_Lung)


####Lung cancer to Circulating inflammtory cytokines

dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_Lung_to_inf")
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_Lung_to_inf/IVs_poor")
inf_name <- list.files("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inflammation cytosines")
getIV_poor <- function(x = inf_name){
  dir <- list.dirs("/Users/llls2012163.com/Lung cancer/CancerTypeExposure")
  SNP <- c()
  for(i in dir[2:12]){
    Files <- list.files(i)
    path <- paste0(i,'/',Files)
    data_promt <- read.delim(path)
    iv <- data_promt$SNP
    if(length(SNP)==0){
        SNP = iv
      }else{
        SNP = c(SNP,iv)
      }
    SNP <- SNP[!duplicated(SNP)]
  }
    for(j in x){
      path_outcome <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/inflammation cytosines/",j)
      print(j)
      data_promt <- read.delim(path_outcome)
      data_promt <- data_promt[data_promt$MarkerName %in% SNP,]
      output_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_Lung_to_inf/IVs_poor/",j,'.txt')
      write.table(data_promt,file=output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
    }
}
getIV_poor()

dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_Lung_to_inf/outcome")

MapIVs_Lung_inf <- function(x = inf_name){
  exposure_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/CancerTypeExposure")[2:12]
  for(i in exposure_dir){
    File <- list.files(i)
    exposure_path <- paste0(i,'/',File)
    exposure_data <- read.delim(exposure_path)
    Lung_name <- limma::strsplit2(i,'/')[,6]
    output_path_pre <- paste0("~/Lung cancer/inf_Lung/outcome_Lung_to_inf/outcome/",Lung_name)
    dir.create(output_path_pre)
    for(j in x){
      IVs_poor_path <- paste0("~/Lung cancer/inf_Lung/outcome_Lung_to_inf/IVs_poor/",j,'.txt')
      IVs_poor <- read.delim(IVs_poor_path)
      colnames(IVs_poor) <- c('SNP',"effect_allele.outcome","other_allele.outcome","eaf.outcome",'seEaf.outcome',
                              'MIAF.outcome','MAF.outcome','beta.outcome','se.outcome','pval_origin.outcome',
                              'direction.outcome','het.outcome','hetchi.outcome','hetdf.outcome','hetpvalue.outcome',
                              'samplesize.outcome')
      output_path <- paste0(output_path_pre,'/',j)
      data_promt <- IVs_poor[IVs_poor$SNP %in% exposure_data$SNP,]
      write.table(data_promt,file=output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
      
    }
  }
}
MapIVs_Lung_inf()  

# harmnoise data
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_Lung_to_inf/Harmnoise_data")

Harmonise_Lung_inf <- function(x = inf_name){
  exposure_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/CancerTypeExposure")[2:12]
  for(i in exposure_dir){
    File <- list.files(i)
    exposure_path <- paste0(i,'/',File)
    exposure_data <- read.delim(exposure_path)
    Lung_name <- limma::strsplit2(i,'/')[,6]
    outcome_path_pre <- paste0("~/Lung cancer/inf_Lung/outcome_Lung_to_inf/outcome/",Lung_name)
    output_path_pre <- paste0("~/Lung cancer/inf_Lung/outcome_Lung_to_inf/Harmnoise_data/",Lung_name)
    dir.create(output_path_pre)
    for(j in x){
      outcome_path <- paste0(outcome_path_pre,'/',j)
      outcome_data <- read.delim(outcome_path)
      inf <- limma::strsplit2(j,'\\.')[1]
      outcome_data$id.outcome <- rep(inf,dim(outcome_data)[1])
      outcome_data$outcome <- rep(inf,dim(outcome_data)[1])
      data_promt <- TwoSampleMR::harmonise_data(exposure_dat = exposure_data,outcome_dat = outcome_data)
      output_path <- paste0(output_path_pre,'/',j)
      write.table(data_promt,file=output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
      
    }
  }
}
Harmonise_Lung_inf()

# Calculate F value
Lung_dirs <- list.dirs("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_Lung_to_inf/Harmnoise_data")[2:12]
Lung_dirs <- limma::strsplit2(Lung_dirs,'/')[,8]

calculate_F_lager_than10 <- function(x = Lung_dirs){
  for(i in x){
    print(i)
    file_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_Lung_to_inf/Harmnoise_data/",i)
    File <- list.files(file_path)
    if(i == Lung_dirs[9]){
      for(j in File){
        path_promt <- paste0(file_path,'/',j)
        if(!file.size(path_promt) ==0){
          data_promt <- read.delim(path_promt)
          data_promt$samplesize.exposure <- rep(218792,dim(data_promt)[1])
          R_squre <- TwoSampleMR::get_r_from_bsen(data_promt$beta.exposure,data_promt$se.exposure,data_promt$samplesize.exposure)
          F_value <- as.numeric(R_squre)*(as.numeric(data_promt$samplesize.exposure) - 2)/(1 - as.numeric(R_squre))
          data_promt <- cbind(data_promt, R2 = as.numeric(R_squre), F_value = F_value)
          data_finall <- data_promt[as.numeric(data_promt$F_value) > 10,]
          write.table(data_finall,file = path_promt,row.names = F,sep = '\t',quote = F)
          print(j)
        }else{
          message('No IVs were mapped')
        }
        
      }
    }else{
      for(j in File){
        path_promt <- paste0(file_path,'/',j)
        if(!file.size(path_promt) ==0){
          data_promt <- read.delim(path_promt)
          R_squre <- TwoSampleMR::get_r_from_bsen(data_promt$beta.exposure,data_promt$se.exposure,data_promt$samplesize.exposure)
          F_value <- as.numeric(R_squre)*(as.numeric(data_promt$samplesize.exposure) - 2)/(1 - as.numeric(R_squre))
          data_promt <- cbind(data_promt, R2 = as.numeric(R_squre), F_value = F_value)
          data_finall <- data_promt[as.numeric(data_promt$F_value) > 10,]
          write.table(data_finall,file = path_promt,row.names = F,sep = '\t',quote = F)
          print(j)
        }else{
          message('No IVs were mapped')
        }
        
      }
    }
    
  }
}
calculate_F_lager_than10(Lung_dirs)


# MR analysis
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_resultLungToinf")
MR_analysis_LungToinf <- function(x = Lung_dirs){
  for(i in x){
    print(i)
    file_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_Lung_to_inf/Harmnoise_data/",i)
    File <- list.files(file_path)
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_resultLungToinf/",i,'.txt')
    for(j in File){
        path_promt <- paste0(file_path,'/',j)
        if(!file.size(path_promt) ==0){
          data_promt <- read.delim(path_promt)
          result <- TwoSampleMR::mr(data_promt,method_list = c('mr_ivw_fe','mr_ivw_mre','mr_egger_regression','mr_weighted_median','mr_weighted_mode','mr_simple_mode'))
          if(!file.exists(output_path)){
            write.table(result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
          }else{
            write.table(result,file = output_path,row.names = F,sep = '\t',quote = F,col.names = F,append = TRUE)
          }
            print(j)
        }else{
          message('No IVs were mapped')
        }
        
      }
    }
}
MR_analysis_LungToinf(Lung_dirs[1])

library(snow)
c1 <- makeCluster(8,type = 'SOCK')
clusterEvalQ(c1,library(snowfall))
clusterApply(c1,Lung_dirs[2:11],MR_analysis_LungToinf)

# heterogeneity
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_HeteroLungToinf")
MR_heterogenity_Lung_inf <- function(x = Lung_dirs){
  for(i in x){
    print(i)
    file_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_Lung_to_inf/Harmnoise_data/",i)
    File <- list.files(file_path)
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_HeteroLungToinf/",i,'.txt')
    for(j in File){
      path_promt <- paste0(file_path,'/',j)
      if(!file.size(path_promt) ==0){
        data_promt <- read.delim(path_promt)
        result <- TwoSampleMR::mr_heterogeneity(data_promt)
        if(!file.exists(output_path)){
          write.table(result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
        }else{
          write.table(result,file = output_path,row.names = F,sep = '\t',quote = F,col.names = F,append = TRUE)
        }
        print(j)
      }else{
        message('No IVs were mapped')
      }
      
    }
  }
}
MR_heterogenity_Lung_inf()

# MR pleiotropy analysis
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_pleiotropyLungToinf")
MR_pleiotropy_Lung_inf <- function(x = Lung_dirs){
  for(i in x){
    print(i)
    file_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_Lung_to_inf/Harmnoise_data/",i)
    File <- list.files(file_path)
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_pleiotropyLungToinf/",i,'.txt')
    for(j in File){
      path_promt <- paste0(file_path,'/',j)
      if(!file.size(path_promt) ==0){
        data_promt <- read.delim(path_promt)
        result <- TwoSampleMR::mr_pleiotropy_test(data_promt)
        if(!file.exists(output_path)){
          write.table(result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
        }else{
          write.table(result,file = output_path,row.names = F,sep = '\t',quote = F,col.names = F,append = TRUE)
        }
        print(j)
      }else{
        message('No IVs were mapped')
      }
      
    }
  }
}
MR_pleiotropy_Lung_inf()

# MR presso
dir.create("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_pressoLungToinf")

MR_presso_Lung_inf <- function(x = Lung_dirs){
  for(i in x){
    print(i)
    file_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/outcome_Lung_to_inf/Harmnoise_data/",i)
    File <- list.files(file_path)
    output_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/result/MR_pressoLungToinf/",i,'.txt')
    for(j in File){
      path_promt <- paste0(file_path,'/',j)
        data_promt <- read.delim(path_promt)
        print(j)
        inf <- limma::strsplit2(j,'\\.')[,1]
        if(dim(data_promt)[1] >=4 & dim(data_promt)[1] < 1000){
          het_adb_ada <- MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome',BetaExposure = "beta.exposure",SdOutcome = 'se.outcome',SdExposure = "se.exposure",
                                             OUTLIERtest = TRUE,DISTORTIONtest = TRUE,data = data_promt,NbDistribution = 1000,SignifThreshold = 0.05)
          result <- cbind(id.exposure = data_promt$id.exposure[2],het_adb_ada$`Main MR results`)
          result$global_pval <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$Pvalue,'NA')
          result$global_RSSobs <- c(het_adb_ada$`MR-PRESSO results`$`Global Test`$RSSobs,'NA')
          result$outcome <- rep(inf,2)
          if(!file.exists(output_path)){
            write.table(result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
          }else{
            write.table(result,file = output_path,row.names = F,sep = '\t',quote = F,col.names = F,append = TRUE)
          }
          
        }else{
          message('nSNPs out of the range')
      
    }
    }
  }
}
MR_presso_Lung_inf(Lung_dirs[1])
clusterApply(c1,Lung_dirs,MR_presso_Lung_inf)



###Data summary

####1. 731 immune cell traits--------->Lung cancer
Result_storge_dir <- "/Users/llls2012163.com/Lung cancer/result_file_immTocancerType"
MR_result_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType")[5]
MR_het_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType")[2]
MR_file <- list.files(MR_result_dir)

FDR_adjust <- function(x = MR_file){
  for(i in x){
    print(i)
    path_prmot <- paste0(list.dirs("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType")[5],'/',i)
    data_prmot <- read.delim(path_prmot)
    data_prmot <- na.omit(data_prmot)
    promt <- data.frame()
    for(j in data_prmot$method[!duplicated(data_prmot$method)]){
      print(j)
      promt1 <- data_prmot[data_prmot$method %in% j,]
      promt1$FDR <- p.adjust(promt1$pval,method = 'BH',length(promt1$pval))
      
      if(dim(promt)[1]==0){
        promt <- promt1
      }else{
        promt <- rbind(promt,promt1)
      }
    }
  }
  is <- as.vector(as.numeric(rownames(promt)))
  promt <- cbind(is = is,promt)
  promt <- promt[order(is),]
  return(promt)
}


res <- FDR_adjust(MR_file[1])

MR_result_sig_select <- function(x = MR_file){
  MR_sig <- data.frame()
  for(i in x){
    print(i)
    path_prmot <- paste0(list.dirs("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType")[5],'/',i)
    data_prmot <- read.delim(path_prmot)
    
    #data_prmot <- FDR_adjust(i)
    
    MR_result_sig_fix <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (fixed effects)',]
    MR_result_sig_fix <- MR_result_sig_fix[MR_result_sig_fix$pval < 5e-02,]
    MR_result_sig_random <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (multiplicative random effects)',]
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$pval < 5e-02,]
    
    path_prmot_het <- paste0("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType/MR_heterogeneity/",i)
    data_prmot_het <- read.delim(path_prmot_het)
    data_prmot_het <- data_prmot_het[data_prmot_het$method == 'Inverse variance weighted',]
    data_prmot_het <- data_prmot_het[data_prmot_het$Q_pval < 0.05,]
    
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$id.exposure %in% data_prmot_het$id.exposure & MR_result_sig_random$id.outcome %in% data_prmot_het$id.outcome,]
    MR_result_sig_fix <- MR_result_sig_fix[! MR_result_sig_fix$id.exposure %in% data_prmot_het$id.exposure & MR_result_sig_fix$id.outcome %in% data_prmot_het$id.outcome,]
    
    MR_sig_IVW <- rbind(MR_result_sig_fix,MR_result_sig_random)
    
    MR_result_sing <- data_prmot[data_prmot$id.exposure %in% MR_sig_IVW$id.exposure & data_prmot$id.outcome %in% MR_sig_IVW$id.outcome,]
    
    MR_result_sing1 <- MR_result_sing[MR_result_sing$id.exposure %in% data_prmot_het$id.exposure & MR_result_sing$id.outcome %in% data_prmot_het$id.outcome & MR_result_sing$method !='Inverse variance weighted (fixed effects)',]
    MR_result_sing2 <- MR_result_sing[!MR_result_sing$id.exposure %in% data_prmot_het$id.exposure & MR_result_sing$id.outcome %in% data_prmot_het$id.outcome & MR_result_sing$method !='Inverse variance weighted (multiplicative random effects)',]
    
    MR_result <- rbind(MR_result_sing1,MR_result_sing2)
    
    if(dim(MR_sig)[1]==0){
      MR_sig <- MR_result
    }else{
      MR_sig <- rbind(MR_sig,MR_result)
    }
    
  }
  return(MR_sig)
}
res <- MR_result_sig_select(MR_file)

write.table(res,file = "/Users/llls2012163.com/Lung cancer/result_file_immTocancerType/MR_sig.txt",row.names = F,sep = '\t',quote = F,append = FALSE,col.names = TRUE)



# 1455 has the heterogenity
MR_pleiotropy_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType")[3]
MR_file <- list.files(MR_pleiotropy_dir)
MR_pleiotropy_phe <- function(x = MR_file){
  ple_phe <- data.frame()
  No_ple_phe <- data.frame()
  for(i in x){
    print(i)
    ple_path <- paste0("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType/MR_pleiotropy/",i)
    data_promt <- read.delim(ple_path)
    data_ple <- data_promt[data_promt$pval < 0.05,]
    data_no_ple <- data_promt[data_promt$pval >=0.05,]
    
    if(dim(ple_phe)[1]==0){
      ple_phe <- data_ple
    }else{
      ple_phe <- rbind(ple_phe,data_ple)
    }
  }
  
  ple_phe <- na.omit(ple_phe)
  if(dim(No_ple_phe)[1]==0){
    No_ple_phe <- data_no_ple
  }else{
    No_ple_phe <- rbind(No_ple_phe,data_no_ple)
  }
  
  ple_distinc <- list(ple = ple_phe,unple = No_ple_phe)
  
  return(ple_distinc)
}
MR_ple <- MR_pleiotropy_phe()

# 479 has pleiotropy
# MR_presso
MR_presso_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType")[4]
MR_file <- list.files(MR_presso_dir)
MR_presso_pleiotropy_phe <- function(x = MR_file){
  id_data_frame <- data.frame(id.outcome = c('ebi-a-GCST004747','ebi-a-GCST004749','ukb-a-54','ieu-a-984','ebi-a-GCST900013972',
                                        'ukb-b-14521','ukb-b-20176','ukb-b-15826','finn-b-C3_LUNG_NOSMALL_EXALLC','ieu-a-988',
                                        'ieu-a-989'),MR_file = x)
  ple_phe <- data.frame()
  No_ple_phe <- data.frame()
  for(i in c(1:11)){
    print(id_data_frame$MR_file[i])
    presso_path <- paste0(MR_presso_dir,'/',id_data_frame$MR_file[i])
    presso_data <- read.delim(presso_path)
    presso_data <- cbind(presso_data,id.outcome = rep(id_data_frame$id.outcome[i],dim(presso_data)[1]))
    presso_ple <- presso_data[presso_data$global_pval < 0.05,]
    presso_ple <- na.omit(presso_ple)
    presso_ple <- presso_data[presso_data$id.exposure %in% presso_ple$id.exposure,]
    presso_no_ple <- presso_data[presso_data$global_pval >= 0.05,]
    if(dim(ple_phe)[1]==0){
      ple_phe <- presso_ple
    }else{
      ple_phe <- rbind(ple_phe,presso_ple)
    }
    
    if(dim(No_ple_phe)[1]==0){
      No_ple_phe <- presso_no_ple
    }else{
      No_ple_phe <- rbind(No_ple_phe,presso_no_ple)
    }
    
    
  }
  
  Presso_ple <- list(ple_phe = ple_phe,No_ple_phe=No_ple_phe)
  return(Presso_ple)
}
Presso_ple <- MR_presso_pleiotropy_phe()

MR_sig <- read.delim("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType/MR_sig.txt")
MR_ple_ple <- MR_ple$ple[MR_ple$ple$id.exposure %in% MR_sig$id.exposure & MR_ple$ple$id.outcome %in% MR_sig$id.outcome,]
MR_presso_ple <- Presso_ple$ple_phe[Presso_ple$ple_phe$id.exposure %in% MR_sig$id.exposure & Presso_ple$ple_phe$id.outcome %in% MR_sig$id.outcome,]
MR_presso_no_ple <- Presso_ple$No_ple_phe[Presso_ple$No_ple_phe$id.exposure %in% MR_sig$id.exposure & Presso_ple$No_ple_phe$id.outcome %in% MR_sig$id.outcome,]

#Stableest result
MR_stable_result <- MR_sig[MR_sig$id.exposure %in% MR_ple$unple$id.exposure & MR_sig$id.outcome %in% MR_ple$ple$id.outcome,]
MR_stable_result <- MR_stable_result[MR_stable_result$id.exposure %in% Presso_ple$No_ple_phe$id.exposure & MR_stable_result$id.outcome %in% Presso_ple$No_ple_phe$id.outcome,]
MR_stable_result <- na.omit(MR_stable_result)

load("/Users/llls2012163.com/GWAS/TraitID.Rdata")
ID_to_immTrait <- ID_Trait$ImmID
res <- data.frame()
for(i in 1:731){
  data_promt <- MR_stable_result[MR_stable_result$id.exposure == ID_to_immTrait$id[i],]
  data_promt <- cbind(Trait = rep(ID_to_immTrait$X[i],dim(data_promt)[1]),data_promt)
  if(dim(res)[1]==0){
    res <- data_promt
  }else{
    res <- rbind(res,data_promt)
  }
}

res <- res[order(res$id.outcome),]


res <- cbind(res,OR=exp(res$b),'Lower 95% CI'=exp(res$b - 1.96*res$se),'Upper 95% CI' = exp(res$b + 1.96*res$se))
write.table(res,file = "/Users/llls2012163.com/Lung cancer/result_file_immTocancerType/MR_sig_stable.txt",sep = '\t',quote = FALSE,row.names = FALSE)
rm(list = ls())

#data visualization
MR_sig_stable <- read.delim("~/Lung cancer/result_file_immTocancerType/MR_sig_stable.txt")
MR_sig_stable$outcome <- limma::strsplit2(MR_sig_stable$outcome," \\|\\| ")[,1]
MR_sig_stable <- MR_sig_stable[,c(-3,-5)]

load("/Users/llls2012163.com/GWAS/TraitID.RData")
id_trait <- ID_Trait$ImmID

traitpanel <- read_excel("Desktop/traitpanel.xlsx")
q <- data.frame()
for(i in 1:dim(MR_sig_stable)[1]){
  trait <- MR_sig_stable$id.exposure[i]
  ids <- limma::strsplit2(trait,'-')[,3]
  a <- traitpanel[traitpanel$`GWAS Catalog Accession Number`==ids,]
  if(length(q)==0){
    q <- a
  }else{
    q <- rbind(q,a)
  }
}
MR_sig_stable <- cbind(q,MR_sig_stable)
write.table(MR_sig_stable,file = "/Users/llls2012163.com/add4.txt",sep = '\t',quote = FALSE,row.names = FALSE)

MR_sig_stable1 <- MR_sig_stable[MR_sig_stable$method %in% c('Inverse variance weighted (fixed effects)','Inverse variance weighted (multiplicative random effects)'),]
colnames(MR_sig_stable1)[1] <- 'TraitType'

library(webr)
library(tidyverse)
p1<-  PieDonut(MR_sig_stable1,aes(outcome, Panel),
               r0=0.2,
               r1=0.6,
               r2=0.9,
               pieLabelSize = 2.5,
               donutLabelSize = 2.5,labelposition = 0)


p2 <- PieDonut(MR_sig_stable1,aes(outcome, TraitType),
               r0=0.2,
               r1=0.6,
               r2=0.9,
               pieLabelSize = 2.5,
               donutLabelSize = 2.5,labelposition = 0)

library(ggplot2)

p3 <- ggplot(MR_sig_stable1,aes(OR,-log10(pval)))+geom_point()



####2. Lung cancer--------------------->731 immunce cell traits
Result_storge_dir <- "/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm"

MR_result_dir <- list.dirs(Result_storge_dir)[5]
MR_het_dir <- list.dirs(Result_storge_dir)[2]
MR_file <- list.files(MR_result_dir)

FDR_adjust <- function(x = MR_file){
  for(i in x){
    print(i)
    path_prmot <- paste0(list.dirs("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType")[5],'/',i)
    data_prmot <- read.delim(path_prmot)
    data_prmot <- na.omit(data_prmot)
    promt <- data.frame()
    for(j in data_prmot$method[!duplicated(data_prmot$method)]){
      print(j)
      promt1 <- data_prmot[data_prmot$method %in% j,]
      promt1$FDR <- p.adjust(promt1$pval,method = 'BH',length(promt1$pval))
      
      if(dim(promt)[1]==0){
        promt <- promt1
      }else{
        promt <- rbind(promt,promt1)
      }
    }
  }
  is <- as.vector(as.numeric(rownames(promt)))
  promt <- cbind(is = is,promt)
  promt <- promt[order(is),]
  return(promt)
}
MR_result_sig_select <- function(x = MR_file){
  output_path <- paste0("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_sig.txt")
  MR_sig <- data.frame()
  for(i in MR_file){
    print(i)
    path_prmot <- paste0(list.dirs("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm")[5],'/',i)
    data_prmot <- read.delim(path_prmot)
    #data_prmot <- FDR_adjust(i)
    MR_result_sig_fix <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (fixed effects)',]
    MR_result_sig_fix <- MR_result_sig_fix[MR_result_sig_fix$FDR < 5e-02,]
    MR_result_sig_random <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (multiplicative random effects)',]
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$FDR < 5e-02,]
    
    path_prmot_het <- paste0("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_heterogeneity/",i)
    data_prmot_het <- read.delim(path_prmot_het)
    data_prmot_het <- data_prmot_het[data_prmot_het$method == 'Inverse variance weighted',]
    data_prmot_het <- data_prmot_het[data_prmot_het$Q_pval < 0.05,]
    
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$id.exposure %in% data_prmot_het$id.exposure & MR_result_sig_random$id.outcome %in% data_prmot_het$id.outcome,]
    MR_result_sig_fix <- MR_result_sig_fix[! MR_result_sig_fix$id.exposure %in% data_prmot_het$id.exposure & MR_result_sig_fix$id.outcome %in% data_prmot_het$id.outcome,]
    
    MR_sig_IVW <- rbind(MR_result_sig_fix,MR_result_sig_random)
    MR_sig_IVW <- MR_sig_IVW[MR_sig_IVW$pval < 0.05,]
    MR_result_sing <- data_prmot[data_prmot$id.exposure %in% MR_sig_IVW$id.exposure & data_prmot$id.outcome %in% MR_sig_IVW$id.outcome,]
    
    MR_result_sing1 <- MR_result_sing[MR_result_sing$id.exposure %in% data_prmot_het$id.exposure & MR_result_sing$id.outcome %in% data_prmot_het$id.outcome & MR_result_sing$method !='Inverse variance weighted (fixed effects)',]
    MR_result_sing2 <- MR_result_sing[!MR_result_sing$id.exposure %in% data_prmot_het$id.exposure & MR_result_sing$id.outcome %in% data_prmot_het$id.outcome & MR_result_sing$method !='Inverse variance weighted (multiplicative random effects)',]
    
    MR_result <- rbind(MR_result_sing1,MR_result_sing2)
    
    if(dim(MR_sig)[1]==0){
      MR_sig <- MR_result
    }else{
      MR_sig <- rbind(MR_sig,MR_result)
    }
    
    
    if(!file.exists(output_path)){
      write.table(MR_result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
    }else{
      write.table(MR_result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = F)
    }
    
    
  }
  return(MR_sig)
}


res <- MR_result_sig_select(MR_file)

MR_pleiotropy_dir <- list.dirs(Result_storge_dir)[3]
MR_file <- list.files(MR_pleiotropy_dir)
MR_pleiotropy_phe <- function(x = MR_file){
  ple_phe <- data.frame()
  No_ple_phe <- data.frame()
  for(i in x){
    print(i)
    ple_path <- paste0("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_pleiotropy/",i)
    data_promt <- read.delim(ple_path)
    data_ple <- data_promt[data_promt$pval < 0.05,]
    data_no_ple <- data_promt[data_promt$pval >=0.05,]
    
    if(dim(ple_phe)[1]==0){
      ple_phe <- data_ple
    }else{
      ple_phe <- rbind(ple_phe,data_ple)
    }
  }
  
  ple_phe <- na.omit(ple_phe)
  if(dim(No_ple_phe)[1]==0){
    No_ple_phe <- data_no_ple
  }else{
    No_ple_phe <- rbind(No_ple_phe,data_no_ple)
  }
  
  ple_distinc <- list(ple = ple_phe,unple = No_ple_phe)
  
  return(ple_distinc)
}
MR_ple <- MR_pleiotropy_phe()

MR_presso_dir <- list.dirs(Result_storge_dir)[4]
MR_file <- list.files(MR_presso_dir)
MR_presso_pleiotropy_phe <- function(x = MR_file){
  id_data_frame <- data.frame(id.outcome = c('ebi-a-GCST004747','ebi-a-GCST004749','ukb-a-54','ieu-a-984','ebi-a-GCST900013972',
                                             'ukb-b-14521','ukb-b-20176','ukb-b-15826','finn-b-C3_LUNG_NOSMALL_EXALLC','ieu-a-988',
                                             'ieu-a-989'),MR_file = x)
  ple_phe <- data.frame()
  No_ple_phe <- data.frame()
  for(i in c(1:11)){
    print(id_data_frame$MR_file[i])
    presso_path <- paste0(MR_presso_dir,'/',id_data_frame$MR_file[i])
    presso_data <- read.delim(presso_path)
    presso_data <- cbind(presso_data,id.outcome = rep(id_data_frame$id.outcome[i],dim(presso_data)[1]))
    presso_ple <- presso_data[presso_data$global_pval < 0.05,]
    presso_ple <- na.omit(presso_ple)
    presso_ple <- presso_data[presso_data$id.exposure %in% presso_ple$id.exposure,]
    presso_no_ple <- presso_data[presso_data$global_pval >= 0.05,]
    if(dim(ple_phe)[1]==0){
      ple_phe <- presso_ple
    }else{
      ple_phe <- rbind(ple_phe,presso_ple)
    }
    
    if(dim(No_ple_phe)[1]==0){
      No_ple_phe <- presso_no_ple
    }else{
      No_ple_phe <- rbind(No_ple_phe,presso_no_ple)
    }
    
    
  }
  
  Presso_ple <- list(ple_phe = ple_phe,No_ple_phe=No_ple_phe)
  return(Presso_ple)
}
Presso_ple <- MR_presso_pleiotropy_phe()

MR_sig <- read.delim("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_sig.txt")
MR_ple_ple <- MR_ple$ple[MR_ple$ple$id.exposure %in% MR_sig$id.exposure & MR_ple$ple$id.outcome %in% MR_sig$id.outcome,]
MR_presso_ple <- Presso_ple$ple_phe[Presso_ple$ple_phe$id.exposure %in% MR_sig$id.exposure & Presso_ple$ple_phe$id.outcome %in% MR_sig$id.outcome,]
MR_presso_no_ple <- Presso_ple$No_ple_phe[Presso_ple$No_ple_phe$id.exposure %in% MR_sig$id.exposure & Presso_ple$No_ple_phe$id.outcome %in% MR_sig$id.outcome,]

#Stbleest result
MR_stable_result <- MR_sig[MR_sig$id.exposure %in% MR_ple$unple$id.exposure & MR_sig$id.outcome %in% MR_ple$ple$id.outcome,]
MR_stable_result <- MR_stable_result[MR_stable_result$id.exposure %in% Presso_ple$No_ple_phe$id.exposure & MR_stable_result$id.outcome %in% Presso_ple$No_ple_phe$id.outcome,]
MR_stable_result <- na.omit(MR_stable_result)

load("/Users/llls2012163.com/GWAS/TraitID.Rdata")
ID_to_immTrait <- ID_Trait$ImmID
res <- data.frame()
for(i in 1:731){
  data_promt <- MR_stable_result[MR_stable_result$id.outcome == ID_to_immTrait$id[i],]
  data_promt <- cbind(Trait = rep(ID_to_immTrait$X[i],dim(data_promt)[1]),data_promt)
  if(dim(res)[1]==0){
    res <- data_promt
  }else{
    res <- rbind(res,data_promt)
  }
}

res <- res[order(res$id.outcome),]
res <- cbind(res[,1:9],OR=exp(res$b),'Lower 95% CI'=exp(res$b - 1.96*res$se),'Upper 95% CI' = exp(res$b + 1.96*res$se),pval = res[,10])
write.table(res,file = "/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_sig_stable.txt",sep = '\t',quote = FALSE,row.names = FALSE)

rm(list = ls())

# data visualization


MR_sig_stable <- read.delim("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_sig_stable.txt")
MR_sig_stable$outcome <- limma::strsplit2(MR_sig_stable$outcome," \\|\\| ")[,1]
MR_sig_stable <- MR_sig_stable[,c(1,2,4,6:13)]

load("/Users/llls2012163.com/GWAS/TraitID.RData")
id_trait <- ID_Trait$ImmID

traitpanel <- read_excel("Desktop/traitpanel.xlsx")
q <- data.frame()
for(i in 1:dim(MR_sig_stable)[1]){
  trait <- MR_sig_stable$id.exposure[i]
  ids <- limma::strsplit2(trait,'-')[,3]
  a <- traitpanel[traitpanel$`GWAS Catalog Accession Number`==ids,]
  if(length(q)==0){
    q <- a
  }else{
    q <- rbind(q,a)
  }
}
MR_sig_stable <- cbind(q,MR_sig_stable)
MR_sig_stable1 <- MR_sig_stable[MR_sig_stable$method %in% c('Inverse variance weighted (fixed effects)','Inverse variance weighted (multiplicative random effects)'),]
colnames(MR_sig_stable1)[1] <- 'TraitType'

library(webr)
library(tidyverse)
p1<-  PieDonut(MR_sig_stable1,aes(Panel, outcome),
               r0=0.2,
               r1=0.6,
               r2=0.9,
               pieLabelSize = 2.5,
               donutLabelSize = 2.5,labelposition = 0)


p2 <- PieDonut(MR_sig_stable1,aes(TraitType, outcome),
               r0=0.2,
               r1=0.6,
               r2=0.9,
               pieLabelSize = 2.5,
               donutLabelSize = 2.5,labelposition = 0)

library(cowplot)
library(patchwork)
plot_grid(p1,p2,ncol = 2)
p1+p2+plot_layout(ncol = 2)

###### NULL result

####3. 731 immune cell traits--------->41 circulating inflammatory cytokines

Result_storge_dir <- "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf"

MR_result_dir <- list.dirs(Result_storge_dir)[5]
MR_het_dir <- list.dirs(Result_storge_dir)[2]
MR_file <- list.files(MR_result_dir)

FDR_adjust <- function(x = MR_file){
  for(i in x){
    print(i)
    path_prmot <- paste0(list.dirs("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf")[5],'/',i)
    data_prmot <- read.delim(path_prmot)
    data_prmot <- na.omit(data_prmot)
    promt <- data.frame()
    for(j in data_prmot$method[!duplicated(data_prmot$method)]){
      print(j)
      promt1 <- data_prmot[data_prmot$method %in% j,]
      promt1$FDR <- p.adjust(promt1$pval,method = 'BH',length(promt1$pval))
      
      if(dim(promt)[1]==0){
        promt <- promt1
      }else{
        promt <- rbind(promt,promt1)
      }
    }
  }
  is <- as.vector(as.numeric(rownames(promt)))
  promt <- cbind(is = is,promt)
  promt <- promt[order(is),]
  return(promt)
}


MR_result_sig_select <- function(x = MR_file){
  output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf/MR_sig.txt")
  MR_sig <- data.frame()
  for(i in x){
    print(i)
    path_prmot <- paste0(list.dirs("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf")[5],'/',i)
    data_prmot <- read.delim(path_prmot)
    #data_prmot <- FDR_adjust(i)
    MR_result_sig_fix <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (fixed effects)',]
    MR_result_sig_fix <- MR_result_sig_fix[MR_result_sig_fix$pval < 5e-02,]
    MR_result_sig_random <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (multiplicative random effects)',]
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$pval < 5e-02,]
    
    path_prmot_het <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf/MR_HeterogeneityimmToinf/",i)
    data_prmot_het <- read.delim(path_prmot_het)
    data_prmot_het <- data_prmot_het[data_prmot_het$method == 'Inverse variance weighted',]
    data_prmot_het <- data_prmot_het[data_prmot_het$Q_pval < 0.05,]
    
    
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$id.exposure %in% data_prmot_het$id.exposure & MR_result_sig_random$id.outcome %in% data_prmot_het$id.outcome,]
    MR_result_sig_fix <- MR_result_sig_fix[! MR_result_sig_fix$id.exposure %in% data_prmot_het$id.exposure & MR_result_sig_fix$id.outcome %in% data_prmot_het$id.outcome,]
    
    MR_sig_IVW <- rbind(MR_result_sig_fix,MR_result_sig_random)
    MR_sig_IVW <- MR_sig_IVW[MR_sig_IVW$pval < 0.05,]
    MR_result_sing <- data_prmot[data_prmot$id.exposure %in% MR_sig_IVW$id.exposure & data_prmot$id.outcome %in% MR_sig_IVW$id.outcome,]
    
    MR_result_sing1 <- MR_result_sing[MR_result_sing$id.exposure %in% data_prmot_het$id.exposure & MR_result_sing$id.outcome %in% data_prmot_het$id.outcome & MR_result_sing$method !='Inverse variance weighted (fixed effects)',]
    MR_result_sing2 <- MR_result_sing[!MR_result_sing$id.exposure %in% data_prmot_het$id.exposure & MR_result_sing$id.outcome %in% data_prmot_het$id.outcome & MR_result_sing$method !='Inverse variance weighted (multiplicative random effects)',]
    
    MR_result <- rbind(MR_result_sing1,MR_result_sing2)
    
    if(dim(MR_sig)[1]==0){
      MR_sig <- MR_result
    }else{
      MR_sig <- rbind(MR_sig,MR_result)
    }
    
    
    if(!file.exists(output_path)){
      write.table(MR_result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
    }else{
      write.table(MR_result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = F)
    }
    
    
  }
  return(MR_sig)
  
}
res <- MR_result_sig_select(MR_file)



MR_pleiotropy_dir <- list.dirs(Result_storge_dir)[3]
MR_file <- list.files(MR_pleiotropy_dir)
MR_pleiotropy_phe <- function(x = MR_file){
  ple_phe <- data.frame()
  No_ple_phe <- data.frame()
  for(i in MR_file){
    print(i)
    ple_path <- paste0(MR_pleiotropy_dir,'/',i)
    data_promt <- read.delim(ple_path)
    data_ple <- data_promt[data_promt$pval < 0.05,]
    data_no_ple <- data_promt[data_promt$pval >=0.05,]
    
    if(dim(ple_phe)[1]==0){
      ple_phe <- data_ple
    }else{
      ple_phe <- rbind(ple_phe,data_ple)
    }
    
    ple_phe <- na.omit(ple_phe)
    if(dim(No_ple_phe)[1]==0){
      No_ple_phe <- data_no_ple
    }else{
      No_ple_phe <- rbind(No_ple_phe,data_no_ple)
    }
  }
  
  
  
  ple_distinc <- list(ple = ple_phe,unple = No_ple_phe)
  
  return(ple_distinc)
}
MR_ple <- MR_pleiotropy_phe()



MR_presso_dir <- list.dirs(Result_storge_dir)[4]
MR_file <- list.files(MR_presso_dir)

MR_presso_pleiotropy_phe <- function(x = MR_file){
  id_data_frame <- data.frame(id.outcome = c('bNGF','CTACK','EOTAXIN','FGF-BASIC','G-CSF','GROA','HGF','IFN-G','IL_10','IL_17',"IL-12p70","IL-13",
                                             "IL-16","IL-18","IL-1a","IL-1ra","IL-2","IL-2a","IL-4","IL-5","IL-6","IL-7","IL-8","IL-9",
                                             "IP-10","M-CSF","MCP-1:MCAF","MCP-3","MIF","MIG","MIP-1a","MIP-1b","PDGF-BB","RANTES","SCF","SCGF-B",
                                             "SDF-1a","TNF-a","TNF-B","TRAIL","VEGF"),MR_file = MR_file)
  ple_phe <- data.frame()
  No_ple_phe <- data.frame()
  for(i in c(1:41)){
    print(id_data_frame$MR_file[i])
    presso_path <- paste0(MR_presso_dir,'/',id_data_frame$MR_file[i])
    presso_data <- read.delim(presso_path)
    presso_data <- cbind(presso_data[,-10],id.outcome = rep(id_data_frame$id.outcome[i],dim(presso_data)[1]))
    presso_ple <- presso_data[presso_data$global_pval < 0.05,]
    presso_ple <- na.omit(presso_ple)
    presso_ple <- presso_data[presso_data$id.exposure %in% presso_ple$id.exposure,]
    presso_no_ple <- presso_data[presso_data$global_pval >= 0.05,]
    if(dim(ple_phe)[1]==0){
      ple_phe <- presso_ple
    }else{
      ple_phe <- rbind(ple_phe,presso_ple)
    }
    
    if(dim(No_ple_phe)[1]==0){
      No_ple_phe <- presso_no_ple
    }else{
      No_ple_phe <- rbind(No_ple_phe,presso_no_ple)
    }
    
    
  }
  
  Presso_ple <- list(ple_phe = ple_phe,No_ple_phe=No_ple_phe)
  return(Presso_ple)
}

Presso_ple <- MR_presso_pleiotropy_phe()

MR_sig <- read.delim("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf/MR_sig.txt")
MR_ple_ple <- MR_ple$ple[MR_ple$ple$id.exposure %in% MR_sig$id.exposure & MR_ple$ple$id.outcome %in% MR_sig$id.outcome,]
MR_presso_ple <- Presso_ple$ple_phe[Presso_ple$ple_phe$id.exposure %in% MR_sig$id.exposure & Presso_ple$ple_phe$id.outcome %in% MR_sig$id.outcome,]
MR_presso_no_ple <- Presso_ple$No_ple_phe[Presso_ple$No_ple_phe$id.exposure %in% MR_sig$id.exposure & Presso_ple$No_ple_phe$id.outcome %in% MR_sig$id.outcome,]


MR_stable_result <- MR_sig[MR_sig$id.exposure %in% MR_ple$unple$id.exposure & MR_sig$id.outcome %in% MR_ple$ple$id.outcome,]
MR_stable_result <- MR_stable_result[MR_stable_result$id.exposure %in% Presso_ple$No_ple_phe$id.exposure & MR_stable_result$id.outcome %in% Presso_ple$No_ple_phe$id.outcome,]
MR_stable_result <- na.omit(MR_stable_result)


load("/Users/llls2012163.com/GWAS/TraitID.Rdata")
ID_to_immTrait <- ID_Trait$ImmID
res <- data.frame()
for(i in 1:731){
  data_promt <- MR_stable_result[MR_stable_result$id.exposure == ID_to_immTrait$id[i],]
  data_promt <- cbind(Trait = rep(ID_to_immTrait$X[i],dim(data_promt)[1]),data_promt)
  if(dim(res)[1]==0){
    res <- data_promt
  }else{
    res <- rbind(res,data_promt)
  }
}

res <- res[order(res$id.outcome),]
res <- cbind(res[,1:9],OR=exp(res$b),'Lower 95% CI'=exp(res$b - 1.96*res$se),'Upper 95% CI' = exp(res$b + 1.96*res$se),pval=res[,10])
write.table(res,file = "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf/MR_sig_stable.txt",sep = '\t',quote = FALSE,row.names = FALSE)
rm(list = ls())

###### data visualization---------------------------------------------------

MR_sig_stable <- read.delim("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf/MR_sig_stable.txt")

MR_sig_stable <- MR_sig_stable[,c(1,2,3,5:13)]

load("/Users/llls2012163.com/GWAS/TraitID.RData")
id_trait <- ID_Trait$ImmID

traitpanel <- read_excel("Desktop/traitpanel.xlsx")
q <- data.frame()
for(i in 1:dim(MR_sig_stable)[1]){
  trait <- MR_sig_stable$id.exposure[i]
  ids <- limma::strsplit2(trait,'-')[,3]
  a <- traitpanel[traitpanel$`GWAS Catalog Accession Number`==ids,]
  if(length(q)==0){
    q <- a
  }else{
    q <- rbind(q,a)
  }
}
MR_sig_stable <- cbind(q,MR_sig_stable)

write.table(MR_sig_stable,file = "/Users/llls2012163.com/add5.txt",sep = '\t',quote = FALSE,row.names = FALSE)

MR_sig_stable1 <- MR_sig_stable[MR_sig_stable$method %in% c('Inverse variance weighted (fixed effects)','Inverse variance weighted (multiplicative random effects)'),]
colnames(MR_sig_stable1)[1] <- 'TraitType'

seq <- as.data.frame(table(MR_sig_stable1$id.outcome))
seq <- seq[order(seq$Freq,decreasing = TRUE),]
colnames(seq) <- c('exposure','seq')
seq$exposure <- levels(seq$exposure)

library(ggplot2)
library(ggprism)
p <- ggplot(seq,aes(exposure,seq))+geom_bar(stat = 'identity',fill = '#40E0D0')+theme_prism()+
  theme(axis.text.x = element_text(angle = 30,size= 10))+geom_text(aes(label = seq),vjust = -0.2)+
  labs(x= '',y='Number of immune cell traits')
p

library(webr)
library(tidyverse)
p1<-  PieDonut(MR_sig_stable1,aes(Panel, TraitType),
               r0=0.2,
               r1=0.6,
               r2=0.9,
               pieLabelSize = 2.5,
               donutLabelSize = 2.5,labelposition = 0)


####4. 41 circulating inflammatory cytokines-------------------->731 immune cell traits
Result_storge_dir <- "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/inf_imm"

MR_result_dir <- list.dirs(Result_storge_dir)[5]
MR_het_dir <- list.dirs(Result_storge_dir)[2]
MR_file <- list.files(MR_result_dir)
MR_result_sig_select <- function(x = MR_file){
  output_path <- paste0("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/inf_imm/MR_sig.txt")
  MR_sig <- data.frame()
  for(i in x){
    print(i)
    path_prmot <- paste0(MR_result_dir,'/',i)
    data_prmot <- read.delim(path_prmot)
    MR_result_sig_fix <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (fixed effects)',]
    MR_result_sig_fix <- MR_result_sig_fix[MR_result_sig_fix$pval < 5e-02,]
    MR_result_sig_random <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (multiplicative random effects)',]
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$pval < 5e-02,]
    
    path_prmot_het <- paste0(MR_het_dir,'/',i)
    data_prmot_het <- read.delim(path_prmot_het)
    data_prmot_het <- data_prmot_het[data_prmot_het$method == 'Inverse variance weighted',]
    data_prmot_het <- data_prmot_het[data_prmot_het$Q_pval < 0.05,]
    
    
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$id.outcome %in% data_prmot_het$id.outcome,]
    MR_result_sig_fix <- MR_result_sig_fix[!MR_result_sig_fix$id.outcome %in% data_prmot_het$id.outcome,]
    
    MR_sig_IVW <- rbind(MR_result_sig_fix,MR_result_sig_random)
    MR_sig_IVW <- MR_sig_IVW[MR_sig_IVW$pval < 0.05,]
    MR_result_sing <- data_prmot[data_prmot$id.outcome %in% MR_sig_IVW$id.outcome,]
    
    MR_result_sing1 <- MR_result_sing[ MR_result_sing$id.outcome %in% data_prmot_het$id.outcome & MR_result_sing$method !='Inverse variance weighted (fixed effects)',]
    MR_result_sing2 <- MR_result_sing[!MR_result_sing$id.outcome %in% data_prmot_het$id.outcome & MR_result_sing$method !='Inverse variance weighted (multiplicative random effects)',]
    
    MR_result <- rbind(MR_result_sing1,MR_result_sing2)
    
    if(dim(MR_sig)[1]==0){
      MR_sig <- MR_result
    }else{
      MR_sig <- rbind(MR_sig,MR_result)
    }
    
    
    if(!file.exists(output_path)){
      write.table(MR_result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
    }else{
      write.table(MR_result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = F)
    }
    
    
  }
  return(MR_sig)
  
}
res <- MR_result_sig_select(MR_file)

MR_pleiotropy_dir <- list.dirs(Result_storge_dir)[3]
MR_file <- list.files(MR_pleiotropy_dir)

MR_pleiotropy_phe <- function(x = MR_file){
  ple_phe <- data.frame()
  No_ple_phe <- data.frame()
  for(i in MR_file){
    print(i)
    ple_path <- paste0(MR_pleiotropy_dir,'/',i)
    data_promt <- read.delim(ple_path)
    data_ple <- data_promt[data_promt$pval < 0.05,]
    data_no_ple <- data_promt[data_promt$pval >=0.05,]
    
    if(dim(ple_phe)[1]==0){
      ple_phe <- data_ple
    }else{
      ple_phe <- rbind(ple_phe,data_ple)
    }
    
    ple_phe <- na.omit(ple_phe)
    if(dim(No_ple_phe)[1]==0){
      No_ple_phe <- data_no_ple
    }else{
      No_ple_phe <- rbind(No_ple_phe,data_no_ple)
    }
  }
  
  
  
  ple_distinc <- list(ple = ple_phe,unple = No_ple_phe)
  
  return(ple_distinc)
}
MR_ple <- MR_pleiotropy_phe()


MR_presso_dir <- list.dirs(Result_storge_dir)[4]
MR_file <- list.files(MR_presso_dir)

MR_presso_pleiotropy_phe <- function(x = MR_file){
  id_data_frame <- data.frame(id.exposure = c('CTACK','EOTAXIN','FGF-BASIC','GROA','IFN-G','IL_10','IL_17',"IL-12p70",
                                             "IL-16","IL-18","IL-1a","IL-2","IL-2a","L-4","IL-5","IL-6","IL-7","IL-8",
                                             "IP-10","M-CSF","MCP-1:MCAF","MIF","MIG","MIP-1b","PDGF-BB","SCF","SCGF-B",
                                             "SDF-1a","TNF-B","TRAIL","VEGF"),MR_file = MR_file)
  ple_phe <- data.frame()
  No_ple_phe <- data.frame()
  for(i in c(1:31)){
    print(id_data_frame$MR_file[i])
    presso_path <- paste0(MR_presso_dir,'/',id_data_frame$MR_file[i])
    presso_data <- read.delim(presso_path)
    presso_data <- cbind(presso_data[,-1],id.exposure = rep(id_data_frame$id.exposure[i],dim(presso_data)[1]))
    presso_ple <- presso_data[presso_data$global_pval < 0.05,]
    presso_ple <- na.omit(presso_ple)
    presso_ple <- presso_data[presso_data$id.outcome %in% presso_ple$id.outcome,]
    presso_no_ple <- presso_data[presso_data$global_pval >= 0.05,]
    presso_no_ple <- na.omit(presso_no_ple)
    if(dim(ple_phe)[1]==0){
      ple_phe <- presso_ple
    }else{
      ple_phe <- rbind(ple_phe,presso_ple)
    }
    
    if(dim(No_ple_phe)[1]==0){
      No_ple_phe <- presso_no_ple
    }else{
      No_ple_phe <- rbind(No_ple_phe,presso_no_ple)
    }
    
    
  }
  
  Presso_ple <- list(ple_phe = ple_phe,No_ple_phe=No_ple_phe)
  return(Presso_ple)
}

Presso_ple <- MR_presso_pleiotropy_phe()


MR_sig <- read.delim("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/inf_imm/MR_sig.txt")
MR_ple_ple <- MR_ple$ple[MR_ple$ple$exposure %in% MR_sig$exposure & MR_ple$ple$id.outcome %in% MR_sig$id.outcome,]
MR_presso_ple <- Presso_ple$ple_phe[Presso_ple$ple_phe$id.exposure %in% MR_sig$exposure & Presso_ple$ple_phe$id.outcome %in% MR_sig$id.outcome,]
MR_presso_no_ple <- Presso_ple$No_ple_phe[Presso_ple$No_ple_phe$id.exposure %in% MR_sig$exposure & Presso_ple$No_ple_phe$id.outcome %in% MR_sig$id.outcome,]


MR_stable_result <- MR_sig[MR_sig$id.exposure %in% MR_ple$unple$id.exposure & MR_sig$id.outcome %in% MR_ple$ple$id.outcome,]
MR_stable_result <- MR_stable_result[MR_stable_result$exposure %in% Presso_ple$No_ple_phe$id.exposure & MR_stable_result$id.outcome %in% Presso_ple$No_ple_phe$id.outcome,]
MR_stable_result <- na.omit(MR_stable_result)


load("/Users/llls2012163.com/GWAS/TraitID.Rdata")
ID_to_immTrait <- ID_Trait$ImmID
res <- data.frame()
for(i in 1:731){
  data_promt <- MR_stable_result[MR_stable_result$id.outcome == ID_to_immTrait$id[i],]
  data_promt <- cbind(Trait = rep(ID_to_immTrait$X[i],dim(data_promt)[1]),data_promt)
  if(dim(res)[1]==0){
    res <- data_promt
  }else{
    res <- rbind(res,data_promt)
  }
}

res <- res[order(res$id.exposure),]
res <- cbind(res[,1:9],OR=exp(res$b),'Lower 95% CI'=exp(res$b - 1.96*res$se),'Upper 95% CI' = exp(res$b + 1.96*res$se),pval = res[,10])
write.table(res,file = "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/inf_imm/MR_sig_stable.txt",sep = '\t',quote = FALSE,row.names = FALSE)
rm(list = ls())

####5. 41 circulating inflammatory cytokines-------------------->Lung cancer
Result_storge_dir <- "/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung"

MR_result_dir <- list.dirs(Result_storge_dir)[5]
MR_het_dir <- list.dirs(Result_storge_dir)[2]
MR_file <- list.files(MR_result_dir)

FDR_adjust <- function(x = MR_file){
  for(i in x){
    print(i)
    path_prmot <- paste0(list.dirs("/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung")[5],'/',i)
    data_prmot <- read.delim(path_prmot)
    data_prmot <- na.omit(data_prmot)
    promt <- data.frame()
    for(j in data_prmot$method[!duplicated(data_prmot$method)]){
      print(j)
      promt1 <- data_prmot[data_prmot$method %in% j,]
      promt1$FDR <- p.adjust(promt1$pval,method = "bonferroni",length(promt1$pval))
      
      if(dim(promt)[1]==0){
        promt <- promt1
      }else{
        promt <- rbind(promt,promt1)
      }
    }
  }
  is <- as.vector(as.numeric(rownames(promt)))
  promt <- cbind(is = is,promt)
  promt <- promt[order(is),]
  return(promt)
}
MR_result_sig_select <- function(x = MR_file){
  output_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung/MR_sig.txt")
  MR_sig <- data.frame()
  for(i in x){
    print(i)
    path_prmot <- paste0(MR_result_dir,'/',i)
    data_prmot <- read.delim(path_prmot)
    #data_prmot <- FDR_adjust(i)
    MR_result_sig_fix <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (fixed effects)',]
    MR_result_sig_fix <- MR_result_sig_fix[MR_result_sig_fix$pval < 5e-02,]
    MR_result_sig_random <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (multiplicative random effects)',]
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$pval < 5e-02,]
    
    path_prmot_het <- paste0(MR_het_dir,'/',i)
    data_prmot_het <- read.delim(path_prmot_het)
    data_prmot_het <- data_prmot_het[data_prmot_het$method == 'Inverse variance weighted',]
    data_prmot_het <- data_prmot_het[data_prmot_het$Q_pval < 0.05,]
    
    
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$exposure_inf %in% data_prmot_het$exposure_inf,]
    MR_result_sig_fix <- MR_result_sig_fix[!MR_result_sig_fix$exposure_inf %in% data_prmot_het$exposure_inf,]
    
    MR_sig_IVW <- rbind(MR_result_sig_fix,MR_result_sig_random)
    MR_sig_IVW <- MR_sig_IVW[MR_sig_IVW$pval < 0.05,]
    MR_result_sing <- data_prmot[data_prmot$exposure_inf %in% MR_sig_IVW$exposure_inf,]
    
    MR_result_sing1 <- MR_result_sing[ MR_result_sing$exposure_inf %in% data_prmot_het$exposure_inf & MR_result_sing$method !='Inverse variance weighted (fixed effects)',]
    MR_result_sing2 <- MR_result_sing[!MR_result_sing$exposure_inf %in% data_prmot_het$exposure_inf & MR_result_sing$method !='Inverse variance weighted (multiplicative random effects)',]
    
    MR_result <- rbind(MR_result_sing1,MR_result_sing2)
    cancerType <- limma::strsplit2(i,'\\.')[,1]
    MR_result <- cbind(id.outcome = rep(cancerType,dim(MR_result)[1]),MR_result[,-3])
    
    if(dim(MR_sig)[1]==0){
      MR_sig <- MR_result
    }else{
      MR_sig <- rbind(MR_sig,MR_result)
    }
    
    
    if(!file.exists(output_path)){
      write.table(MR_result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
    }else{
      write.table(MR_result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = F)
    }
    
    
  }
  return(MR_sig)
  
}
res <- MR_result_sig_select(MR_file)

res <- read.delim("/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung/MR_sig.txt")



MR_pleiotropy_dir <- list.dirs(Result_storge_dir)[3]
MR_file <- list.files(MR_pleiotropy_dir)
MR_pleiotropy_phe <- function(x = MR_file){
  ple_phe <- data.frame()
  No_ple_phe <- data.frame()
  for(i in MR_file){
    print(i)
    ple_path <- paste0(MR_pleiotropy_dir,'/',i)
    data_promt <- read.delim(ple_path)
    cancerType <- limma::strsplit2(i,'\\.')[,1]
    data_promt <- cbind(id.outcome = rep(cancerType,dim(data_promt)[1]),data_promt[,-3])
    
    data_ple <- data_promt[data_promt$pval < 0.05,]
    data_no_ple <- data_promt[data_promt$pval >=0.05,]
    
    if(dim(ple_phe)[1]==0){
      ple_phe <- data_ple
    }else{
      ple_phe <- rbind(ple_phe,data_ple)
    }
    
    ple_phe <- na.omit(ple_phe)
    if(dim(No_ple_phe)[1]==0){
      No_ple_phe <- data_no_ple
    }else{
      No_ple_phe <- rbind(No_ple_phe,data_no_ple)
    }
  }
  
  
  
  ple_distinc <- list(ple = ple_phe,unple = No_ple_phe)
  
  return(ple_distinc)
}
MR_ple <- MR_pleiotropy_phe()


MR_presso_dir <- list.dirs(Result_storge_dir)[4]
MR_file <- list.files(MR_presso_dir)

MR_presso_pleiotropy_phe <- function(x = MR_file){
  No_ple_phe <- data.frame()
  ple_phe <- data.frame()
  for(i in x){
    print(i)
    presso_path <- paste0(MR_presso_dir,'/',i)
    presso_data <- read.delim(presso_path)
    cancerType <- limma::strsplit2(i,'\\.')[,1]
    presso_data <- cbind(id.outcome = rep(cancerType,dim(presso_data)[1]),presso_data)
    presso_ple <- presso_data[presso_data$global_pval < 0.05,]
    presso_ple <- na.omit(presso_ple)
    presso_ple <- presso_data[presso_data$id.exposure %in% presso_ple$id.exposure,]
    presso_no_ple <- presso_data[presso_data$global_pval >= 0.05,]
    presso_no_ple <- na.omit(presso_no_ple)
    if(dim(ple_phe)[1]==0){
      ple_phe <- presso_ple
    }else{
      ple_phe <- rbind(ple_phe,presso_ple)
    }
    
    if(dim(No_ple_phe)[1]==0){
      No_ple_phe <- presso_no_ple
    }else{
      No_ple_phe <- rbind(No_ple_phe,presso_no_ple)
    }
    
    
  }
  
  Presso_ple <- list(ple_phe = ple_phe,No_ple_phe=No_ple_phe)
  return(Presso_ple)
}
Presso_ple <- MR_presso_pleiotropy_phe()

MR_sig <- read.delim("/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung/MR_sig.txt")
MR_ple_ple <- MR_ple$ple[MR_ple$ple$id.exposure %in% MR_sig$id.exposure & MR_ple$ple$id.outcome %in% MR_sig$id.outcome,]
MR_presso_ple <- Presso_ple$ple_phe[Presso_ple$ple_phe$id.exposure %in% MR_sig$id.exposure & Presso_ple$ple_phe$id.outcome %in% MR_sig$id.outcome,]
MR_presso_no_ple <- Presso_ple$No_ple_phe[Presso_ple$No_ple_phe$id.exposure %in% MR_sig$id.exposure & Presso_ple$No_ple_phe$id.outcome %in% MR_sig$id.outcome,]


MR_stable_result <- MR_sig[MR_sig$id.exposure %in% MR_ple$unple$id.exposure & MR_sig$id.outcome %in% MR_ple$ple$id.outcome,]
MR_stable_result <- MR_stable_result[MR_stable_result$id.exposure %in% Presso_ple$No_ple_phe$id.exposure & MR_stable_result$id.outcome %in% Presso_ple$No_ple_phe$id.outcome,]
MR_stable_result <- na.omit(MR_stable_result)

res <- MR_stable_result
res <- res[order(res$id.exposure),]
res <- cbind(res[,1:9],OR=exp(res$b),'Lower 95% CI'=exp(res$b - 1.96*res$se),'Upper 95% CI' = exp(res$b + 1.96*res$se),pval = res[,10])
write.table(res,file = "/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung/MR_sig_stable.txt",sep = '\t',quote = FALSE,row.names = FALSE)
rm(list = ls())

# data visulation

MR_sig_stable <- read.delim("/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung/MR_sig_stable.txt")

write.table(MR_sig_stable,file = "/Users/llls2012163.com/MR_sig_stable.txt",sep = '\t',quote = FALSE,row.names = FALSE)

MR_sig_stable1 <- MR_sig_stable[MR_sig_stable$method %in% c('Inverse variance weighted (fixed effects)','Inverse variance weighted (multiplicative random effects)'),]

library(forestplot)
result_inf_Lung_med <- MR_sig_stable1

result_inf_Lung_med <- result_inf_Lung_med[order(result_inf_Lung_med$id.outcome),]
result_inf_Lung_med <- result_inf_Lung_med[,c(1,2,10,11,12,13)]


CI <- paste0(round(result_inf_Lung_med$OR,4),'(',round(result_inf_Lung_med$Lower.95..CI,4),'-',round(result_inf_Lung_med$Upper.95..CI,4),')')
result_inf_Lung_med$pval <- round(result_inf_Lung_med$pval,4)
result_inf_Lung_med <- cbind(result_inf_Lung_med,CI)
result_inf_Lung_med$id.outcome <- c('Lung adenocarcinoma',rep(NA,6),'Small cell lung carcinoma',rep(NA,9),'Squamous cell lung cancer',rep(NA,3))
data_forest <- rbind(c('outcome','exposure',NA,NA,NA,'P.value','OR(95%CI)'),result_inf_Lung_med)


forestplot(
  data_forest[,c(1,2,7,6)],        # 需要显示在森林图中的列
  mean = data_forest[, 3],             # 均值列（HR），它将显示为森林图的小方块或其他形状哈哈哈哈哈
  lower = data_forest[, 4],            # 95%置信区间的下限数据列
  upper = data_forest[, 5],            # 95%置信区间的上限数据列
  zero = 1,                        # 均值为1时的参考线，也就是零线
  boxsize = 0.2,                   # 方框的大小
  graph.pos = "right",             # 森林图在右侧显示
  hrzl_lines = list(               # 水平线样式的设置
    "1" = gpar(lty = 1, lwd = 1),  # 均值线
    "2" = gpar(lty = 2),           # 下限和上限之间的虚线
    "14" = gpar(lwd = 2, lty = 0)# 下限和上限线
  ),xticks = c(0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
)



####6. Lung cancer---------------------> 41 circulating inflammatory

Result_storge_dir <- "/Users/llls2012163.com/Lung cancer/inf_Lung/result/LungToinf"

MR_result_dir <- list.dirs(Result_storge_dir)[5]
MR_het_dir <- list.dirs(Result_storge_dir)[2]
MR_file <- list.files(MR_result_dir)
MR_result_sig_select <- function(x = MR_file){
  output_path <- paste0("/Users/llls2012163.com/Lung cancer/inf_Lung/result/LungToinf/MR_sig.txt")
  MR_sig <- data.frame()
  for(i in x){
    print(i)
    path_prmot <- paste0(MR_result_dir,'/',i)
    data_prmot <- read.delim(path_prmot)
    MR_result_sig_fix <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (fixed effects)',]
    MR_result_sig_fix <- MR_result_sig_fix[MR_result_sig_fix$pval < 5e-02,]
    MR_result_sig_random <- data_prmot[data_prmot$method %in% 'Inverse variance weighted (multiplicative random effects)',]
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$pval < 5e-02,]
    
    path_prmot_het <- paste0(MR_het_dir,'/',i)
    data_prmot_het <- read.delim(path_prmot_het)
    data_prmot_het <- data_prmot_het[data_prmot_het$method == 'Inverse variance weighted',]
    data_prmot_het <- data_prmot_het[data_prmot_het$Q_pval < 0.05,]
    
    
    MR_result_sig_random <- MR_result_sig_random[MR_result_sig_random$id.outcome %in% data_prmot_het$id.outcome,]
    MR_result_sig_fix <- MR_result_sig_fix[!MR_result_sig_fix$id.outcome %in% data_prmot_het$id.outcome,]
    
    MR_sig_IVW <- rbind(MR_result_sig_fix,MR_result_sig_random)
    MR_sig_IVW <- MR_sig_IVW[MR_sig_IVW$pval < 0.05,]
    MR_result_sing <- data_prmot[data_prmot$id.outcome %in% MR_sig_IVW$id.outcome,]
    
    MR_result_sing1 <- MR_result_sing[ MR_result_sing$id.outcome %in% data_prmot_het$id.outcome & MR_result_sing$method !='Inverse variance weighted (fixed effects)',]
    MR_result_sing2 <- MR_result_sing[!MR_result_sing$id.outcome %in% data_prmot_het$id.outcome & MR_result_sing$method !='Inverse variance weighted (multiplicative random effects)',]
    
    MR_result <- rbind(MR_result_sing1,MR_result_sing2)
    cancerType <- limma::strsplit2(i,'\\.')[,1]
    MR_result <- cbind(exposure = rep(cancerType,dim(MR_result)[1]),MR_result[,-4])
    
    if(dim(MR_sig)[1]==0){
      MR_sig <- MR_result
    }else{
      MR_sig <- rbind(MR_sig,MR_result)
    }
    
    
    if(!file.exists(output_path)){
      write.table(MR_result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE)
    }else{
      write.table(MR_result,file = output_path,row.names = F,sep = '\t',quote = F,append = TRUE,col.names = F)
    }
    
    
  }
  return(MR_sig)
  
}
res <- MR_result_sig_select(MR_file)


MR_pleiotropy_dir <- list.dirs(Result_storge_dir)[3]
MR_file <- list.files(MR_pleiotropy_dir)
MR_pleiotropy_phe <- function(x = MR_file){
  ple_phe <- data.frame()
  No_ple_phe <- data.frame()
  for(i in MR_file){
    print(i)
    ple_path <- paste0(MR_pleiotropy_dir,'/',i)
    data_promt <- read.delim(ple_path)
    cancerType <- limma::strsplit2(i,'\\.')[,1]
    data_promt <- cbind(exposure = rep(cancerType,dim(data_promt)[1]),data_promt[,-4])
    
    data_ple <- data_promt[data_promt$pval < 0.05,]
    data_no_ple <- data_promt[data_promt$pval >=0.05,]
    
    if(dim(ple_phe)[1]==0){
      ple_phe <- data_ple
    }else{
      ple_phe <- rbind(ple_phe,data_ple)
    }
    
    ple_phe <- na.omit(ple_phe)
    if(dim(No_ple_phe)[1]==0){
      No_ple_phe <- data_no_ple
    }else{
      No_ple_phe <- rbind(No_ple_phe,data_no_ple)
    }
  }
  
  
  
  ple_distinc <- list(ple = ple_phe,unple = No_ple_phe)
  
  return(ple_distinc)
}
MR_ple <- MR_pleiotropy_phe()


MR_presso_dir <- list.dirs(Result_storge_dir)[4]
MR_file <- list.files(MR_presso_dir)

MR_presso_pleiotropy_phe <- function(x = MR_file){
  No_ple_phe <- data.frame()
  ple_phe <- data.frame()
  for(i in x){
    print(i)
    presso_path <- paste0(MR_presso_dir,'/',i)
    presso_data <- read.delim(presso_path)
    cancerType <- limma::strsplit2(i,'\\.')[,1]
    presso_data <- cbind(exposure = rep(cancerType,dim(presso_data)[1]),presso_data)
    presso_ple <- presso_data[presso_data$global_pval < 0.05,]
    presso_ple <- na.omit(presso_ple)
    presso_ple <- presso_data[presso_data$id.exposure %in% presso_ple$id.exposure,]
    presso_no_ple <- presso_data[presso_data$global_pval >= 0.05,]
    presso_no_ple <- na.omit(presso_no_ple)
    if(dim(ple_phe)[1]==0){
      ple_phe <- presso_ple
    }else{
      ple_phe <- rbind(ple_phe,presso_ple)
    }
    
    if(dim(No_ple_phe)[1]==0){
      No_ple_phe <- presso_no_ple
    }else{
      No_ple_phe <- rbind(No_ple_phe,presso_no_ple)
    }
    
    
  }
  
  Presso_ple <- list(ple_phe = ple_phe,No_ple_phe=No_ple_phe)
  return(Presso_ple)
}
Presso_ple <- MR_presso_pleiotropy_phe()

MR_sig <- read.delim("/Users/llls2012163.com/Lung cancer/inf_Lung/result/LungToinf/MR_sig.txt")
MR_ple_ple <- MR_ple$ple[MR_ple$ple$id.exposure %in% MR_sig$id.exposure & MR_ple$ple$id.outcome %in% MR_sig$id.outcome,]
MR_presso_ple <- Presso_ple$ple_phe[Presso_ple$ple_phe$id.exposure %in% MR_sig$id.exposure & Presso_ple$ple_phe$id.outcome %in% MR_sig$id.outcome,]
MR_presso_no_ple <- Presso_ple$No_ple_phe[Presso_ple$No_ple_phe$id.exposure %in% MR_sig$id.exposure & Presso_ple$No_ple_phe$id.outcome %in% MR_sig$id.outcome,]


MR_stable_result <- MR_sig[MR_sig$id.exposure %in% MR_ple$unple$id.exposure & MR_sig$id.outcome %in% MR_ple$ple$id.outcome,]
MR_stable_result <- MR_stable_result[MR_stable_result$exposure %in% Presso_ple$No_ple_phe$exposure & MR_stable_result$outcome %in% Presso_ple$No_ple_phe$outcome,]
MR_stable_result <- na.omit(MR_stable_result)

res <- MR_stable_result
res <- res[order(res$id.exposure),]
res <- cbind(res[,1:8],OR=exp(res$b),'Lower 95% CI'=exp(res$b - 1.96*res$se),'Upper 95% CI' = exp(res$b + 1.96*res$se),pval = res[,9])
write.table(res,file = "/Users/llls2012163.com/Lung cancer/inf_Lung/result/LungToinf/MR_sig_stable.txt",sep = '\t',quote = FALSE,row.names = FALSE)


dir.create("/Users/llls2012163.com/Lung cancer/Finall_result")



###### 
#####
#######Screen mediations
library(ggplot2)
library(forestplot)

## median filter Immune cell traits ----------> inflammation cytokines ---------> Lung cancer

result_Imm_inf <- read.delim("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf/MR_sig_stable.txt")
result_inf_Lung <- read.delim("/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung/MR_sig_stable.txt")
result_Imm_Lung <- read.delim("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType/MR_sig_stable.txt")


# reverse MR result
result_Lung_imm <- read.delim("/Users/llls2012163.com/Lung cancer/result_file_CancerTypeToimm/MR_sig_stable.txt")
result_Lung_inf <- read.delim("/Users/llls2012163.com/Lung cancer/inf_Lung/result/LungtoInf/MR_sig_stable.txt")
result_inf_Imm <- read.delim("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/inf_imm/MR_sig_stable.txt")

table(result_inf_Lung$id.outcome)

result_inf_Lung <- result_inf_Lung[order(result_inf_Lung$id.outcome),]
mediators <- result_inf_Lung[result_inf_Lung$method %in% c('Inverse variance weighted (fixed effects)','Inverse variance weighted (multiplicative random effects)'),]
mediators_la <- mediators$exposure_inf[1:7]
mediators_sc <- mediators$exposure_inf[8:17]
mediators_sq <- mediators$exposure_inf[18:21]

result_Imm_Lung $outcome <- limma::strsplit2(result_Imm_Lung$outcome," \\|\\| ")[,1]
result_Imm_Lung <- result_Imm_Lung[result_Imm_Lung$method %in% c('Inverse variance weighted (fixed effects)','Inverse variance weighted (multiplicative random effects)'),]
overall_la <- result_Imm_Lung[result_Imm_Lung$outcome == 'Lung adenocarcinoma',]
overall_sc <- result_Imm_Lung[result_Imm_Lung$outcome == 'Small cell lung carcinoma',]
overall_sq <- result_Imm_Lung[result_Imm_Lung$outcome == 'Squamous cell lung cancer',]

result_Imm_inf <- result_Imm_inf[result_Imm_inf$method %in% c('Inverse variance weighted (fixed effects)','Inverse variance weighted (multiplicative random effects)'),]
imm_inf_la <- result_Imm_inf[result_Imm_inf$id.outcome %in% mediators_la,]
imm_inf_sc <- result_Imm_inf[result_Imm_inf$id.outcome %in% mediators_sc,]
imm_inf_sq <- result_Imm_inf[result_Imm_inf$id.outcome %in% mediators_sq,]

imm_inf_la <- imm_inf_la[imm_inf_la$Trait %in% intersect(imm_inf_la$Trait,overall_la$Trait),]
imm_inf_sc <- imm_inf_sc[imm_inf_sc$Trait %in% intersect(imm_inf_sc$Trait,overall_sc$Trait),]
imm_inf_sq <- imm_inf_sq[imm_inf_sq$Trait %in% intersect(imm_inf_sq$Trait,overall_sq$Trait),]

mediators_la <- intersect(mediators_la,imm_inf_la$outcome)
mediators_sc <- intersect(mediators_sc,imm_inf_sc$outcome)
mediators_sq <- intersect(mediators_sq,imm_inf_sq$outcome)

overall_la <- overall_la[overall_la$Trait %in% intersect(overall_la$Trait,imm_inf_la$Trait),]
overall_sc <- overall_sc[overall_sc$Trait %in% intersect(overall_sc$Trait,imm_inf_sc$Trait),]
overall_sq <- overall_sq[overall_sq$Trait %in% intersect(overall_sq$Trait,imm_inf_sq$Trait),]


imm_inf_la <- imm_inf_la[imm_inf_la$Trait %in% intersect(overall_la$Trait,imm_inf_la$Trait),]
imm_inf_sc <- imm_inf_sc[imm_inf_sc$Trait %in% intersect(overall_sc$Trait,imm_inf_sc$Trait),]
imm_inf_sq <- imm_inf_sq[imm_inf_sq$Trait %in% intersect(overall_sq$Trait,imm_inf_sq$Trait),]

mediators <- mediators[-5,]

overall_la <- overall_la[order(overall_la$Trait),]
overall_sc <- overall_sc[order(overall_sc$Trait),]
overall_sq <- overall_sq[order(overall_sq$Trait),]

imm_inf_la <- imm_inf_la[order(imm_inf_la$Trait),]
imm_inf_sc <- imm_inf_sc[order(imm_inf_sc$Trait),]
imm_inf_sq <- imm_inf_sq[order(imm_inf_sq$Trait),]

Prepare_forest_overall <- function(x){
  data_promt <- x[,c(1,4,11,12,13,10)]
  CI <- paste0(round(x$OR,4),'(',round(x$Lower.95..CI,4),'-',round(x$Upper.95..CI,4),')')
  data_promt$pval <- as.numeric(data_promt$pval)
  pval <- vector()
  
  for(i in length(data_promt$pval)){
    if(data_promt$pval[i] > 1e-04){pval1 = round(data_promt$pval,4)}else{
      pval1 <- format(data_promt$pval[i],scientific = TRUE,digits=100)}
    if(length(pval)==0){
      pval <- pval1
    }else{
      pval <- c(pval,pval1)
    }
  }
  
  data_promt$pval <- pval
  data_promt <- cbind(data_promt,CI)
  data_promt <- rbind(c('exposure','outcome',NA,NA,NA,'P.value','OR(95%CI)'),data_promt)
  return(data_promt)
}

overall_la_forest <- Prepare_forest_overall(overall_la)

forestplot(
  overall_la_forest[,c(1,7,6)],        # 需要显示在森林图中的列
  mean = overall_la_forest[, 3],             # 均值列（HR），它将显示为森林图的小方块或其他形状哈哈哈哈哈
  lower = overall_la_forest[, 4],            # 95%置信区间的下限数据列
  upper = overall_la_forest[, 5],            # 95%置信区间的上限数据列
  zero = 1,                        # 均值为1时的参考线，也就是零线
  boxsize = 0.2,                   # 方框的大小
  graph.pos = "right",             # 森林图在右侧显示
  hrzl_lines = list(               # 水平线样式的设置
    "1" = gpar(lty = 1, lwd = 1),  # 均值线
    "2" = gpar(lty = 2),           # 下限和上限之间的虚线
    "14" = gpar(lwd = 2, lty = 0)# 下限和上限线
  ),xticks = c(0.6,0.8,1.0,1.2,1.4)
)

overall_sc_forest <- Prepare_forest_overall(overall_sc)

forestplot(
  overall_sc_forest[,c(1,7,6)],        # 需要显示在森林图中的列
  mean = overall_sc_forest[, 3],             # 均值列（HR），它将显示为森林图的小方块或其他形状哈哈哈哈哈
  lower = overall_sc_forest[, 4],            # 95%置信区间的下限数据列
  upper = overall_sc_forest[, 5],            # 95%置信区间的上限数据列
  zero = 1,                        # 均值为1时的参考线，也就是零线
  boxsize = 0.2,                   # 方框的大小
  graph.pos = "right",             # 森林图在右侧显示
  hrzl_lines = list(               # 水平线样式的设置
    "1" = gpar(lty = 1, lwd = 1),  # 均值线
    "2" = gpar(lty = 2),           # 下限和上限之间的虚线
    "14" = gpar(lwd = 2, lty = 0)# 下限和上限线
  ),xticks = c(0.6,0.8,1.0,1.2,1.4,1.6)
)

overall_sq_forest <- Prepare_forest_overall(overall_sq)

forestplot(
  overall_sq_forest[,c(1,7,6)],        # 需要显示在森林图中的列
  mean = overall_sq_forest[, 3],             # 均值列（HR），它将显示为森林图的小方块或其他形状哈哈哈哈哈
  lower = overall_sq_forest[, 4],            # 95%置信区间的下限数据列
  upper = overall_sq_forest[, 5],            # 95%置信区间的上限数据列
  zero = 1,                        # 均值为1时的参考线，也就是零线
  boxsize = 0.2,                   # 方框的大小
  graph.pos = "right",             # 森林图在右侧显示
  hrzl_lines = list(               # 水平线样式的设置
    "1" = gpar(lty = 1, lwd = 1),  # 均值线
    "2" = gpar(lty = 2),           # 下限和上限之间的虚线
    "14" = gpar(lwd = 2, lty = 0)# 下限和上限线
  ),xticks = c(0.8,0.9,1.0,1.1,1.2)
)

Prepare_forest_mediators <- function(x){
  data_promt <- x[,c(5,1,10:13)]
  CI <-  paste0(round(x$OR,4),'(',round(x$Lower.95..CI,4),'-',round(x$Upper.95..CI,4),')')
  data_promt$pval <- round(data_promt$pval,4)
  data_promt <- cbind(data_promt,CI)
  data_promt <- rbind(c('exposure','outcome',NA,NA,NA,'P.value','OR(95%CI)'),data_promt)
  return(data_promt)
}

mediators_forest <- Prepare_forest_mediators(mediators)

forestplot(
  mediators_forest[,c(1,2,7,6)],        # 需要显示在森林图中的列
  mean = mediators_forest[, 3],             # 均值列（HR），它将显示为森林图的小方块或其他形状哈哈哈哈哈
  lower = mediators_forest[, 4],            # 95%置信区间的下限数据列
  upper = mediators_forest[, 5],            # 95%置信区间的上限数据列
  zero = 1,                        # 均值为1时的参考线，也就是零线
  boxsize = 0.2,                   # 方框的大小
  graph.pos = "right",             # 森林图在右侧显示
  hrzl_lines = list(               # 水平线样式的设置
    "1" = gpar(lty = 1, lwd = 1),  # 均值线
    "2" = gpar(lty = 2),           # 下限和上限之间的虚线
    "14" = gpar(lwd = 2, lty = 0)# 下限和上限线
  ),xticks = c(0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
)

Prepare_forest_inf_lung <- function(x){
  data_promt <- x[,c(1,4,10:13)]
  CI <-  paste0(round(x$OR,4),'(',round(x$Lower.95..CI,4),'-',round(x$Upper.95..CI,4),')')
  data_promt$pval <- round(data_promt$pval,4)
  data_promt <- cbind(data_promt,CI)
  data_promt <- rbind(c('exposure','outcome',NA,NA,NA,'P.value','OR(95%CI)'),data_promt)
  return(data_promt)
}

imm_inf_la_forest <- Prepare_forest_inf_lung(imm_inf_la)

forestplot(
  imm_inf_la_forest[,c(1,2,7,6)],        # 需要显示在森林图中的列
  mean = imm_inf_la_forest[, 3],             # 均值列（HR），它将显示为森林图的小方块或其他形状哈哈哈哈哈
  lower = imm_inf_la_forest[, 4],            # 95%置信区间的下限数据列
  upper = imm_inf_la_forest[, 5],            # 95%置信区间的上限数据列
  zero = 1,                        # 均值为1时的参考线，也就是零线
  boxsize = 0.2,                   # 方框的大小
  graph.pos = "right",             # 森林图在右侧显示
  hrzl_lines = list(               # 水平线样式的设置
    "1" = gpar(lty = 1, lwd = 1),  # 均值线
    "2" = gpar(lty = 2),           # 下限和上限之间的虚线
    "14" = gpar(lwd = 2, lty = 0)# 下限和上限线
  ),xticks = c(0.6,0.8,1.0,1.2,1.4,1.6)
)

imm_inf_sc_forest <- Prepare_forest_inf_lung(imm_inf_sc)

forestplot(
  imm_inf_sc[,c(1,2,7,6)],        # 需要显示在森林图中的列
  mean = imm_inf_sc[, 3],             # 均值列（HR），它将显示为森林图的小方块或其他形状哈哈哈哈哈
  lower = imm_inf_sc[, 4],            # 95%置信区间的下限数据列
  upper = imm_inf_sc[, 5],            # 95%置信区间的上限数据列
  zero = 1,                        # 均值为1时的参考线，也就是零线
  boxsize = 0.2,                   # 方框的大小
  graph.pos = "right",             # 森林图在右侧显示
  hrzl_lines = list(               # 水平线样式的设置
    "1" = gpar(lty = 1, lwd = 1),  # 均值线
    "2" = gpar(lty = 2),           # 下限和上限之间的虚线
    "14" = gpar(lwd = 2, lty = 0)# 下限和上限线
  ),xticks = c(0.6,0.8,1.0,1.2,1.4,1.6)
)


imm_inf_sq_forest <- Prepare_forest_inf_lung(imm_inf_sq)
forestplot(
  imm_inf_sq_forest[,c(1,2,7,6)],        # 需要显示在森林图中的列
  mean = imm_inf_sq_forest[, 3],             # 均值列（HR），它将显示为森林图的小方块或其他形状哈哈哈哈哈
  lower = imm_inf_sq_forest[, 4],            # 95%置信区间的下限数据列
  upper = imm_inf_sq_forest[, 5],            # 95%置信区间的上限数据列
  zero = 1,                        # 均值为1时的参考线，也就是零线
  boxsize = 0.2,                   # 方框的大小
  graph.pos = "right",             # 森林图在右侧显示
  hrzl_lines = list(               # 水平线样式的设置
    "1" = gpar(lty = 1, lwd = 1),  # 均值线
    "2" = gpar(lty = 2),           # 下限和上限之间的虚线
    "14" = gpar(lwd = 2, lty = 0)# 下限和上限线
  ),xticks = c(0.6,0.8,1.0,1.2,1.4)
)


# calculate propotion
# LUAD
b1 <- imm_inf_la$b


b2 <- data.frame()
for(i in mediators[mediators$id.outcome == 'Lung adenocarcinoma',]$exposure){
  step1 <- imm_inf_la[imm_inf_la$id.outcome ==i,][,c('b','se')]
  if(dim(b2)[1]==0){
    b2 <- step1
  }else{
    b2 <- rbind(b2,step1)
  }
}

colnames(b2) <- c('b2','se2')


btotal <- data.frame()
for(i in imm_inf_la$Trait){
  step1 <- overall_la[overall_la$Trait ==i,][,c('b','se')]
  if(dim(btotal)[1]==0){
    btotal <- step1
  }else{
    btotal <- rbind(btotal,step1)
  }
}

colnames(btotal) <- c('b0','se0')

med_imm_inf_la <- cbind(imm_inf_la[,c(1,4)],btotal,imm_inf_la[,c(8,9)],b2)

colnames(med_imm_inf_la) <- c("Trait", "outcome","b0","se0","b1", "se1" , "b2","se2" )

med_effect <- b1*b2$b2
direct_effect <- btotal$b0 - med_effect
med_pro <- med_effect/btotal$b0

med_imm_inf_la$direct_effect = direct_effect
med_imm_inf_la$mediation_effect = med_effect
med_imm_inf_la$'med_proportion(%)' = med_pro*100
med_imm_inf_la$'med_proportion(%)' = ifelse(med_imm_inf_la$'med_proportion(%)' >0, med_imm_inf_la$'med_proportion(%)',NA)

result <- data.frame()
for(i in 1:59){
  res <- RMediation::medci(mu.x = med_imm_inf_la$b1[i],mu.y = med_imm_inf_la$b2[i],type = 'asymp',
                           se.x = med_imm_inf_la$se1[i],se.y = med_imm_inf_la$se2[i],
                           rho=0, alpha=.05,plot=FALSE, plotCI=FALSE)
  `Lower95% CI` <- res$`95% CI`[1]
  `Upper95% CI` <- res$`95% CI`[2]
  q <- data.frame(`Lower95% CI` = `Lower95% CI`, `Upper95% CI` = `Upper95% CI`, Estimate = res$Estimate, se = res$SE)
  if(dim(result)[1]== 0){
    result <- q
  }else{
    result <- rbind(result,q)
  }
}

med_imm_inf_la <- cbind(med_imm_inf_la,result[,c(1,2)])

library(ggplot2)
library(ggsankey)
library(cols4all)
library(tidyverse)
library(dittoSeq)
data <- med_imm_inf_la[,c(1,2)]
data2 <- data %>% make_long(Trait,outcome)
data2$node <- factor(data2$node,levels = c(rev(unique(data$Trait)),rev(unique(data$outcome))))

ggplot(data2, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                fill = node, label = node)) +
  #设置
  geom_sankey(flow.fill="#DFDFDF",#连线颜色
              flow.alpha = 0.5, ## 条带透明度
              flow.color="grey60",#连线边框颜色
              #node.fill=dittoColors()[1:36],#节点颜色，[1:36]数值需要根据自己的数据进行修改
              width=0.2) + #node的宽度
  #设置桑葚图文字
  geom_sankey_text(size = 3,#文字大小
                   color= "black",#文字颜色
                   hjust=1) + #文字位置，右对齐
  theme_void()+
  #隐藏图例
  theme(legend.position = 'none') 


# SC 

b1 <- imm_inf_sc$b


b2 <- data.frame()
for(i in mediators[mediators$id.outcome == 'small cell lung carcinoma',]$exposure){
  step1 <- imm_inf_sc[imm_inf_sc$id.outcome ==i,][,c('b','se')]
  if(dim(b2)[1]==0){
    b2 <- step1
  }else{
    b2 <- rbind(b2,step1)
  }
}

colnames(b2) <- c('b2','se2')


btotal <- data.frame()
for(i in imm_inf_sc$Trait){
  step1 <- overall_sc[overall_sc$Trait ==i,][,c('b','se')]
  if(dim(btotal)[1]==0){
    btotal <- step1
  }else{
    btotal <- rbind(btotal,step1)
  }
}

colnames(btotal) <- c('b0','se0')

med_imm_inf_sc <- cbind(imm_inf_sc[,c(1,4)],btotal,imm_inf_sc[,c(8,9)],b2)

colnames(med_imm_inf_sc) <- c("Trait", "outcome","b0","se0","b1", "se1" , "b2","se2" )

med_effect <- b1*b2$b2
direct_effect <- btotal$b0 - med_effect
med_pro <- med_effect/btotal$b0

med_imm_inf_sc$direct_effect = direct_effect
med_imm_inf_sc$mediation_effect = med_effect
med_imm_inf_sc$'med_proportion(%)' = med_pro*100
med_imm_inf_sc$'med_proportion(%)' = ifelse(med_imm_inf_sc$'med_proportion(%)' >0, med_imm_inf_sc$'med_proportion(%)',NA)

result <- data.frame()
for(i in 1:68){
  res <- RMediation::medci(mu.x = med_imm_inf_sc$b1[i],mu.y = med_imm_inf_sc$b2[i],type = 'asymp',
                           se.x = med_imm_inf_sc$se1[i],se.y = med_imm_inf_sc$se2[i],
                           rho=0, alpha=.05,plot=FALSE, plotCI=FALSE)
  `Lower95% CI` <- res$`95% CI`[1]
  `Upper95% CI` <- res$`95% CI`[2]
  q <- data.frame(`Lower95% CI` = `Lower95% CI`, `Upper95% CI` = `Upper95% CI`, Estimate = res$Estimate, se = res$SE)
  if(dim(result)[1]== 0){
    result <- q
  }else{
    result <- rbind(result,q)
  }
}

med_imm_inf_sc <- cbind(med_imm_inf_sc,result[,c(1,2)])

library(ggplot2)
library(ggsankey)
library(cols4all)
library(tidyverse)
library(dittoSeq)
data <- med_imm_inf_sc[,c(1,2)]
data2 <- data %>% make_long(Trait,outcome)
data2$node <- factor(data2$node,levels = c(rev(unique(data$Trait)),rev(unique(data$outcome))))

ggplot(data2, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                  fill = node, label = node)) +
  #设置
  geom_sankey(flow.fill="#DFDFDF",#连线颜色
              flow.alpha = 0.5, ## 条带透明度
              flow.color="grey60",#连线边框颜色
              #node.fill=dittoColors()[1:36],#节点颜色，[1:36]数值需要根据自己的数据进行修改
              width=0.2) + #node的宽度
  #设置桑葚图文字
  geom_sankey_text(size = 3,#文字大小
                   color= "black",#文字颜色
                   hjust=1) + #文字位置，右对齐
  theme_void()+
  #隐藏图例
  theme(legend.position = 'none') 


# sq

b1 <- imm_inf_sq$b


b2 <- data.frame()
for(i in mediators[mediators$id.outcome == 'squamous cell lung cancer',]$exposure){
  step1 <- imm_inf_sq[imm_inf_sq$id.outcome ==i,][,c('b','se')]
  if(dim(b2)[1]==0){
    b2 <- step1
  }else{
    b2 <- rbind(b2,step1)
  }
}

colnames(b2) <- c('b2','se2')

sq <- data.frame()
for(i in imm_inf_sq$id.exposure){
  step1 <- overall_sq[overall_sq$id.exposure == i,]
  if(dim(sq)[1]==0){
    sq = step1
  }else{
    sq = rbind(sq,step1)
  }
}




btotal <- sq[,c(8,9)]


colnames(btotal) <- c('b0','se0')

med_imm_inf_sq <- cbind(imm_inf_sq[,c(1,4)],btotal,imm_inf_sq[,c(8,9)],b2)

colnames(med_imm_inf_sq) <- c("Trait", "outcome","b0","se0","b1", "se1" , "b2","se2" )

med_effect <- b1*b2$b2
direct_effect <- btotal$b0 - med_effect
med_pro <- med_effect/btotal$b0

med_imm_inf_sq$direct_effect = direct_effect
med_imm_inf_sq$mediation_effect = med_effect
med_imm_inf_sq$'med_proportion(%)' = med_pro*100
med_imm_inf_sq$'med_proportion(%)' = ifelse(med_imm_inf_sq$'med_proportion(%)' >0, med_imm_inf_sq$'med_proportion(%)',NA)

result <- data.frame()
for(i in 1:36){
  res <- RMediation::medci(mu.x = med_imm_inf_sq$b1[i],mu.y = med_imm_inf_sq$b2[i],type = 'asymp',
                           se.x = med_imm_inf_sq$se1[i],se.y = med_imm_inf_sq$se2[i],
                           rho=0, alpha=.05,plot=FALSE, plotCI=FALSE)
  `Lower95% CI` <- res$`95% CI`[1]
  `Upper95% CI` <- res$`95% CI`[2]
  q <- data.frame(`Lower95% CI` = `Lower95% CI`, `Upper95% CI` = `Upper95% CI`, Estimate = res$Estimate, se = res$SE)
  if(dim(result)[1]== 0){
    result <- q
  }else{
    result <- rbind(result,q)
  }
}

med_imm_inf_sq <- cbind(med_imm_inf_sq,result[,c(1,2)])

library(ggplot2)
library(ggsankey)
library(cols4all)
library(tidyverse)
library(dittoSeq)
data <- med_imm_inf_sq[,c(1,2)]
data2 <- data %>% make_long(Trait,outcome)
data2$node <- factor(data2$node,levels = c(rev(unique(data$Trait)),rev(unique(data$outcome))))

ggplot(data2, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                  fill = node, label = node)) +
  #设置
  geom_sankey(flow.fill="#DFDFDF",#连线颜色
              flow.alpha = 0.5, ## 条带透明度
              flow.color="grey60",#连线边框颜色
              #node.fill=dittoColors()[1:36],#节点颜色，[1:36]数值需要根据自己的数据进行修改
              width=0.2) + #node的宽度
  #设置桑葚图文字
  geom_sankey_text(size = 3,#文字大小
                   color= "black",#文字颜色
                   hjust=1) + #文字位置，右对齐
  theme_void()+
  #隐藏图例
  theme(legend.position = 'none') 

med_analysis_result <- list(med_imm_inf_la = med_imm_inf_la,med_imm_inf_sc = med_imm_inf_sc, med_imm_inf_sq = med_imm_inf_sq)
save(med_analysis_result,file = 'med_analysis_result.Rdata')

#### mvMR

do_mvMR <- function(x = 'med' ,y = 'overall'){
  med <- as.data.frame(table(x$outcome))[,1]
  for(i in med){
    exposure <- x[x$outcome == i,]$Trait
    id.expoure <- y[y$Trait %in% exposure,]$id.exposure
    for(j in id.expoure){
      path_promt <- paste0("/Users/llls2012163.com/Lung cancer/Harmnoise_immToCancer/ImmToLungAdenocarinoma/",j,'.txt')
    }
  }
  
  
}

write.table(med_imm_inf_la,file = "/Users/llls2012163.com/med_imm_inf_la.txt",row.names = F,sep = '\t',quote = F,append = FALSE,col.names = TRUE)
write.table(med_imm_inf_sc,file = "/Users/llls2012163.com/med_imm_inf_sc.txt",row.names = F,sep = '\t',quote = F,append = FALSE,col.names = TRUE)
write.table(med_imm_inf_sq,file = "/Users/llls2012163.com/med_imm_inf_sq.txt",row.names = F,sep = '\t',quote = F,append = FALSE,col.names = TRUE)


### sensitive analysis data summary
# 
Result_storge_dir <- "/Users/llls2012163.com/Lung cancer/result_file_immTocancerType"
MR_result_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType")[5]
MR_het_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType")[2]
MR_presso_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType")[4]
MR_ple_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/result_file_immTocancerType")[3]
MR_file <- list.files(MR_result_dir)


Sensitive_analysis <- function(x = MR_file){
  MR_sig <- data.frame()
  for(i in x){
    print(i)
    #MR_result_path <- paste0(MR_result_dir,'/',i)
    #MR_result_promt <- read.delim(MR_result_path)
    
    path_prmot_het <- paste0(MR_het_dir,'/',i)
    data_prmot_het <- read.delim(path_prmot_het)
    data_prmot_het <- data_prmot_het[data_prmot_het$method == 'Inverse variance weighted',]
    
    a <- data_prmot_het
    # dplyr::inner_join(MR_result_promt,data_prmot_het,by='id.exposure')
    
    
    path_promt_pleiotropy <- paste0(MR_ple_dir,'/',i)
    data_promt_pleiotropy <- read.delim(path_promt_pleiotropy)
    
    b <-  dplyr::inner_join(a,data_promt_pleiotropy,by='id.exposure')
    
    path_promt_presso <- paste0(MR_presso_dir,'/',i)
    data_promt_presso <- read.delim(path_promt_presso)
    #data_promt_presso <- data_promt_presso[1,]
    
    c <- dplyr::inner_join(b,data_promt_presso,by='id.exposure')
    if(dim(MR_sig)[1]==0){
      MR_sig <- c
    }else{
      MR_sig <- rbind(MR_sig,c)
    }
    
  }
  MR_sig <- MR_sig[,c(1:3,6:8,12:14,16:23)]
  MR_sig$outcome.x <- limma::strsplit2(MR_sig$outcome.x,'\\|\\|')[,1]
  return(MR_sig)
}



a <- Sensitive_analysis(MR_file)

write.table(a,file="/Users/llls2012163.com/a.txt",sep = '\t',quote = FALSE,row.names = FALSE)


Result_storge_dir <- "/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung"
MR_result_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung")[5]
MR_het_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung")[2]
MR_presso_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung")[4]
MR_ple_dir <- list.dirs("/Users/llls2012163.com/Lung cancer/inf_Lung/result/infToLung")[3]
MR_file <- list.files(MR_result_dir)

Sensitive_analysis <- function(x = MR_file){
  MR_sig <- data.frame()
  for(i in x){
    print(i)
    #MR_result_path <- paste0(MR_result_dir,'/',i)
    #MR_result_promt <- read.delim(MR_result_path)
    
    path_prmot_het <- paste0(MR_het_dir,'/',i)
    data_prmot_het <- read.delim(path_prmot_het)
    data_prmot_het <- data_prmot_het[data_prmot_het$method == 'Inverse variance weighted',]
    
    a <- data_prmot_het
    # dplyr::inner_join(MR_result_promt,data_prmot_het,by='id.exposure')
    
    
    path_promt_pleiotropy <- paste0(MR_ple_dir,'/',i)
    data_promt_pleiotropy <- read.delim(path_promt_pleiotropy)
    
    b <-  dplyr::inner_join(a,data_promt_pleiotropy,by='id.exposure')
    
    path_promt_presso <- paste0(MR_presso_dir,'/',i)
    data_promt_presso <- read.delim(path_promt_presso)
    #data_promt_presso <- data_promt_presso[1,]
    
    c <- dplyr::inner_join(b,data_promt_presso,by='id.exposure')
    if(dim(MR_sig)[1]==0){
      MR_sig <- c
    }else{
      MR_sig <- rbind(MR_sig,c)
    }
    
  }
  # MR_sig$outcome.x <- limma::strsplit2(MR_sig$outcome.x,'\\|\\|')[,1]
  MR_sig <- MR_sig[,c(1:4,7:9,14:16,19:25)]
  return(MR_sig)
}

b <- Sensitive_analysis(MR_file)
write.table(b,file="/Users/llls2012163.com/b.txt",sep = '\t',quote = FALSE,row.names = FALSE)


Result_storge_dir <- "/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf"
MR_result_dir <- list.dirs("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf")[5]
MR_het_dir <- list.dirs("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf")[2]
MR_presso_dir <- list.dirs("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf")[4]
MR_ple_dir <- list.dirs("/Users/llls2012163.com/Circulating inflammatory cytokines GWAS/Results/imm_inf")[3]
MR_file <- list.files(MR_result_dir)

Sensitive_analysis <- function(x = MR_file){
  MR_sig <- data.frame()
  for(i in x){
    print(i)
    #MR_result_path <- paste0(MR_result_dir,'/',i)
    #MR_result_promt <- read.delim(MR_result_path)
    
    path_prmot_het <- paste0(MR_het_dir,'/',i)
    data_prmot_het <- read.delim(path_prmot_het)
    data_prmot_het <- data_prmot_het[data_prmot_het$method == 'Inverse variance weighted',]
    
    a <- data_prmot_het
    # dplyr::inner_join(MR_result_promt,data_prmot_het,by='id.exposure')
    
    
    path_promt_pleiotropy <- paste0(MR_ple_dir,'/',i)
    data_promt_pleiotropy <- read.delim(path_promt_pleiotropy)
    
    b <-  dplyr::inner_join(a,data_promt_pleiotropy,by='id.exposure')
    
    path_promt_presso <- paste0(MR_presso_dir,'/',i)
    data_promt_presso <- read.delim(path_promt_presso)
    #data_promt_presso <- data_promt_presso[1,]
    
    c <- dplyr::inner_join(b,data_promt_presso,by='id.exposure')
    if(dim(MR_sig)[1]==0){
      MR_sig <- c
    }else{
      MR_sig <- rbind(MR_sig,c)
    }
    
  }
  #MR_sig$outcome.x <- limma::strsplit2(MR_sig$outcome.x,'\\|\\|')[,1]
  MR_sig <- MR_sig[,c(1,2,6:8,12:14,16:22)]
  return(MR_sig)
}
c <- Sensitive_analysis(MR_file)
write.table(c,file="/Users/llls2012163.com/c.txt",sep = '\t',quote = FALSE,row.names = FALSE)






