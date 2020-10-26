# CID-S model free analysis , pre-processing of data
# created NL - 10-01-2019



# Set up ----------------------------------------------------------------------
rm(list=ls()) # Clears everything from environment

# Below: function to load necessary packages and install if necessary
init <- function(need) {
  ip <- .packages(all.available = T)
  if (any((need %in% ip) == F)) {
    install.packages(need[!(need %in% ip)])
  }
  ok <- sapply(1:length(need), function(p) require(need[[p]], 
                                                   character.only = T))
}

init(c("plyr","ggplot2","reshape2","car","afex","lme4","nloptr","grid","gtable",
       "cowplot","nortest","tidyverse"))


# Experiment to pre-process ---------------------------------------------------
Experiment <- "CIDRS22" 
# CIDRS14: Experiment 1 in the paper
# CIDRS19: Experiment 2 in the paper
# CIDRS21: Experiment 3 in the paper
# CIDRS22: Experiment 4 in the paper


setwd("./..")
getwd()

# Folder where individual data files are stored
datapath<-paste0("Data/",Experiment)


# Load test files from all participants --------------------------------------

if (Experiment %in% c("CIDRS19","CIDRS21","CIDRS14")){
  
  load_memory <- function(path) { 
    files <- dir(path, pattern = 'dataTEST',full.names = TRUE,recursive = TRUE) 
    # search for test data (dataTEST) file in all subdirectories (recursive) of path
    tables <- lapply(files, read.csv,header=FALSE)
    do.call(rbind, tables)
    # collate all of them into 1 dataframe
  }
  
  load_study <- function(path) { 
    files <- dir(path, pattern = '_ID_CHECK',full.names = TRUE,recursive = TRUE) 
    # search for test data (dataTEST) file in all subdirectories (recursive) of path
    tables <- lapply(files, read.csv,header=FALSE)
    do.call(rbind, tables)
    # collate all of them into 1 dataframe
  }
  
  
  # PREPROCESSING/EXCLUSION --------------------------------------------------
  # Exclude demo (subject 0) and invalid trials (testtrial -99)
  
  memdata<-load_memory(datapath) 
  colnames(memdata)<-c("TestTrial","Block","StimType","StimID","CIDCorrect","idRT",
                       "RecJud","recRT","SDTclass","RecConf","SourceJud","SourceRT",
                       "SourceCorr","SourceConf","Bottom1Top2New0","SubjID","SubjAge",
                       "M1F2","TotExpDur")
  
  # pre-analysis exclusion
  # Exclude demo (subject 0) and invalid trials (testtrial -99)
  
  # Set up Exclusion Column
  memdata$ExclReason<-"valid"
  memdata<-memdata[!memdata$SubjID==0,] #demo
  memdata[memdata$CIDCorrect=="NaN",]$ExclReason<-"NoResponse" # messed up trial
  memdata[memdata$CIDCorrect=="0",]$ExclReason<-"CIDincorrect" # Incorrect identification
  
  #Exclude too fast (< 200 ms) and too slow (M + 3SD) identification responses
  memdata_cleaned<-NULL
  for(i in unique(memdata$SubjID)){
    
    rawdata_sub<-memdata[memdata$SubjID==i,]
    mintime<-200
    maxtime<-mean(rawdata_sub[rawdata_sub$CIDCorrect=="1" & rawdata_sub$idRT>mintime,]$idRT)+
      3*sd(rawdata_sub[rawdata_sub$CIDCorrect=="1" & rawdata_sub$idRT>mintime,]$idRT)
    
    if (nrow(rawdata_sub[rawdata_sub$idRT<mintime & 
                         rawdata_sub$idRT!="-99" & 
                         rawdata_sub$CIDCorrect=="1",])>=1) {
      rawdata_sub[rawdata_sub$idRT<mintime & rawdata_sub$idRT!="-99" & rawdata_sub$CIDCorrect=="1",]$ExclReason<-"TooFast"
    }
    if(nrow(rawdata_sub[rawdata_sub$idRT>maxtime & rawdata_sub$CIDCorrect=="1",]))
    {
      rawdata_sub[rawdata_sub$idRT>maxtime & rawdata_sub$CIDCorrect=="1",]$ExclReason<-"TooSlow"
    }
    
    memdata_cleaned<-rbind(memdata_cleaned,rawdata_sub)
  }
  
  
  idwords<-load_study(datapath)
  
  # add stimulus word and participant responses from study file
  memdata_cleaned$item<-idwords$V1
  memdata_cleaned$item_presponse<-idwords$V2
} else {
  
  
  load_memory <- function(path) { 
    files <- dir(path, pattern = 'dataTESTJUD',full.names = TRUE,recursive = TRUE) 
    # search for test data (dataTEST) file in all subdirectories (recursive) of path
    tables <- lapply(files, read.csv,header=FALSE)
    do.call(rbind, tables)
    # collate all of them into 1 dataframe
  }
  
  load_idwords <- function(path) { 
    files <- dir(path, pattern = '_ID_CHECK',full.names = TRUE,recursive = TRUE) 
    # search for test data (dataTEST) file in all subdirectories (recursive) of path
    tables <- lapply(files, read.csv,header=FALSE)
    do.call(rbind, tables)
    # collate all of them into 1 dataframe
  }
  
  load_idnumbers <- function(path) { 
    files <- dir(path, pattern = 'dataTESTID',full.names = TRUE,recursive = TRUE) 
    # search for test data (dataTEST) file in all subdirectories (recursive) of path
    tables <- lapply(files, read.csv,header=FALSE)
    do.call(rbind, tables)
    # collate all of them into 1 dataframe
  }
  
  
  
  
  
  # PREPROCESSING/EXCLUSION --------------------------------------------------
  # Exclude demo (subject 0) and invalid trials (testtrial -99)
  
  memdata<-load_memory(datapath) 
  
  headerpath <- paste0(datapath,"/", unique(memdata$V10)[[1]],"/TestJUDHeads.csv")
  cat("\n", file = headerpath, append = TRUE)
  colnames(memdata)<-unname(unlist(read.csv(paste0(datapath,"/", unique(memdata$V10)[[1]],"/TestJUDHeads.csv"),header = FALSE))) 
  
  
  idnumbers<-load_idnumbers(datapath)
  idwords<-load_idwords(datapath)
  idwords$V4<-idnumbers$V4
  idwords$SubjID<-rep(unique(idnumbers$V8),each=nrow(idwords)/length(unique(idnumbers$V8)))
  
  mem_idmem<-NULL
  for (i in unique(memdata$SubjID)){
    
    onepart<-idwords[idwords$SubjID==i,]
    idRTtimes <- idnumbers[idnumbers$V8 == i,]
    mempart<-memdata[memdata$SubjID==i,]
    mempart_ordered<-mempart[order(mempart$StimID),]
    onepart_ordered<-onepart[order(onepart$V4),]
    idRTtimes_ordered<-idRTtimes[order(idRTtimes$V4),]
    mempart_ordered$CIDCorrect<-onepart_ordered$V3
    mempart_ordered$idRT <- idRTtimes_ordered$V6
    mempart_ordered$item<-onepart_ordered$V1
    mempart_ordered$item_presponse<-onepart_ordered$V2
    mem_idmem<-rbind(mem_idmem,mempart_ordered)
    
  }
  
  memdata <- mem_idmem
  # pre-analysis exclusion
  # Exclude demo (subject 0) and invalid trials (testtrial -99)
  memdata<-memdata[!memdata$SubjID %in% c(0,99),] 
  memdata$ExclReason<-"valid"
  if (nrow(memdata[memdata$CIDCorrect=="NaN",])>0){
    memdata[memdata$CIDCorrect=="NaN",]$ExclReason<-"NoResponse"
  }
  memdata[memdata$CIDCorrect=="0",]$ExclReason<-"CIDincorrect"
  
  memdata_cleaned<-NULL
  for(i in unique(memdata$SubjID)){
    
    rawdata_sub<-memdata[memdata$SubjID==i,]
    mintime<-200
    maxtime<-mean(rawdata_sub[rawdata_sub$CIDCorrect=="1" & rawdata_sub$idRT>mintime,]$idRT)+
      3*sd(rawdata_sub[rawdata_sub$CIDCorrect=="1" & rawdata_sub$idRT>mintime,]$idRT)
    
    if (nrow(rawdata_sub[rawdata_sub$idRT<mintime & 
                         rawdata_sub$idRT!="-99" & 
                         rawdata_sub$CIDCorrect=="1",])>=1) {
      rawdata_sub[rawdata_sub$idRT<mintime & rawdata_sub$idRT!="-99" & rawdata_sub$CIDCorrect=="1",]$ExclReason<-"TooFast"
    }
    if(nrow(rawdata_sub[rawdata_sub$idRT>maxtime & rawdata_sub$CIDCorrect=="1",]))
    {
      rawdata_sub[rawdata_sub$idRT>maxtime & rawdata_sub$CIDCorrect=="1",]$ExclReason<-"TooSlow"
    }
    
    memdata_cleaned<-rbind(memdata_cleaned,rawdata_sub)
  }
  
  
  memdata_cleaned<-memdata_cleaned[order(memdata_cleaned$SubjID,memdata_cleaned$TestTrial),]
  
}

# Demographics ----
# Demographics ----------------------------------------------------------------------


memdata_cleaned %>%
  group_by(SubjID) %>%
  summarize(age = unique(SubjAge),
            gender=unique(M1F2)) %>%
  summarize(NumberParticipants = length(gender),
            Male = length(which(gender=="1")),
            mAge = mean(age),
            sdAge = sd(age))

#A priori exclude participants


if (Experiment %in% "CIDRS19"){
  excludeppt <- c()
} else if (Experiment %in% "CIDRS21"){
  excludeppt <- c()
} else if (Experiment %in% "CIDRS22"){
  excludeppt <- c(20) #texting during experiment
} else if (Experiment %in% "CIDRS14"){
  excludeppt <- c(33) #block1 no identification attempt
}

memdata_cleaned <- memdata_cleaned[!memdata_cleaned$SubjID %in% excludeppt,]


# Export single file with exclusion column and stimulus/responses --------------------
setwd("./Standard Analysis/")
write.csv(memdata_cleaned,paste0(Experiment, "_fulldatascreened.csv"),row.names=FALSE)


# Preprocessing report in paper -------------------------------------------------------


# Number a-priori excluded participants

length(excludeppt)


# Percentage correct
# Screening
memdata_cleaned$SubjID <- as.factor(as.character(memdata_cleaned$SubjID))
memdata_cleaned$ExclReason<- as.factor(as.character(memdata_cleaned$ExclReason))

screening<-ddply(memdata_cleaned,.(SubjID,ExclReason),summarize,
                 number = length(ExclReason)/(nrow(memdata_cleaned)/length(unique(memdata_cleaned$SubjID))))

emptydf<-with(screening,expand.grid(SubjID=levels(SubjID),
                                    ExclReason=levels(ExclReason)))

screeningf<-merge(screening,emptydf,all.y=T)
screeningf$number[is.na(screeningf$number)]<-0
total<-ddply(screeningf,.(ExclReason),summarize,
             means = mean(number),
             sds = sd(number))
total
min(screening[screening$ExclReason == "valid",]$number)

  # Percentage too fast/ too slow
  memdata_cleaned$additional<-2
  memdata_cleaned[memdata_cleaned$ExclReason=="TooFast" | memdata_cleaned$ExclReason=="TooSlow",]$additional<-1
  
  screening<-ddply(memdata_cleaned,.(SubjID,additional),summarize,
                   number = length(ExclReason)/(nrow(memdata_cleaned)/length(unique(memdata_cleaned$SubjID))))
  
  total<-ddply(screening,.(additional),summarize,
               means = mean(number),
               sds = sd(number))
