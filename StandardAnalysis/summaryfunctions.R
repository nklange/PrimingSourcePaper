init <- function(need) {
  ip <- .packages(all.available = T)
  if (any((need %in% ip) == F)) {
    install.packages(need[!(need %in% ip)])
  }
  ok <- sapply(1:length(need), function(p) require(need[[p]], 
                                                   character.only = T))
}

serror <- function(x) sd(x,na.rm=T)/sqrt(length(x[!is.na(x)]))
serror_ind <- function(x,y) sqrt((sd(x)^2/length(x)) + (sd(y)^2/length(y)))

cohensdz<-function(ttest,n){
  dz<-ttest$statistic/sqrt(n)
  return(dz)
}

cohensd<-function(ttest,n1,n2){
  d<-ttest$statistic*(n1+n2)/(sqrt(ttest$parameter)*sqrt(n1*n2))
  return(d)
}

ttestAgainstZero <- function(data){
  
  result<-t.test(data,mu=0) 
  BFt <- ttestBF(data,mu=0)
  BFt <- extractBF(BFt)
  
  data.frame(Mean = result$estimate,
             SE = serror(data),
             df = result$parameter,
             tval = result$statistic,
             pval = result$p.value,
             cohensdz = cohensdz(result,result$parameter[[1]] + 1)[[1]],
             BF = BFt$bf,
             order =  (nchar(trunc(BFt$bf)) - 1))
}


ttestTwoCond<- function(data1,data2, paired = TRUE){
  
  result<-t.test(data1,data2,paired = paired) 
  BFt <- ttestBF(data1,data2,paired = paired)
  BFt <- extractBF(BFt)
  
  if (paired == TRUE) {
    SEs <- serror(data1 - data2)
    CD <- cohensdz(result,result$parameter[[1]] + 1)[[1]]
  } else {
    SEs <- serror_ind(data1,data2)
    CD <- cohensd(result,nrow(data1),nrow(data2))
  }
  data.frame(Mean = result$estimate,
             SE = SEs,
             df = result$parameter,
             tval = result$statistic,
             pval = result$p.value,
             cohensd = CD,
             BF = BFt$bf,
             order =  (nchar(trunc(BFt$bf)) - 1))
}


# Summary Functions ------------------------------------------------------------------

# Summary for graphs with within-subjects confidence intervals from R cookbook
## Summarizes data.
## Gives count, SUM, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE2 <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                       conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                      # mean = mean   (xx[,col], na.rm=na.rm), 
                      sum = sum (xx[,col], na.rm=na.rm), #change mean to sum in this case
                      sd   = sd     (xx[,col], na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("sum"=measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}


## Summarizes data.
## Gives count, MEAN, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                      mean = mean   (xx[,col], na.rm=na.rm), 
                      sd   = sd     (xx[,col], na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean"=measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

normDataWithin <- function(data=NULL, idvar, measurevar,
                           betweenvars=NULL, na.rm=FALSE, .drop=TRUE) {
  
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], 
                                             na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}


## Summarizes data, handling within-subjects variables by removing 
## inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same
## between-group mean),
##   standard deviation, standard error of the mean, and confidence
## interval.
## If there are within-subject variables, calculate adjusted values 
## using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to 
## be summariezed
##   betweenvars: a vector containing names of columns that are 
## between-subjects variables
##   withinvars: a vector containing names of columns that are
## within-subjects variables
##   idvar: the name of a column that identifies each subject 
## (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval 
## (default is 95%)

summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL,
                            withinvars=NULL,idvar=NULL, na.rm=FALSE, conf.interval=.95, 
                            .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), 
                            drop=FALSE], FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following 
            non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  
  datac <- summarySE2(data, measurevar, groupvars=c(betweenvars,
                                                    withinvars),na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated 
  # with normed data)
  
  
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  
  ndata <- normDataWithin(data, idvar, measurevar, 
                          betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and 
  # nwithin vars the same
  
  ndatac <- summarySE2(ndata, measurevar_n, 
                       groupvars=c(betweenvars, withinvars),
                       na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error 
  # and confidence interval
  #  Get the product of the number of conditions of 
  # within-S variables
  
  nWithinGroups <- prod(vapply(ndatac[,withinvars, drop=FALSE],
                               FUN=nlevels, FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

# ROC curves ---------------------------------------------------------------------------
# source ROC

makesourceROC<-function(alldata){
  fullrocrating<-NULL
  
  for (recrating in sort(unique(alldata$RecConf),decreasing=T)){
    
    forroc<-alldata[alldata$RecConf==recrating,]
    
    rocrating<-NULL
    for (rating in sort(unique(alldata$SourceConf),decreasing=T)){  
      
      hit<-sum(forroc[forroc$SourceInput=="Top" &
                        forroc$SourceConf>=rating,]$items)/
        sum(forroc[forroc$SourceInput=="Top",]$items)
      
      fa<-sum(forroc[forroc$SourceInput=="Bottom" & 
                       forroc$SourceConf>=rating,]$items)/
        sum(forroc[forroc$SourceInput=="Bottom",]$items)
      
      adjhit<-sum(forroc[forroc$SourceInput=="Top" &
                           forroc$SourceConf>=rating,]$adjnumber)/
        sum(forroc[forroc$SourceInput=="Top",]$adjnumber)
      
      adjfa<-sum(forroc[forroc$SourceInput=="Bottom" & 
                          forroc$SourceConf>=rating,]$adjnumber)/
        sum(forroc[forroc$SourceInput=="Bottom",]$adjnumber)
      
      pair<-cbind(hit,fa,adjhit,adjfa)
      rocrating<-rbind(rocrating,pair)
    }
    
    pretty<-cbind(recrating,rocrating)
    fullrocrating<-rbind(fullrocrating,pretty)
    fullrocrating<-as.data.frame(fullrocrating)
  }
  
  
  acrossrocrating<-NULL
  acrosssource<-ddply(alldata,.(SourceInput,SourceConf),summarise,
                      sumitems=sum(items))
  acrosssource$adjnumber<-acrosssource$sumitems+0.5
  
  rocrating<-NULL
  for (rating in sort(unique(acrosssource$SourceConf),decreasing=T)){  
    
    hit<-sum(acrosssource[acrosssource$SourceInput=="Top" &
                            acrosssource$SourceConf>=rating,]$sumitems)/
      sum(acrosssource[acrosssource$SourceInput=="Top",]$sumitems)
    
    fa<-sum(acrosssource[acrosssource$SourceInput=="Bottom" & 
                           acrosssource$SourceConf>=rating,]$sumitems)/
      sum(acrosssource[acrosssource$SourceInput=="Bottom",]$sumitems)
    
    adjhit<-sum(acrosssource[acrosssource$SourceInput=="Top" &
                               acrosssource$SourceConf>=rating,]$adjnumber)/
      sum(acrosssource[acrosssource$SourceInput=="Top",]$adjnumber)
    
    adjfa<-sum(acrosssource[acrosssource$SourceInput=="Bottom" & 
                              acrosssource$SourceConf>=rating,]$adjnumber)/
      sum(acrosssource[acrosssource$SourceInput=="Bottom",]$adjnumber)
    
    
    pair<-cbind(hit,fa,adjhit,adjfa)
    rocrating<-rbind(rocrating,pair)
  }
  
  rocrating<-as.data.frame(rocrating)
  rocrating$zhit<-qnorm(rocrating$adjhit)
  rocrating$zfa<-qnorm(rocrating$adjfa)
  rocrating$recrating<-"across"
  
  fullrocrating$zhit<-qnorm(fullrocrating$adjhit)
  fullrocrating$zfa<-qnorm(fullrocrating$adjfa)
  
  fullrocrating<-rbind(fullrocrating,rocrating)
  
  return(fullrocrating)
}


ggplot_source_split<-function(alldata){
  source_roc<-ggplot(alldata[alldata$recrating!="across",],aes(x=fa,y=hit)) +
    geom_point(size=5,aes(shape=recrating),fill="black") +
    theme_bw(base_size = 18)+ 
    theme(legend.position=c(0.75,0.26),
          legend.direction="vertical") +
    scale_shape_manual(name="Recognition",labels=c("Certain New","Probably New","Guess New","Guess Old",
                                                   "Probably Old","Certain Old"),values=c(0,1,2,5,6,22)) +
    scale_x_continuous(name="FA",limits=c(0,1))+
    scale_y_continuous(name="Hit",limits=c(0,1))+
    geom_abline(intercept=0,slope=1)+
    theme(plot.title=element_text(size=18, face="bold"),
          legend.key = element_blank(),
          axis.text.x  = element_text(size=12,colour="black"),
          axis.title.x = element_text(size=12),
          axis.text.y  = element_text(size=12,colour="black"),
          axis.title.y = element_text(size=12,vjust=1.2),
          legend.text = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title=element_text(size=14,colour="black"),
          panel.border = element_rect(colour="black"),
          strip.text.x = element_text(colour = "white", size = 14, face="bold"),
          strip.background = element_rect(fill="black",colour="black"),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA))
  
  source_zroc<-ggplot(alldata[alldata$recrating!="across",],aes(x=zfa,y=zhit)) +
    geom_point(size=5,aes(shape=recrating),fill="black") +
    theme_bw(base_size = 18)+ 
    theme(legend.position=c(0.8,0.26),
          legend.direction="vertical") +
    scale_x_continuous(name="zFA",limits=c(-3,3))+
    scale_y_continuous(name="zHit",limits=c(-3,3))+
    scale_shape_manual(name="Recognition",labels=c("Certain New","Probably New","Guess New","Guess Old",
                                                   "Probably Old","Certain Old"),values=c(0,1,2,5,6,22)) +
    theme_bw(base_size = 18)+ 
    theme(plot.title=element_text(size=18, face="bold"),
          axis.text.x  = element_text(size=12,colour="black"),
          axis.title.x = element_text(size=12),
          axis.text.y  = element_text(size=12,colour="black"),
          axis.title.y = element_text(size=12,vjust=1.2),
          legend.text = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title=element_text(size=14,colour="black"),
          legend.key=element_blank(),
          legend.background = element_rect(fill = "transparent"),
          panel.border = element_rect(colour="black"),
          strip.text.x = element_text(colour = "white", size = 14, face="bold"),
          strip.background = element_rect(fill="black",colour="black"),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA))
  
  plot_grid(source_roc, source_zroc, labels=c("A", "B"), ncol = 2, nrow = 1)
}


# Recog ROC


makerecogROC<-function(alldata){
  
  fullrocrating_rec<-NULL
  
  
  for (sourcerating in sort(unique(alldata$SourceConf),decreasing=T)){
    
    forroc<-alldata[alldata$SourceConf==sourcerating,]
    
    Trocrating<-NULL
    Brocrating<-NULL
    for (rating in sort(unique(alldata$RecConf),decreasing=T)){  
      
      Thit<-sum(forroc[forroc$SourceInput=="Top" &
                         forroc$RecConf>=rating,]$items)/
        sum(forroc[forroc$SourceInput=="Top",]$items)
      
      Tfa<-sum(forroc[forroc$SourceInput=="New" & 
                        forroc$RecConf>=rating,]$items)/
        sum(forroc[forroc$SourceInput=="New",]$items)
      
      Bhit<-sum(forroc[forroc$SourceInput=="Bottom" &
                         forroc$RecConf>=rating,]$items)/
        sum(forroc[forroc$SourceInput=="Bottom",]$items)
      
      Bfa<-sum(forroc[forroc$SourceInput=="New" & 
                        forroc$RecConf>=rating,]$items)/
        sum(forroc[forroc$SourceInput=="New",]$items)
      
      Tpair<-cbind("Top",Thit,Tfa)
      Bpair<-cbind("Bottom",Bhit,Bfa)
      Trocrating<-rbind(Trocrating,Tpair)
      Brocrating<-rbind(Brocrating,Bpair)
    }
    rocrating<-rbind(Trocrating,Brocrating)
    pretty<-cbind(sourcerating,rocrating)
    fullrocrating_rec<-rbind(fullrocrating_rec,pretty)
    fullrocrating_rec<-as.data.frame(fullrocrating_rec)
  }
  
  
  
  acrosssource<-ddply(alldata,.(SourceInput,RecConf),summarise,
                      sumitems=sum(items))
  emptydf <- with(acrosssource, expand.grid(SourceInput = levels(SourceInput), 
                                            RecConf = levels(RecConf)))
  Trocrating<-NULL
  Brocrating<-NULL
  
  for (rating in sort(unique(acrosssource$RecConf),decreasing=T)){  
    
    Thit<-sum(acrosssource[acrosssource$SourceInput=="Top" &
                             acrosssource$RecConf>=rating,]$sumitems)/
      sum(acrosssource[acrosssource$SourceInput=="Top",]$sumitems)
    
    Bhit<-sum(acrosssource[acrosssource$SourceInput=="Bottom" & 
                             acrosssource$RecConf>=rating,]$sumitems)/
      sum(acrosssource[acrosssource$SourceInput=="Bottom",]$sumitems)
    
    Tfa<-sum(acrosssource[acrosssource$SourceInput=="New" & 
                            acrosssource$RecConf>=rating,]$sumitems)/
      sum(acrosssource[acrosssource$SourceInput=="New",]$sumitems)
    
    Bfa<-sum(acrosssource[acrosssource$SourceInput=="New" & 
                            acrosssource$RecConf>=rating,]$sumitems)/
      sum(acrosssource[acrosssource$SourceInput=="New",]$sumitems)
    
    Tpair<-cbind("Top",Thit,Tfa)
    Bpair<-cbind("Bottom",Bhit,Bfa)
    Trocrating<-rbind(Trocrating,Tpair)
    Brocrating<-rbind(Brocrating,Bpair)
  }
  rocrating<-rbind(Trocrating,Brocrating)
  rocrating<-as.data.frame(rocrating)
  rocrating$sourcerating<-"across"
  names(rocrating)[1] <- "V2"
  
  
  fullrocrating_rec<-rbind(fullrocrating_rec,rocrating)
  
  fullrocrating_rec$Thit<-as.numeric(as.character(fullrocrating_rec$Thit))
  fullrocrating_rec$Tzhit<-qnorm(fullrocrating_rec$Thit)
  fullrocrating_rec$Tfa<-as.numeric(as.character(fullrocrating_rec$Tfa))
  fullrocrating_rec$Tzfa<-qnorm(fullrocrating_rec$Tfa)
  
  return(fullrocrating_rec)
}

ggplot_recog_split<-function(alldata){
  
  
  recogbysource.roc<-ggplot(alldata[alldata$sourcerating!="across",],
                            aes(x=Tfa,y=Thit,shape=sourcerating)) +
    geom_point(size=5)+
    scale_shape_manual(labels=c("Certain Bottom","Probably Bottom","Guess Bottom","Guess Top",
                                "Probably Top","Certain Top"),name="Source",
                       values=c(16,17,15,3,7,8)) +
    scale_x_continuous(name="FA",limits=c(0,1))+
    scale_y_continuous(name="Hit",limits=c(0,1))+
    geom_abline(intercept = 0,slope=1)+
    #geom_line()+
    facet_grid(.~V2)+
    theme_bw(base_size = 18)+ 
    theme(legend.position=c(0.85,0.26),
          legend.direction="vertical") +
    theme(plot.title=element_text(size=18, face="bold"),
          axis.text.x  = element_text(size=12, colour="black"),
          axis.title.x = element_text(size=12, colour="black"),
          axis.text.y  = element_text(size=12, colour="black"),
          axis.title.y = element_text(size=12,vjust=1.2, colour="black"),
          legend.text = element_text(size=12, colour="black"),
          legend.key = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title=element_text(colour="black",size=14),
          legend.background=element_rect(fill="transparent"),
          panel.border = element_rect(colour="black"),
          strip.text.x = element_text(colour = "white", size = 14, face="bold"),
          strip.background = element_rect(fill="black",colour="black"),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA))
  
  recogbysource.zroc<-ggplot(alldata[alldata$sourcerating!="across",],
                             aes(x=Tzfa,y=Tzhit,shape=sourcerating)) +
    geom_point(size=5) +
    #geom_line()+
    scale_shape_discrete(name="Source",labels=c("Certain Bottom","Maybe Bottom","Guess Bottom","Guess Top",
                                                "Maybe Top","Certain Top")) +
    scale_x_continuous(name="zFA")+
    scale_y_continuous(name="zHit")+
    facet_grid(.~V2)+
    theme_bw(base_size = 18)+ 
    theme(legend.position=c(0.9,0.3),
          legend.direction="vertical") +
    theme(plot.title=element_text(size=18, face="bold"),
          axis.text.x  = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y  = element_text(size=12),
          axis.title.y = element_text(size=12,vjust=1.2),
          legend.text = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title=element_text(colour="black"),
          panel.border = element_rect(colour="black"),
          strip.text.x = element_text(colour = "white", size = 14, face="bold"),
          strip.background = element_rect(fill="black",colour="black"),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA))
  
  plot_grid(recogbysource.roc, recogbysource.zroc, labels=c("A", "B"), ncol = 1, nrow = 2)
  
}


makerecogcollapseROC<-function(alldata){
  
  
  
  alldata<-ddply(alldata,.(StimType,RecConf,SourceConfrecoded),summarise,
                 item=sum(items))
  
  fullrocrating_rec<-NULL
  
  
  for (sourcerating in sort(unique(alldata$SourceConfrecoded),decreasing=T)){
    
    forroc<-alldata[alldata$SourceConfrecoded==sourcerating,]
    
    Trocrating<-NULL
    
    for (rating in sort(unique(alldata$RecConf),decreasing=T)){  
      
      Thit<-sum(forroc[forroc$StimType=="1" &
                         forroc$RecConf>=rating,]$item)/
        sum(forroc[forroc$StimType=="1",]$item)
      
      Tfa<-sum(forroc[forroc$StimType=="0" & 
                        forroc$RecConf>=rating,]$item)/
        sum(forroc[forroc$StimType=="0",]$item)
      
      
      Tpair<-cbind("OldNew",Thit,Tfa)
      Trocrating<-rbind(Trocrating,Tpair)
    }
    pretty<-cbind(sourcerating,Trocrating)
    fullrocrating_rec<-rbind(fullrocrating_rec,pretty)
    fullrocrating_rec<-as.data.frame(fullrocrating_rec)
  }
  
  
  
  acrosssource<-ddply(alldata,.(StimType,RecConf),summarise,
                      sumitems=sum(item))
  
  Trocrating<-NULL
  
  for (rating in sort(unique(acrosssource$RecConf),decreasing=T)){  
    
    Thit<-sum(acrosssource[acrosssource$StimType=="1" &
                             acrosssource$RecConf>=rating,]$sumitems)/
      sum(acrosssource[acrosssource$StimType=="1",]$sumitems)
    
    Tfa<-sum(acrosssource[acrosssource$StimType=="0" & 
                            acrosssource$RecConf>=rating,]$sumitems)/
      sum(acrosssource[acrosssource$StimType=="0",]$sumitems)
    
    
    Tpair<-cbind("OldNew",Thit,Tfa)
    
    Trocrating<-rbind(Trocrating,Tpair)
  }
  rocrating<-as.data.frame(Trocrating)
  rocrating$sourcerating<-"across"
  names(rocrating)[1] <- "V2"
  
  
  fullrocrating_rec<-rbind(fullrocrating_rec,rocrating)
  
  fullrocrating_rec$Thit<-as.numeric(as.character(fullrocrating_rec$Thit))
  fullrocrating_rec$Tzhit<-qnorm(fullrocrating_rec$Thit)
  fullrocrating_rec$Tfa<-as.numeric(as.character(fullrocrating_rec$Tfa))
  fullrocrating_rec$Tzfa<-qnorm(fullrocrating_rec$Tfa)
  
  return(fullrocrating_rec)
}

gridplotidnumber<-function(RTplot,Numberplot){
  grid.newpage()
  
  # two plots
  
  # extract gtable
  g1 <- ggplot_gtable(ggplot_build(RTplot))
  g2 <- ggplot_gtable(ggplot_build(Numberplot))
  
  # overlap the panel of 2nd plot on that of 1st plot
  pp <- c(subset(g1$layout, name == "panel", se = t:r))
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t,
                       4, pp$b, 4)
  
  
  # axis tweaks
  ia <- which(g2$layout$name == "axis-l")
  ga <- g2$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
  
  g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths)-1)
  
  
  g <- gtable_add_grob(g, ax, pp$t, length(g$widths)-1, pp$b)
  
  g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], 9)
  g <- gtable_add_grob(g, g2$grobs[[which(g2$layout$name == "ylab-l")]], pp$t, 9, pp$b)
  
  g$layout$clip[grep("layout",g$layout$name)] <- "off"
  
  return(grid.draw(g))}



set_panel_size <- function(p=NULL, g=ggplotGrob(p), file=NULL, 
                           margin = unit(1,"mm"),
                           width=unit(7.2, "cm"), 
                           height=unit(7, "cm")){
  
  panels <- grep("panel", g$layout$name)
  panel_index_w<- unique(g$layout$l[panels])
  panel_index_h<- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  
  if(getRversion() < "3.3.0"){
    
    # the following conversion is necessary
    # because there is no `[<-`.unit method
    # so promoting to unit.list allows standard list indexing
    g$widths <- grid:::unit.list(g$widths)
    g$heights <- grid:::unit.list(g$heights)
    
    g$widths[panel_index_w] <-  rep(list(width),  nw)
    g$heights[panel_index_h] <- rep(list(height), nh)
    
  } else {
    
    g$widths[panel_index_w] <-  rep(width,  nw)
    g$heights[panel_index_h] <- rep(height, nh)
    
  }
  
  if(!is.null(file))
    ggsave(file, g, 
           width = convertWidth(sum(g$widths) + margin, 
                                unitTo = "in", valueOnly = TRUE),
           height = convertHeight(sum(g$heights) + margin,  
                                  unitTo = "in", valueOnly = TRUE))
  
  g
}


set_panel_size2 <- function(p=NULL, g=ggplotGrob(p), file=NULL, 
                            margin = unit(1,"mm"),
                            width=unit(2.2, "cm"), 
                            height=unit(7, "cm")){
  
  panels <- grep("panel", g$layout$name)
  panel_index_w<- unique(g$layout$l[panels])
  panel_index_h<- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  
  if(getRversion() < "3.3.0"){
    
    # the following conversion is necessary
    # because there is no `[<-`.unit method
    # so promoting to unit.list allows standard list indexing
    g$widths <- grid:::unit.list(g$widths)
    g$heights <- grid:::unit.list(g$heights)
    
    g$widths[panel_index_w] <-  rep(list(width),  nw)
    g$heights[panel_index_h] <- rep(list(height), nh)
    
  } else {
    
    g$widths[panel_index_w] <-  rep(width,  nw)
    g$heights[panel_index_h] <- rep(height, nh)
    
  }
  
  if(!is.null(file))
    ggsave(file, g, 
           width = convertWidth(sum(g$widths) + margin, 
                                unitTo = "in", valueOnly = TRUE),
           height = convertHeight(sum(g$heights) + margin,  
                                  unitTo = "in", valueOnly = TRUE))
  
  g
}


set_panel_size3 <- function(p=NULL, g=ggplotGrob(p), file=NULL, 
                            margin = unit(1,"mm"),
                            width=unit(8, "cm"), 
                            height=unit(8, "cm")){
  
  panels <- grep("panel", g$layout$name)
  panel_index_w<- unique(g$layout$l[panels])
  panel_index_h<- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  
  if(getRversion() < "3.3.0"){
    
    # the following conversion is necessary
    # because there is no `[<-`.unit method
    # so promoting to unit.list allows standard list indexing
    g$widths <- grid:::unit.list(g$widths)
    g$heights <- grid:::unit.list(g$heights)
    
    g$widths[panel_index_w] <-  rep(list(width),  nw)
    g$heights[panel_index_h] <- rep(list(height), nh)
    
  } else {
    
    g$widths[panel_index_w] <-  rep(width,  nw)
    g$heights[panel_index_h] <- rep(height, nh)
    
  }
  
  if(!is.null(file))
    ggsave(file, g, 
           width = convertWidth(sum(g$widths) + margin, 
                                unitTo = "in", valueOnly = TRUE),
           height = convertHeight(sum(g$heights) + margin,  
                                  unitTo = "in", valueOnly = TRUE))
  
  g
}

