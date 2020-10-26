# CID-RS model free analysis (continuous identification, recognition, source memory)
# created NL - 07-Dec-2016
# amended NL - 02-Feb-2017
# amended for cidrs22 NL - 19-12-2019

# Set up ----------------------------------------------------------------------
rm(list=ls()) # Clears everything from environment

source("summaryfunctions.R")
init(c("plyr","ggplot2","reshape2","car","afex","grid","emmeans",
       "nortest","BayesFactor","dplyr","tidyr","magrittr","tibble"))

Experiment <- "CIDRS22"
# Import preprocessed data ----------------------------------------------------

memdata<-read.csv(paste0(Experiment,"_fulldatascreened.csv"))

# Only analyse valid trials
memdata_cleaned<-memdata[memdata$ExclReason=="valid",]


# Set columntypes -------------------------------------------------------------

memdata_cleaned[colnames(memdata_cleaned)] <- lapply(memdata_cleaned[colnames(memdata_cleaned)], as.character) #all columns as character

factorcols<-c("SubjID","Block", "StimType","SourceJud",
              "SourceCorr","SourceConf","Bottom1Top2New0","M1F2") #find factor columns
numbercols<-c("TestTrial","idRT","SourceRT","SubjAge","TotExpDur") #find numeric columns
memdata_cleaned[factorcols] <- lapply(memdata_cleaned[factorcols], as.factor) #subset columns as factors
memdata_cleaned[numbercols] <- lapply(memdata_cleaned[numbercols], as.numeric) #subset columns as numeric

# change labels of Source
memdata_cleaned$Bottom1Top2New0<-revalue(memdata_cleaned$Bottom1Top2New0,c("0" ="New", "1"="bBottom", "2"="aTop")) 

# add SplitHalf column
memdata_cleaned$SplitHalf<-"Odd"
memdata_cleaned[memdata_cleaned$TestTrial %% 2 ==0,]$SplitHalf<-"Even"
memdata_cleaned$SplitHalf<-as.factor(as.character(memdata_cleaned$SplitHalf))


# Add SourceConfrecoded and SourceConfLMH column ---------------------------------------

memdata_cleaned$SourceConfrecoded<-memdata_cleaned$SourceConf
memdata_cleaned$SourceConfrecoded<-revalue(memdata_cleaned$SourceConfrecoded,
                                           c("1"="S1","2"="S2","3"="S3",
                                             "4"="S4","5"="S5","6"="S6"))
memdata_cleaned$SourceConfrecoded<-as.character(memdata_cleaned$SourceConfrecoded)
memdata_cleaned[memdata_cleaned$Bottom1Top2New0=="bBottom" & 
                  memdata_cleaned$SourceConfrecoded=="S1",]$SourceConfrecoded<-6
memdata_cleaned[memdata_cleaned$Bottom1Top2New0=="bBottom" & 
                  memdata_cleaned$SourceConfrecoded=="S2",]$SourceConfrecoded<-5
memdata_cleaned[memdata_cleaned$Bottom1Top2New0=="bBottom" & 
                  memdata_cleaned$SourceConfrecoded=="S3",]$SourceConfrecoded<-4
memdata_cleaned[memdata_cleaned$Bottom1Top2New0=="bBottom" & 
                  memdata_cleaned$SourceConfrecoded=="S4",]$SourceConfrecoded<-3
memdata_cleaned[memdata_cleaned$Bottom1Top2New0=="bBottom" & 
                  memdata_cleaned$SourceConfrecoded=="S5",]$SourceConfrecoded<-2
memdata_cleaned[memdata_cleaned$Bottom1Top2New0=="bBottom" & 
                  memdata_cleaned$SourceConfrecoded=="S6",]$SourceConfrecoded<-1
memdata_cleaned[memdata_cleaned$SourceConfrecoded=="1",]$SourceConfrecoded<-"S1"
memdata_cleaned[memdata_cleaned$SourceConfrecoded=="2",]$SourceConfrecoded<-"S2"
memdata_cleaned[memdata_cleaned$SourceConfrecoded=="3",]$SourceConfrecoded<-"S3"
memdata_cleaned[memdata_cleaned$SourceConfrecoded=="4",]$SourceConfrecoded<-"S4"
memdata_cleaned[memdata_cleaned$SourceConfrecoded=="5",]$SourceConfrecoded<-"S5"
memdata_cleaned[memdata_cleaned$SourceConfrecoded=="6",]$SourceConfrecoded<-"S6"
memdata_cleaned$SourceConfrecoded<-as.factor(as.character(memdata_cleaned$SourceConfrecoded))

# SourceConf ID LMH coding 
memdata_cleaned$SourceConfLMH<-"aL"
memdata_cleaned[memdata_cleaned$SourceConfrecoded=="S1" | 
                  memdata_cleaned$SourceConfrecoded=="S6",]$SourceConfLMH<-"cH"
memdata_cleaned[memdata_cleaned$SourceConfrecoded=="S2" | 
                  memdata_cleaned$SourceConfrecoded=="S5",]$SourceConfLMH<-"bM"
memdata_cleaned$SourceConfLMH<-as.factor(as.character(memdata_cleaned$SourceConfLMH))


# Normality of RT -----------------------------------------------------


normalities<-function(dataparticipant){
  lillie.raw<-dlply(dataparticipant,.(StimType),
                    function(x) lillie.test(x$idRT))
  ad.raw<-dlply(dataparticipant,.(StimType),
                function(x) ad.test(x$idRT))
  shapiro.raw<-dlply(dataparticipant,.(StimType),
                     function(x) shapiro.test(x$idRT))
  
  lillie.log<-dlply(dataparticipant,.(StimType),
                    function(x) lillie.test(log(x$idRT)))
  ad.log<-dlply(dataparticipant,.(StimType),
                function(x) ad.test(log(x$idRT)))
  shapiro.log<-dlply(dataparticipant,.(StimType),
                     function(x) shapiro.test(log(x$idRT)))
  
  collect<-data.frame(Participant = rep(unique(dataparticipant$SubjID),12),
                      StimType = rep(c(1,2),6),
                      transformation = rep(c("raw","logged"),each=6),
                      tests = rep(rep(c("lilliefor","ad","shapiro"),each=2),2),
                      pval = c(lillie.raw$`1`$p.value,lillie.raw$`2`$p.value,
                               ad.raw$`1`$p.value,ad.raw$`2`$p.value,
                               shapiro.raw$`1`$p.value,shapiro.log$`2`$p.value,
                               lillie.log$`1`$p.value,lillie.log$`2`$p.value,
                               ad.log$`1`$p.value,ad.log$`2`$p.value,
                               shapiro.log$`1`$p.value,shapiro.log$`2`$p.value)
  )
  collect$normal<-0
  collect[collect$pval>0.050,]$normal<-1
  return(collect)
}
normalityRT<-NULL
for (i in unique(memdata_cleaned$SubjID)){
  
  results<-normalities(memdata_cleaned[memdata_cleaned$SubjID==i,])
  normalityRT<-rbind(normalityRT,results)
}

ddply(normalityRT,.(transformation,tests,StimType),summarize,
      ppnormal = sum(normal),
      percentage = ppnormal/length(unique(memdata_cleaned$SubjID)))



# SOURCE MEMORY ---------------------------------------------------------
# main source memory 


source_all_est <- memdata_cleaned %>%
  group_by(SubjID, Bottom1Top2New0, SourceCorr) %>%
  summarize(number = length(SourceCorr)) %>%
  complete_(c("SubjID","Bottom1Top2New0","SourceCorr")) %>%
  mutate_each(funs(replace(., is.na(.), 0)), number) %>%
  filter((Bottom1Top2New0 == "aTop"  & SourceCorr == "1") |
           (Bottom1Top2New0 == "bBottom" & SourceCorr == "2") |
           (Bottom1Top2New0 == "bBottom"  & SourceCorr == "1") |
           (Bottom1Top2New0 == "aTop" & SourceCorr == "2")) %>%
  distinct() %>%
  arrange(SubjID,Bottom1Top2New0,SourceCorr) %>%
  mutate(adjnumber = number + 0.5) %>%
  mutate(rates = number/sum(number),
         adjrates = adjnumber/sum(adjnumber)) %>%
  filter((Bottom1Top2New0=="aTop" & SourceCorr=="1") | 
           (Bottom1Top2New0=="bBottom" & SourceCorr=="2")) %>%
  group_by(SubjID) %>%
  summarize(pr_recog = diff(rates),
            dprime = diff(qnorm(adjrates)),
            cbias = -0.5 * sum(qnorm(adjrates)))

## test against chance: dprime
ttestAgainstZero(source_all_est %>% .$dprime)

## test against chance: c-bias
ttestAgainstZero(source_all_est %>% .$cbias)


# source memory split
source_split_est <- memdata_cleaned %>%
  group_by(SubjID, SplitHalf, Bottom1Top2New0, SourceCorr) %>%
  summarize(number = length(SourceCorr)) %>%
  complete_(c("SubjID","SplitHalf","Bottom1Top2New0","SourceCorr")) %>%
  mutate_each(funs(replace(., is.na(.), 0)), number) %>%
  filter((Bottom1Top2New0 == "aTop"  & SourceCorr == "1") |
           (Bottom1Top2New0 == "bBottom" & SourceCorr == "2") |
           (Bottom1Top2New0 == "bBottom"  & SourceCorr == "1") |
           (Bottom1Top2New0 == "aTop" & SourceCorr == "2")) %>%
  distinct() %>%
  arrange(SubjID, SplitHalf,Bottom1Top2New0,SourceCorr) %>%
  mutate(adjnumber = number + 0.5) %>%
  mutate(rates = number/sum(number),
         adjrates = adjnumber/sum(adjnumber)) %>%
  filter((Bottom1Top2New0=="aTop" & SourceCorr=="1") | 
           (Bottom1Top2New0=="bBottom" & SourceCorr=="2")) %>%
  group_by(SubjID,SplitHalf) %>%
  summarize(pr_recog = diff(rates),
            dprime = diff(qnorm(adjrates)),
            cbias = -0.5 * sum(qnorm(adjrates)))


cor.test(source_split_est %>% filter(SplitHalf == "Even") %>% .$dprime,
         source_split_est %>% filter(SplitHalf == "Odd" ) %>% .$dprime)
correlationBF(source_split_est %>% filter(SplitHalf == "Even") %>% .$dprime,
              source_split_est %>% filter(SplitHalf == "Odd") %>% .$dprime)


# Mean measures for table ----------------

memdata_cleaned %>%
  group_by(SubjID,StimType) %>%
  summarize(ONRT = mean(idRT)) %>%
  group_by(StimType) %>% 
  summarize(meanRT = mean(ONRT),
            seRT = serror(ONRT))


source_all<-ddply(memdata_cleaned[memdata_cleaned$Bottom1Top2New0!="New",],.(SubjID,Bottom1Top2New0,SourceCorr),summarize,
                  number=length(SourceCorr))
source_all_rates<-ddply(source_all, c("SubjID","Bottom1Top2New0"), mutate,
                        rates = number/sum(number))


ddply(source_all_rates,.(SourceCorr),summarize,
      means = mean(rates),
      ses = serror(rates))

# SOURCE MEMORY and ID -------------------------------------------------------------
#1. RT for correct vs incorrect source 

RTbysourcecorrectd <- memdata_cleaned %>%
  filter(SourceCorr != 0) %>%
  group_by(SubjID,SourceCorr) %>%
  summarize(sourceID=mean(idRT),
            numberID = length(idRT)) %>%
  complete_("SourceCorr") %>% 
  filter(SourceCorr != 0) %>% 
  group_by(SubjID) %>% 
  filter(all(!is.na(sourceID))) %>%
  group_by(SubjID) %>% 
  summarize(sourceiddiff = diff(sourceID))


ttestAgainstZero(RTbysourcecorrectd %>% .$sourceiddiff)


#2. Percentage participants source correct RT < source incorrect RT

RTbysourcecorrectd %>%
  summarize(percentageRT = sum(sourceiddiff > 0)/length(sourceiddiff))


#3. ID across source confidence 

source_all_idf <- memdata_cleaned %>%
  group_by(SubjID,SourceCorr,SourceConfLMH) %>%
  summarize(sourceRT = mean(idRT),
            numberRT = length(idRT)) %>%
  complete_("SourceConfLMH") %>% 
  filter(SourceCorr != 0)

source_all_lw <- source_all_idf %>%
  group_by(SubjID) %>% 
  filter(all(!is.na(sourceRT)))


source_all_idf %>%
  group_by(SourceCorr,SourceConfLMH) %>%
  summarize(alln =length(numberRT[!is.na(numberRT)]),
            means = mean(numberRT,na.rm=T),
            ses = serror(numberRT),
            meanRT = mean(sourceRT,na.rm=T),
            seRT = serror(sourceRT))

source_all_lw %>%
  group_by(SourceCorr,SourceConfLMH) %>%
  summarize(alln =length(numberRT[!is.na(numberRT)]),
            means = mean(numberRT,na.rm=T),
            ses = serror(numberRT),
            meanRT = mean(sourceRT,na.rm=T),
            seRT = serror(sourceRT)) 


IDLMH_anova<-aov_car(sourceRT ~ SourceConfLMH*SourceCorr+
                       Error(SubjID/SourceConfLMH*SourceCorr),
                     data=source_all_lw, anova_table = list(correction = "GG", es = "pes"))
IDLMH_anova$anova_table$`Pr(>F)`
summary(IDLMH_anova)
IDLMH_anova



IDLMH_conftrend <- emmeans(IDLMH_anova, ~SourceConfLMH)
contrast(IDLMH_conftrend, "poly")
contrast(IDLMH_conftrend, "pairwise", adjust="bonferroni")


IDLMH_BFanova <- anovaBF(sourceRT ~ SourceConfLMH*SourceCorr + SubjID, 
                         data = source_all_lw, 
                         whichRandom = "SubjID",
                         whichModels = "withmain",
                         progress=FALSE)
IDLMH_BFanova
IDLMH_BFanova[4]/IDLMH_BFanova[3]


# Graph


source_all_id_data<-ddply(memdata_cleaned,
                          .(SubjID,StimType,SourceCorr,SourceConfLMH),summarize,
                          sourceRT=mean(idRT),
                          numberSID=length(idRT))

source_all_id_data<-source_all_id_data[source_all_id_data$StimType!="2" & source_all_id_data$SourceCorr!="0",]
source_all_id_data$SourceCorr<-as.factor(as.character(source_all_id_data$SourceCorr))
emptydf<-with(source_all_id_data, expand.grid(
  SubjID=levels(SubjID),
  SourceCorr = levels(SourceCorr),
  SourceConfLMH = levels(SourceConfLMH)))
source_all_idf<-merge(source_all_id_data,emptydf,all.y=TRUE)

source_all_lw<-source_all_idf[!(source_all_idf$SubjID %in% 
                                  unique(source_all_idf[is.na(source_all_idf$sourceRT),]$SubjID)),]

SourceID<-source_all_lw %>% 
  group_by(SourceCorr,SourceConfLMH) %>% 
  summarize(sourceRTs = mean(sourceRT),
            ci= serror(sourceRT) * qt(.975/2 + .5,length(unique(SubjID)) - 1))


source_all_id<-ddply(memdata_cleaned,
                     .(SubjID,StimType,SourceCorr),summarize,
                     sourceRT=mean(idRT),
                     number=length(idRT))
SourceID_all<-source_all_id %>% 
                          filter(StimType !=2) %>% 
                        group_by(SourceCorr) %>% 
                          summarize(sourceRTs = mean(sourceRT),
                                    ci= serror(sourceRT) * qt(.975/2 + .5,length(unique(SubjID)) - 1))
SourceID_all$SourceConfLMH<-"aAll"
AllSource<-bind_rows(SourceID,SourceID_all)
AllSource$grouping<-c(1,1,1,2,2,2,3,4)
AllSource$grouping<-as.factor(as.character(AllSource$grouping))
AllSource$SourceCorr<-as.factor(as.character(AllSource$SourceCorr))
AllSource$SourceConfLMH<-as.factor(as.character(AllSource$SourceConfLMH))


emp_cidsource<-ggplot(AllSource, aes(x=SourceConfLMH, y=sourceRTs, group=grouping, fill = SourceCorr)) +
  geom_line(data=AllSource, aes(x=SourceConfLMH, y=sourceRTs),size=0.5,linetype="solid",position=position_dodge(0.4)) + 
  
  geom_errorbar(data=AllSource, aes(x=SourceConfLMH,ymin=sourceRTs-ci, ymax=sourceRTs+ci), 
                width=0,position=position_dodge(0.4))+
  geom_point(size=4,data=AllSource, aes(x=SourceConfLMH, y=sourceRTs,fill=SourceCorr),shape=22,
             position=position_dodge(0.3), stroke = 0.5) +# draw bar graphs
  
  geom_segment(x = 1.5, xend = 1.5, y = 1300, yend=2600,linetype="dashed")+
  scale_fill_manual(name='Item type',values=c("#000000","#FFFFFF"), labels=c("Correct","Incorrect"))+
  guides(fill = guide_legend(override.aes = list(fill=c("#000000","#FFFFFF"),
                                                 size = 4))) +
  scale_x_discrete(name='Source confidence',labels=c('Overall\n','Guess','Probably','Sure'))+
  scale_y_continuous(name='Identification RT (ms)',limits=c(1400,2600),breaks=seq(1400,2600,200)) +
  #facet_wrap(~type,labeller=labeller(type = idlabels))+
  # ggtitle(label="B")+
  theme_bw(base_size = 18)+ 
  theme(legend.position = c(1,1),
        #legend.position = "none",
        plot.margin = unit(c(1,1,1,2), "cm"),
        legend.justification = c(1, 1), 
        legend.direction="vertical",
        plot.title=element_text(size=18, face="bold"),
        aspect.ratio=1,
        plot.background = element_rect(fill="transparent",colour="transparent"),
        axis.text.x  = element_text(size=12,colour="black"),
        axis.title.x = element_blank(),
        axis.text.y  = element_text(size=12,colour="black"),
        axis.title.y = element_text(size=12,colour="black"),
        legend.text = element_text(size = 12),
        legend.key=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour="black", size=0.5),
        strip.text.x = element_text(colour = "black", size = 12, face="bold"),
        strip.background = element_rect(fill="transparent",colour="transparent"))
emp_cidsource

cidsource3B<-set_panel_size3(emp_cidsource)
lay4<-rbind(c(1))
grid.arrange(cidsource3B,
             widths = c(6)) #14.5x4; 1400 x 400

