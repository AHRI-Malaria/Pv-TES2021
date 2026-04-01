### Pv_2021

## code used for generating figures in the manuscript entitled "Efficacy of chloroquine plus low-dose primaquine and pyronaridine–artesunate plus low-dose primaquine against Plasmodium vivax malaria in adults and transmission to mosquitoes in Ethiopia: a phase 2/3, observer-masked, randomised, parallel-group trial"

library(readxl)
library(tidyverse)
library(survival)
library(survminer)
library(ggprism)
library(ggpubr)
library(pammtools)
library(ggrepel)
library(mgcv)
library(epitools)
library(blme)
library(flextable)
library(dplyr)
library(ggplot2)

parmic.w <- read_excel("data/PvTES_data.xlsx",sheet="Pv_TES_All")
parmic.w = parmic.w %>% 
  mutate(id = `study id`, parmic0 = pv_asec_dens_d0,density_0 = pv_asec_dens_d0,density_1 = as.numeric(pv_asec_dens_d1),density_2=as.numeric(pv_asec_dens_d2),density_3 = as.numeric(pv_asec_dens_d3),density_7 = as.numeric(pv_asec_dens_d7),density_14 = as.numeric(pv_asec_dens_d14),density_21=as.numeric(pv_asec_dens_d21),density_28 = as.numeric(pv_asec_dens_d28),density_35 = as.numeric(pv_asec_dens_d35),density_42 = as.numeric(pv_asec_dens_d42)) %>%
  dplyr:: select(id, `studyarm_d0`, parmic0, density_0, density_1, density_2, density_3, density_7, density_14, density_21, density_28, density_35, density_42)


parmic <- parmic.w %>%
  pivot_longer(cols = starts_with(c("density")),
               names_to = c(".value", "visit"),
               names_sep = "_") %>%
  mutate(visit = as.numeric(visit),
         id = as.factor(id)) %>% 
  na.omit()


parmic.s <- parmic %>% 
  group_by(id) %>%
  mutate(cleared = ifelse(lag(density,default=0)>0 & lead(density,default = 0)==0 & density==0,1,0),
         everclear = sum(cleared),rev=n():1,keep= (cleared==1 & cumsum(cleared)==1) | (everclear==0 & rev==1), reinfected = ifelse(lag(density)==0 & lag(density, n=2)==0 & density>0, 1, 0 ), reinfected=sum(reinfected,na.rm=T)) %>%
  filter(keep) %>%
  dplyr::select(-keep,-rev,-everclear,-density)
### microscopy gametocyte density
# pvgam_dens_d0= microscopy gametocyte density at day 0

gammic.w <- read_excel("data/PvTES_data.xlsx",sheet="Pv_TES_All")

gammic.w = gammic.w %>% 
  mutate(id = `study id`, gammic0 = pvgam_dens_d0,density_0 = pvgam_dens_d0,density_1 = as.numeric(pvgam_dens_d1),density_2=as.numeric(pvgam_dens_d2),density_3 = as.numeric(pvgam_dens_d3),density_7 = as.numeric(pvgam_dens_d7),density_14 = as.numeric(pvgam_dens_d14),density_21=as.numeric(pvgam_dens_d21),density_28 = as.numeric(pvgam_dens_d28),density_35 = as.numeric(pvgam_dens_d35),density_42 = as.numeric(pvgam_dens_d42)) %>%
  dplyr::select(id, `studyarm_d0`, gammic0, density_0, density_1, density_2, density_3, density_7, density_14, density_21, density_28, density_35, density_42)

gammic <- gammic.w %>%
  pivot_longer(cols = starts_with(c("density")),
               names_to = c(".value", "visit"),
               names_sep = "_") %>%
  mutate(visit = as.numeric(visit),
         id = as.factor(id)) %>% 
  na.omit()


gammic.s <- gammic  %>% filter(gammic0>0)%>% 
  group_by(id) %>%
  mutate(cleared = ifelse(lag(density,default=0)>0 & lead(density,default = 0)==0 & density==0,1,0),
         everclear = sum(cleared),rev=n():1,keep= (cleared==1 & cumsum(cleared)==1) | (everclear==0 & rev==1), reinfected = ifelse(lag(density)==0 & lag(density, n=2)==0 & density>0, 1, 0 ), reinfected=sum(reinfected,na.rm=T)) %>%
  filter(keep) %>%
  dplyr::select(-keep,-rev,-everclear,-density)


##### qPCR parasite density
##PvparapuLbdd0 = qPCR parasite density at day 0
parpcr.w <- read_excel("data/PvTES_data.xlsx",sheet="Pv_TES_All")
parpcr.w = parpcr.w %>% 
  mutate(id = `study id`, parpcr0 = PvparapuLbdd0,density_0 = PvparapuLbdd0,density_1 = as.numeric(PvparapuLbdd1),density_2=as.numeric(PvparapuLbdd2),density_3 = as.numeric(PvparapuLbdd3),density_7 = as.numeric(PvparapuLbdd7),density_14 = as.numeric(PvparapuLbdd14)) %>%
  dplyr::select(id, `studyarm_d0`, parpcr0, density_0, density_1, density_2, density_3, density_7, density_14)

parpcr <- parpcr.w %>%
  pivot_longer(cols = starts_with(c("density")),
               names_to = c(".value", "visit"),
               names_sep = "_") %>%
  mutate(visit = as.numeric(visit),
         id = as.factor(id),
         density = ifelse(density<3,0,density)) %>% 
  na.omit()

parpcr.s <- parpcr %>% 
  group_by(id) %>%
  mutate(cleared = ifelse(lag(density,default=0)>0 & lead(density,default = 0)==0 & density==0,1,0),
         everclear = sum(cleared),rev=n():1,keep= (cleared==1 & cumsum(cleared)==1) | (everclear==0 & rev==1), reinfected = ifelse(lag(density)==0 & lag(density, n=2)==0 & density>0, 1, 0 ), reinfected=sum(reinfected,na.rm=T)) %>%
  filter(keep) %>%
  dplyr::select(-keep,-rev,-everclear,-density)


### qPCR gametocyte density
# pvgampulbdd0 = qPCR gametocyte density at day 0 

gampcr.w <- read_excel("data/PvTES_data.xlsx",sheet="Pv_TES_All")
gampcr.w = gampcr.w %>% 
  mutate(id = `study id`, gampcr0 = pvgampulbdd0,density_0 = pvgampulbdd0,density_1 = as.numeric(pvgampulbdd1),density_2=as.numeric(pvgampulbdd2),density_3 = as.numeric(pvgampulbdd3),density_7 = as.numeric(pvgampulbdd7),density_14 = as.numeric(pvgampulbdd14)) %>%
  dplyr::select(id, `studyarm_d0`, gampcr0, density_0, density_1, density_2, density_3, density_7, density_14)

gampcr <- gampcr.w %>%
  pivot_longer(cols = starts_with(c("density")),
               names_to = c(".value", "visit"),
               names_sep = "_") %>%
  mutate(visit = as.numeric(visit),
         id = as.factor(id),
         density = ifelse(density<3,0,density)) %>% 
  na.omit()


gampcr.s <- gampcr %>% filter(gampcr0>0) %>% 
  group_by(id) %>%
  mutate(cleared = ifelse(lag(density,default=0)>0 & lead(density,default = 0)==0 & density==0,1,0),
         everclear = sum(cleared),rev=n():1,keep= (cleared==1 & cumsum(cleared)==1) | (everclear==0 & rev==1), reinfected = ifelse(lag(density)==0 & lag(density, n=2)==0 & density>0, 1, 0 ), reinfected=sum(reinfected,na.rm=T)) %>%
  filter(keep) %>%
  dplyr::select(-keep,-rev,-everclear,-density)
## long form

longdt <- parmic %>% 
  full_join(gammic, by = c("id", "studyarm_d0", "visit")) %>%
  full_join(parpcr, by = c("id", "studyarm_d0", "visit")) %>%
  full_join(gampcr, by = c("id", "studyarm_d0", "visit"))

dt0 = longdt[longdt$visit==0,]
longdt_0 = longdt[longdt$visit>0,]


# Survival analysis -------------------------------------------------------


km6a = survfit(Surv(visit,cleared)~`studyarm_d0`, data= parmic.s)
survout6a = survminer::ggsurvplot(km6a, conf.int = TRUE, risk.table = TRUE,xlim=c(0,14),risk.table.fontsize = 3, break.time.by=1, ris.table.height=0.2)
dt6a = survout6a$data.survplot


km6b = survfit(Surv(visit,cleared)~`studyarm_d0`, data= gammic.s)
survout6b = survminer::ggsurvplot(km6b, conf.int = TRUE, risk.table = TRUE,xlim=c(0,14),risk.table.fontsize = 3, break.time.by=1, ris.table.height=0.2)
dt6b = survout6b$data.survplot


km6c = survfit(Surv(visit,cleared)~`studyarm_d0`, data= parpcr.s)
survout6c = survminer::ggsurvplot(km6c, conf.int = TRUE, risk.table = TRUE,xlim=c(0,14),risk.table.fontsize = 3, break.time.by=1, ris.table.height=0.2)
dt6c = survout6c$data.survplot

km6d = survfit(Surv(visit,cleared)~`studyarm_d0`, data= gampcr.s)
survout6d = survminer::ggsurvplot(km6d, conf.int = TRUE, risk.table = TRUE,xlim=c(0,14),risk.table.fontsize = 3, break.time.by=1, ris.table.height=0.2)
dt6d = survout6d$data.survplot


dt6a$type = "Asexual parasites (microscopy)"
dt6b$type = "Gametocytes (microscopy)"
dt6c$type = "Asexual parasites (pcr)"
dt6d$type = "Gametocytes (pcr)"

dt6 = rbind(dt6a, dt6b, dt6c, dt6d) %>%
  dplyr::select(time, surv, n.censor, upper, lower, `studyarm_d0`, type)

dt6$type = factor(dt6$type, levels = c("Asexual parasites (microscopy)","Gametocytes (microscopy)","Asexual parasites (pcr)","Gametocytes (pcr)"))

dt6sub = expand.grid(time=0, surv=1, upper=1, lower=1, `studyarm_d0`=sort(unique(dt6$`studyarm_d0`)),type=sort(unique(dt6$type)))
# Define the columns in dt6sub to match dt6
dt6sub$n.censor = 0

# Reorder the columns in dt6sub to match the order in dt6
dt6sub = dt6sub[, names(dt6)]
dt6 = unique(rbind(dt6,dt6sub))
###


cox1a = coxph(Surv(visit,cleared)~`studyarm_d0`+log10(parmic0), data= parmic.s)
sm1a=summary(cox1a)
cf = confint(cox1a)
p1a = sm1a$coefficients[1,5]
p1a = ifelse(p1a<0.0001, "<0.0001",ifelse(p1a>0.9999,">0.9999",paste0("=",format(round(p1a,4),nsmall=4))))
(HR1a = paste0("HR = ",format(round(exp(cox1a$coefficients[1]),2),nsmall=2)," (95% CI: ",format(round(exp(cf[1,1]),2),nsmall=2),", ",format(round(exp(cf[1,2]),2),nsmall=2),"), p",p1a))

cox1b = coxph(Surv(visit,cleared)~`studyarm_d0`+log10(gammic0), data= gammic.s)
sm1b=summary(cox1b)
cf = confint(cox1b)
p1b = sm1b$coefficients[1,5]
p1b = ifelse(p1b<0.0001, "<0.0001",ifelse(p1b>0.9999,">0.9999",paste0("=",format(round(p1b,4),nsmall=4))))
(HR1b = paste0("HR = ",format(round(exp(cox1b$coefficients[1]),2),nsmall=2)," (95% CI: ",format(round(exp(cf[1,1]),2),nsmall=2),", ",format(round(exp(cf[1,2]),2),nsmall=2),"), p",p1b))



cox1c = coxph(Surv(visit,cleared)~`studyarm_d0`+log10(parpcr0), data= parpcr.s)
sm1c=summary(cox1c)
cf = confint(cox1c)
p1c = sm1c$coefficients[1,5]
p1c = ifelse(p1c< 0.0001, "<0.0001",ifelse(p1c>0.9999,">0.9999",paste0("=",format(round(p1c,4),nsmall=4))))
(HR1c = paste0("HR = ",format(round(exp(cox1c$coefficients[1]),2),nsmall=2)," (95% CI: ",format(round(exp(cf[1,1]),2),nsmall=2),", ",format(round(exp(cf[1,2]),2),nsmall=2),"), p",p1c))

cox1d = coxph(Surv(visit,cleared)~`studyarm_d0`+log10(gampcr0), data= gampcr.s)
sm1d=summary(cox1d)
cf = confint(cox1d)
p1d = sm1d$coefficients[1,5]
p1d = ifelse(p1d<0.0001, "<0.0001",ifelse(p1d>0.9999,">0.9999",paste0("=",format(round(p1d,4),nsmall=4))))
(HR1d = paste0("HR = ",format(round(exp(cox1d$coefficients[1]),2),nsmall=2)," (95% CI: ",format(round(exp(cf[1,1]),2),nsmall=2),", ",format(round(exp(cf[1,2]),2),nsmall=2),"), p",p1d))

#### Km plot
### For Appendix 2.1 : Parasite and gametocyte positivity.
Fig_S21 = ggarrange(
  
  ggarrange(
   
    
    survout6a$plot+
      xlab("Time (days)")+
      #geom_text(data=survout6a$data.survtable[survout6a$data.survtable$time %in% c(0:3,7,14) & survout6a$data.survtable$n.censor>0,] %>% mutate(y=ifelse(strata=="studyarm=CQ",1.1,1.05)), aes(x=time, y=y,label=paste0("c=",n.censor),col=strata))+
      #scale_y_continuous(labels=scales::percent)+
      #scale_y_continuous(labels = function(x) format(x * 100, nsmall = 0), breaks = seq(0, 100, by = 20), limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100"))+
      ylab("Proportion parasite \npositive (microscopy)")+
      scale_x_continuous(limits=c(0,14),breaks=c(0:3,7,14))+
      theme_classic()+
      theme(axis.line = element_line(size = 0.5))+
      scale_color_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
      scale_fill_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
      geom_text(aes(x=10,y=0.9),label=HR1a),
    
    survout6a$data.survtable[survout6a$data.survtable$time %in% c(0:3,7,14),] %>% ggplot(aes(x=time, y=strata,label=n.risk)) +
      scale_y_discrete(labels=c("CQ+PQ","PA+PQ"))+
      theme_classic()+
      theme(
        axis.title.y = element_blank()
      )+
      xlab("Time (days)")+
      scale_x_continuous(limits=c(0,14),breaks=c(0:3,7,14))+
      geom_text(data = survout6a$data.survtable, aes(x = time, y = strata, label = paste0(n.risk, "(", n.censor , ")")), vjust = -0.5, size = 3) +
      theme(plot.title = element_blank())+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.text.y = element_text(colour = c("red","darkblue"),size = 10,vjust = -0.5)),
    
    nrow = 2,
    heights =c(0.2,0.1),
    common.legend = TRUE,
    align="hv",
    legend="none"
  )
  ,
  
  ggarrange(
    
    survout6c$plot+
      #geom_text(data=survout6c$data.survtable[survout6c$data.survtable$time %in% c(0:3,7,14) & survout6c$data.survtable$n.censor>0,]%>% mutate(y=ifelse(strata=="studyarm=CQ",1.1,1.05)), aes(x=time, y=y,label=paste0("c=",n.censor),col=strata))+
      xlab("Time (days)")+
     # scale_y_continuous(labels=scales::percent)+
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100"))+
      ylab("Proportion Pv18S \npositive")+
      scale_x_continuous(limits=c(0,14),breaks=c(0:3,7,14))+
      theme_classic()+
      theme(axis.line = element_line(size = 0.5))+
      scale_color_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
      scale_fill_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
      geom_text(aes(x=10,y=0.9),label=HR1c),
    
    survout6c$data.survtable[survout6c$data.survtable$time %in% c(0:3,7,14),] %>% ggplot(aes(x=time, y=strata,label=n.risk)) + 
      geom_text(data = survout6c$data.survtable, aes(x = time, y = strata, label = paste0(n.risk, "(", n.censor , ")")), vjust = -0.5, size = 3) +
      scale_y_discrete(labels=c("CQ+PQ","PA+PQ"))+
      theme_classic()+
      theme(
        axis.title.y = element_blank()
      )+
      xlab("Time (days)")+
      scale_x_continuous(limits=c(0,14),breaks=c(0:3,7,14))+
      theme(plot.title = element_blank())+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            
            axis.text.y = element_text(colour = c("red","darkblue"),size = 10,vjust = -0.5)),
      
    
    nrow = 2,
    heights =c(0.2,0.1),
    common.legend = TRUE,
    align="hv",
    legend="none"
  )
  ,
  
  ggarrange(
    
    survout6b$plot+
      xlab("Time (days)")+
      #geom_text(data=survout6b$data.survtable[survout6b$data.survtable$time %in% c(0:3,7,14) & survout6b$data.survtable$n.censor>0,]%>% mutate(y=ifelse(strata=="studyarm=CQ",1.1,1.05)), aes(x=time, y=y,label=paste0("c=",n.censor),col=strata))+
      #scale_y_continuous(labels=scales::percent)+
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100")) +
      ylab("Proportion gametocyte \npositive (microscopy)")+
      scale_x_continuous(limits=c(0,14),breaks=c(0:3,7,14))+
      theme_classic()+
      theme(axis.line = element_line(size = 0.5))+
      scale_color_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
      scale_fill_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
      geom_text(aes(x=10,y=0.9),label=HR1b),
    
    survout6b$data.survtable[survout6b$data.survtable$time %in% c(0:3,7,14),] %>% ggplot(aes(x=time, y=strata,label=n.risk)) + 
      geom_text(data = survout6b$data.survtable, aes(x = time, y = strata, label = paste0(n.risk, "(", n.censor , ")")), vjust = -0.5, size = 3) +
      scale_y_discrete(labels=c("CQ+PQ","PA+PQ"))+
      theme_classic()+
      theme(
        axis.title.y = element_blank()
      )+
      xlab("Time (days)")+
      scale_x_continuous(limits=c(0,14),breaks=c(0:3,7,14))+
      theme(plot.title = element_blank())+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            
            axis.text.y = element_text(colour = c("red","darkblue"),size = 10,vjust = -0.5)),
    
    
    nrow = 2,
    heights = c(0.2,0.1),
    common.legend = TRUE,
    align="hv",
    legend="none"
  )
  
  ,

  
  ggarrange(
    
    survout6d$plot+
      xlab("Time (days)")+
      #geom_text(data=survout6d$data.survtable[survout6d$data.survtable$time %in% c(0:3,7,14) & survout6d$data.survtable$n.censor>0,]%>% mutate(y=ifelse(strata=="studyarm=CQ",1.1,1.05)), aes(x=time, y=y,label=paste0("c=",n.censor),col=strata))+
      #scale_y_continuous(labels=scales::percent)+
      #scale_y_continuous(labels = function(x) format(x * 100, nsmall = 0), breaks = seq(0, 100, by = 20),limits = c(0, 100)/100) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100")) +
      ylab("Proportion Pvs25 \npositive")+
      scale_x_continuous(limits=c(0,14),breaks=c(0:3,7,14))+
      theme_classic()+
      theme(axis.line = element_line(size = 0.5))+
      scale_color_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
      scale_fill_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
      geom_text(aes(x=10,y=0.9),label=HR1d),
    survout6d$data.survtable[survout6d$data.survtable$time %in% c(0:3,7,14),] %>% ggplot(aes(x=time, y=strata,label=n.risk)) +
      geom_text(data = survout6d$data.survtable, aes(x = time, y = strata, label = paste0(n.risk, "(", n.censor , ")")), vjust = -0.5, size = 3) +
      scale_y_discrete(labels=c("CQ+PQ","PA+PQ"))+
      theme_classic()+
      theme(
        axis.title.y = element_blank()
      )+
      xlab("Time (days)")+
      scale_x_continuous(limits=c(0,14),breaks=c(0:3,7,14))+
      theme(plot.title = element_blank())+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            
            axis.text.y = element_text(colour = c("red","darkblue"),size = 10,vjust = -0.5)),
    
    
    nrow = 2,
    heights = c(0.2,0.1),
    common.legend = TRUE,
    align="hv",
    legend="none"
  )
  ,
  
  nrow=2, ncol=2,
  labels=c("(A)","(B)","(C)","(D)"),
  common.legend = TRUE,
  legend="none",
  vjust= 2,
  hjust=-1)
Fig_S21
## Fig_S14
##For Appendix 1. 4: Parasite and gametocyte densities during follow-up.
# Longitudinal analysis redone --------------------------------------------

longdt_complete <- longdt[longdt$visit <= 14,] 

# Convert visit to a factor with explicit levels
longdt_complete$visit <- factor(longdt_complete$visit, levels = sort(unique(longdt$visit)))


options(scipen=999)
lplot1a <- ggplot(data=longdt_complete, aes(x=as.factor(visit), y=density.x, col=`studyarm_d0`), size=2) +
  geom_jitter(alpha=0.2, size=2, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  geom_boxplot(aes(fill=`studyarm_d0`), alpha=0.5, position=position_dodge(width=0.75)) +
  theme_classic() +
  ylab("Asexual parasite density/µL (mic)") +
  xlab("Days ") +
  scale_color_manual(values=c("red", "darkblue"), labels=c("CQ+PQ","PA+PQ")) +
  scale_fill_manual(values=c("red", "darkblue"), labels=c("CQ+PQ","PA+PQ")) +
  theme(legend.position = "top") +
  theme(axis.line = element_line(size = 0.5))+
  scale_y_log10(limits=c(0.5, 10000000), breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000)) +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))

lplot1a
  

lplot1b = ggplot(data=longdt_complete, aes(x=as.factor(visit), y=density.y, col=`studyarm_d0`),  size=2, width=0.2)+
  geom_jitter(alpha=0.2, size=2, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  geom_boxplot(aes(fill=`studyarm_d0`), alpha=0.5, position=position_dodge(width=0.75)) +
  theme_classic()+
  #coord_cartesian(xlim=c(0,14))+
  ylab("Gametocyte density/µL (mic)")+
  xlab("Days")+
  #scale_x_continuous(breaks = c(0:3,7,14))+
  scale_color_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  scale_fill_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  theme(legend.position = "top")+
  theme(axis.line = element_line(size = 0.5))+
  scale_y_log10(lim=c(0.5,10000000), breaks=c(1,10,100,1000,10000,100000,1000000)) +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))
lplot1b


lplot1c = ggplot(data=longdt_complete, aes(x=as.factor(visit), y=density.x.x, col=`studyarm_d0`),  size=2, width=0.2)+
  geom_jitter(alpha=0.2, size=2, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  geom_boxplot(aes(fill=`studyarm_d0`), alpha=0.5, position=position_dodge(width=0.75)) +
  theme_classic()+
  #coord_cartesian(xlim=c(0,14))+
  ylab("Pv18S copies/µL")+
  xlab("Days ")+
  #scale_x_continuous(breaks = c(0:3,7,14))+
  scale_color_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  scale_fill_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  theme(legend.position = "top")+
  theme(axis.line = element_line(size = 0.5))+
  scale_y_log10(lim=c(0.5,10000000), breaks=c(1,10,100,1000,10000,100000,1000000))+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))
lplot1c


lplot1d = ggplot(data=longdt_complete, aes(x=as.factor(visit), y=density.y.y, col=`studyarm_d0`),  size=2, width=0.2)+
  geom_jitter(alpha=0.2, size=2)+
  geom_boxplot(aes(fill=`studyarm_d0`), alpha=0.5)+
  theme_classic()+
  #coord_cartesian(xlim=c(0,14))+
  ylab("Pvs25 copies/µL")+
  xlab("Days")+
  #scale_x_continuous(breaks = c(0:3,7,14))+
  scale_color_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  scale_fill_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  theme(legend.position = "top")+
  theme(axis.line = element_line(size = 0.5))+
  scale_y_log10(lim=c(0.5,10000000), breaks=c(1,10,100,1000,10000,100000,1000000)) +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))
lplot1d


Fig_S14 =ggarrange(
  lplot1a,
  lplot1c,
  lplot1b,
  lplot1d,
  nrow=2, ncol=2,
  align="hv",
  common.legend = TRUE,
  legend="top",
  labels=c("(A)","(B)","(C)","(D)"),
  vjust=1.3,
  hjust=-1.8) 
Fig_S14

## Fig_S13
## For Appendix 1. 3: Raw prevalences and relative decrease in asexual parasite and gametocyte prevalence.
# Longitudinal prevalence -------------------------------------------------

longdt3 = longdt[longdt$visit<=3 & longdt$visit>0,]

mod2a = gam(I(density.x>0) ~ `studyarm_d0`+as.factor(visit) + log(parmic0) + 
              s(id, bs="re"),
            family=poisson(link="log"),
            data = longdt3)
(sm2a=summary(mod2a))
cf = cbind(sm2a$p.coeff-1.96*sm2a$se[1:5],sm2a$p.coeff+1.96*sm2a$se[1:5])
p2a = sm2a$p.pv[[2]]
p2a = ifelse(p2a<0.0001, "<0.0001",ifelse(p2a>0.9999,">0.9999",paste0("=",format(round(p2a,4),nsmall=4))))
(DR2a = paste0("RR = ",format(round(exp(mod2a$coefficients[[2]]),2),nsmall=2)," (95% CI: ",format(round(exp(cf[2,1]),2),nsmall=2),", ",format(round(exp(cf[2,2]),2),nsmall=2),"), p",p2a))

newd = longdt3 %>% 
  group_by(`studyarm_d0`, visit) %>%
  mutate(prev = I(density.x>0)) %>%
  summarise(n=n(),count=sum(prev),prev = mean(prev)) %>% 
  mutate(lwr = binom.exact(count, n, conf.level = 0.95)[[4]], upr=binom.exact(count, n, conf.level = 0.95)[[5]])
  


options(scipen=999)
lplot2a = ggplot(data=newd, aes(x=visit, y=prev,fill=as.factor(`studyarm_d0`)))+
  geom_bar(stat="identity",alpha=0.9, position=position_dodge())+
  theme_classic()+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,by=0.2), limits=c(0,1))+
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.1, linewidth=1.2, position=position_dodge(0.9))+
  ylab("Asexual parasite prevalence (mic)")+
  xlab("Days after treatment")+
  scale_x_continuous(breaks = 1:3)+
  scale_color_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  scale_fill_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  geom_text(aes(x=2,y=1),label=DR2a)+
  theme(axis.line = element_line(size = 0.1))+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme(legend.position = "top")
lplot2a



mod2b = gam(I(density.y>0) ~ `studyarm_d0`+as.factor(visit) + log(gammic0) + 
              s(id, bs="re"),
            family=poisson(link="log"),
            data = longdt3[longdt3$gammic0>0,])
(sm2b=summary(mod2b))
cf = cbind(sm2b$p.coeff-1.96*sm2b$se[1:5],sm2b$p.coeff+1.96*sm2b$se[1:5])
p2b = sm2b$p.pv[[2]]
p2b = ifelse(p2b<0.0001, "<0.0001",ifelse(p2b>0.9999,">0.9999",paste0("=",format(round(p2b,4),nsmall=4))))
(DR2b = paste0("RR = ",format(round(exp(mod2b$coefficients[[2]]),2),nsmall=2)," (95% CI: ",format(round(exp(cf[2,1]),2),nsmall=2),", ",format(round(exp(cf[2,2]),2),nsmall=2),"), p",p2b))

######## .....
newd <- longdt3[longdt3$gammic0 > 0,] %>%
  group_by(studyarm_d0, visit) %>%
  mutate(prev = as.integer(density.y > 0)) %>%  # Convert logical to integer
  summarise(n = n(), count = sum(prev), prev = mean(prev)) %>%
  filter(count >= 0) %>%  # Filter out rows with negative counts
  mutate(
    lwr = binom.exact(count, n, conf.level = 0.95)[[4]],
    upr = binom.exact(count, n, conf.level = 0.95)[[5]]
  )


options(scipen=999)
lplot2b = ggplot(data=newd, aes(x=visit, y=prev,fill=as.factor(`studyarm_d0`)))+
  geom_bar(stat="identity",alpha=0.9, position=position_dodge())+
  theme_classic()+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,by=0.2), limits=c(0,1))+
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.1, linewidth=1.2, position=position_dodge(0.9))+
  ylab("Gametocyte prevalence (mic)")+
  xlab("Days after treatment")+
  scale_x_continuous(breaks = 1:3)+
  scale_color_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  scale_fill_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  geom_text(aes(x=2,y=1),label=DR2b)+
  theme(axis.line = element_line(size = 0.1))+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme(legend.position = "top")
lplot2b



mod2c = gam(I(density.x.x>0) ~ `studyarm_d0`+as.factor(visit) + log(parpcr0) + 
              s(id, bs="re"),
            family=poisson(link="log"),
            data = longdt3)
(sm2c=summary(mod2c))
cf = cbind(sm2c$p.coeff-1.96*sm2c$se[1:5],sm2c$p.coeff+1.96*sm2c$se[1:5])
p2c = sm2c$p.pv[[2]]
p2c = ifelse(p2c<0.0001, "<0.0001",ifelse(p2c>0.9999,">0.9999",paste0("=",format(round(p2c,4),nsmall=4))))
(DR2c = paste0("RR = ",format(round(exp(mod2c$coefficients[[2]]),2),nsmall=2)," (95% CI: ",format(round(exp(cf[2,1]),2),nsmall=2),", ",format(round(exp(cf[2,2]),2),nsmall=2),"), p",p2c))

newd = longdt3 %>% 
  group_by(`studyarm_d0`, visit) %>%
  mutate(prev = I(density.x.x>0)) %>%
  summarise(n=n(),count=sum(prev,na.rm=T),prev = mean(prev, na.rm=T)) %>% 
  mutate(lwr = binom.exact(count, n, conf.level = 0.95)[[4]], upr=binom.exact(count, n, conf.level = 0.95)[[5]])



options(scipen=999)
lplot2c = ggplot(data=newd, aes(x=visit, y=prev,fill=as.factor(`studyarm_d0`)))+
  geom_bar(stat="identity",alpha=0.9, position=position_dodge())+
  theme_classic()+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,by=0.2), limits=c(0,1))+
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.1, linewidth=1.2, position=position_dodge(0.9))+
  ylab("Pv18S prevalence")+
  xlab("Days after treatment ")+
  scale_x_continuous(breaks = 1:3)+
  scale_color_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  scale_fill_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  geom_text(aes(x=2,y=1),label=DR2c)+
  theme(axis.line = element_line(size = 0.1))+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme(legend.position = "top")
lplot2c




mod2d = gam(I(density.y.y>0) ~ `studyarm_d0`+as.factor(visit) + log(gampcr0) + 
              s(id, bs="re"),
            family=poisson(link="log"),
            data = longdt3[longdt3$gampcr0>0,])
(sm2d=summary(mod2d))
cf = cbind(sm2d$p.coeff-1.96*sm2d$se[1:5],sm2d$p.coeff+1.96*sm2d$se[1:5])
p2d = sm2d$p.pv[[2]]
p2d = ifelse(p2d<0.0001, "<0.0001",ifelse(p2d>0.9999,">0.9999",paste0("=",format(round(p2d,4),nsmall=4))))
(DR2d = paste0("RR = ",format(round(exp(mod2d$coefficients[[2]]),2),nsmall=2)," (95% CI: ",format(round(exp(cf[2,1]),2),nsmall=2),", ",format(round(exp(cf[2,2]),2),nsmall=2),"), p",p2d))

newd = longdt3[longdt3$gampcr0>0,] %>% 
  group_by(`studyarm_d0`, visit) %>%
  mutate(prev = I(density.y.y>0)) %>%
  summarise(n=n(),count=sum(prev,na.rm=T),prev = mean(prev,na.rm=T)) %>% 
  mutate(lwr = binom.exact(count, n, conf.level = 0.95)[[4]], upr=binom.exact(count, n, conf.level = 0.95)[[5]])



options(scipen=999)
lplot2d = ggplot(data=newd, aes(x=visit, y=prev,fill=as.factor(`studyarm_d0`)))+
  geom_bar(stat="identity",alpha=0.9, position=position_dodge())+
  theme_classic()+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,by=0.2), limits=c(0,1))+
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.1, linewidth=1.2, position=position_dodge(0.9))+
  ylab("Pvs25 prevalence")+
  xlab("Days after treatment")+
  scale_x_continuous(breaks = 1:3)+
  scale_color_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  scale_fill_manual(values=c("red","darkblue"),labels=c("CQ+PQ","PA+PQ"))+
  geom_text(aes(x=2,y=1),label=DR2d)+
  theme(axis.line = element_line(size = 0.1))+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme(legend.position = "top")
lplot2d





Fig_S13=ggarrange(
  lplot2a,
  lplot2c,
  lplot2b,
  lplot2d,
  nrow=2, ncol=2,
  align="hv",
  common.legend = TRUE,
  legend="top",labels=c("(A)","(B)","(C)","(D)"),
  vjust=0.6,
  hjust=-1)

Fig_S13
#### Figure 2: Infectiousness of P. vivax patients to mosquitoes before and after treatment using schizonticides (CQ and PA) with or without primaquine. 
# Add mosquito infection data ---------------------------------------------

mi= read_excel("data/Mosqito infection summary.xlsx")
colnames(mi) <- c("id",colnames(mi)[-1])
mi <- mi %>% rename(visit=visit_day)

feed_pm <- read_excel("data/Parasite density by microscopy.xlsx") %>% 
  mutate(id = `study id`, parmic0 = pv_para_dens_d0,density_0 = pv_para_dens_d0,density_1 = as.numeric(pv_para_dens_d1),density_2=as.numeric(pv_para_dens_d2),density_3 = as.numeric(pv_para_dens_d3),density_7 = as.numeric(pv_para_dens_d7),density_14 = as.numeric(pv_para_dens_d14),density_21=as.numeric(pv_para_dens_d21),density_28 = as.numeric(pv_para_dens_d28),density_35 = as.numeric(pv_para_dens_d35),density_42 = as.numeric(pv_para_dens_d42)) %>%
  dplyr::select(id, studyarm, parmic0, density_0, density_1, density_2, density_3, density_7, density_14, density_21, density_28, density_35, density_42)%>%
  pivot_longer(cols = starts_with(c("density")),
               names_to = c(".value", "visit"),
               names_sep = "_") %>%
  mutate(visit = as.numeric(visit),
         id = as.factor(id)) %>% 
  na.omit()%>% rename(dens.pm = density, dens.pm0 = parmic0)


feed_gm <- read_excel("data/Gametocyte density microscopy.xlsx")%>% 
  mutate(id = `study id`, studyarm=studyarm_d0, gammic0 = pvgam_dens_d0,density_0 = pvgam_dens_d0,density_1 = as.numeric(pvgam_dens_d1),density_2=as.numeric(pvgam_dens_d2),density_3 = as.numeric(pvgam_dens_d3),density_7 = as.numeric(pvgam_dens_d7),density_14 = as.numeric(pvgam_dens_d14),density_21=as.numeric(pvgam_dens_d21),density_28 = as.numeric(pvgam_dens_d28),density_35 = as.numeric(pvgam_dens_d35),density_42 = as.numeric(pvgam_dens_d42)) %>%
  dplyr::select(id, studyarm, gammic0, density_0, density_1, density_2, density_3, density_7, density_14, density_21, density_28, density_35, density_42)%>%
  pivot_longer(cols = starts_with(c("density")),
               names_to = c(".value", "visit"),
               names_sep = "_") %>%
  mutate(visit = as.numeric(visit),
         id = as.factor(id)) %>% 
  na.omit()%>% rename(dens.gm = density, dens.gm0 = gammic0)


feed_pp <- read_excel("data/Parasite density pcr.xlsx")%>% 
  mutate(id = `study id`, parpcr0 = PvparapuLbd_d0,density_0 = PvparapuLbd_d0,density_1 = as.numeric(PvparapuLbd_d1),density_2=as.numeric(PvparapuLbd_d2),density_3 = as.numeric(PvparapuLbd_d3),density_7 = as.numeric(PvparapuLbd_d7),density_14 = as.numeric(PvparapuLbd_d14)) %>%
  dplyr::select(id, studyarm, parpcr0, density_0, density_1, density_2, density_3, density_7, density_14)%>%
  pivot_longer(cols = starts_with(c("density")),
               names_to = c(".value", "visit"),
               names_sep = "_") %>%
  mutate(visit = as.numeric(visit),
         id = as.factor(id),
         density = ifelse(density<3,0,density)) %>% 
  na.omit()%>% rename(dens.pp = density, dens.pp0 = parpcr0)
#above qpcr densities <3 were classified as negative



feed_gp <- read_excel("data/Gametocyte density pcr.xlsx") %>% 
  mutate(id = `study id`, gampcr0 = pvs25copypul_d0,density_0 = pvs25copypul_d0,density_1 = as.numeric(pvs25copypul_d1),density_2=as.numeric(pvs25copypul_d2),density_3 = as.numeric(pvs25copypul_d3),density_7 = as.numeric(pvs25copypul_d7),density_14 = as.numeric(pvs25copypul_d14)) %>%
  dplyr::select(id, studyarm, gampcr0, density_0, density_1, density_2, density_3, density_7, density_14)%>%
  pivot_longer(cols = starts_with(c("density")),
               names_to = c(".value", "visit"),
               names_sep = "_") %>%
  mutate(visit = as.numeric(visit),
         id = as.factor(id),
         density = ifelse(density<3,0,density)) %>% 
  na.omit()%>% rename(dens.gp = density, dens.gp0 = gampcr0)
#above qpcr densities <3 were classified as negative


longfdt = feed_pm %>% 
  full_join(feed_gm, by=c("id", "studyarm", "visit")) %>%
  full_join(feed_pp, by=c("id", "studyarm", "visit")) %>%
  full_join(feed_gp, by=c("id", "studyarm", "visit")) 

fdt0 = longfdt[longfdt$visit==0,]
longfdt_0 = longfdt[longfdt$visit>0,]




mi <- mi %>% left_join(longfdt %>% 
                         dplyr::select(id,visit, studyarm, dens.pm, dens.gm, dens.pp, dens.gp), by=c("id","visit"))

mi0 = mi[mi$visit==0,]


#adjusting for any density might be hard, from those who didn't get PQ at baseline we only know one individual's densities. Need this data.

# Plots of infectious likelihood, infectiousness and intensity of infection ---------


inf = mi %>% 
  filter(visit<=2) %>%
  group_by(Arm, visit) %>%
  mutate(prev = I(Infected>0)) %>%
  summarise(n=n(),count=sum(prev),prev = mean(prev)) %>% 
  mutate(lwr = binom.exact(count, n, conf.level = 0.95)[[4]], upr=binom.exact(count, n, conf.level = 0.95)[[5]])

inf

options(scipen=999)
inf$prev <- inf$prev * 100
inf$lwr= inf$lwr*100
inf$upr =inf$upr*100
lplot5a = ggplot(data=inf, aes(x=visit, y=prev,fill=Arm, alpha=Arm))+
  geom_bar(stat="identity", position=position_dodge(),width=0.5,size=0.1)+
  theme_classic()+
  #scale_y_continuous(labels = scales::percent, breaks=seq(0,1,by=0.2), limits=c(0,1))+
  #scale_y_continuous(labels = scales*100, breaks=seq(0,100,by=20), limits=c(0,100))+
  #scale_y_continuous(breaks = seq(0, 100, by = 0.2)*100, limits = c(0, 100)) +
  scale_y_continuous(labels = function(x) format(round(x), nsmall = 0), breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.25, linewidth=0.5, position=position_dodge(0.9), alpha=0.25)+
  ylab("Proportion infectious individuals")+
  xlab("Days after treatment")+
  scale_x_continuous(breaks = 0:2)+
  scale_color_manual(values=c("red","red","darkblue","darkblue"))+
  scale_fill_manual(values=c("red","red","darkblue","darkblue"))+
  scale_alpha_manual(values=c(1,0.25,1,0.25))+
  theme(legend.title = element_blank(),legend.position = "top",strip.background = element_blank(),  # This line removes the facet box
        strip.text = element_text(size = 12))+
  facet_wrap(~Arm, nrow=1) +
  theme(axis.line = element_line(size = 0.5)) 
lplot5a


inf2 = mi %>% 
  filter(visit<=2) %>%
  ungroup() %>%
  mutate(prev = (Infected)/(Dissected)) 



options(scipen=999)
inf2$prev=inf2$prev*100
inf2$lwr= inf2$lwr*100
inf2$upr =inf2$upr*100
lplot5b = ggplot(data=inf2, aes(x=visit, y=prev,fill=Arm, alpha=Arm, group=visit))+
  geom_boxplot( position=position_dodge(),width=0.5)+
  geom_jitter(width=0.25, alpha=0.4)+
  theme_classic()+
  #scale_y_continuous(labels = scales::percent, breaks=seq(0,1,by=0.2), limits=c(0,1))+
  #scale_y_continuous(breaks = seq(0, 1, by = 0.2) * 100, limits = c(0, 100)) +
  scale_y_continuous(labels = function(x) format(round(x), nsmall = 0), breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  ylab("Proportion infected mosquitos")+
  xlab("Days after treatment")+
  scale_x_continuous(breaks = 0:2)+
  scale_color_manual(values=c("red","red","darkblue","darkblue"))+
  scale_fill_manual(values=c("red","red","darkblue","darkblue"))+
  scale_alpha_manual(values=c(1,0.25,1,0.25))+
  theme(legend.title = element_blank(),legend.position = "top",strip.background = element_blank(),  # This line removes the facet box
        strip.text = element_text(size = 12))+
  facet_wrap(~Arm, nrow=1) +
  theme(axis.line = element_line(size = 0.5))
lplot5b




inf3 = mi %>% 
  filter(visit<=2 & Infected>0) %>%
  mutate(meano = `Total Oocyst Density`/Infected) 



options(scipen=999)
lplot5c = ggplot(data=inf3, aes(x=visit, y=meano,fill=Arm, alpha=Arm, group=visit))+
  geom_boxplot( position=position_dodge(),width=0.5)+
  geom_jitter(width=0.25, alpha=0.4)+
  theme_classic()+
 # ylab("Mean oocyst density per infected mosquito")+
  ylab("Mean oocyst density")+
  xlab("Days after treatment")+
  scale_x_continuous(limits=c(-0.3,2),breaks=0:2)+
  scale_color_manual(values=c("red","red","darkblue","darkblue"))+
  scale_fill_manual(values=c("red","red","darkblue","darkblue"))+
  scale_alpha_manual(values=c(1,0.25,1,0.25))+
  theme(legend.title = element_blank(),legend.position = "top",strip.background = element_blank(),  # This line removes the facet box
        strip.text = element_text(size = 12))+
  facet_wrap(~Arm, nrow=1)+
  scale_y_log10(breaks=c(1,10,100), limits=c(1,250)) +
  theme(axis.line = element_line(size = 0.5)) 
 
lplot5c


fig_2=ggarrange(
  lplot5a,
  lplot5b,
  lplot5c,
  nrow=3,
  align="hv",
  common.legend = TRUE,
  legend="top",labels=c("(A)","(B)","(C)","(D)"),
  vjust=0.6,
  hjust=-1
)

fig_2
## fig_S15
# for Appendix 1. 5: The relative reduction in infectivity to mosquitoes before and after treatment in P. vivax. 

# Plots and tables for relative reduction in densities, infectious --------

##### Prop Infected mosquitos #####
mi$prop = mi$Infected/mi$Dissected
mi$arm = factor(mi$Arm)
ref = levels(mi$arm)
mi$arm2 = mi$arm 
mi$arm = relevel(mi$arm,ref[1])
mi$visit = as.factor(mi$visit)
mi$id = as.factor(mi$id)
mi = mi %>% 
  group_by(id) %>% 
  arrange(visit) %>% 
  mutate(dens.pm0 = first(dens.pm),
         dens.gm0 = first(dens.gm),
         dens.pp0 = first(dens.pp),
         dens.gp0 = first(dens.gp))


mod_control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

mod1 = bglmer(cbind(Infected, Dissected) ~ visit*arm  + (1 | id), 
       data = mi, 
       family = binomial(link = "log"),
       fixef.prior = normal(2), 
       cov.prior = wishart(2),
       control=mod_control)

sb1=summary(mod1)
sb1$p.coeff = sb1$coefficients[,1]
sb1$se = sb1$coefficients[,2]
sb1$p.pv = sb1$coefficients[,4]

newd0 = data.frame(arm=levels(mi$arm)) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(mi$arm)[1],"visit1", paste0("visit1:arm", arm)),
    beta = sb1$p.coeff[[values]],
    se = sb1$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


day2newd0 = data.frame(arm=levels(mi$arm)) 
day2newd0 = day2newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(mi$arm)[1],"visit2", paste0("visit2:arm", arm)),
    beta = sb1$p.coeff[[values]],
    se = sb1$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


day3newd0 = data.frame(arm=levels(mi$arm)) 
day3newd0 = day3newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(mi$arm)[1],"visit3", paste0("visit3:arm", arm)),
    beta = sb1$p.coeff[[values]],
    se = sb1$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^3+se^3)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


pdat = rbind(
  newd0 %>% mutate(day=1),
  day2newd0 %>% mutate(day=2),
  day3newd0 %>% mutate(day=3)
)  

pdat$ci_lower = ifelse(pdat$ci_lower<0, 0, pdat$ci_lower)
pdat$arm = factor(pdat$arm,levels = c("CQ","PA","CQ+PQ","PA+PQ"))
#custom_colors <- c("#B3CDE3FF","#FBB4AEFF")  # Define custom colors for the slices
pdat$TRA=pdat$TRA*100
pdat$ci_lower=pdat$ci_lower*100
  pdat$ci_upper= pdat$ci_upper*100
plot6a = ggplot(data=pdat, aes(x=day, y=TRA, fill=arm, group=arm))+
  geom_bar(stat="identity", position = position_dodge(), width=0.75, aes(alpha=arm))+
  geom_point(position = position_dodge(0.8),alpha=1)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(0.8), width=0.25, linewidth=0.5,alpha=0.25)+
  theme_classic()+
  #theme(legend.position = "top")+
  #scale_y_continuous(labels=scales::percent, limits=c(0,1),breaks=seq(0,1,0.2))+
  scale_y_continuous(labels = function(x) format(round(x), nsmall = 0), breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  scale_x_continuous(limits=c(0,4),breaks=1:3)+
  ylab("Relative reduction(%) in \nproportion infected mosquitos")+
  facet_wrap(~arm, nrow=1)+
  xlab("Days after treatment")+
  #scale_x_continuous(limits=c(-0.3,2),breaks=0:2)+
 ## scale_color_manual(values=c("#FBB4AEFF","#FBB4AEFF","#B3CDE3FF","#B3CDE3FF"))+
  #scale_fill_manual(values=c("#FBB4AEFF","#FBB4AEFF","#B3CDE3FF","#B3CDE3FF"))+
  scale_color_manual(values=c("red","red","darkblue","darkblue"))+
  scale_fill_manual(values=c("red","red","darkblue","darkblue"))+
  theme(axis.line = element_line(size = 0.5)) +
  theme(legend.title = element_blank(),legend.position = "top",strip.background = element_blank(),  # This line removes the facet box
     strip.text = element_text(size = 12))+
  #scale_alpha_manual(values=c(1,1,1,1))
  
  scale_alpha_manual(values=c(1,0.25,1,0.25))
  
plot6a



newd0 = newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day2newd0 = day2newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day3newd0 = day3newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

current = ref[1]

newd =  data.frame(arm=levels(mi$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visit1:arm", arm), reference = current,
    pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day2newd =  data.frame(arm=levels(mi$arm)) 

day2newd = day2newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visit2:arm", arm), reference = current,
    pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day3newd =  data.frame(arm=levels(mi$arm)) 

day3newd = day3newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visit3:arm", arm), reference = current,
    pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)



for (i in 2:length(ref)){
  
  current = ref[i]
  
  mi$arm = relevel(mi$arm2,current)
  
  mod1 = bglmer(cbind(Infected, Dissected-Infected) ~ visit*arm + (1 | id), 
                data = mi, 
                family = binomial(link = "log"),
                fixef.prior = normal(2), 
                cov.prior = wishart(2),
                control=mod_control)
  
  sb1=summary(mod1)
  sb1$p.coeff = sb1$coefficients[,1]
  sb1$se = sb1$coefficients[,2]
  sb1$p.pv = sb1$coefficients[,4]
  
  
  newd1 =  data.frame(arm=levels(mi$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visit1:arm", arm), reference = current,
      pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
  day2newd1 =  data.frame(arm=levels(mi$arm)) 
  
  day2newd1 = day2newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visit2:arm", arm), reference = current,
      pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day2newd = rbind(day2newd, day2newd1)
  
  
  day3newd1 =  data.frame(arm=levels(mi$arm)) 
  
  day3newd1 = day3newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visit3:arm", arm), reference = current,
      pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day3newd = rbind(day3newd, day3newd1)
  
  
  print(i)
}



newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>%
               dplyr::select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix1 <- pivot_wider(newd, names_from = arm, values_from = result)

result_matrix1[1,3]=""
result_matrix1[1,4]=""
result_matrix1[1,5]=""
result_matrix1[2,4]=""
result_matrix1[2,5]=""
result_matrix1[3,5]=""
result_matrix1
#flextable(result_matrix1) %>% autofit() %>% save_as_docx( path = "output/Prop_infected_reduction_day1.docx")



day2newd = rbind(day2newd,
                 day2newd0 %>% mutate(reference = arm, result=RRG) %>% 
                   dplyr::select(arm, reference, result)) %>% arrange(reference, arm)

day2newd = unique(as.data.frame(day2newd))

result_matrix12 <- pivot_wider(day2newd, names_from = arm, values_from = result)

result_matrix12[1,3]=""
result_matrix12[1,4]=""
result_matrix12[1,5]=""
result_matrix12[2,4]=""
result_matrix12[2,5]=""
result_matrix12[3,5]=""
result_matrix12
#flextable(result_matrix12) %>% autofit() %>% save_as_docx( path = "output/Prop_infected_reduction_day2.docx")



day3newd = rbind(day3newd,
                 day3newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day3newd = unique(as.data.frame(day3newd))

result_matrix13 <- pivot_wider(day3newd, names_from = arm, values_from = result)
result_matrix13[1,3]=""
result_matrix13[1,4]=""
result_matrix13[1,5]=""
result_matrix13[2,4]=""
result_matrix13[2,5]=""
result_matrix13[3,5]=""
#flextable(result_matrix13) %>% autofit() %>% save_as_docx( path = "output/Prop_infected_reduction_day3.docx")
result_matrix13



##### Reduction in individuals infectiousness #####

mi$prop = as.numeric(mi$Infected>1)
mi$arm = factor(mi$Arm)
ref = levels(mi$arm)
mi$arm2 = mi$arm 
mi$arm = relevel(mi$arm,ref[1])
mi$visit = as.factor(mi$visit)
mi$id = as.factor(mi$id)



mod_control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

mod1 = bglmer(prop ~ visit*arm + (1 | id), 
              data = mi, 
              family = poisson(link = "log"),
              fixef.prior = normal(2), 
              cov.prior = wishart(2),
              control=mod_control)

sb1=summary(mod1)
sb1$p.coeff = sb1$coefficients[,1]
sb1$se = sb1$coefficients[,2]
sb1$p.pv = sb1$coefficients[,4]

newd0 = data.frame(arm=levels(mi$arm)) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(mi$arm)[1],"visit1", paste0("visit1:arm", arm)),
    beta = sb1$p.coeff[[values]],
    se = sb1$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


day2newd0 = data.frame(arm=levels(mi$arm)) 
day2newd0 = day2newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(mi$arm)[1],"visit2", paste0("visit2:arm", arm)),
    beta = sb1$p.coeff[[values]],
    se = sb1$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


day3newd0 = data.frame(arm=levels(mi$arm)) 
day3newd0 = day3newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(mi$arm)[1],"visit3", paste0("visit3:arm", arm)),
    beta = sb1$p.coeff[[values]],
    se = sb1$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^3+se^3)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


pdat = rbind(
  newd0 %>% mutate(day=1),
  day2newd0 %>% mutate(day=2),
  day3newd0 %>% mutate(day=3)
)  

pdat$ci_lower = ifelse(pdat$ci_lower<0, 0, pdat$ci_lower)
long$day <- factor(long$day, levels = c("Day 0", "Day 1", "Day 2", "Day 3", "Day 7", "Day 14", "Day 21", "Day 28"))
table(pdat$arm)
pdat$arm = factor(pdat$arm,levels = c("CQ","PA","CQ+PQ","PA+PQ"))
pdat$TRA=pdat$TRA*100
pdat$ci_lower=pdat$ci_lower*100
pdat$ci_upper=pdat$ci_upper*100
plot6b = ggplot(data=pdat, aes(x=day, y=TRA, fill=arm, group=arm))+
  geom_bar(stat="identity", position = position_dodge(), width=0.75, aes(alpha=arm))+
  geom_point(position = position_dodge(0.8),alpha=1)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(0.8), width=0.25, linewidth=0.5,alpha=0.25)+
  theme_classic()+
  theme(legend.position = "top")+
  #scale_y_continuous(labels=scales::percent, limits=c(0,1),breaks=seq(0,1,0.2))+
  scale_y_continuous(labels = function(x) format(round(x), nsmall = 0), breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  scale_x_continuous(limits=c(0,4),breaks=1:3)+  ylab("Relative reduction (%) in individual's\nlikelihood of infecting mosquitos")+
  facet_wrap(~arm, nrow=1)+
  xlab("Days after treatment")+
  #scale_color_manual(values=c("#FBB4AEFF","#FBB4AEFF","#B3CDE3FF","#B3CDE3FF"))+
  #scale_fill_manual(values=c("#FBB4AEFF","#FBB4AEFF","#B3CDE3FF","#B3CDE3FF"))+
  #scale_color_manual(values=c("red","red","darkblue","darkblue"))+
  #scale_fill_manual(values=c("red","red","darkblue","darkblue"))+
  #scale_alpha_manual(values=c(1,1,1,1))
  #scale_x_continuous(limits=c(-0.3,2),breaks=0:2)+
  scale_color_manual(values=c("red","red","darkblue","darkblue"))+
  scale_fill_manual(values=c("red","red","darkblue","darkblue"))+
  theme(axis.line = element_line(size = 0.5)) +
  theme(legend.title = element_blank(),legend.position = "top",strip.background = element_blank(),  # This line removes the facet box
        strip.text = element_text(size = 12))+
  scale_alpha_manual(values=c(1,0.25,1,0.25))

plot6b


newd0 = newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day2newd0 = day2newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day3newd0 = day3newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

current = ref[1]

newd =  data.frame(arm=levels(mi$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visit1:arm", arm), reference = current,
    pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day2newd =  data.frame(arm=levels(mi$arm)) 

day2newd = day2newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visit2:arm", arm), reference = current,
    pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day3newd =  data.frame(arm=levels(mi$arm)) 

day3newd = day3newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visit3:arm", arm), reference = current,
    pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)



for (i in 2:length(ref)){
  
  current = ref[i]
  
  mi$arm = relevel(mi$arm2,current)
  
  mod1 = bglmer(prop ~ visit*arm + (1 | id), 
                data = mi, 
                family = poisson(link = "log"),
                fixef.prior = normal(2), 
                cov.prior = wishart(2),
                control=mod_control)
  
  sb1=summary(mod1)
  sb1$p.coeff = sb1$coefficients[,1]
  sb1$se = sb1$coefficients[,2]
  sb1$p.pv = sb1$coefficients[,4]
  
  
  newd1 =  data.frame(arm=levels(mi$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visit1:arm", arm), reference = current,
      pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
  day2newd1 =  data.frame(arm=levels(mi$arm)) 
  
  day2newd1 = day2newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visit2:arm", arm), reference = current,
      pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day2newd = rbind(day2newd, day2newd1)
  
  
  day3newd1 =  data.frame(arm=levels(mi$arm)) 
  
  day3newd1 = day3newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visit3:arm", arm), reference = current,
      pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day3newd = rbind(day3newd, day3newd1)
  
  
  print(i)
}



newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>% dplyr:: select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix1 <- pivot_wider(newd, names_from = arm, values_from = result)

result_matrix1[1,3]=""
result_matrix1[1,4]=""
result_matrix1[1,5]=""
result_matrix1[2,4]=""
result_matrix1[2,5]=""
result_matrix1[3,5]=""

#flextable(result_matrix1) %>% autofit() %>% save_as_docx( path = "output/Likelihood_individual_infecting_reduction_day1.docx")
result_matrix1


day2newd = rbind(day2newd,
                 day2newd0 %>% mutate(reference = arm, result=RRG) %>% dplyr::select(arm, reference, result)) %>% arrange(reference, arm)

day2newd = unique(as.data.frame(day2newd))

result_matrix12 <- pivot_wider(day2newd, names_from = arm, values_from = result)

result_matrix12[1,3]=""
result_matrix12[1,4]=""
result_matrix12[1,5]=""
result_matrix12[2,4]=""
result_matrix12[2,5]=""
result_matrix12[3,5]=""

#flextable(result_matrix12) %>% autofit() %>% save_as_docx( path = "output/Likelihood_individual_infecting_reduction_day2.docx")
result_matrix12


day3newd = rbind(day3newd,
                 day3newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day3newd = unique(as.data.frame(day3newd))

result_matrix13 <- pivot_wider(day3newd, names_from = arm, values_from = result)
result_matrix13[1,3]=""
result_matrix13[1,4]=""
result_matrix13[1,5]=""
result_matrix13[2,4]=""
result_matrix13[2,5]=""
result_matrix13[3,5]=""
#flextable(result_matrix13) %>% autofit() %>% save_as_docx( path = "output/Likelihood_individual_infecting_day3.docx")
result_matrix13
### Appendix 1. 5: The relative reduction in infectivity to mosquitoes before and after treatment in P. vivax. 
#### Reduction in mean oocysts amongst those infected 
misub = mi 
misub$prop = misub$`Total Oocyst Density`
misub$arm = factor(misub$Arm)
ref = levels(misub$arm)
misub$arm2 = misub$arm 
misub$arm = relevel(misub$arm,ref[1])
misub$visit = as.factor(misub$visit)
misub$id = as.factor(misub$id)
str(misub$visit)  
str(misub$arm)    
misub$visit <- as.factor(misub$visit)
misub$arm <- as.factor(misub$arm)
misub$visit <- relevel(misub$visit, ref = "0")  # Set "0" as the reference level
misub$arm <- relevel(misub$arm, ref = "CQ")     

mod_control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

# Create an interaction term manually
misub$visit_arm <- interaction(misub$visit, misub$arm)
#####
# Set sum-to-zero contrasts
options(contrasts = c("contr.sum", "contr.poly"))

#######
mod1 = bglmer(I(prop/Infected) ~ visit*arm + (1 | id), 
              data = misub, 
              family = gaussian(link="log"),
              fixef.prior = normal(2), 
              cov.prior = wishart(2),
              control=mod_control)


sb1=summary(mod1)
sb1$p.coeff = sb1$coefficients[,1]
sb1$se = sb1$coefficients[,2]
sb1$p.pv = sb1$coefficients[,4]
newd0 = data.frame(arm=levels(misub$arm)) 
newd0 = data.frame(arm=levels(misub$arm)[1:4]) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(misub$arm)[1],"visit1", paste0("visit1:arm", arm)),
    beta = if (values %in% names(sb1$p.coeff)) sb1$p.coeff[[values]] else NA,  # Check if values exists
    se = if (values %in% names(sb1$se)) sb1$se[[values]] else NA  # Check if values exists
    #beta = sb1$p.coeff[[values]],
    #se = sb1$se[[values]]
    ) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))



pdat = rbind(
  newd0 %>% mutate(day=1)
  )

pdat[4,"values"]="visit1"
pdat[4,"day"]= 1
pdat$arm[4]="PA+PQ"
pdat$TRA[4]=1
pdat$ci_lower[4]=1
pdat$ci_upper[4]=1
pdat$pval[4]="<0.0001"
########
pdat[5,"values"]="visit2"
pdat[5,"day"]= 2
pdat$arm[5]="PA+PQ"
pdat$TRA[5]=1
pdat$ci_lower[5]=1
pdat$ci_upper[5]=1
pdat$pval[5]="<0.0001"
######
pdat[6,"values"]="visit3"
pdat[6,"day"]= 3
pdat$arm[6]="PA+PQ"
pdat$TRA[6]=1
pdat$ci_lower[6]=1
pdat$ci_upper[6]=1
pdat$pval[6]="<0.0001"
####
######## CQ
pdat[7,"values"]="visit2"
pdat[7,"day"]= 2
pdat$arm[7]="CQ"
pdat$TRA[7]=1
pdat$ci_lower[7]=1
pdat$ci_upper[7]=1
pdat$pval[7]="<0.0001"
######
pdat[8,"values"]="visit3"
pdat[8,"day"]= 3
pdat$arm[8]="CQ"
pdat$TRA[8]=1
pdat$ci_lower[8]=1
pdat$ci_upper[8]=1
pdat$pval[8]="<0.0001"
##### PA
########
pdat[9,"values"]="visit2"
pdat[9,"day"]= 2
pdat$arm[9]="PA"
pdat$TRA[9]=1
pdat$ci_lower[9]=1
pdat$ci_upper[9]=1
pdat$pval[9]="<0.0001"
######
pdat[10,"values"]="visit3"
pdat[10,"day"]= 3
pdat$arm[10]="PA"
pdat$TRA[10]=1
pdat$ci_lower[10]=1
pdat$ci_upper[10]=1
pdat$pval[10]="<0.0001"
#### CQ+PQ
########
pdat[11,"values"]="visit2"
pdat[11,"day"]= 2
pdat$arm[11]="CQ+PQ"
pdat$TRA[11]=1
pdat$ci_lower[11]=1
pdat$ci_upper[11]=1
pdat$pval[11]="<0.0001"
######
pdat[12,"values"]="visit3"
pdat[12,"day"]= 3
pdat$arm[12]="CQ+PQ"
pdat$TRA[12]=1
pdat$ci_lower[12]=1
pdat$ci_upper[12]=1
pdat$pval[12]="<0.0001"


pdat$ci_lower = ifelse(pdat$ci_lower<0, 0, pdat$ci_lower)

pdat$TRA=pdat$TRA*100
pdat$ci_lower=pdat$ci_lower*100
pdat$ci_upper=pdat$ci_upper*100
plot6c = ggplot(data=pdat, aes(x=day, y=TRA, fill=arm, group=arm))+
  geom_bar(stat="identity", position = position_dodge(), width=0.75, aes(alpha=arm))+
  geom_point(position = position_dodge(0.8),alpha=1)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(0.8), width=0.25, linewidth=0.5,alpha=0.25)+
  theme_classic()+
  theme(legend.position = "top")+
  #scale_y_continuous(labels=scales::percent, limits=c(0,1),breaks=seq(0,1,0.2))+
  scale_y_continuous(labels = function(x) format(round(x), nsmall = 0), breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  scale_x_continuous(limits=c(0,4),breaks=1:3)+
  ylab("Relative reduction (%) in oocysts \n per infected mosquito")+
  facet_wrap(~arm, nrow=1)+
  xlab("Days after treatment")+
  #scale_x_continuous(limits=c(-0.3,2),breaks=0:2)+
  scale_color_manual(values=c("red","red","darkblue","darkblue","grey"))+
  scale_fill_manual(values=c("red","red","darkblue","darkblue","grey"))+
  scale_alpha_manual(values=c(1,0.25,1,0.25,1))+
  theme(axis.line = element_line(size = 0.5)) +
  theme(legend.title = element_blank(),legend.position = "top",strip.background = element_blank(),  # This line removes the facet box
        strip.text = element_text(size = 12))+
  scale_alpha_manual(values=c(1,0.25,1,0.25))

plot6c

#ggsave("output/reduction_prop_infected.png",device="png", width=12, height=8, dpi=600)

newd0 = data.frame(arm = levels(misub$arm)[1:4]) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm == levels(misub$arm)[1], "visit1", paste0("visit1:arm", arm)),
    beta = if (values %in% names(sb1$p.coeff)) sb1$p.coeff[[values]] else NA,  # Check if values exists
    se = if (values %in% names(sb1$se)) sb1$se[[values]] else NA  # Check if values exists
  ) %>%
  ungroup() %>%
  mutate(
    n = 1:n(), 
    beta = ifelse(n == 1, beta, first(beta) + beta), 
    se = ifelse(n == 1, se, sqrt(first(se)^2 + se^2)),
    TRA = ifelse(arm == "PA+PQ", 1, 1 - exp(beta)),  # Set TRA to 100 for PA+PQ
    ci_lower = ifelse(arm == "PA+PQ", 1, 1 - exp(beta + 1.96 * se)),  # Set ci_lower to 100 for PA+PQ
    ci_upper = ifelse(arm == "PA+PQ", 1, 1 - exp(beta - 1.96 * se)),  # Set ci_upper to 100 for PA+PQ
    p.pv = ifelse(arm == "PA+PQ", 0.00001, pnorm(abs(beta / se), lower.tail = F) * 2),  # Set p.pv to 0.00001 for PA+PQ
    pval = ifelse(arm == "PA+PQ", "<0.0001",  # Set pval to "<0.0001" for PA+PQ
                  ifelse(p.pv < 0.0001, "<0.0001", 
                         ifelse(p.pv > 0.9999, ">0.9999", 
                                paste0("=", format(round(p.pv, 4), nsmall = 4)))))
  )

newd0 = newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)


current = ref[1]

newd =  data.frame(arm=levels(misub$arm)[1:3]) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visit1:arm", arm), reference = current,
    #values = paste0("visit1:arm", arm), reference = current,
    pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


##################
day2newd0 = data.frame(arm=levels(misub$arm)) 
day2newd0 = day2newd0 %>% 
  group_by(arm) %>%
  mutate(
 values = ifelse(arm==levels(misub$arm)[1],"visit2", paste0("visit2:arm", arm)),
   beta = sb1$p.coeff[[values]],
 se = sb1$se[[values]]) %>%
 ungroup() %>%
 mutate(n = 1:n(), 
 beta=ifelse(n==1,beta,first(beta)+beta), 
   se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
    TRA = 1-exp(beta),
       ci_lower = 1-exp(beta+1.96*se),
        ci_upper = 1-exp(beta-1.96*se),
          p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

day2newd0 = day2newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")) %>%
  dplyr::select(arm, RRG)

day2newd =  data.frame(arm=levels(misub$arm)) 
 day2newd = day2newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
     values = paste0("visit2:arm", arm), reference = current,
   pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
dplyr::select(arm, reference, result)




for (i in 2:length(ref)){
  
  current = ref[i]
  
  misub$arm = relevel(misub$arm2,current)
  
  mod_control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  
  mod1 = bglmer(I(prop/Infected) ~ visit*arm + (1 | id), 
                offset=log(Infected),
                data = misub, 
                family = gaussian(link="log"),
                fixef.prior = normal(2), 
                cov.prior = wishart(2),
                control=mod_control)
  
  sb1=summary(mod1)
  sb1$p.coeff = sb1$coefficients[,1]
  sb1$se = sb1$coefficients[,2]
  sb1$p.pv = sb1$coefficients[,4]  
  
  newd1 =  data.frame(arm=levels(misub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visit1:arm", arm), reference = current,
      pval = sb1$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
 
  print(i)
}



newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>%dplyr:: select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix1 <- pivot_wider(newd, names_from = arm, values_from = result)
result_matrix1$`PA+PQ`=NA

result_matrix1[1,3]=""
result_matrix1[1,4]=""
result_matrix1[1,5]=""
result_matrix1[2,4]=""
result_matrix1[2,5]=""
result_matrix1[3,5]=""

result_matrix1[4,]=""
result_matrix1$reference[4]="PA+PQ"
result_matrix1$`PA+PQ`[4]="No infected mosquitos"

#flextable(result_matrix1) %>% autofit() %>% save_as_docx( path = "output/mean_oocysts_reduction.docx")


fig_S15=ggarrange(
  plot6b,
  plot6a,
  plot6c,
  nrow=3,
  align="hv",
  common.legend = TRUE,
  legend="top",labels=c("(A)","(B)","(C)","(D)"),
  vjust=0.6,
  hjust=-1
)

fig_S15
## for faver clearance
# Plot the identity bar chart with facets and specified colors
fever_clear=read_xlsx("data/fever_clearance.xlsx")
fig_S11= ggplot(fever_clear, aes(x = as.factor(day), y = prop, color = Arm, fill = Arm)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7,show.legend = FALSE) +  # Adjust width for space
  labs(title = "",
       x = "Follow up days",
       y = "Fever clearance (%)") +
  facet_wrap(~ Arm) +
  scale_y_continuous(breaks = seq(0,100, by = 20)) + 
  scale_fill_manual(values = c("red", "blue")) +  
  scale_color_manual(values = c("red", "blue")) +  
  theme_pubr() +
  theme(axis.ticks = element_line(size = 0.1),  
        axis.ticks.length = unit(0.05, "cm"),  
        axis.line = element_line(size = 0.15),  
        strip.background = element_blank(),  
        #strip.text = element_text(size = 8, face = "bold"), 
        strip.text = element_text(size = 6),
        axis.title = element_text(size=8),
        axis.text.x = element_text(size = 6),  
        axis.text.y = element_text(size = 6))   
fig_S11
