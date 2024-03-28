#设定读取地址
setwd("C:/Users/Eliane andar/Desktop")
library(readxl)
library(pROC)
library(reportROC)
library(gmodels)
library(MASS)
library(tableone)
library(MatchIt)
library(mice)
library(autoReg)
library(rms)
library(caret)
library(rmda)
library(dplyr)
library(rrtable)
library(Hmisc)
library(broom)
data1 <- data.frame(read_xlsx("BCR TIME-1year.xlsx", sheet=1))
str(data1)
shapiro.test(as.numeric(data1$age))#
shapiro.test(as.numeric(data1$PSA))
shapiro.test(as.numeric(data1$PSAD))
shapiro.test(as.numeric(data1$F.T))
shapiro.test(as.numeric(data1$prostate_volume))
shapiro.test(as.numeric(data1$MTD))
shapiro.test(as.numeric(data1$UD))
shapiro.test(as.numeric(data1$depth))#
shapiro.test(as.numeric(data1$CCL))#
shapiro.test(as.numeric(data1$posi.core.percent))
shapiro.test(as.numeric(data1$surgical.experience))
shapiro.test(as.numeric(data1$API))#
shapiro.test(as.numeric(data1$PD))#
shapiro.test(as.numeric(data1$PCI))
shapiro.test(as.numeric(data1$symphysis.angle))
shapiro.test(as.numeric(data1$ISD))
shapiro.test(as.numeric(data1$BFW))#
shapiro.test(as.numeric(data1$maximal.prostate.width))
shapiro.test(as.numeric(data1$lower.conjugate.of.pelvic.midplane))#
shapiro.test(as.numeric(data1$urethral.width))#

vars <- c("PCI","age","PSA","EPE","cT","ISUP","prostate_volume","PIRADS","MTD","UD","Apex.depth","CCL","posi.core.percent","pSM","pT","surgical.experience","API","PD","BFW","maximal.prostate.width","ISD","symphysis.angle","lower.conjugate.of.pelvic.midplane","urethral.width","PSAD")
tableone <- CreateTableOne(vars = vars, data = data1)
table1<-print(tableone,nonnormal = c("PSA","MTD","UD","CCL","surgical.experience","posi.core.percent","PSAD","PCI","prostate_volume","symphysis.angle","ISD","maximal.prostate.width"),showAllLevels = TRUE)
write.csv(table1,file="table1.csv")

library(survival)
library(forgein)
library(survminer)
attach(melanom)
names(melanom)


dataframe<- data.frame(read_xlsx("BCR TIME-1year.xlsx", sheet=1))
str(dataframe)
res.cut <- surv_cutpoint(dataframe, time = "BCR.TIME2", event = "BCR.1year",
                         variables = c("PSA", "MTD", "UD", "Apex.depth" ,"prostate_volume",
                                       "CCL","posi.core.percent","surgical.experience",
                                       "API","PD","PCI","symphysis.angle","ISD","BFW","IFI",
                                       "maximal.prostate.width","lower.conjugate.of.pelvic.midplane",
                                       "urethral.width","age")) 
summary(res.cut)
res.cut2 <- surv_categorize(res.cut)
head(res.cut2)

#age
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ age, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 


p <-ggsurvplot(fit, data = res.cut2,  
               conf.int = TRUE,  pval = TRUE, 
               surv.median.line = "hv",
               xlab = "Time in months",   
               ylab = "BCR-free survival",
               break.time.by = 10,
               ggtheme = theme_bw(),
               legend.title = "Group",  
               legend.labs = c("age ≤ 73 years", "age > 73 years"), 
               legend = c(0.8, 0.1))

p$plot <- p$plot +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))

p$plot <- p$plot +
  theme(plot.margin = margin(2, 0.2, 1, 1, "cm"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14))

print(p)   
#PSA
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ PSA, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 


#MTD
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ MTD, data = data1)
ggsurvplot(fit, data = data1,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 

           
p <-ggsurvplot(fit, data = data1,  
               conf.int = TRUE,  pval = TRUE, 
               surv.median.line = "hv",
               xlab = "Time in months",   
               ylab = "BCR-free survival",
               break.time.by = 10,
               ggtheme = theme_bw(),
               legend.title = "Group",   
               legend.labs = c("Maximum tumor diameter ≤ 23.9mm", "Maximum tumor diameter > 23.9mm"), 
               legend = c(0.6, 0.1))

p$plot <- p$plot +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))

p$plot <- p$plot +
  theme(plot.margin = margin(2, 0.2, 1, 1, "cm"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16))

print(p)  

#maximal.prostate.width 
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ maximal.prostate.width, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 
     
#UD
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ UD, data = data1)
ggsurvplot(fit, data = data1,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 

p <-ggsurvplot(fit, data = data1,  
               conf.int = TRUE,  pval = TRUE, 
               surv.median.line = "hv",
               xlab = "Time in months",   
               ylab = "BCR-free survival",
               break.time.by = 10,
               ggtheme = theme_bw(),
               legend.title = "Group",  
               legend.labs = c("UD ≤ 4mm", "UD > 4mm"), 
               legend = c(0.8, 0.1))

p$plot <- p$plot +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))

p$plot <- p$plot +
  theme(plot.margin = margin(2, 0.2, 1, 1, "cm"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16))

print(p) 
#Apex depth 0.12  median not
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ Apex.depth, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 

#CCL  median not
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ CCL, data = data1)
ggsurvplot(fit, data = data1,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 

p <-ggsurvplot(fit, data = data1,  
               conf.int = TRUE,  pval = TRUE, 
               surv.median.line = "hv",
               xlab = "Time in months",  
               ylab = "BCR-free survival",
               break.time.by = 10,
               ggtheme = theme_bw(),
               legend.title = "Group",   
               legend.labs = c("CCL ≤ 19.9mm", "CCL > 19.9mm"), 
               legend = c(0.8, 0.1))

p$plot <- p$plot +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))

p$plot <- p$plot +
  theme(plot.margin = margin(2, 0.2, 1, 1, "cm"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16))

print(p) 

fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ API, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 

fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ PD, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 
#PCI  median not 
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ PCI, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 
#symphysis.angle
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ symphysis.angle, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 
#ISD  median not
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ ISD, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 

#BFW p0.098 median not
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ BFW, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 
#IFI median not
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ IFI, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 

#prostate_volume 0.0034
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ prostate.volume, data = data1)
str(data1)
ggsurvplot(fit, data = data1,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 

p <-ggsurvplot(fit, data = data1,  
               conf.int = TRUE,  pval = TRUE, 
               surv.median.line = "hv",
               xlab = "Time in months",   
               ylab = "BCR-free survival",
               break.time.by = 10,
               ggtheme = theme_bw(),
               legend.title = "Group",   
               legend.labs = c("prostate volume ≤ 31mm", "prostate volume > 31mm"), 
               legend = c(0.7, 0.1))

p$plot <- p$plot +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))

p$plot <- p$plot +
  theme(plot.margin = margin(2, 0.2, 1, 1, "cm"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16))

print(p) 

#posi.core.percent 0.081 median not
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ posi.core.percent, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue")

#Surgical experience 0.00 median not
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ Surgical experience, data = data1)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue")


#plower.conjugate.of.pelvic.midplane 0.13
fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ lower.conjugate.of.pelvic.midplane, data = res.cut2)
ggsurvplot(fit, data = res.cut2,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 

#urethral.width p0.011

fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ urethral.width, data = data1)
p <-ggsurvplot(fit, data = data1,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv",
           xlab = "Time in months",   
           ylab = "BCR-free survival",
           break.time.by = 10,
           ggtheme = theme_bw(),
           legend.title = "Group",   
           legend.labs = c("PSA", "PSA2"), 
           legend = c(0.7, 0.1),
           risk.table = TRUE, palette = "hue") 

p <-ggsurvplot(fit, data = data1,  
               conf.int = TRUE,  pval = TRUE, 
               surv.median.line = "hv",
               xlab = "Time in months",   
               ylab = "BCR-free survival",
               break.time.by = 10,
               ggtheme = theme_bw(),
               legend.title = "Group",   
               legend.labs = c("urethral width ≤ 1.08mm", "urethral width > 1.08mm"), 
               legend = c(0.7, 0.1))

p$plot <- p$plot +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))

p$plot <- p$plot +
  theme(plot.margin = margin(2, 0.2, 1, 1, "cm"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16))

print(p) 



fit <- survfit(Surv(BCR.TIME2, BCR.1year) ~ 1, data = data1)
p <-ggsurvplot(fit, data = data1,  
               conf.int = TRUE,  pval = TRUE, 
               surv.median.line = "hv",
               xlab = "Time in months",   
               ylab = "BCR-free survival",
               break.time.by = 10,
               ggtheme = theme_bw(),
               legend.title = "Group", 
               legend.labs = c("Overall. Suvival probability at 100 months"), 
               legend = c(0.6, 0.1)) 

p$plot <- p$plot +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))

p$plot <- p$plot +
  theme(plot.margin = margin(2, 0.2, 1, 1, "cm"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16))

print(p) 


#COX
data1$age.73 <- factor(data1$age.73, levels = c('0','1'))
data1$pT <- factor(data1$pT, levels = c('2','3'))
data1$cT <- factor(data1$cT, levels = c('2','3'))
data1$ISUP <- factor(data1$ISUP, levels = c('1','2','3','4'))
data1$prostate.volume <- factor(data1$prostate.volume, levels = c('0','1'))
data1$MTD <- factor(data1$MTD, levels = c('0','1'))
data1$UD <- factor(data1$UD, levels = c('0','1'))
data1$CCL <- factor(data1$CCL, levels = c('0','1'))
data1$Surgical.experience <- factor(data1$Surgical.experience, levels = c('0','1'))
data1$urethral.width <- factor(data1$urethral.width, levels = c('0','1'))
data1$EPE <- factor(data1$EPE, levels = c('0','1'))
data1$PIRADS <- factor(data1$PIRADS, levels = c('3','4','5'))


coxmod<-coxph(Surv(BCR.TIME2,BCR.1year)~age + PSA + PSAD + cT + 
                prostate.volume + PIRADS + MTD + UD + maximal.prostate.width  
              + CCL + EPE + ISUP + posi.core.percent +Surgical.experience
              +urethral.width,
              data=data1)
summary(coxmod)

ft3<-autoReg(coxmod,uni=TRUE,threshold=0.05, final= TRUE)
myft(ft3) 

mycox3<-coxph(Surv(BCR.TIME2,BCR.1year)~Maximum.tumor.diameter+UD+ISUP.grade.of.biopsy+work.list+urethral.width,data=data1)
summary(mycox3)

ft4<-autoReg(mycox3,uni=TRUE,threshold=0.05, final= TRUE)
myft(ft4)
c_index <- summary(mycox3)$concordance
c_index#0.843  (se = 0.028 )


data1$age.73 <- factor(data1$age.73, levels = c('0','1'))
data1$PSA.CAPRA <- factor(data1$PSA.CAPRA, levels = c('1','2','3','4'))
data1$ISUP.CAPRA <- factor(data1$ISUP.CAPRA, levels = c('1','2','3'))
data1$posi.core.percent.CAPRA <- factor(data1$posi.core.percent.CAPRA, levels = c('0','1'))
data1$cT <- factor(data1$cT_summary, levels = c('2','3'))
#CAPRA
mycox2<-coxph(Surv(BCR.TIME2,BCR.1year)~age.73+PSA.CAPRA+ISUP.CAPRA+cT+posi.core.percent.CAPRA,data=data1)
c_index <- summary(mycox2)$concordance
c_index# C 0.68063212 se(C)  0.06706341 
ft3<-autoReg(mycox2,uni=TRUE,threshold=0.05, final= TRUE)
myft(ft3)
#Partin Risk Model
mycox4<-coxph(Surv(BCR.TIME2,BCR.1year)~PSA+ISUP+cT,
              data=data1)
summary(mycox4)
ft4<-autoReg(mycox4,uni=TRUE,threshold=0.05, final= TRUE)
myft(ft4)#0.665  (se = 0.069 )

anova(mycox2,mycox3) 
anova(mycox3,mycox4) 

##ROC
library(timeROC)
data1$lp <- predict(mycox3, newdata = data1,type ="lp")
ROC1 <- timeROC(
  T =data1$BCR.TIME2,
  delta = data1$BCR.1year,
  marker = data1$lp,
  cause = 1,
  weighting="marginal",
  times = 18,
  ROC = TRUE,
  iid = TRUE
)
ROC1$AUC
confint(ROC1,level = 0.95)$CI_AUC 

data1$lp <- predict(mycox2, newdata = data1,type ="lp")
ROC2 <- timeROC(
  T =data1$BCR.TIME2,
  delta = data1$BCR.1year,
  marker = data1$lp,
  cause = 1,
  weighting="marginal",
  times = 18,
  ROC = TRUE,
  iid = TRUE
)
ROC2$AUC
confint(ROC2,level = 0.95)$CI_AUC #0.6478821 46.24-83.34

compare(ROC1, ROC2)
compare(ROC1, ROC3)
pdf('NC ROCs in 18 year OS.pdf',width = 6,height = 6)
plot(ROC1, time = 18, col='#F94343', lwd=3, title = "")
plot(ROC2, time = 18, col='#00A5E0', lwd=3, add = T)

legend("bottomright",
       c(paste0("AUC of the preoperative nomogram:0.89 (0.83-0.95)"), 
         paste0("AUC of the UCSF CAPRA risk model:0.65 (0.462-0.83)")),
       col=c('#F94343', '#00A5E0'),
       lty=1, lwd=3,bty = "n")
dev.off()

library(timeROC)
data1$lp <- predict(mycox3, newdata = data1,type ="lp")
ROC1 <- timeROC(
  T =data1$BCR.TIME2,
  delta = data1$BCR.1year,
  marker = data1$lp,
  cause = 1,
  weighting="marginal",
  times = 18,
  ROC = TRUE,
  iid = TRUE
)
ROC1$AUC
confint(ROC1,level = 0.95)$CI_AUC #0.8851219 82.45-94.57

data1$lp <- predict(mycox2, newdata = data1,type ="lp")
ROC2 <- timeROC(
  T =data1$BCR.TIME2,
  delta = data1$BCR.1year,
  marker = data1$lp,
  cause = 1,
  weighting="marginal",
  times = 18,
  ROC = TRUE,
  iid = TRUE
)
ROC2$AUC
confint(ROC2,level = 0.95)$CI_AUC 

data1$lp <- predict(mycox4, newdata = data1,type ="lp")
ROC3 <- timeROC(
  T =data1$BCR.TIME2,
  delta = data1$BCR.1year,
  marker = data1$lp,
  cause = 1,
  weighting="marginal",
  times = 18,
  ROC = TRUE,
  iid = TRUE
)
ROC3$AUC
confint(ROC3,level = 0.95)$CI_AUC ##AUC0.6273127 43.93-81.53

pdf('NC ROCs in 18 year OS3.pdf',width = 6,height = 6)
plot(ROC1, time = 18, col='#F94343', lwd=3, title = "")
plot(ROC2, time = 18, col='#00A5E0', lwd=3, add = T)
plot(ROC3, time = 18, col='#458B00', lwd=3, add = T)
legend("bottomright",
       c(paste0("AUC of the preoperative nomogram:0.89 (0.83-0.95)"), 
         paste0("AUC of the UCSF CAPRA risk model:0.65 (0.46-0.83)"),
         paste0("AUC of the Partin tables:0.627 (0.44-0.82)")),
       col=c('#F94343', '#00A5E0','#458B00'),
       lty=1, lwd=3,bty = "n")
dev.off()

library(foreign)


#nomogram 
data1$MTD <- factor(data1$MTD, levels = c(0,1), labels = c("≤ 23.9mm", "> 23.9mm"))
data1$UD <- factor(data1$UD, levels = c(0,1), labels = c("≤ 4mm", "> 4mm"))
data1$Surgical.experience <- factor(data1$Surgical.experience, levels = c(0,1), labels = c("≤ 100 cases", "> 100 cases"))
data1$urethral.width <- factor(data1$urethral.width, levels = c(0,1), labels = c("≤ 1.08mm", "> 1.08mm"))
#另一种计算列线图方法
coxm <- cph(Surv(BCR.TIME2,BCR.1year)~MTD+UD+ISUP+Surgical.experience+urethral.width,x=T,y=T,data=data1,surv=T)
nom <- nomogram(coxm,fun=list(function(x)surv(1*18,x)), 
                                      funlabel = c('BCR-free survival at 18 months'),
                        
                lp=F,
                fun.at=c('0.9','0.8','0.70','0.6','0.5','0.4','0.3','0.2'))
par(mar=c(2,5,3,2),cex=1.1)##mar 图形空白边界  cex 文本和符号大小
plot(nom,xfrac=0.5)


legend('topright', legend=col.terms, lwd=2, col=1:length(col.terms), cex=0.8)
terms <- c('MTD'= 'Maximum tumor diameter', 'UD'= 'UD', 'ISUP'= 'ISUP grade of biopsy', 'Surgical.experience'= 'Surgical.experience', 'urethral.width'= 'Urethral width')

plot(nom, xfrac=0.6, lwd=2, col.terms=list('MTD'= 'Maximum tumor diameter', 'UD'= 'UD', 'ISUP'= 'ISUP grade of biopsy', 'Surgical.experience'= 'Surgical.experience', 'urethral.width'= 'Urethral width'))


#calibrate curve
coxm_1 <- cph(Surv(BCR.TIME2,BCR.1year)~MTD+UD+ISUP+Surgical.experience+urethral.width,data=data1,surv=T,x=T,y=T,time.inc = 1*18)
cal_1<-calibrate(coxm_1,u=1*18,cmethod='KM',m=40,B=1000)
par(mar=c(7,4,4,3),cex=1.0)
plot(cal_1,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue = 255))
     xlab='Nomogram Predicted Probability of 18 months BCR-free survival',
     ylab='Actualproportion of 18 months BCR-free survival',
     col=c(rgb(192,98,83,maxColorValue = 255)),
     xlim = c(0.7,1),ylim = c(0.2,1)) 


#DCA
library(dcurves)
library(foreign)
library(survival)
library(dplyr)
str(data1)
f1<-coxph(Surv(BCR.TIME2,BCR.1year)~MTD+UD+ISUP+Surgical.experience+urethral.width,data=data1)
f2<-coxph(Surv(BCR.TIME2,BCR.1year)~age.73+ISUP.CAPRA+cT,data=data1)
f3<-coxph(Surv(BCR.TIME2,BCR.1year)~PSA+cT+ISUP,data=data1)
summary(f1)
summary(f2)
summary(f3)

data1$my = c(1- (summary(survfit(f1, newdata=data1), times=18)$surv))
data1$CAPRA = c(1- (summary(survfit(f2, newdata=data1), times=18)$surv))
data1$ab = c(1- (summary(survfit(f3, newdata=data1), times=18)$surv))

dca(Surv(BCR.TIME2,BCR.1year) ~ my, 
    data = data1,
    time = 18,
    thresholds = 1:50 / 100) %>%
  plot(smooth = T)
dca(Surv(BCR.TIME2,BCR.1year) ~CAPRA, 
    data = data1,
    time = 18,
    thresholds = 1:50 / 100) %>%
  plot(smooth = T)
dca(Surv(BCR.TIME2,BCR.1year) ~Partin, 
    data = data1,
    time = 18,
    thresholds = 1:50 / 100) %>%
  plot(smooth = T)

dca(Surv(BCR.TIME2,BCR.1year) ~ my+CAPRA+Partin, 
    data = data1,
    time = 18,
    thresholds = 1:50 / 100) %>%
  plot(smooth = T)

#DCA2
library("dcurves")
library("ggplot2")
data1$my = c(1- (summary(survfit(ft1, newdata=data1), times=18)$surv))
data1$CAPRA = c(1- (summary(survfit(ft2, newdata=data1), times=18)$surv))
data1$Partin = c(1- (summary(survfit(ft3, newdata=data1), times=18)$surv))

#Three risk categories models
library(openxlsx)
library(stringr)
library(dplyr)

data = read.xlsx("./COX.xlsx")

multi_cox <- function(data,status_colname,time_colname){
  
  library(dplyr)
  library(purrr)
  library(survival)
  
  # data = survdata
  # status_colname = "status"
  # time_colname = "os.time"
  
  
  
  target = colnames(data)[!colnames(data) %in% c(status_colname,time_colname)]
  
  FML <- as.formula(paste0("Surv(",time_colname,",",status_colname,")~",paste(paste0("`",target,"`"),collapse = "+")))
  cox <- coxph(FML, data = data)
  sum_cox <- summary(cox)
  P <- round(sum_cox$coefficients[,5],3) 
  HR <- round(sum_cox$coefficients[,2],3) 
  CI_down <- round(sum_cox$conf.int[,3],3) 
  CI_up <- round(sum_cox$conf.int[,4],3) 
  CI <- paste0(HR," (",CI_down," - ",CI_up,")")
  
  names(HR) = names(HR) %>% str_extract(pattern = paste0(target,collapse = "|"))
  
  # i = 9
  # i = 4
  
  class_num_sum = 0                                             # 
  for (i in 1:length(target)) {
    print(i)
    pos = i
    if (!data[[target[i]]] %>% is.numeric()) {
      
      class_num = data[[target[i]]] %>% unique() %>% na.omit() %>%  length()     
      if (class_num > 2) {
        P = append(P,values = c(NA),after = pos+class_num_sum-1)
        
        HR = append(HR,values = c(NA),after = pos+class_num_sum-1)
        names(HR)[pos+class_num_sum] = target[i]
        names(HR)[c((pos+class_num_sum+1):(pos+class_num_sum+1+class_num-2))] = as.factor(data[[target[i]]]) %>% levels() %>% .[-1]
        
        CI_down = append(CI_down,values = c(NA),after = pos+class_num_sum-1)
        
        CI_up = append(CI_up,values = c(NA),after = pos+class_num_sum-1)
        
        CI = append(CI,values = c(NA),after = pos+class_num_sum-1)
        
        
        class_num_sum = class_num_sum + class_num -1
      }
    }
    
    
  }
  
  
  cox_res <- data.frame(
    "Symbol" =  names(HR),
    "HR" = HR,
    "lower_95" = CI_down,    
    "upper_95" = CI_up ,       
    # "HR_95_CI" = paste0(HR," (",CI_down,"-",CI_up,")"),
    "HR_95_CI" = CI,
    p.value = P,
    
    stringsAsFactors =F
  )
  
  
  cox_res <- cox_res %>% dplyr::mutate(Type = ifelse(.data$p.value < 0.05 & .data$HR > 
                                                       1, "Risky", ifelse(.data$p.value < 0.05 & .data$HR < 
                                                                            1, "Protective", "NS"))) %>% 
    dplyr::mutate(Type = factor(Type,levels = c("NS", "Risky", "Protective")))
  
  return(list(cox_res = cox_res,cox = cox))
}


colnames(data) = str_replace(colnames(data),"\\-","\\.")

multi_cox_res = multi_cox(data = data,status_colname = "BCR.1year",time_colname = "BCR.TIME2")

multi_cox_res$cox_res

cox = multi_cox_res$cox

rs = predict(cox,type ="lp")
data$Risk_score = rs


data = data %>% mutate(Risk_type1 = ifelse(Risk_score > median(Risk_score),"High","Low"),
                       Risk_type2 = case_when(Risk_score < quantile(Risk_score,0.25) ~"Low",
                                              Risk_score < quantile(Risk_score,0.75) ~"Median",
                                              T~"High"
                       )
)


str(data)

#' Title KM 
#'
#' @param survdata 
#' @param var_colname 
#' @param cut_method 
#' @param status_colname 
#' @param time_colname 
#' @param conf.int 
#'
#' @return
#' @export
#'
#' @examples
draw_surv_plot <- function(survdata,var_colname,cut_method,status_colname,time_colname,conf.int = F){
  
  
  library(survival)
  library(survminer)
  
  ##' Extract hazard ratio and confidence intervals from a coxph object
  ##'
  ##' Convenience function to extract hazard ratio and confidence intervals from a coxph object,
  ##' as generated by a call to coxph(), comparing survival times between two groups
  ##'
  ##' @param coxphData coxph object; coxphData <- coxph(Surv(TTE, Cens) ~ group, data=df)
  ##'
  ##' @author Dorothee Nickles
  ##' @export
  ##' @import survival
  ##'
  getHRandCIfromCoxph <- function(coxphData) {
    
    stopifnot(is(coxphData, "coxph"))
    
    tmp <- cbind(
      summary(coxphData)$coef[,
                              c("Pr(>|z|)", "exp(coef)"),
                              drop=FALSE],
      summary(coxphData)$conf[,
                              c("lower .95", "upper .95"),
                              drop=FALSE]
    )
    colnames(tmp) <- c("P","HR","CI_low_0.95","CI_up_0.95")
    tmp[,2:4]<-round(tmp[,2:4],digits=4)
    return(tmp)
  }
  
  
  
  if (is.numeric(survdata[[var_colname]])) {
    
    
    if (cut_method == "best") {
      
      
      tryCatch(expr = {
        cut_point <- as.numeric(surv_cutpoint(survdata,  
                                              time = time_colname,            
                                              event = status_colname,          
                                              variables =gene            
        )$cutpoint[1])
      },
      error = function(e){
        # surv_pvalue <- data.frame(gene = gene, pvalue = paste0("error: ",e))
        surv_pvalue <- data.frame(gene = gene, pvalue = NA)
        return(surv_pvalue)
      })
      
      survdata$Group <- ifelse(survdata[,gene,drop = T] >cut_point,"High","Low")
      survdata$cut_point <- cut_point
    }
    
    
    if (cut_method == "median") {
      survdata$Group <- ifelse(survdata[,gene,drop = T] >median(survdata[,gene,drop = T]),"High","Low")
      survdata$cut_point <- median(survdata[,gene,drop = T])
    }
    
    
  }else{
    
    survdata$Group = survdata[[var_colname]]
  }
  
  

  if (survdata$Group %>% unique() %>% length() == 2) {
    
   
    pvalue<-getHRandCIfromCoxph(coxph(Surv(survdata[,time_colname],survdata[,status_colname])~survdata[,"Group"],data = survdata))
    HR <- paste("Hazard Ratio = ", round(pvalue[,2],2), sep = "")
    CI <- paste("95% CI: ", paste(round(pvalue[,3],2), round(pvalue[,4],2), sep = " - "), sep = "")
  }
  
  
  
  
  curv_str = as.formula(paste0("Surv(",time_colname,",",status_colname,"==1)~","Group"))
  sfit <- surv_fit(curv_str, data = survdata)
  pvalue_surv <- surv_pvalue(sfit)$pval        # Compute p-value from survfit objects  compared using the log-rank test
  print(pvalue_surv)
  
  xlimit = max(survdata[[time_colname]])
  my_plot = ggsurvplot(sfit,
                       # Add p-value and tervals
                       risk.table = TRUE,         
                       tables.height = 0.3,
                       # pval = TRUE,
                       # pval = paste(pval = ifelse(pvalue[,1] < 0.0001, "P < 0.0001",
                       #                            paste("P = ",round(pvalue[,1],4), sep = "")),
                       #              HR, CI,sep = "\n"),
                       
                       
                       # pval = paste(pval = ifelse(pvalue_surv < 0.0001, "P < 0.0001",
                       #                            paste("P = ",round(pvalue_surv,4), sep = "")),
                       #              HR, CI,sep = "\n"),
                       
                       
                      
                       # pval = ifelse(pvalue_surv < 0.0001, "P < 0.0001",
                       #                            paste("P = ",round(pvalue_surv,4), sep = "")),
                       
                       
                      
                       pval = ifelse(survdata$Group %>% unique() %>% length() > 2,
                                     
                                     ifelse(pvalue_surv < 0.0001, "P < 0.0001",
                                            paste("P = ",round(pvalue_surv,4), sep = "")),
                                     
                                     paste(pval = ifelse(pvalue_surv < 0.0001, "P < 0.0001",
                                                         paste("P = ",round(pvalue_surv,4), sep = "")),
                                           HR, CI,sep = "\n")),
                       
                       
                       
                       pval.size = 5,
                       pval.coord = c(0.01*xlimit,0.22),
                       conf.int = conf.int,             
                       size = 0.8, 
                       censor = T,
                       # censor.shape = "|",
                       break.time.by = 12,                  
                       linetype = 1,
                       palette = "nejm",
                       xlim = c(0,xlimit),
                       legend = c(0.77, 0.85),
                       # legend.labs = c(paste0('Low ',mini_sig),paste0("High ",mini_sig)),
                       # legend.labs = c("Low","High"),
                       legend.title= var_colname,
                       font.legend = c(15),
                       risk.table.fontsize = c(3.5),    # table内字体大小
                       font.y = c(12,"bold"),
                       font.x = c(12,"bold"),
                       
                       xlab= "Survival Time (Months)",
                       ggtheme = theme(
                         
                         # aspect.ratio = 1,
                         # line = element_line(size = 0.25),
                         legend.position = "top",
                         legend.background=element_blank(),
                         legend.key = element_blank(),
                         axis.text.x=element_text(size=rel(1),face="bold",vjust =0.4), 
                         axis.text.y=element_text(size=rel(1),face="bold",vjust =0.4), 
                         axis.title.x = element_text(size=rel(1),face="bold"),
                         axis.title.y = element_text(size=rel(1),face="bold"),
                         axis.line.x = element_line(size = 0.5, colour = "black"),
                         axis.line.y = element_line(size = 0.5, colour = "black"),
                         legend.text= element_blank(),
                         legend.title=element_blank(),
                         panel.border = element_blank(),
                         panel.grid = element_blank(),
                         panel.background = element_blank()
                         
                       ), # Change ggplot2 theme
                       
                       tables.theme =  theme(
                         # aspect.ratio = 1,
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         axis.text.y=element_text(size=10,face="bold")
                       )
  ) #%++%
  # scale_x_continuous(breaks = scales::breaks_width(width = 24))      
  
  
  
  my_plot_l <- arrange_ggsurvplots(list(my_plot), print = FALSE, ncol = 1, nrow = 1)
  # ggsave(filename = "./km_curv_add_HR_CI_pvalue.pdf",plot = my_plot_l,width = 6,height = 5)
  return(my_plot_l)
  
  
}


draw_surv_plot(survdata = data,var_colname = "Risk_type2",status_colname = "BCR.1year",time_colname = "BCR.TIME2")

export::graph2pdf(file = "plot1",width = 7,height = 6)


draw_surv_plot(survdata = data,var_colname = "Risk_type1",status_colname = "BCR.1year",time_colname = "BCR.TIME2")

export::graph2pdf(file = "plot2",width = 7,height = 6)

