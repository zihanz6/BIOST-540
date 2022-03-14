library(ggplot2)
library(tidyverse)
library(nlme)
library(joineR)
library(MASS)
library(lattice)
library(dplyr)

cd4<-cd4[,-1]
cd4$age_category<-cut(cd4$age, c(14,35,55,75))
cd4$period<-cut(cd4$week, c(0,8,16,24,32,40), include.lowest=TRUE)

#Compare baseline characteristics between treatment groups
with(cd4, table(group, sex))

group1<-subset(cd4, group==1)
mean(group1$age)
sd(group1$age)
median(group1$age)

group2<-subset(cd4, group==2)
mean(group2$age)
sd(group2$age)
median(group2$age)

group3<-subset(cd4, group==3)
mean(group3$age)
sd(group3$age)
median(group3$age)

group4<-subset(cd4, group==4)
mean(group4$age)
sd(group4$age)
median(group4$age)

#Compare by treatment groups
#Plot
cd4$group<-factor(cd4$group)
p<-ggplot(data=cd4, aes(x=week, y=logcd4, group=group, col=group))
p+geom_line(alpha=0.2)+facet_grid(.~group)+geom_line(data=cd4%>%
                                                       group_by(group, week)%>%summarise(logcd4=mean(logcd4)),
                                                     aes(x=week, y=logcd4, group=group))+theme_bw()+ylab("log(CD4 count + 1)")

#Descriptive Statistics
stat<-cd4%>%group_by(period, group)%>%summarise(number_of_subjects=sum(!is.na(logcd4))%>%round(2),
                                                logcd4_mean=mean(logcd4, na.rm=TRUE)%>%round(2), logcd4_sd=sd(logcd4, na.rm=TRUE)%>%round(2))
knitr::kable(stat, col.names = c("Week","Treatment Group","Number of Subjects",
                                 "Mean of log-CD4","SD of log-CD4"), align="rrrrr")
#Compare by baseline age
#Plot
p<-ggplot(data=cd4, aes(x=week, y=logcd4, group=age_category, col=age_category))
p+geom_line(alpha=0.2)+facet_grid(.~age_category)+geom_line(data=cd4%>%
                                                              group_by(age_category, week)%>%summarise(logcd4=mean(logcd4)),
                                                            aes(x=week, y=logcd4, group=age_category))+theme_bw()+ylab("log(CD4 count +1)")

#Descriptive Statistics
stat1<-cd4%>%group_by(period, age_category)%>%summarise(number_of_subjects=sum(!is.na(logcd4)),
                                                        logcd4_mean=mean(logcd4, na.rm=TRUE)%>%round(2), logcd4_sd=sd(logcd4, na.rm=TRUE)%>%round(2))
knitr::kable(stat1, col.names = c("Week","Baseline Age","Number of Subjects",
                                  "Mean of log-CD4","SD of log-CD4"), align="rrrrr")

#Compare by sex
#Plot
cd4$sex<-factor(cd4$sex)
p<-ggplot(data=cd4, aes(x=week, y=logcd4, group=sex, col=sex))
p+geom_line(alpha=0.2)+facet_grid(.~sex)+geom_line(data=cd4%>%
                                                     group_by(sex, week)%>%summarise(logcd4=mean(logcd4)),
                                                   aes(x=week, y=logcd4, group=sex))+theme_bw()+ylab("log(CD4 count + 1)")

#Descriptive Statistics
stat2<-cd4%>%group_by(period, sex)%>%summarise(number_of_subjects=sum(!is.na(logcd4)),
                                               logcd4_mean=mean(logcd4, na.rm=TRUE)%>%round(2), logcd4_sd=sd(logcd4, na.rm=TRUE)%>%round(2))
knitr::kable(stat2, col.names = c("Week","Sex","Number of Subjects",
                                  "Mean of log-CD4","SD of log-CD4"), align="rrrrr")

#Missing Data
visit_counts<-cd4%>%group_by(id)%>%summarise(n=n(), nobs=sum(!is.na(logcd4)))
cd4$visit_times<-visit_counts$nobs[match(cd4$id, visit_counts$id)]
cd4_baseline<-cd4[,c("id", "sex", "age", "visit_times")]
cd4_baseline<-cd4_baseline[!duplicated(cd4_baseline),]
complete_case<-cd4_baseline[cd4_baseline$visit_times>=6,]
missing_case<-cd4_baseline[cd4_baseline$visit_times<6,]
table(complete_case$sex)
table(missing_case$sex)
summary(complete_case$age)
summary(missing_case$age)


Part 2
cd4.long <- read.csv("~/Downloads/cd4.csv")
cd4.long$sex <- relevel(factor(cd4.long$sex), ref=1) # making males the reference 
cd4.long$group <- relevel(factor(cd4.long$group), ref=1) # make 1 = zidovudine alternating monthly with 400 mg didanosine the reference group
cd4.long$time <- as.numeric(cut(cd4.long$week, c(0,8,16,24,32,40), include.lowest=TRUE)) 
head(cd4.long)

### unadjusted models
# add knot at week 16
cd4.long$weekSpline16 <- (cd4.long$week - 16)*(cd4.long$week > 16)
##### LMM Random intercepts + slopes uncorrelated with linear spline 
mod1 <- lme( logcd4 ~ group*week + weekSpline16, 
             method = "ML", data = cd4.long, 
             random = reStruct( ~ 1 + week | id, pdClass="pdDiag", REML=F))

summary(mod1)

##### LMM Random intercepts + slopes correlated with linear spline 
mod2 <- lme( logcd4 ~ group*week + weekSpline16, 
             method = "ML", data = cd4.long, 
             random = reStruct( ~ 1 + week | id, pdClass="pdSymm", REML=F))
mod2_result<-summary(mod2)
mod2_result

####### LMM Random intercepts with linear spline 
mod3 <- lme( logcd4 ~ group*week + weekSpline16, 
             method = "ML", data = cd4.long, 
             random = reStruct( ~ 1 | id, REML=F))
summary(mod3)
anova(mod1,mod2,mod3) # pick model2


interval1 <- intervals(mod2,which = "fixed")

#table of result
table1 <- data.frame(coefficients = round(mod2_result$coefficients$fixed,3),
                     lower = round(interval1$fixed[,1],3),
                     upper = round(interval1$fixed[,3],3),
                     pvalue = round(mod2_result$tTable[,5],3))

knitr::kable(table1, col.names = gsub("[.]"," ", names(table1)), align = "lccrr",caption = "Table4 Results from model2")


### Hypothesis Testing
reduced.mod2 <- lme( logcd4 ~ group + week + weekSpline16, 
                     method = "ML", data = cd4.long, 
                     random = reStruct( ~ 1 + week | id, pdClass="pdSymm", REML=F))
anova(reduced.mod2,mod2)
anova(mod2)


### adjusted models
# add knot at week 16
cd4.long$weekSpline16 <- (cd4.long$week - 16)*(cd4.long$week > 16)
##### LMM Random intercepts + slopes uncorrelated with linear spline ######
mod4 <- lme( logcd4 ~ group*week + weekSpline16, 
             method = "ML", data = cd4.long, 
             random = reStruct( ~ 1 + week | id, pdClass="pdDiag", REML=F))

summary(mod4)

##### LMM Random intercepts + slopes correlated with linear spline ######
mod5 <- lme( logcd4 ~ age + sex + group*week + weekSpline16, 
             method = "ML", data = cd4.long, 
             random = reStruct( ~ 1 + week | id, pdClass="pdSymm", REML=F))
mod5_result <- summary(mod5)
mod5_result

####### LMM Random intercepts with linear spline 
mod6 <- lme( logcd4 ~ age + sex + group*week + weekSpline16, 
             method = "ML", data = cd4.long, 
             random = reStruct( ~ 1 | id, REML=F))
summary(mod6)
anova(mod4,mod5,mod6) # pick model5
interval2 <- intervals(mod5,which = "fixed")
#table of results
table2 <- data.frame(coefficients = round(mod5_result$coefficients$fixed,3),
                     lower = round(interval2$fixed[,1],3),
                     upper = round(interval2$fixed[,3],3),
                     pvalue = round(mod5_result$tTable[,5],3))

knitr::kable(table2, col.names = gsub("[.]", " ", names(table2)), align = "lccrr",caption = "Table5 Results from model5")

### Hypothesis Testing
reduced.mod5 <- lme(logcd4 ~ age + sex + group+week + weekSpline16, 
                    method = "ML", data = cd4.long, 
                    random = reStruct( ~ 1 + week | id, pdClass="pdSymm", REML=F))
anova(reduced.mod5,mod5)
anova(mod5)

#mod2
cluster.res <- resid(mod2, level=1, type="normalized")
pop.res <- resid(mod2, level=0, type="normalized")
bi <- mod2$coefficients$random$id
plot(cd4.long$week, cluster.res, xlab="week")
par(mfrow=c(2,2))
qqnorm(bi[,1], main="Random Intercepts")
qqline(bi[,1])
qqnorm(bi[,2], main="Random Slopes")
qqline(bi[,2])
qqnorm(pop.res, main="Population Residuals")
qqline(pop.res)
qqnorm(cluster.res, main="Cluster Residuals")
qqline(cluster.res)

#mod5
cluster.res <- resid(mod5, level=1, type="normalized")
pop.res <- resid(mod5, level=0, type="normalized")
bi <- mod5$coefficients$random$id
plot(cd4.long$week, cluster.res, xlab="week")
par(mfrow=c(2,2))
qqnorm(bi[,1], main="Random Intercepts")
qqline(bi[,1])
qqnorm(bi[,2], main="Random Slopes")
qqline(bi[,2])
qqnorm(pop.res, main="Population Residuals")
qqline(pop.res)
qqnorm(cluster.res, main="Cluster Residuals")
qqline(cluster.res)
