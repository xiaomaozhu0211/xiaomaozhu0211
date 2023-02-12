library(plyr)
library(jskm)
library(survival)
library(survminer)
library(tableone)
library(condsurv)
library(ggprism)
library(ggsci)
library(ggalt)
library(riskRegression)
library(ggDCA)
library(forestplot)
library(stringr)
library(patchwork)
library(rms)
library(glmnet)

rm(list = ls())
data0<-read.csv("tnbc.csv")
data<-filter(data0, data0$time !='0')#n=32836

aa <- data.frame(data)
for (i in names(aa)[c(5:13)]){aa[,i] <- as.factor(aa[,i])}

#7:3
set.seed(1)
sub<-sample(1:nrow(aa),round(nrow(aa)*7/10))
length(sub) 
train<-aa[sub,];test<-aa[-sub,]
group <- c(rep(0,times=dim(train)[1]))
traindata <- cbind(group,train)
group<- c(rep(1,times=dim(test)[1]))
testdata <- cbind(group,test)
aa <- rbind(traindata,testdata)
names(traindata)

dput(names(aa))
myVars <- c("age", "age1", "marry", "race", "t", "n","tnm", "g", 
            "surg", "che", "rt", "bcss", "status", "time")
catVars <- c("age1", "marry", "race", "t", "n","tnm", "g", 
             "surg", "che", "rt", "bcss", "status")
table <- CreateTableOne(vars = myVars, 
                        factorVars = catVars,
                        strata = "group", 
                        data = aa,
                        addOverall = TRUE)
table1<- print(table, 
               catDigits =1,contDigits =1,pDigits = 3, 
               showAllLevels=TRUE, 
               quote = FALSE, 
               printToggle = TRUE)



#--------------------------------------------------------------------------KM
aa <- data.frame(data)
f <- survfit(Surv(time/12, status==1) ~ 1, data =aa)
cond_times <- seq(0,9.92)
p1 <- gg_conditional_surv(basekm = f, 
                          at = cond_times, 
                          main = "",
                          xlab = "Time Since Diagnosis (years)",
                          ylab = "Overall Survival (%)", 
                          lwd =1)+
  theme_prism(base_size =15)+
  theme(legend.position=c(0.1,0.5)) +
  guides(fill=guide_legend(nrow=1))+
  theme(panel.grid=element_blank(),
        panel.grid.major.y = element_line(colour = "grey80"))+
  labs(color = "Given")+
  scale_colour_d3()+
  scale_y_continuous( breaks = seq(0.65,1,by=0.05),expand = c(0,0),
                      limits = c(0.65,1))+
  scale_x_continuous( breaks = seq(0,10,by=1),expand = c(0,0),
                      limits = c(0,10.1));p1

#----------------------------------------------------------------------------AHR
library(ggalt)
aa$time <- aa$time/12
aa$time1 <- ceiling(aa$time)
km1<-survfit(Surv(time1,status==1)~1,data=aa)
res1<-summary(km1)
df <- cbind(res1$time,res1$n.event/res1$n.risk)
df <- data.frame(df)
p3 <-ggplot(df)+
  geom_point(aes(x=X1,y=X2),size=3)+
  geom_xspline(aes(x=X1,y=X2),spline_shape = 0.4,size=1,color="red")+
  theme_classic()+
  theme_prism(base_size =15)+
  scale_x_continuous(limits=c(1,10),breaks=seq(0,10,1),expand = c(0.02,0))+
  scale_y_continuous(limits=c(0,0.1),breaks=seq(0,0.1,0.02),expand = c(0.02,0))+
  labs(x="Time Since Diagnosis (years)", y="Annual Hazard Rate (%)")

#--------------------------------------------------------------------------LASSO
aa <- data.frame(data)
for (i in names(aa)[c(5:13)]){aa[,i] <- as.factor(aa[,i])}

x <- data.matrix(traindata[,c(6:14)])
y <- data.matrix(Surv(traindata$time,traindata$status))
fit<- glmnet(x, y,family="cox",alpha=1)
plot(fit, xvar="lambda", label=T)# 
set.seed(100)
cv.fit <- cv.glmnet(x, y,family="cox",nfolds=10,
                    type.measure="mse", 
                    alpha = 1);plot(cv.fit)
log(cv.fit$lambda.1se)#-4.094898
ridge.coef1 <- predict(fit,s=cv.fit$lambda.1se,type = "coefficients");ridge.coef1
cv.fit$lambda.1se
lasso<-coxph(Surv(time,status==1)~age1+marry+race+t+n+che+surg+rt,data=traindata,x=T);summary(lasso)

#------------------------------------------------------------------------nomogram
aa <- traindata  
nomo<-datadist(aa)
options(datadist='nomo')
#
nomo1 <- cph(Surv(time,status==1)~age1+marry+race+t+n+che+surg+rt,x=T,y=T,
             data=aa,surv=T,time.inc=12*5);nomo1
#
surv <- Survival(nomo1)
surv1 <- function(x)surv(60,lp=x)
surv2 <- function(x)surv(119,lp=x)

#3
nomo2<-nomogram(nomo1,
                fun=list(surv1,
                         surv2),
                funlabel=c('5-year OS probability',
                           '10-year OS probability'),
                lp =F, 
                maxscale=100,
                fun.at=c("0.99","0.95",'0.9','0.8',
                         '0.7','0.6','0.5','0.4',
                         '0.3','0.2','0.1',"0.05","0.01"));plot(nomo2)
#-------------
aa <- traindata  
nomo<-datadist(aa)
options(datadist='nomo')
nomo2 <- cph(Surv(time,status==1)~age1+marry+race+t+n+che+surg+rt,
             x=T,y=T,data=aa,surv=T,time.inc=60)
output2<- calibrate(nomo2, cmethod='KM',method='boot', u=60,m=1900, B=1000)
plot(output2,add=F,conf.int=T,cex.subtitles=0.5, riskdist=TRUE,lwd=2,lty=1,
     errbar.col=c("#1F77B4FF"),xlim=c(0,1), ylim=c(0,1),
     xlab="Predicted OS",
     ylab="Actual OS ",
     col=c("#1F77B4FF"))

#------------------------------------------------------------------------ ROC --
aa <- traindata
f1 <- coxph(Surv(time/12,status==1)~age1+marry+race+t+n+che+surg+rt,
            data=aa,x=TRUE,y=TRUE)
fig1<-Score(list('model1'=f1),
            formula=Hist(time/12, status)~1,
            data = aa,
            se.fit=1L,
            times=c(1:9),
            plots="ROC",
            metrics ="auc")
cs <- read.csv("auc.csv")
pauc <-ggplot(cs)+
  geom_point(aes(x=time,y=auc,group = name,color=name),size=3)+
  geom_xspline(aes(x=time,y=auc,group = name,color=name),size=1)+
  geom_errorbar(aes(x=time,
                    ymin =ci1,
                    ymax =ci2,
                    color=name),width =0.2,size=1)+
  scale_colour_d3()+
  theme_prism(base_size =15)+
  theme(legend.position=c(0.1,0.2))+
  scale_y_continuous( breaks = seq(70,90,by=2),expand = c(0.01,0),
                      limits = c(70,90.01))+
  scale_x_continuous(  breaks = seq(1,10,by=1),expand = c(0.01,0),
                       limits = c(1,10.1))+
  labs(x="Time Since Diagnosis (years)", y="Time-dependent AUC")

#------------------------------------------------------------------------DCA----
f3 <- coxph(Surv(time,status==1)~age1+marry+race+t+n+che+surg+rt,
            data=traindata,x=TRUE,y=TRUE)
f4 <- coxph(Surv(time,status==1)~age1+marry+race+t+n+che+surg+rt,
            data=testdata,x=TRUE,y=TRUE)
fig3<- dca(f3,times=c(60,119))
pdca1 <- ggplot(fig3,
                color = c("#1F77B4FF","#FF7F0EFF","grey","grey","#374E54"),
                linetype =F,lwd = 1)+
  theme_classic()+   
  theme_prism(base_size =10)+
  theme(legend.position="top");pdca1
fig4<- dca(f4,times=c(60,119))
pdca2 <-ggplot(fig4,
               color = c("#1F77B4FF","#FF7F0EFF","grey","grey","#374E54"),
               linetype =F,lwd = 1)+
  theme_classic()+   
  theme_prism(base_size =10)+
  theme(legend.position="top")


#cs-nomo
nomo1 <- cph(Surv(time,status==1)~age1+marry+race+t+n+che+surg+rt,x=T,y=T,
             data=aa,surv=T,time.inc=12*5);nomo1
#
surv <- Survival(nomo1)
surv1 <- function(x)surv(12*1,lp=x)
surv2 <- function(x)surv(12*2,lp=x)
surv3 <- function(x)surv(12*3,lp=x)
surv4 <- function(x)surv(12*4,lp=x)
surv5 <- function(x)surv(12*5,lp=x)
surv6 <- function(x)surv(12*6,lp=x)
surv7 <- function(x)surv(12*7,lp=x)
surv8 <- function(x)surv(12*8,lp=x)
surv9 <- function(x)surv(12*9,lp=x)
surv10 <- function(x)surv(119,lp=x)

#3-nomogram模型建立
nomo2<-nomogram(nomo1,
                fun=list(surv1,surv2,surv3,surv4,
                         surv5,surv6,surv7,surv8,
                         surv9,surv10),
                funlabel=c('1-year OS probability',
                           '2-year OS probability',
                           '3-year OS probability',
                           '4-year OS probability',
                           '5-year OS probability',
                           '6-year OS probability',
                           '7-year OS probability',
                           '8-year OS probability',
                           '9-year OS probability',
                           '10-year OS probability'),
                lp =F, 
                maxscale=100,
                fun.at=c('0.995','0.99','0.985','0.98','0.975','0.97','0.965','0.96','0.955','0.95','0.945','0.94','0.935','0.93','0.925','0.92','0.915',"0.91",'0.905','0.9',
                         '0.895','0.89','0.885','0.88','0.875','0.87','0.865','0.86','0.855','0.85','0.845','0.84','0.835','0.83','0.825','0.82','0.815',"0.18",'0.805','0.8',
                         '0.795','0.79','0.785','0.78','0.775','0.77','0.765','0.76','0.755','0.75','0.745','0.74','0.735','0.73','0.725','0.72','0.715',"0.71",'0.705','0.7',
                         '0.695','0.69','0.685','0.68','0.675','0.67','0.665','0.66','0.655','0.65','0.645','0.64','0.635','0.63','0.625','0.62','0.615',"0.61",'0.605','0.6',
                         '0.595','0.59','0.585','0.58','0.575','0.57','0.565','0.56','0.555','0.55','0.545','0.54','0.535','0.53','0.525','0.52','0.515',"0.51",'0.505','0.5',
                         '0.495','0.49','0.485','0.48','0.475','0.47','0.465','0.46','0.455','0.45','0.445','0.44','0.435','0.43','0.425','0.42','0.415',"0.41",'0.405','0.4',
                         '0.395','0.39','0.385','0.38','0.375','0.37','0.365','0.36','0.355','0.35','0.345','0.34','0.335','0.33','0.325','0.32','0.315',"0.31",'0.305','0.3',
                         '0.295','0.29','0.285','0.28','0.275','0.27','0.265','0.26','0.255','0.25','0.245','0.24','0.235','0.23','0.225','0.22','0.215',"0.21",'0.205','0.2',
                         '0.195','0.19','0.185','0.18','0.175','0.17','0.165','0.16','0.155','0.15','0.145','0.14','0.135','0.13','0.125','0.12','0.115',"0.11",'0.105','0.1',
                         '0.095','0.09','0.085','0.08','0.075','0.07','0.065','0.06','0.055','0.05','0.045','0.04','0.035','0.03','0.025','0.02','0.015',"0.01",'0.005','0.0'));nomo2

cs1<- data.frame(nomo2$`1-year OS probability`$x,nomo2$`1-year OS probability`$x.real)
cs2<- data.frame(nomo2$`2-year OS probability`$x,nomo2$`2-year OS probability`$x.real)
cs3<- data.frame(nomo2$`3-year OS probability`$x,nomo2$`3-year OS probability`$x.real)
cs4<- data.frame(nomo2$`4-year OS probability`$x,nomo2$`4-year OS probability`$x.real)
cs5<- data.frame(nomo2$`5-year OS probability`$x,nomo2$`5-year OS probability`$x.real)
cs6<- data.frame(nomo2$`6-year OS probability`$x,nomo2$`6-year OS probability`$x.real)
cs7<- data.frame(nomo2$`7-year OS probability`$x,nomo2$`7-year OS probability`$x.real)
cs8<- data.frame(nomo2$`8-year OS probability`$x,nomo2$`8-year OS probability`$x.real)
cs9<- data.frame(nomo2$`9-year OS probability`$x,nomo2$`9-year OS probability`$x.real)
cs10<-data.frame(nomo2$`10-year OS probability`$x,nomo2$`10-year OS probability`$x.real)