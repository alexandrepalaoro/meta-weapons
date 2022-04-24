library(metafor)
library(scales)
library(ape)
library(rotl)
library(tidyverse)


rm(list=ls())

full<-read.csv("full_dataset_weapons2.csv",h=T,sep=';')

### BUILDING TREE ####

spp<-tnrs_match_names(full$sp, context_name = "Animals")

#View(spp)

## It found a couple of problems with some beetle species. So, we're gonna 
## substitute those entries with close relationships. 

spp[11,4]<-1029823 #Librodor being substituted by the sub-family
spp[48,4]<-495641 #Chaging Onthophagus taurus for the family 
spp[51,4]<-765934 #And changing acuminatus for the sub-family Dynastinae

my_tree <- tol_induced_subtree(spp$ott_id)
x11()
plot(my_tree,no.margin = T)
dev.off()

## computing the branch lengths of the tree

my_tree.ult<-compute.brlen(my_tree, method = "Grafen")
is.ultrametric(my_tree.ult)
plot(my_tree.ult,no.margin = T)

# Changing the tip labels
my_tree.ult$tip.label<-c("Rangifer tarandus","Capreolus capreolus",
                         "Hipposideros armiger","Ficedula hypoleuca",
                         "Ctenophorus maculosus","Anolis valencienni",
                         "Anolis lineatopus","Anolis sagrei",
                         "Anolis cristatellus","Anolis gundlachi",
                         "Anolis evermanni","Anolis distichus",
                         "Anolis carolinensis","Anolis angusticeps",
                         "Gallotia galloti","Carinascincus microlepidotus",
                         "Gonatodes albogularis","Gnatocerus cornutus",
                         "Sagra femorata","Librodor japonicus",
                         "Heterochelus chiragricus","Dicronocephalus wallichii",
                         "Trypoxylus dichotomus","Onthophagus acuminatus",
                         "Onthophagus taurus", "Aegus chelifer",
                         "Narnia fermorata","Riptortus pedestris",
                         "Anisolabis maritima","Gryllus firmus",
                         "Gryllus pennsylvanicus","Teleogryllus commodus",
                         "Loxoblemmus doenitzi","Acheta domesticus",
                         "Melanotes ornata","Pachyrhamma waitomoensis",
                         "Hemideina crassidens","Kosciuscola tristis",
                         "Scopimera globosa","Austruca annulipes",
                         "Carcinus maenas","Cherax dispar",
                         "Parastacus pilimanus","Parastacus brasiliensis",
                         "Faxonius obscurus","Cambarus robustus",
                         "Cambarus carinirostris","Aegla longirostri",
                         "Diogenes nitidimanus","Neogonodactylus bredini",
                         "Tetranychus urticae","Cambridgea foliata")



#### PLOTTING FIGURE S2 #### 

#First step is figuring out the numbers of the edges and nodes.

plot(my_tree.ult,no.margin = T,cex=1.2,label.offset = 0.05,
     edge.width=2)
nodelabels()
tiplabels()
edgelabels()
dev.off()

# beetle = darkgreen; orthopt = darkorange; crab = blue; mammal = deeppink; bird = steel blue/
# lizard = brown4; arach = azure4

#pdf("my_tree.pdf",h=14,w=12)
tree.colors<-c("darkgreen",rep("darkorange",2),rep("blue",2),"black",rep("deeppink",5),
               "black","steelblue",rep("brown4",25),rep("black",3),rep("darkgreen",16),
               rep("darkslateblue",3),rep("darkorange",14),rep("blue",21),rep("azure4",3))
plot(my_tree.ult,no.margin = T,cex=1.2,label.offset = 0.05,
     edge.width=4,edge.color=tree.colors)
dev.off()


# Making sure the tree is ultrametric

is.ultrametric(my_tree.ult)

# Making a variance-covariance matrix to input in metafor

cov.matrix<-vcv(my_tree.ult,corr=T)

full$phylo<-full$sp

## Getting the values of body size out, keeping just the weapon sizes

weapons<-full[full$trait=='weapon',]

weap <- weapons %>%
  filter(measure %in% c("linear","performance"))

pruned<-drop.tip(my_tree.ult,my_tree.ult$tip.label[-match(weap$phylo,my_tree.ult$tip.label)])
pruned.vcv<-vcv.phylo(pruned)

### Loading the matrix that contains the correlations between effect sizes
## i.e., the allometries between wepaon traits

study.mat<-read.csv("correl_matrix_weap_LIN-PERF.csv",h=T)
head(study.mat)
study.mat<-as.matrix(study.mat)
colnames(study.mat)<-weap$study
rownames(study.mat)<-weap$study
isSymmetric.matrix(study.mat)

#### PERFORMANCE VERSUS LINEAR ####

meta.t<-rma.mv(yi~measure,V=vi,data=weap,method="REML",
               random=list(~1|phylo, ~1|study, ~1|genus/dyad.type),
               R = list(phylo = pruned.vcv, study = study.mat))
summary(meta.t)


# I removed the intercept to make it simpler to plot & calculate the confidence intervals

meta.t2<-rma.mv(yi~measure-1,V=vi,data=weap,method="REML",
                random=list(~1|phylo, ~1|study, ~1|genus/dyad.type),
                R = list(phylo = pruned.vcv, study = study.mat))
summary(meta.t2)
confint(meta.t2)

### FIGURE 4 ####

#tiff("plot-components.tiff",w=160,h=140,units='mm',res=600,compression='lzw')
png("plot-components.png",w=160,h=140,units='mm',res=600)

plot(y=meta.t2$beta[1:2], x = 0:1, bty='l',las=1,cex=2,pch=21,
     bg="grey50",xlim=c(-0.5,1.6),ylim=c(-0.5,1.2),xlab="Measure",
     xaxt='n',ylab="Effect size (Hedges' g)",col=F)
axis(side=1, at=0:1, labels=c("Linear","Performance"),
     las=1)
segments(x0=0:1,y0=meta.t2$ci.lb[1:2],y1=meta.t2$ci.ub[1:2],lwd=5,
         col=c("grey50"))
abline(h=0,lwd=3,col='red',lty=2)
points(y=meta.t2$beta[1:2],x=0:1,pch=21,cex=2.3,bg="grey50")
weap %>%
  count(measure)
####
# Linear : N = 74; Performance = 23
####

text(0,1,expression(italic("74")))
text(1,1,expression(italic("23")))

dev.off()


#### CALCULATING I-SQUARED (i.e., heterogeneity) ####
##I-squared calculations ## Special thanks to Eduardo S. A. Santos!

weap$wi <- 1/sqrt(weap$vi) # precision = 1 / standard error of effect size (Equation 20; Nakagawa & Santos 2012)

s2m.0 <- sum(weap$wi*(length(weap$wi)-1))/(sum(weap$wi)^2-sum(weap$wi^2)) # Equation 22

I2.total <- ((meta.t$sigma2[1]+meta.t$sigma2[2]+meta.t$sigma2[3]+
                meta.t$sigma2[4])/(
                  meta.t$sigma2[1]+meta.t$sigma2[2]+meta.t$sigma2[3]+
                    meta.t$sigma2[4]+s2m.0)) * 100
I2.total

## and 95% CI for I2.total:
I2.total - qchisq(.95, df=1)/2; I2.total + qchisq(.95, df=1)/2

#--- phylo ID ---

I2.phylo <- ((meta.t$sigma2[1])/(
  meta.t$sigma2[1]+meta.t$sigma2[2]+meta.t$sigma2[3]+
    meta.t$sigma2[4]+s2m.0)) * 100

I2.phylo
## and 95% CI for I2.phylo:
I2.phylo - qchisq(.95, df=1)/2; I2.phylo + qchisq(.95, df=1)/2

#--- study ---

I2.std <- ((meta.t$sigma2[2])/(
  meta.t$sigma2[1]+meta.t$sigma2[2]+meta.t$sigma2[3]+
    meta.t$sigma2[4]+s2m.0)) * 100


round(I2.std,digits=5)
## and 95% CI for I2.spp:
I2.std - qchisq(.95, df=1)/2; I2.std + qchisq(.95, df=1)/2

#--- species ID ----

I2.spp <- ((meta.t$sigma2[3])/(
  meta.t$sigma2[1]+meta.t$sigma2[2]+meta.t$sigma2[3]+
    meta.t$sigma2[4]+s2m.0)) * 100

I2.spp

## and 95% CI for I2.phylo:
I2.spp - qchisq(.95, df=1)/2; I2.spp + qchisq(.95, df=1)/2

#--- dyad.type ---

I2.dyad <- ((meta.t$sigma2[4])/(
  meta.t$sigma2[1]+meta.t$sigma2[2]+meta.t$sigma2[3]+
    meta.t$sigma2[4]+s2m.0)) * 100

I2.dyad

## and 95% CI for I2.phylo:
I2.dyad - qchisq(.95, df=1)/2; I2.dyad + qchisq(.95, df=1)/2


##--- H2-squared ----

H2.ini <- ((meta.t$sigma2[1])/(meta.t$sigma2[1]+meta.t$sigma2[2]+meta.t$sigma2[3]+
                                 meta.t$sigma2[4]))
round(H2.ini,digits=5)

H2.ini - qchisq(.95, df=1)/2; H2.ini + qchisq(.95, df=1)/2

#### EGGER'S TEST ####

reg.test.t<-lm(residuals.rma(meta.t)~sqrt(vi),data=weap)
summary(reg.test.t)
confint(reg.test.t)


########
#### DISPLAY ANALYSIS ####
#######

weap %>%
  count(measure,display,disp)

display<- weap %>%
  filter(measure == "linear")

pruned.disp<-drop.tip(my_tree.ult,my_tree.ult$tip.label[-match(display$phylo,my_tree.ult$tip.label)])
disp.vcv<-vcv.phylo(pruned.disp)

#write.csv(display,"display_data.csv")

is.ultrametric(pruned.disp)

## Loading the correlation matrix

disp.mat<-read.csv("correl_matrix_display.csv",h=T)
head(disp.mat)
disp.mat<-as.matrix(disp.mat)
colnames(disp.mat)<-display$study
rownames(disp.mat)<-display$study
isSymmetric.matrix(study.mat)


### MODEL TO TEST DIFFERENCES IN DISPLAYS ####

meta.dis<-rma.mv(yi~display,V=vi,data=display,method="REML",
                 random=list(~1|phylo,~1|study,~1|sp/dyad.type),
                 R = list(phylo = disp.vcv, study = disp.mat))
summary(meta.dis)
anova(meta.dis)

## taking the intercept out just to make it easier to plot & 
## calculate the confidence intervals

meta.dis2<-rma.mv(yi~display-1,V=vi,data=display,method="REML",
                  random=list(~1|phylo,~1|study,~1|sp/dyad.type),
                  R = list(phylo = disp.vcv, study = disp.mat))
summary(meta.dis2)

### FIGURE 5 ####

png("display-plot.png", res=600, units = 'mm', w=140,h=120)
plot(y=meta.dis2$beta[1:3],x=0:2,bty='l',las=1,cex=1.5,pch=21,
     bg="grey50",ylab="Effect size (Hedges' g)",ylim=c(-0.2,1.2),xlim=c(-0.5,2.5),
     xlab="Type of display",xaxt='n')
axis(side=1, at=0:2, labels=c("Body size","Absent","Weapon"))
segments(x0=0:2,y0=meta.dis2$ci.lb[1:3],y1=meta.dis2$ci.ub[1:3],lwd=3,
         col="grey50")
points(y=meta.dis2$beta[1:3],x=0:2,cex=2.3,pch=21,bg='grey50')
abline(h=0,lwd=2,col='red',lty=2)

display %>%
  count(display)

####
# No N = 6; Body = 27; Weapon N = 41
####

text(0,1.1,expression(italic("27")))
text(1,1.1,expression(italic("6")))
text(2,1.1,expression(italic("41")))

dev.off()


#### Display - I-squared calculations  ####

weap$wi <- 1/sqrt(weap$vi) # precision = 1 / standard error of effect size (Equation 20; Nakagawa & Santos 2012)

s2m.0 <- sum(weap$wi*(length(weap$wi)-1))/(sum(weap$wi)^2-sum(weap$wi^2)) # Equation 22

I2.total <- ((meta.dis$sigma2[1]+meta.dis$sigma2[2]+meta.dis$sigma2[3]+
                meta.dis$sigma2[4])/(
                  meta.dis$sigma2[1]+meta.dis$sigma2[2]+meta.dis$sigma2[3]+
                    meta.dis$sigma2[4]+s2m.0)) * 100
I2.total

## and 95% CI for I2.total:
I2.total - qchisq(.95, df=1)/2; I2.total + qchisq(.95, df=1)/2

#--- phylo ID ---

I2.phylo <- ((meta.dis$sigma2[1])/(
  meta.dis$sigma2[1]+meta.dis$sigma2[2]+meta.dis$sigma2[3]+
    meta.dis$sigma2[4]+s2m.0)) * 100

I2.phylo
## and 95% CI for I2.phylo:
I2.phylo - qchisq(.95, df=1)/2; I2.phylo + qchisq(.95, df=1)/2

#--- study ---

I2.std <- ((meta.dis$sigma2[2])/(
  meta.dis$sigma2[1]+meta.dis$sigma2[2]+meta.dis$sigma2[3]+
    meta.dis$sigma2[4]+s2m.0)) * 100


round(I2.std,digits=5)
## and 95% CI for I2.spp:
I2.std - qchisq(.95, df=1)/2; I2.std + qchisq(.95, df=1)/2

#--- species ID ----

I2.spp <- ((meta.dis$sigma2[3])/(
  meta.dis$sigma2[1]+meta.dis$sigma2[2]+meta.dis$sigma2[3]+
    meta.dis$sigma2[4]+s2m.0)) * 100

round(I2.spp,digits=7)

## and 95% CI for I2.phylo:
I2.spp - qchisq(.95, df=1)/2; I2.spp + qchisq(.95, df=1)/2

#--- dyad.type ---

I2.dyad <- ((meta.dis$sigma2[4])/(
  meta.dis$sigma2[1]+meta.dis$sigma2[2]+meta.dis$sigma2[3]+
    meta.dis$sigma2[4]+s2m.0)) * 100

I2.dyad

## and 95% CI for I2.phylo:
I2.dyad - qchisq(.95, df=1)/2; I2.dyad + qchisq(.95, df=1)/2


##--- H2-squared ----

H2.ini <- ((meta.dis$sigma2[1])/(meta.dis$sigma2[1]+meta.dis$sigma2[2]+meta.dis$sigma2[3]+
                                   meta.dis$sigma2[4]))
round(H2.ini,digits=5)

H2.ini - qchisq(.95, df=1)/2; H2.ini + qchisq(.95, df=1)/2


confint(meta.dis)

#### EGGER'S TEST ####


reg.test.t<-lm(residuals.rma(meta.dis)~sqrt(vi),data=display)
summary(reg.test.t)
confint(reg.test.t)



####----------------------###
##### FIGHTING STYLE ANALYSES #####
####----------------------###

weap %>%
  count(measure,fight.cat)

fight.ana <- weap %>%
  filter(measure %in% "linear")

fight.ana %>%
  count(measure,fight.cat)
fight.ana %>% count(sp)

fight.ana

## Loading the correlation matrix

lin.mat<-read.csv("correl_matrix_weap_LIN.csv",h=T,sep=';')
head(lin.mat)
lin.mat<-as.matrix(lin.mat)
length(lin.mat[,1])
length(fight.ana$art)
colnames(lin.mat)<-fight.ana$study
rownames(lin.mat)<-fight.ana$study
isSymmetric.matrix(lin.mat)

pruned.funct<-drop.tip(my_tree.ult,my_tree.ult$tip.label[-match(fight.ana$phylo,my_tree.ult$tip.label)])
funct.vcv<-vcv.phylo(pruned.funct)

colnames(funct.vcv) %in% fight.ana$phylo 

pruned.funct$tip.label


### TESTING THE DIFFERENCE BETWEEN FUNCTIONS ####

meta.fight <- rma.mv(yi~fight.cat,V=vi,data=funct.ana,method="REML",
                    random=list(~1|phylo, ~1|study, ~1|genus/dyad.type),
                    R = list(phylo = funct.vcv, study = lin.mat))

summary(meta.fight)

## Testing the differences between the levels within function manually

anova(meta.fight,L=c(1,0,-1))
anova(meta.fight,L=c(0,1,-1))


## taking the intercept out just to make it easier to plot & 
## calculate the confidence intervals

meta.fight2<-rma.mv(yi~fight.cat-1,V=vi,data=funct.ana,method="REML",
                   random=list(~1|phylo, ~1|study, ~1|genus/dyad.type),
                   R = list(phylo = funct.vcv, study = lin.mat))
summary(meta.fight2)

### FIGURE 6 ####

png("function-plot.png",res=600,units='mm',w=180,h=140)

plot(y=meta.fight2$beta[1:3],x=0:2,bty='l',las=1,cex=1.5,pch=21,
     bg=c("grey50","grey50","grey50","black","black"),
     ylab="Effect size (Hedges' g)",ylim=c(-0.2,1.2),xlim=c(-0.5,2.5),
     xlab="Fighting style",xaxt='n')

axis(side=1, at=c(0:2), labels=c("Size-emphasis","Intermediate","Performance-emphasis"))

segments(x0=c(0,1,2),y0=meta.fight2$ci.lb[1:3],y1=meta.fight2$ci.ub[1:3],lwd=3,
         col="grey50")
points(y=meta.fight2$beta[1:3],x=c(0:2),cex=2.3,pch=21,bg="grey50")
abline(h=0,lwd=2,col='red',lty=2)

fight.ana %>%
  count(fight.cat)

text(0,1.2,expression(italic("10")))
text(1,1.2,expression(italic("47")))
text(2,1.2,expression(italic("17")))


dev.off()


#### Function - I-squared calculations ####

weap$wi <- 1/sqrt(weap$vi) # precision = 1 / standard error of effect size (Equation 20; Nakagawa & Santos 2012)

s2m.0 <- sum(weap$wi*(length(weap$wi)-1))/(sum(weap$wi)^2-sum(weap$wi^2)) # Equation 22

I2.total <- ((meta.fight$sigma2[1]+meta.fight$sigma2[2]+meta.fight$sigma2[3]+
                meta.fight$sigma2[4])/(
                  meta.fight$sigma2[1]+meta.fight$sigma2[2]+meta.fight$sigma2[3]+
                    meta.fight$sigma2[4]+s2m.0)) * 100
I2.total

## and 95% CI for I2.total:
I2.total - qchisq(.95, df=1)/2; I2.total + qchisq(.95, df=1)/2

#--- phylo ID ---

I2.phylo <- ((meta.fight$sigma2[1])/(
  meta.fight$sigma2[1]+meta.fight$sigma2[2]+meta.fight$sigma2[3]+
    meta.fight$sigma2[4]+s2m.0)) * 100

I2.phylo
## and 95% CI for I2.phylo:
I2.phylo - qchisq(.95, df=1)/2; I2.phylo + qchisq(.95, df=1)/2

#--- study ---

I2.std <- ((meta.fight$sigma2[2])/(
  meta.fight$sigma2[1]+meta.fight$sigma2[2]+meta.fight$sigma2[3]+
    meta.fight$sigma2[4]+s2m.0)) * 100


round(I2.std,digits=5)
## and 95% CI for I2.spp:
I2.std - qchisq(.95, df=1)/2; I2.std + qchisq(.95, df=1)/2

#--- species ID ----

I2.spp <- ((meta.fight$sigma2[3])/(
  meta.fight$sigma2[1]+meta.fight$sigma2[2]+meta.fight$sigma2[3]+
    meta.fight$sigma2[4]+s2m.0)) * 100

round(I2.spp,digits=7)

## and 95% CI for I2.phylo:
I2.spp - qchisq(.95, df=1)/2; I2.spp + qchisq(.95, df=1)/2

#--- dyad.type ---

I2.dyad <- ((meta.fight$sigma2[4])/(
  meta.fight$sigma2[1]+meta.fight$sigma2[2]+meta.fight$sigma2[3]+
    meta.fight$sigma2[4]+s2m.0)) * 100

I2.dyad

## and 95% CI for I2.phylo:
I2.dyad - qchisq(.95, df=1)/2; I2.dyad + qchisq(.95, df=1)/2


##--- H2-squared ----

H2.ini <- ((meta.fight$sigma2[1])/(meta.fight$sigma2[1]+meta.fight$sigma2[2]+meta.fight$sigma2[3]+
                                     meta.fight$sigma2[4]))
round(H2.ini,digits=5)

H2.ini - qchisq(.95, df=1)/2; H2.ini + qchisq(.95, df=1)/2


confint(meta.fight)

#### EGGER'S TEST ####


reg.test.t<-lm(residuals.rma(meta.fight)~sqrt(vi),data=weap)
summary(reg.test.t)
confint(reg.test.t)


#### TRAIT PERFORMANCE VS LINEAR COMPARISON WITHIN DECAPODS AND LIZARDS ####

weap2 <- weap %>% filter(order %in% c("Decapoda","Lizard"))

pruned.post<-drop.tip(my_tree.ult,my_tree.ult$tip.label[-match(weap2$phylo,my_tree.ult$tip.label)])
pruned.postvcv<-vcv.phylo(pruned.post)

meta.post<-rma.mv(yi~measure,V=vi,data=weap2,method="REML",
                  random=list(~1|phylo, ~1|study, ~1|genus/dyad.type),
                  R = list(phylo = pruned.postvcv))

summary(meta.post)

#### FIGURE S3 ####

png("plot-performanceONLY.png",w=160,h=140,units='mm',res=600)

plot(y=meta.post$beta[1:2], x = 0:1, bty='l',las=1,cex=2,pch=21,
     bg="grey50",xlim=c(-0.5,1.6),ylim=c(-0.5,1.2),xlab="Measure",
     xaxt='n',ylab="Effect size (Hedges' g)",col=F)
axis(side=1, at=0:1, labels=c("Linear","Performance"),
     las=1)
segments(x0=0:1,y0=meta.post$ci.lb[1:2],y1=meta.post$ci.ub[1:2],lwd=5,
         col=c("grey50"))
abline(h=0,lwd=3,col='red',lty=2)
points(y=meta.post$beta[1:2],x=0:1,pch=21,cex=2.3,bg="grey50")
weap2 %>%
  count(measure)
####
# Linear : N = 31; Performance = 21
####

text(0,1,expression(italic("31")))
text(1,1,expression(italic("21")))

dev.off()
