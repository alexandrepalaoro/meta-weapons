###Simulation model to evaluate whether the importance of a
###trait in determining victory affects the expression 
###its mean values in winners and losers

library(faux)
library(sciplot)
library(boot)

#Div1 and div2 defines the maximum probability of inversion
#in which the weaker individual was able to win the contest.
#Since div1 and div2 are equal to 1/div, the values 2 and 4
#determine maximum probabilities of 0.5 and 0.25 respectively.

#Creating data frames and vectors to insert the values of the
#simulations
data=data.frame(list())
final.mean1=c()
final.mean2=c()
final.mean3=c()
final.sd1w=c()
final.sd2w=c()
final.sd3w=c()
final.sd1l=c()
final.sd2l=c()
final.sd3l=c()

#Determining the maximum probability of inversion (when the 
#weaker rival wins the contest)
div1=4
div2=2

#Performing the simulation 1000 times 
for(i in 1:1000){

#Creating two groups of 100 individuals with equal mean trait
#values. Individuals in the same row will fight in the simulation
  traitind1=rnorm(100, 50, 20)
    
  traitind2=rnorm(100, 50, 20)
    
  
#inserting the trait values in the data frame
  data=data.frame(list(traitind1,traitind2))
    
  colnames(data)=c('traitind1', 'traitind2')
    
#Calculating trait differences between each pair of individuals  
  data$diff=data$traitind1-data$traitind2
  
#Calculating the absolute values of trait differences
#between each pair of individuals  
  data$abs.diff=abs(data$traitind1-data$traitind2)
 
##Calculating the probability function that should determine
##the chances of inversion for each pair of individuals.
##The probability is a negative function of absolute trait
##difference between rivals.
#Prob1 - determine the probability function based on trait differences 
#between rivals. Still a positive function. The pair with the
#greatest trait difference will have probability of inversion
#equal to 1/div1
  data$prob1=(data$abs.diff/max(data$abs.diff))/div1

#Prob2 - Transforming prob 1 into a negative function of trait 
#differences between rivals. Now the pair with the smallest
#trait difference will have 1/div1 chance of showing a 
#victory inversion
  data$prob2=(rep(max(data$prob1),100))-data$prob1

#Prob3 - Transforming negative values (if they occur) of prob2
#in 0 (no chance of inversion)
  data$prob3=ifelse(data$prob2<0,0,data$prob2)
  
#Prob4 - Determining whether each pair should, in fact, show
#a victory inversion based on prob3
  data$prob4=rbinom(length(data$prob3),1,data$prob3) #1 should change the winner to the weaker individual
  
#Repating the same steps of prob1 to prob4, but using div2
#instead of div1
  data$prob5=(data$abs.diff/max(data$abs.diff))/div2
  data$prob6=(rep(max(data$prob5),100))-data$prob5
  data$prob7=ifelse(data$prob6<0,0,data$prob6)
  data$prob8=rbinom(length(data$prob7),1,data$prob7)
    
#Registering the trait value of winners and losers 
#when there is no chance of inversion (i.e. weapon differences
#are the most important in determining victory)
  data$trait.win.ne=ifelse(data$diff<0, data$traitind2, data$traitind1)
  data$trait.los.ne=ifelse(data$diff<0, data$traitind1, data$traitind2)
    
#Registering the trait value of winners and losers 
#when there is a chance of inversion determined by 1/div1
#(i.e. weapon differences are less important in determining
#victory when compared to the prior situation)
  data$trait.win.we1=ifelse(data$prob4==0, data$trait.win.ne, data$trait.los.ne)
  data$trait.los.we1=ifelse(data$prob4==0, data$trait.los.ne, data$trait.win.ne)
    
#Registering the trait value of winners and losers 
#when there is a chance of inversion determined by 1/div1
#(i.e. weapon differences are less important in determining
#victory when compared to the prior situation)
  data$trait.win.we2=ifelse(data$prob8==0, data$trait.win.ne, data$trait.los.ne)
  data$trait.los.we2=ifelse(data$prob8==0, data$trait.los.ne, data$trait.win.ne)
  
#Calculating the differences in mean trait values of winners
#and losers for the three scenarios (no error, 1/div1 and 1/div2)
  final.mean1[i]=mean(data$trait.win.ne)-mean(data$trait.los.ne)
  final.mean2[i]=mean(data$trait.win.we1)-mean(data$trait.los.we1)
  final.mean3[i]=mean(data$trait.win.we2)-mean(data$trait.los.we2)

#Calculating the variation in the trait value differences
#between winners and losers
  final.sd1w[i]=sd(data$trait.win.ne)
  final.sd2w[i]=sd(data$trait.win.we1)
  final.sd3w[i]=sd(data$trait.win.we2)
  final.sd1l[i]=sd(data$trait.los.ne)
  final.sd2l[i]=sd(data$trait.los.we1)
  final.sd3l[i]=sd(data$trait.los.we2)}

#Calculating the mean values for the scenario without
#error for  the 1000 simulation results
a=sort(final.mean1)
meana=mean(final.mean1)
cia=(a[975]-a[25])/2


#Calculating the mean values for the scenario with
#max probs of inversion equal to 1/div1 for the
#1000 simulation results
b=sort(final.mean2)
meanb=mean(final.mean2)
cib=(b[975]-b[25])/2

#Calculating the mean values for the scenario with
#max prob of inversion equal to 1/div2 for the
#1000 simulation results
c=sort(final.mean3)
meanc=mean(final.mean3)
cic=(c[975]-c[25])/2

#Plotting the final figure for mean trait values
par(mar=c(4,5,1,1), mgp=c(2.5,1,0))
plot(c(1,2,3), c(meana, meanb, meanc),ylim=c(meanc-cic,meana+cia),
     ylab='Mean difference between \n winners and losers', 
     xlab='Maximum chance victory by the smaller rival', bty='l', xaxt='n')
axis(1, at=c(1,2,3), labels=c('0','0.25','0.5'))
lines(c(1,1),c(meana-cia,meana+cia))
lines(c(2,2),c(meanb-cib,meanb+cib))
lines(c(3,3),c(meanc-cic,meanc+cic))



#Calculating the mean values of the standard deviation
# of winners and losers for the scenario without
#error for  the 1000 simulation results
##winners
d=sort(final.sd1w)
meand=mean(final.sd1w)
cid=(d[975]-d[25])/2

##losers
g=sort(final.sd1l)
meang=mean(final.sd1l)
cig=(g[975]-g[25])/2


#Calculating the mean values of the standard deviation
#of winners and losers for the scenario with
#max probability of inversion equal to 1/div1 for the
#1000 simulation results

##winners
e=sort(final.sd2w)
meane=mean(final.sd2w)
cie=(e[975]-e[25])/2

##losers
h=sort(final.sd2l)
meanh=mean(final.sd2l)
cih=(h[975]-h[25])/2

#Calculating the mean values of the standard deviation
#of winners and losers for the scenario with
#max probability of inversion equal to 1/div2 for the
#1000 simulation results

##winners
f=sort(final.sd3w)
meanf=mean(final.sd3w)
cif=(f[975]-f[25])/2

##losers
i=sort(final.sd3l)
meani=mean(final.sd3l)
cii=(i[975]-i[25])/2

#Plotting the final figure for the starndard deviation values
par(mar=c(4,5,1,1), mgp=c(2.5,1,0))
plot(c(1,2,4,5,7,8), c(meand, meang, meane, meanh, meanf, meani),
     ylim=c(meand-cid, meani+cii),pch=c(1,16),
     ylab='Mean standard deviation', 
     xlab='Maximum chance victory by the smaller rival', bty='l', xaxt='n')
axis(1, at=c(1.5,4.5,7.5), labels=c('0','0.25','0.5'))
lines(c(1,1),c(meand-cid,meand+cid))
lines(c(2,2),c(meang-cig,meang+cig), lty=2)
lines(c(4,4),c(meane-cie,meane+cie))
lines(c(5,5),c(meanh-cih,meanh+cih), lty=2)
lines(c(7,7),c(meanf-cif,meanf+cif))
lines(c(8,8),c(meani-cii,meani+cii), lty=2)
role.txt=c("Winners", "Losers")
legend(1,23, role.txt, pch = c(1,16), lty=c(1,2), bty='n')

