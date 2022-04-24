rm(list=ls())
library(metafor)
library(scales) 
library(tidyverse)

## This is the code for plotting Figure 3 ####

full<-read.csv("full_dataset_weapons2.csv",h=T,sep=';')
weapons<-full[full$trait=='weapon',]

forest.data <-
  weapons %>%
  select(order,dyad.type,sp,yi,vi,N,study,measure,display,fighting.style)

## Now that the data frame has been re-arranged, we can count the fighting styles and displays ####
## to make Table 3 ####

count(weapons,fighting.style,display,sp) 

### Moving forward to prepare what is going to be plotted

test<-weapons[order(weapons$order,weapons$sp,weapons$dyad.type,weapons$measure,weapons$display,weapons$fighting.style),]

ilab.text<-data.frame(test$dyad.type,test$measure,test$display,test$fighting.style,test$study)
colnames(ilab.text)<-c("Pairing method","Component","Display","Fighting Style", "Study")

ilab.text$`Pairing method`<-ifelse(ilab.text$`Pairing method`=="random","Random","Paired")


ilab.text$Display <- case_when(
  ilab.text$Display == "weapon" ~ "Weapon",
  ilab.text$Display == "body" ~ "Body",
  ilab.text$Display == "no" ~ "No"
)

ilab.text$Component <- case_when(
  ilab.text$Component == "linear" ~ "Linear",
  ilab.text$Component == "performance" ~ "Performance",
  ilab.text$Component == "area" ~ "Area",
  ilab.text$Component == "asymmetry" ~ "Asymmetry",
  ilab.text$Component == "mass" ~ "Mass",
  ilab.text$Component == "index" ~ "Index"
)

unique(ilab.text$`Fighting Style`)

ilab.text$`Fighting Style` <- case_when(
  ilab.text$`Fighting Style` == "push" ~ "Push",
  ilab.text$`Fighting Style` == "pierce.push" ~ "Pierce & Push",
  ilab.text$`Fighting Style` == "squeeze.pierce" ~ "Squeeze & Pierce",
  ilab.text$`Fighting Style` == "push.lift.squeeze" ~ "Push, Lift & Squeeze",
  ilab.text$`Fighting Style` == "lift.push" ~ "Lift & Push",
  ilab.text$`Fighting Style` == "lift.squeeze" ~ "Lift & Squeeze",
  ilab.text$`Fighting Style` == "push.squeeze.pull.lift" ~ "Push, Squeeze, Pull & Lift",
  ilab.text$`Fighting Style` == "squeeze.pull" ~ "Squeeze & Pull",
  ilab.text$`Fighting Style` == "squeeze.push" ~ "Squeeze & Push",
  ilab.text$`Fighting Style` == "push.impact.squeeze" ~ "Push, Impact & Squeeze",
  ilab.text$`Fighting Style` == "push.pull.squeeze" ~ "Push, Pull & Squeeze",
  ilab.text$`Fighting Style` == "squeeze.impact" ~ "Squeeze & Impact",
  ilab.text$`Fighting Style` == "impact" ~ "Impact",
  ilab.text$`Fighting Style` == "squeeze" ~ "Squeeze",
  ilab.text$`Fighting Style` == "push.impact.lift" ~ "Push, Impact & Lift",
  ilab.text$`Fighting Style` == "lift.impact" ~ "Lift & Impact",
  ilab.text$`Fighting Style` == "impact.push" ~ "Impact & Push"
)


count(ilab.text, `Fighting Style`, Display) ### These are the 


## Italicizing all species names

# Creating expression commands using paste()
ordem = paste("expression(italic('", test$sp, "'))", sep="")
ordem

# Using sapply and eval(parse()) to generate a list of expressions to be used in the plot
ital = sapply(ordem, function(x){eval(parse(text=x))})


pdf("mygraph.pdf",w=15,h=25)
forest(x=test$yi,vi=test$vi, xlim=c (-20,10),ylim=c(3,120),lwd = 2,
       ilab=ilab.text,ilab.xpos=c(-14,-11.5,-9,-6,-3),  bg='grey50',
       rows=c(1, 3, 5, 7:23, 25:60, 62:63, 65:74, 76, 78:96, 98:101,
              103:117),
       slab=ital,cex=1,efac=0.2,pch=21, xlab=expression(italic("Hedges' " * g)))
text(-19,118.5,expression(bold("Species")),cex=1.1)
text(-14,118.5,expression(bold("Pairing")),cex=1.1)
text(-11.5,118.5,expression(bold("Measure")),cex=1.1)
text(-9,118.5,expression(bold("Display")),cex=1.1)
text(-6,118.5,expression(bold("Fighting Style")),cex=1.1)
text(-3,118.5,expression(bold("Study")),cex=1.1)
text(8.2,118.5,expression(bold("Hedges' ")~bolditalic("g [95% CI]")),cex=1.1)

dev.off()


#### The rectangles and outlines were added in Photoshop
