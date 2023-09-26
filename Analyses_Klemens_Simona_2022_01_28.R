### TESTING THE INTERACTION BETWEEN A HEATWAVE AND BIOTURBATING FAUNA FOR SEED BURIAL (and Germination)
rm(list=ls())
graphics.off()
library(vegan)
library(plyr)
library(lattice)
library(car)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(lme4)
library(Rmisc)

seed <- read.table("data.txt", header=T, sep="\t")
seed=seed[81:240,]

###### ANALYSES
### SEED DEPTH - FIG 4a and statistical results
hist(seed$seed_depth)
hist(seed$seed_depth^0.5)
leveneTest(seed_depth ~ fauna*temp, data=seed)
leveneTest(seed_depth^0.5 ~ fauna*temp, data=seed)
# non transformed is better for Levenes - but the square rooted on looks much better with histogram.

# set contrast for summary command
seed$fauna=relevel(seed$fauna, ref="NO")

seed_depth<-glm(seed_depth^0.5~fauna*temp + run + (temp/block), data=seed)
summary(seed_depth)
Anova(seed_depth)

datac <- summarySEwithin(seed, measurevar="seed_depth", withinvars=c("fauna","temp"))
datac

g <- ggplot(datac, aes(x=fauna, y=seed_depth, fill=temp)) + xlab("") + ylab("mean depth of seeds (cm)") +
  geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=seed_depth-se, ymax=seed_depth+se)) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_x_discrete(position = "bottom",
                   limits = c("NO", "COCKLE", "SHRIMP", "WORM", "MIX"), 
                   labels = c("No Fauna", "Cockles","Shrimps", "Worms", "Mixed Fauna")) +
  scale_y_continuous(breaks=seq(1:100)) +
  theme_bw() +
  geom_hline(yintercept=38) 
g
# ggplot(seed, aes(x=fauna, y=seed_depth, fill=temp)) + xlab("") + ylab("median depth of seeds (cm)") +
  geom_boxplot(aes(fill=factor(temp))) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_x_discrete(position = "bottom",
                   limits = c("NO", "COCKLE", "SHRIMP", "WORM", "MIX"), 
                   labels = c("No Fauna", "Cockles","Shrimps", "Worms", "Mixed Fauna")) +
  scale_y_continuous(breaks=seq(1:100)) +
  theme_bw() +
  geom_hline(yintercept=38) 

seed$fauna=relevel(seed$fauna, ref="NO")

X0_nb<-glm(X0_nb ~ fauna * temp + run + (temp/block), family=poisson, data=seed)
summary(X0_nb) 
# Strong effects of fauna - cockles and mix decrease the nb of seeds on the surface significantly

library(AER)
dispersiontest(X0_nb) # There is significat overdispersion - I therefore use quasipoisson model

X0_nb<-glm(X0_nb ~ fauna * temp + run + (temp/block), family=quasipoisson, data=seed)
summary(X0_nb) 
# Still strong effects of fauna - cockles and mix decrease the nb of seeds on the surface significantly

# datac <- summarySEwithin(seed, measurevar="X0_nb", withinvars=c("fauna","temp"))
# datac

# ggplot(datac, aes(x=fauna, y=X0_nb, fill=temp)) + xlab("") + ylab("mean number of seeds not buried (cm)") +
  geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=X0_nb-se, ymax=X0_nb+se)) +
  coord_cartesian(ylim=c(0,10)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_x_discrete(position = "bottom",
                   limits = c("NO", "COCKLE", "SHRIMP", "WORM", "MIX"), 
                   labels = c("No Fauna", "Cockles","Shrimps", "Worms", "Mixed Fauna")) +
  scale_y_continuous(breaks=seq(1:100)) +
  theme_bw() +
  geom_hline(yintercept=38) 

ggplot(seed, aes(x=fauna, y=X0_nb, fill=temp)) + xlab("") + ylab("median number of seeds not buried (cm)") +
  geom_boxplot(aes(fill=factor(temp))) +
  coord_cartesian(ylim=c(0,12)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_x_discrete(position = "bottom",
                   limits = c("NO", "COCKLE", "SHRIMP", "WORM", "MIX"), 
                   labels = c("No Fauna", "Cockles","Shrimps", "Worms", "Mixed Fauna")) +
  scale_y_continuous(breaks=seq(1:100)) +
  theme_bw() +
  geom_hline(yintercept=38) 

# Buried seeds (the number of seeds minus those found on the surface) FIG 4b and statistical results

seed$nb_deep=seed$seed_nb-seed$X0_nb
seed$fauna=relevel(seed$fauna, ref="NO")

nb_deep<-glm(nb_deep ~ fauna * temp + run + (temp/block), family=poisson, data=seed)
summary.glm(nb_deep)

library(AER)
dispersiontest(nb_deep)

# There is significant overdispersion of the data  - but not much (more than 1 but less than 2).
# Still, I chose to follow the recommendation and use a quasipoisson model


nb_deep<-glm(nb_deep ~ fauna * temp + run + (temp/block), family=quasipoisson, data=seed)
summary.glm(nb_deep)

datac <- summarySEwithin(seed, measurevar="nb_deep", withinvars=c("fauna","temp"))

# ggplot(datac, aes(x=fauna, y=nb_deep, fill=temp)) + xlab("") + ylab("mean number of buried seeds (#)") +
  geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=nb_deep-se, ymax=nb_deep+se)) +
  coord_cartesian(ylim=c(0,10)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_x_discrete(position = "bottom",
                   limits = c("NO", "COCKLE", "SHRIMP", "WORM", "MIX"), 
                   labels = c("No Fauna", "Cockles","Shrimps", "Worms", "Mixed Fauna")) +
  scale_y_continuous(breaks=seq(1:100)) +
  theme_bw() +
  geom_hline(yintercept=38) 

ggplot(seed, aes(x=fauna, y=nb_deep, fill=temp)) + xlab("") + ylab("median number of buried seeds (#)") +
  geom_boxplot(aes(fill=factor(temp))) +
  coord_cartesian(ylim=c(0,12)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_x_discrete(position = "bottom",
                   limits = c("NO", "COCKLE", "SHRIMP", "WORM", "MIX"), 
                   labels = c("No Fauna", "Cockles","Shrimps", "Worms", "Mixed Fauna")) +
  scale_y_continuous(breaks=seq(1:100)) +
  theme_bw() +
  geom_hline(yintercept=38) 

###################################################################################################

# PERCENT GERMINATED USING SURFACE or BURIED SEEDS AS FACTOR IN ANALYSIS

library(nlme)

depth <- read.table("depth_data.txt", header=T, sep="\t")

depth=depth[-which(depth$run==1),c(1:13)]
depth=na.omit(depth) # cleans the data from NA's

# NO FAUNA, COLD, SHALLOW SEEDS AS CONTRAST
depth$fauna=as.factor(depth$fauna)
depth$temp=as.factor(depth$temp)
depth$Germination1=as.factor(depth$Germination1)

depth$fauna=relevel(depth$fauna, ref="NO")
depth$temp=relevel(depth$temp, ref="COLD")
depth$Germination1=relevel(depth$Germination1, ref="0_germ_perc")

germination<-lme(asin(value1^0.25) ~ Germination1 * fauna * temp + run + (temp/block), 
                 random=~1|Sample.ID, data=depth)
summary(germination)

# NO FAUNA, WARM, SHALLOW SEEDS AS CONTRAST
depth$fauna=relevel(depth$fauna, ref="NO")
depth$temp=relevel(depth$temp, ref="WARM")
depth$Germination1=relevel(depth$Germination1, ref="0_germ_perc")

germination<-lme(asin(value1^0.25) ~ Germination1 * fauna * temp + run + (temp/block), 
                 random=~1|Sample.ID, data=depth)
summary(germination)

# NO FAUNA, COLD, BURIED SEEDS AS CONTRAST
depth$fauna=relevel(depth$fauna, ref="NO")
depth$temp=relevel(depth$temp, ref="COLD")
depth$Germination1=relevel(depth$Germination1, ref="germ_perc_deep")

germination<-lme(asin(value1^0.25) ~ Germination1 * fauna * temp * run + (temp/block), 
                 random=~1|Sample.ID, data=depth)
summary(germination)

# NO FAUNA, WARM, BURIED SEEDS AS CONTRAST
depth$fauna=relevel(depth$fauna, ref="NO")
depth$temp=relevel(depth$temp, ref="WARM")
depth$Germination1=relevel(depth$Germination1, ref="germ_perc_deep")

germination<-lme(asin(value1^0.25) ~ Germination1 * fauna * temp * run + (temp/block), 
                 random=~1|Sample.ID, data=depth)
summary(germination)

#germination between runs:
library(Rmisc)
datarun = summarySEwithin(depth, measurevar = "value1", withinvars = "run")
datarun

### Making the FIGURES ###### Separating cold and warm 

# figure with "COLD" data
cold <- depth[which(depth$temp=="COLD"),]
datac =  summarySEwithin(cold, measurevar="value1", withinvars=c("fauna","Germination1"))
datac

depth$fauna=relevel(depth$fauna, ref="NO")
depth$Germination1=relevel(depth$Germination1, ref="X0_germ_perc")

ggplot(datac, aes(x=fauna, y=value1, fill=Germination1)) + xlab("") + ylab("average germinated seeds (%)") +
  geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=value1-se, ymax=value1+se)) +
  coord_cartesian(ylim=c(0,0.3)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_x_discrete(position = "bottom",
                   limits = c("NO", "COCKLE", "SHRIMP", "WORM", "MIX"), 
                   labels = c("No Fauna", "Cockles","Shrimps", "Worms", "Mixed Fauna")) +
  scale_y_continuous(breaks=seq(1:100)) +
  theme_bw() +
  geom_hline(yintercept=38) 

# figure with "WARM" data

warm <- depth[which(depth$temp=="WARM"),]
datac =  summarySEwithin(warm, measurevar="value1", withinvars=c("fauna","Germination1"))
datac

depth$fauna=relevel(depth$fauna, ref="NO")
depth$Germination1=relevel(depth$Germination1, ref="0_germ_perc")

ggplot(datac, aes(x=fauna, y=value1, fill=Germination1)) + xlab("") + ylab("average germinated seeds (%)") +
  geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=value1-se, ymax=value1+se)) +
  coord_cartesian(ylim=c(0,0.3)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_x_discrete(position = "bottom",
                   limits = c("NO", "COCKLE", "SHRIMP", "WORM", "MIX"), 
                   labels = c("No Fauna", "Cockles","Shrimps", "Worms", "Mixed Fauna")) +
  scale_y_continuous(breaks=seq(1:100)) +
  theme_bw() +
  geom_hline(yintercept=38) 
