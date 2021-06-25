###############################################################

## Practical subterranean conservation
## Mammola S. et al.

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.0.3) and R studio (v. 1.4.1103)
# Authors: Stefano Mammola
# Location: Helsinki, June 2021

###############################################################

# Working directory -------------------------------------------------------

setwd("/Users/stefanomammola/Desktop/PAPERS IN CORSO/Troglo_achille/Analysis") #change with your working directory

# Loading R package -------------------------------------------------------


# clean the workspace -----------------------------------------------------

rm(list=ls())

# Loading useful functions ------------------------------------------------

library("BAT")
library("ggplot2")
library("gridExtra")
library("hypervolume")
library("labdsv")
library("lme4")



library(NbClust)
library(fpc)
library(clusteval)
library(ape)
library(geiger)
library(phytools)


# Useful functions --------------------------------------------------------

## functions for logistic curves

logisticline <- function(z,model) {
  eta <- model$coefficients[1]+model$coefficients[2]*z;
  1 / (1 + exp(-eta))
}

logisticline_min <- function(z,model) {
  eta <- model$coefficients[1]+model$coefficients[2]*z - 1.96*summary(model)$coefficients[2,2];
  1 / (1 + exp(-eta))
}

logisticline_max <- function(z,model) {
  eta <- model$coefficients[1]+model$coefficients[2]*z + 1.96*summary(model)$coefficients[2,2];
  1 / (1 + exp(-eta))
}

# Parameters for plots ----------------------------------------------------

xlab_thermal <- expression(Delta * "T (Critical temperature - Cave temperature) (°C)")
ylab_thermal <- "Survival rate"
color_plot   <- c("darkmagenta","orange","turquoise")

###########################################
#### 1. Morphospace analysis
###########################################

# Loading the dataset ----------------------------------------------------

db_morpho <- read.csv("Morphology.csv", header = TRUE, sep=";", dec=".", as.is = FALSE)
str(db_morpho)

#Selecting females
T.matrix <- db_morpho[db_morpho$Sex=="femmina",]

#Calculating gower distance
gower.mat <- BAT::gower(T.matrix[,4:21])

#Calculating PCoA
euc.pco <- labdsv::pco(gower.mat,k=3)
barplot(euc.pco$eig)
plot(euc.pco)

#Extracting the first 4 PC
HV <- data.frame(species = T.matrix$Species,
                 PC1     = euc.pco$points[,1],
                 PC2     = euc.pco$points[,2],
                 PC3     = euc.pco$points[,3])

n.var=4
hv <- hypervolume_gaussian(subset(HV, species==levels(HV$species)[[max(nlevels(HV$species))]])[,c(2:n.var)],
                          name = levels(HV$species)[[max(nlevels(HV$species))]],
                          kde.bandwidth = estimate_bandwidth(subset(HV,species==levels(HV$species)[[max(nlevels(HV$species))]])[,c(2:n.var)],
                          method = "cross-validation"))

for (i in (nlevels(factor(HV$species))-1):1) {
  
  hv2 = hypervolume_gaussian(subset(HV,species==levels(HV$species)[[i]])[,c(2:n.var)],name=levels(HV$species)[[i]],kde.bandwidth = estimate_bandwidth(subset(HV,species==levels(HV$species)[[i]])[,c(2:n.var)],method="cross-validation"))
  hv  = hypervolume_join(hv2,hv)
  
}

#plotting
plot(hv,
     num.points.max.random = 10000,
     show.data=F,
     pairplot=F,
     show.3d=F,
     showdensity=T,
     showrandom=F,
     contour.lwd=0.5,
     cex.names=0.8,
     cex.random=0.3,
     cex.legend=1,
     cex.axis=1,
     contour.filled=T,
     contour.filled.alpha=0.2,
     show.centroid=TRUE,
     cex.centroid=2,
     colors = color_plot)

#extracting stats
BAT::kernel.alpha(hv)*3 #volume
BAT::kernel.similarity(hv)$Distance_centroids #distance centroids
BAT::kernel.beta(hv) #beta diversity

###########################################
#### 2. Thermal tolerance analysis
###########################################

# Loading the dataset ----------------------------------------------------

db_thermal2 <- read.table("Thermal.csv", header = TRUE, sep=";", dec=",", as.is = FALSE)

## Removing controls and missing data
db_thermal2 <- db_thermal2[db_thermal2$Control=="no",]
db_thermal2 <- subset(db_thermal2, !is.na(Delta_T))

# Is the lethal T different according to the cave? ------------------------

for (i in levels(factor(db_thermal2$Species))) {
  
  specie <- db_thermal2[db_thermal2$Species == i,]
  
  if (nlevels(factor(specie$Cave)) < 2) {
    print("Troglohyphantes")
    print(i)
    print("Not possible to fit ANOVA: species recorded for just one cave")
    print("###################################")
  } 
  
  else {
    print("Troglohyphantes")
    print(i)
    print(summary(aov(Delta_T ~ Cave,data = specie)))
    print(summary(lm(Delta_T ~ Cave,data = specie)))
    print("###################################")
  }} 

#Yes! Tuna from Tana for achillis

#clean
rm(i,specie)

# Plot --------------------------------------------------------------------

#Rename
levels(db_thermal2$Cave)    <- c("Bocetto","Buco_valenza","Camoscere_sup","Buco_del_Nebin","Tana_del_Diavolo","Tornini","Tuna")
levels(db_thermal2$Species) <- c("T. achillis","T. delphinicus", "T. vignai")

#Order
db_thermal2$Cave <- factor(db_thermal2$Cave,
                   levels = c("Buco_valenza","Camoscere_sup","Buco_del_Nebin","Bocetto","Tana_del_Diavolo","Tornini","Tuna"),ordered = TRUE)

#Plot
(thermal_plot2 <- 
  ggplot(db_thermal2, aes(x=Cave, y=Delta_T,fill=Species)) + 
  geom_boxplot(outlier.colour="white",outlier.size=2)+
  coord_flip()+
  ylab(xlab_thermal)+
  xlab("")+
  ylim(0,20)+
  scale_fill_manual(values=c("darkmagenta","orange","turquoise"))+
  scale_x_discrete(labels=c("Buco di\nValenza","Grotta delle\nCamoscere sup.","Buco del\nNebin","Bocetto\nmine","Tana del\nDiavolo", "Tornini\nmine","Tuna do\nDiau"))+
  annotate("text",x=6.91,y=17.9, label="*",size=8)+
  annotate("text",x=1,y=19.5, label="(b)",size=11,fontface="bold")+
  theme_classic()+
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid"),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=16, angle=0),
        axis.ticks.x = element_blank(),
        legend.text = element_text( size = 12, face = "italic"),
        legend.title= element_text(size=0),
        legend.position = c(0.8, 0.6))
)

# Dropping the outlier for subsequent analyses ----------------------------

db_thermal <- droplevels(db_thermal2[-which(db_thermal2$Cave == "Tuna"), ] )

# Summary statistics ------------------------------------------------------

# Lethal Temperature
tapply(db_thermal$Lethal_T,db_thermal$Species,mean)
tapply(db_thermal$Lethal_T,db_thermal$Species,sd)
tapply(db_thermal$Lethal_T,db_thermal$Species,min)
tapply(db_thermal$Lethal_T,db_thermal$Species,max)

# Delta temperature
tapply(db_thermal$Delta_T,db_thermal$Species,mean)
tapply(db_thermal$Delta_T,db_thermal$Species,sd)
tapply(db_thermal$Delta_T,db_thermal$Species,min)
tapply(db_thermal$Delta_T,db_thermal$Species,max)

# Survival rate estimation ------------------------------------------------

x  <- as.data.frame(matrix(NA, ncol = 4, nrow =1)) ; colnames(x) <-c("Lethal","Species","dead","alive")
x2 <- x

for (i in levels(factor(db_thermal$Species))) {
      specie <- db_thermal[db_thermal$Species==i,]
      LT <- seq(from=0,to=max(max(db_thermal$Delta_T)),by=1)
      ALIVE <- c()
      DEAD  <- c()
      
      for(k in 1:length(LT))
      {
        alive <- sum(specie$Delta_T>LT[k])
        dead  <- nrow(specie)-alive
        ALIVE <- append(ALIVE,alive)
        DEAD  <- append(DEAD,dead)
      }
      x3    <- data.frame(Lethal=LT,Species=rep(i,length(LT)),dead=DEAD,alive=ALIVE)
      x2    <- rbind(x2,x3)
      ALIVE <- c()
      DEAD  <- c()
    }

x <- rbind(x,x2) ; x <- x[complete.cases(x), ]
x$Species <- as.factor(x$Species)

#clean the space
rm(alive,ALIVE,dead,DEAD,i,k,LT,x2,x3,specie)

# Is the lethal T different according to the developmental space? ----------------------

levels(db_thermal$Gender) <- list(f = c("f"), m = c("m"), i = c("f_i","m_i","i")) 

m0 <- aov(Delta_T ~ Gender + Error(Species), data = db_thermal)
summary(m0) # no! 

# Is the lethal T different between the species? --------------------------

m1 <- aov(Delta_T ~ Species, data = db_thermal)
summary(m1)

#post-hoc
TukeyHSD(m1,"Species")

# Plotting the lethal curves ----------------------------------------------

# Fitting the LT curves
list_DB <- list() 
model   <- list()

for (i in levels(factor(x$Species))) {
  species      <- x[x$Species==i,]
  temperature  <- species$Lethal 
  alive        <- species$alive
  dead         <- species$dead
  list_DB[[i]] <- species
  model[[i]]   <- glm(cbind(alive,dead)~temperature,binomial)
}  

y <- seq(min(x$Lethal,na.rm=T), max(x$Lethal,na.rm=T)+2,0.01)

rm(alive,dead,i,temperature,list_DB)

## Plot

#Fixing issues.
x[19,3:4] =c(2,10)
x[20,3:4] =c(10,2)
x1 <- x[x$Species=="T. achillis",]
x2 <- x[x$Species=="T. delphinicus",]
x3 <- x[x$Species=="T. vignai",]

(thermal_plot1 <- 
  ggplot()+
  xlab(xlab_thermal)+ 
  ylab(ylab_thermal)+
  xlim(0,20)+
  geom_line(aes(y = logisticline(y,model[[1]]), x = y), colour = color_plot[1],linetype="solid",size=1.1)+
  geom_line(aes(y = logisticline(y,model[[2]]), x = y), colour = color_plot[2],linetype="solid",size=1.1)+
  geom_line(aes(y = logisticline(y,model[[3]]), x = y), colour = color_plot[3],linetype="solid",size=1.1)+
  
  geom_point(aes(y=((x1$alive)/(x1$alive + x1$dead)),x=x1$Lethal),col=color_plot[1],alpha=0.5)+
  geom_point(aes(y=((x2$alive)/(x2$alive + x2$dead)),x=x2$Lethal),col=color_plot[2],alpha=0.5)+
  geom_point(aes(y=((x3$alive)/(x3$alive + x3$dead)),x=x3$Lethal),col=color_plot[3],alpha=0.5)+
  
  geom_rect(aes(xmin=10, xmax=16.5, ymin=0.7, ymax=0.94), fill="white",linetype="solid")+
  
  annotate("segment", x = 10.5, xend = 11, y = .9, yend = .9, colour = color_plot[1],size=1)+
  annotate("segment", x = 10.5, xend = 11, y = .82, yend = .82, colour = color_plot[2],size=1)+
  annotate("segment", x = 10.5, xend = 11, y = .74, yend = .74, colour = color_plot[3],size=1)+
  annotate("text",label="T. achillis",col="black", hjust = 0,x=11.5,y=.9,fontface="italic",size=6)+
  annotate("text",label="T. delphinicus",col="black", hjust = 0,x=11.5,y=.82,fontface="italic",size=6)+
  annotate("text",label="T. vignai",col="black", hjust = 0,x=11.5,y=0.74,fontface="italic",size=6)+
 
  annotate("text",x=20,y=0.07, label="(a)",size=11,fontface="bold")+
    theme_classic()+
    theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid"),
          axis.title.y = element_text(size=16),
          axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=14),
          axis.text.x = element_text(size=16, angle=0),
          axis.ticks.x = element_blank())
)
  
# Arrange in a plate ------------------------------------------------------

pdf('Figure_Thermal.pdf',width=16,height=7)

grid.arrange(thermal_plot1, thermal_plot2, nrow = 1)

dev.off()

###########################################
#### 3. Analyses with genitalia
###########################################

# Loading the datasets ----------------------------------------------------

LAM <- read.table("Lamella.csv", header = TRUE, sep=";", dec=",", as.is = FALSE)
EPI <- read.table("Epigynum.csv", header = TRUE, sep=";", dec=",", as.is = FALSE)

# Plotting ----------------------------------------------------------------

(p1 <- ggplot(LAM, aes(x=Specie, y=A/B)) +
  
  annotate("text",label="(a)",fontface="bold", 
           x= 3.2, 
           y = (0.75*(min(LAM$A/LAM$B)+max(LAM$A/LAM$B))), hjust = 0, vjust = 1,size=12)+
  
  annotate("text",label="(a)",fontface="bold", 
           x= 3.2, 
           y = (0.75*(min(LAM$A/LAM$B)+max(LAM$A/LAM$B))), hjust = 0, vjust = 1,size=12)+
  
  geom_boxplot(fill=c("darkmagenta","orange","turquoise"),outlier.colour="black", outlier.shape=1,outlier.size=1)+
  labs(x ="",y="Lamella - A/B",size=14)+
  scale_x_discrete(labels=c("T. achillis","T. delphinicus","T. vignai"))+
  theme_classic()+
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid"),
        axis.title.y = element_text(face="bold",size=16),
        axis.text.y = element_text(angle = 90,hjust=0.5,size=14),
        axis.text.x = element_text(face="italic",size=14, angle=0),
        axis.ticks.x = element_blank())
)

(p2 <- ggplot(EPI, aes(x=Specie, y=A/B)) +
  geom_boxplot(fill=c("darkmagenta","orange","turquoise"),outlier.colour="black", outlier.shape=1,outlier.size=1)+
  annotate("text",label="(b)",fontface="bold", 
           x= 3.2, 
           y = (0.75*(min(EPI$A/EPI$B)+max(EPI$A/EPI$B))), hjust = 0, vjust = 1,size=12)+
  labs(x ="",y="Epigyne - C/D")+
  scale_x_discrete(labels=c("T. achillis","T. delphinicus","T. vignai"))+
  theme_classic()+
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid"),
        axis.title.y = element_text(face="bold",size=16),
        axis.text.y = element_text(angle = 90,hjust=0.5,size=14),
        axis.text.x = element_text(face="italic",size=14, angle=0),
        axis.ticks.x = element_blank())
)

# Is the difference significant? ----------------------------------------------
m3 <- aov(A/B ~ Specie, data = LAM)
summary(m3)
TukeyHSD(m3)

m4 <- aov(A/B ~ Specie, data = EPI)
summary(m4)
TukeyHSD(m4)

# Arrange in a plate ------------------------------------------------------

pdf('Figure_genitals.pdf',width=12,height=5)

grid.arrange(p1,p2,nrow=1)

dev.off()