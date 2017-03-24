#Ehrharta experiment
#Analyzing community data from October 2014
#created by CAR 15 March 2017

########################
#   Biblioteca         #
                       #   
library(vegan)         #
library(dplyr)         #
library(knitr)         #
library(ggplot2)       #
library(tidyr)         #
library(mvabund)       #
                       #
########################


#go back one folder to get wd to the folder Ehrharta Project

setwd("../")

setwd("Main plots/Data/Tidy")

#### Reading Oct 2014 tidy files


oct14<-read.csv("oct2014_tidy_with_refs.csv", header=T)

# See header names and position
names(oct14)

#Avg Oct 14 data by plot

#Step 1: Select species data and Plot column
oct14avg <- oct14[c(1,7:40)]

#Step 2: Select species data by Plot
oct14avg<-oct14avg %>% group_by(Plot) %>% summarise_each(funs(mean))

#Herbicide 3, 7, and 8 have no plants, and aren't analyzable by Bray Curtis or the other
#distance measurements below. I'll probably have to use something like a Manhattan Distance

#Make "plot" the rowname
#There should be a better way to do this but for now I'm going to save file and reopen it 
#with rownames=True

write.csv(oct14avg, "oct2014_tidy_with_refs_averaged.csv", row.names = FALSE)
oct14avg<-read.csv("oct2014_tidy_with_refs_averaged.csv", header=T, row.names =1) 

View(oct14avg)

###################################
# Community Distance calculations #
###################################

#Using Manhattan Distances because my I have rows of all zeros, even after averaging across subplots

oct14.manhat<-vegdist(oct14avg,method="manhattan")

#Below is a way to compare distances, some of which are not compatible with my data re: zero rows
       #Thoughts from Phil Dixon: When I wonder what distance to use, I often ask whether the 
       #choice matters. If two sets of distances are highly correlated, then the choice doesn’t
       #matter. Here, I look at B-C without and with standardization by site total and 
       #Morisita-Horn distance. Note that MH distance computed from the unstandardized data is 
       #more similar to BC after standardization.

              #oct14.bc <- vegdist(oct14avg)
              #oct14.tot <- decostand(oct14avg, 'total')
              #oct14.bc2 <- vegdist(oct14.tot)
              #oct14.mh <- vegdist(oct14avg, method='horn')
              #pairs_plot<-pairs(cbind(oct14.bc, oct14.bc2, oct14.mh)) 

############################################
# Calulating nMDS with Manhattan Distances #
############################################ 

#Stress is the discrepancy between the Euclidean representation and the original distance 
#matrix. Want small. Here I'm comparing different stresses and number of iterations

#oct14.mds1 <- metaMDS(oct14.manhat, k=1, autotransform=FALSE)
 oct14.mds2 <- metaMDS(oct14.manhat, k=2, trymax=200, autotransform = FALSE)
#oct14.mds3 <- metaMDS(oct14.manhat, k=3, autotransform=FALSE)
#oct14.mds4 <- metaMDS(oct14.manhat, k=4, autotransform=FALSE)
#oct14.allmds <- list(oct14.mds1, oct14.mds2, oct14.mds3, oct14.mds4)
#kable(cbind(Dim=1:4, Stress=sapply(oct14.allmds, function(x){x$stress})),
#            caption="Stress for 1, 2, 3, and 4 D")

#Dim|    Stress|
#  |---:|---------:|
#  |   1| 0.1046711|
#  |   2| 0.0473367|
#  |   3| 0.0305746|
#  |   4| 0.0201499|

#Stress is low (<0.1) for 2 dimensions so I´ll stick with k=2

plot(oct14.mds2)

#################
# Graphing nMDS #
#################

#Using the scores function from vegan to extract the site scores and convert to a data.frame
oct14_scores <- as.data.frame(scores(oct14.mds2))  
# create a column of site names, from the rownames of data.scores
oct14_scores$Plot <- rownames(oct14_scores) 
#look at the data
head(oct14_scores)  

#adding two columns, one with plot number and the other with treatment
oct14_scores <- oct14_scores %>% separate(Plot, c("Site", "Treatment"))

as.data.frame(oct14_scores)

#Make Site a factor
as.factor(oct14_scores$Site)

#look at the data
View(oct14_scores)

#make hull data for treatments
herb <- oct14_scores[oct14_scores$Treatment == "Herbicide", ][chull(oct14_scores[oct14_scores$Treatment == "Herbicide", c("NMDS1", "NMDS2")]), ]
pull <- oct14_scores[oct14_scores$Treatment == "Pull", ][chull(oct14_scores[oct14_scores$Treatment == "Pull", c("NMDS1", "NMDS2")]), ]
cont <- oct14_scores[oct14_scores$Treatment == "Control", ][chull(oct14_scores[oct14_scores$Treatment == "Control", c("NMDS1", "NMDS2")]), ]
ref <- oct14_scores[oct14_scores$Treatment == "Reference", ][chull(oct14_scores[oct14_scores$Treatment == "Reference", c("NMDS1", "NMDS2")]), ]

hull.data <- rbind(herb, pull, cont, ref)  #combine hulls
hull.data

K<-ggplot(oct14_scores,aes(NMDS1,NMDS2)) 
L<-K+ geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Treatment,group=Treatment),alpha=0.30) + geom_jitter(aes(colour=factor(Site), shape=Treatment), size=2) +theme_classic() 

#Remove legend for Site with guides()
L + guides(colour=FALSE)

#save graphic
setwd("../../")
setwd("Graphics")

ggsave("NMDS_Oct14.png")

###################################################
# Analyzing NMDS - Significance Across Treatments #
################################################### 

# Does Treatment method influence spp composition?
# using Manhattan as our choice of distance

#Select only plant data 
oct14data <- oct14[7:40]

#Make plant data a matrix
oct14.m <- as.matrix(oct14data)
oct14.mva <- manyglm(oct14.m ~ Treatment, data=oct14)

anova(oct14.mva)
# warning: slow!

#Time elapsed: 0 hr 6 min 31 sec
#Analysis of Deviance Table

#Model: manyglm(formula = oct14.m ~ Treatment, data = oct14)

#Multivariate test:
#  Res.Df Df.diff   Dev Pr(>Dev)    
#(Intercept)    128                           
#Treatment      125       3 256.9    0.001 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Arguments:
#  Test statistics calculated assuming uncorrelated response (for faster computation) 
#P-value calculated using 999 resampling iterations via PIT-trap resampling (to account for correlation in testing).

summary(oct14.mva)
#  warning: slow
#  coefficient-specific information, also slow because perm. p-values

#Test statistics:
#  wald value Pr(>wald)    
#(Intercept)            17.650  0.000999 ***
#  TreatmentHerbicide      6.629  0.000999 ***
#  TreatmentPull           4.361  0.036963 *  
#  TreatmentReference      8.634  0.000999 ***
#  --- 
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#Test statistic:  12.19, p-value: 0.000999 
#Arguments:
#  Test statistics calculated assuming response assumed to be uncorrelated 
#P-value calculated using 1000 resampling iterations via pit.trap resampling (to account for correlation in testing).

oct14.mva
# prints sp specific information about each fit

######################################################
# Philip's code for 'non-Ehrharta' spp comp analysis #
######################################################

#library(dplyr)
#library(tidyr)


#The tidy data if needed
     #oct14<-read.csv("oct2014_tidy_with_refs.csv", header=T)

temp <- gather(oct2014, 'species', 'cover', 8:41)
tot.cover <- temp %>% group_by(Treatment, Site, species) %>% summarize(cover=sum(cover))
tot.cover2 <- spread(tot.cover, species, cover)

nonE.m <- as.matrix(tot.cover2[,3:36])

#library(mvabund)

nonE.mv <- manyglm(nonE.m ~ Site + Treatment, data=tot.cover2)
anova(nonE.mv)

#Time elapsed: 0 hr 4 min 31 sec
#Analysis of Deviance Table

#Model: manyglm(formula = nonE.m ~ Site + Treatment, data = tot.cover2)

#Multivariate test:
#  Res.Df Df.diff    Dev Pr(>Dev)
#(Intercept)     32                        
#Site            31       1  44.04    0.161
#Treatment       29       2 100.62    0.127
#Arguments:
#Test statistics calculated assuming uncorrelated response (for faster computation) 
#P-value calculated using 999 resampling iterations via PIT-trap resampling
#(to account for correlation in testing).

rownames(nonE.m) <- paste(tot.cover2$Treatment, tot.cover2$Site)

#library(vegan)
nonE.nmds <- metaMDS(nonE.m, k=2, autotransform=F, expand=F, dist='manh')
plot(nonE.nmds, display='sites', type='t')

