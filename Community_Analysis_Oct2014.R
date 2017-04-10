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
# Euclidean is another option but it takes the square of all values to get rid of negatives
# where as Manhattan takes the absolute vale, so it's less driven in large differences
# between abundances

oct14.manhat<-vegdist(oct14avg,method="manhattan")

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
oct14.matrix <- as.matrix(oct14data)
oct14.mva <- manyglm(oct14.matrix ~ Treatment, data=oct14)

anova(oct14.mva)
# warning: slow!

#Time elapsed: 0 hr 6 min 31 sec
#Analysis of Deviance Table

#Model: manyglm(formula = oct14.matrix ~ Treatment, data = oct14)

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
nonE.m <- nonE.m[,-7]


#library(mvabund)

nonE.mv <- manyglm(nonE.m ~ Site.f + Treatment, data=tot.cover2)
anova(nonE.mv)

nonE.mv <- manyglm(nonE.m ~ Site.f + Treatment, data=tot.cover2)
anova(nonE.mv)

#Time elapsed: 0 hr 3 min 9 sec
#Analysis of Deviance Table

#Model: manyglm(formula = nonE.m ~ Site.f + Treatment, data = tot.cover2)

#Multivariate test:
#  Res.Df Df.diff   Dev Pr(>Dev)    
#(Intercept)     32                           
#Site.f          22      10 300.8    0.002 ** 
#  Treatment       20       2 180.3    0.001 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Arguments:
#  Test statistics calculated assuming uncorrelated response (for faster computation) 
#P-value calculated using 999 resampling iterations via PIT-trap resampling (to account for correlation in testing).

nonE.mv2 <- manyglm(nonE.m ~ Site.f, data=tot.cover2)

rownames(nonE.m) <- paste(tot.cover2$Treatment, tot.cover2$Site)


nonE.ancova <- manyglm(nonE.m ~ Site.f + Ehrharta_erecta, data=tot.cover2)
#library(vegan)
nonE.nmds <- metaMDS(nonE.m, k=2, autotransform=F, expand=F, dist='manh')
plot(nonE.nmds, display='sites', type='t')

#> summary(anova(nonE.mv))
#Time elapsed: 0 hr 4 min 32 sec
#Length Class      Mode     
#family         1    -none-     character
#p.uni          1    -none-     character
#test           1    -none-     character
#cor.type       1    -none-     character
#resamp         1    -none-     character
#nBoot          1    -none-     numeric  
#shrink.param   3    -none-     numeric  
#n.bootsdone    1    -none-     numeric  
#table          4    data.frame list     
#uni.p        102    -none-     numeric  
#uni.test     102    -none-     numeric  

##################################################
#                                                #   
# Code from meeting with Philip on 24 March 2017 # 
#                                                #
##################################################

#This data/analysis does not include reference plots, which we want to add

#Convert site to a factor since it's right now numeric

tot.cover2$Site.f <- factor(tot.cover2$Site)

nonE.m <- as.matrix(tot.cover2[,3:36])
#Remove Ehrharta erecta from matrix
nonE.m <- nonE.m[,-7]

#To check that it's really gone
colnames(nonE.m)

#GLM model on 2014 data without Ehrharta
nonE.mv <- manyglm(nonE.m ~ Site.f + Treatment, data=tot.cover2)

str(tot.cover2)

anova(nonE.mv)

#Time elapsed: 0 hr 3 min 9 sec
#Analysis of Deviance Table

#Model: manyglm(formula = nonE.m ~ Site.f + Treatment, data = tot.cover2)

#Multivariate test:
#  Res.Df Df.diff   Dev Pr(>Dev)    
#(Intercept)     32                           
#Site.f          22      10 300.8    0.002 ** 
#Treatment       20       2 180.3    0.001 ***
#  ---
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Arguments:
#  Test statistics calculated assuming uncorrelated response (for faster computation) 
#P-value calculated using 999 resampling iterations via PIT-trap resampling (to account for correlation in testing).

nonE.mv

#Call:  manyglm(formula = nonE.m ~ Site.f + Treatment, data = tot.cover2) 
#[1] "negative.binomial"

#Nuisance Parameter(s) phi estimated by the PHI method.
#Adenalcolen_bicolor                    Caprifoliaceae_sp  
#0.007                                1.030  
#Carex_globosa                Clinopodium_douglasii  
#2.372                                1.099  
#Cotoneaster_pannosus                    Dryopteris_arguta  
#0.007                                0.006  
#Eurybia_radulina                       Fragaria_vesca  
#0.006                                0.831  
#Galium_californicum                     Galium_triflorum  
#0.007                                1.030  
#Geranium_molle                     Iris_douglasiana  
#0.006                                1.812  
#Juncus_patens                   Lonicera_hispidula  
#0.007                                0.631  
#Maianthemum_racemosa                      Melica_subulata  
#0.007                                1.859  
#Polystichum_munitum                Psuedotsuga_menziezii  
#0.006                                0.006  
#Pteridium_aquilinum_var_pubescens                    Quercus_agrifolia  
#1.030                                3.057  
#Rubus_ursinus                      Stachys_bullata  
#0.803                                1.293  
#Symphoricarpos_albus_var_laevigatus                Symphoricarpos_mollis  
#0.007                                1.030  
#Symphoricarpos_sp           Toxicodendron_diversilobum  
#0.006                                1.030  
#Umbellularia_californica                                Unk_1  
#0.006                                0.006  
#Unk_2                                Unk_3  
#0.006                                0.007  
#Unk_Aster_1                          Unk_Aster_2  
#0.006                                0.006  
#Viola_ocellata  
#0.006  

#Degrees of Freedom: 32 Total (i.e. Null); 20 Residual

#Adenalcolen_bicolor  Caprifoliaceae_sp  Carex_globosa
#2*log-likelihood:     -4.992                0.000            -17.438     
#Residual Deviance:     0.000                0.000              3.983     
#AIC:                  32.992               28.000             45.438     
#Clinopodium_douglasii  Cotoneaster_pannosus  Dryopteris_arguta
#2*log-likelihood:    -50.492                -11.771                -3.266         
#Residual Deviance:    11.417                  2.306                 0.001         
#AIC:                  78.492                 39.771                31.266         
#Eurybia_radulina  Fragaria_vesca  Galium_californicum
#2*log-likelihood:     -2.001           -45.285          -5.266           
#Residual Deviance:     0.001             9.616           0.001           
#AIC:                  30.001            73.285          33.266           
#Galium_triflorum  Geranium_molle  Iris_douglasiana  Juncus_patens
#2*log-likelihood:      0.000            -2.001         -15.800            -4.054     
#Residual Deviance:     0.000             0.001           4.140             0.001     
#AIC:                  28.000            30.001          43.800            32.054     
#Lonicera_hispidula  Maianthemum_racemosa  Melica_subulata
#2*log-likelihood:    -47.699              -2.992               -15.947       
#Residual Deviance:     8.641               0.001                 4.118       
#AIC:                  75.699              30.992                43.947       
#Polystichum_munitum  Psuedotsuga_menziezii
#2*log-likelihood:     -3.481               -8.896             
#Residual Deviance:     0.001                0.001             
#AIC:                  31.481               36.896             
#Pteridium_aquilinum_var_pubescens  Quercus_agrifolia  Rubus_ursinus
#2*log-likelihood:      0.000                            -24.077           -161.023     
#Residual Deviance:     0.000                              6.864             28.940     
#AIC:                  28.000                             52.077            189.023     
#Stachys_bullata  Symphoricarpos_albus_var_laevigatus
#2*log-likelihood:   -111.479           -7.408                           
#Residual Deviance:    27.319            0.000                           
#AIC:                 139.479           35.408                           
#Symphoricarpos_mollis  Symphoricarpos_sp  Toxicodendron_diversilobum
#2*log-likelihood:      0.000                 -6.930              0.000                  
#Residual Deviance:     0.000                  0.000              0.000                  
#AIC:                  28.000                 34.930             28.000                  
#Umbellularia_californica  Unk_1     Unk_2     Unk_3     Unk_Aster_1
#2*log-likelihood:     -8.012                    -4.614    -9.687    -4.158    -2.614   
#Residual Deviance:     0.001                     0.000     0.000     0.001     0.001   
#AIC:                  36.012                    32.614    37.687    32.158    30.614   
#Unk_Aster_2  Viola_ocellata
#2*log-likelihood:     -2.614       -2.614      
#Residual Deviance:     0.001        0.001      
#AIC:                  30.614       30.614      

#nonE.mv has treatment as a variable
#nonE.mv <- manyglm(nonE.m ~ Site.f + Treatment, data=tot.cover2)

#Does not have treatment as a variable
nonE.mv2 <- manyglm(nonE.m ~ Site.f, data=tot.cover2)

#Compare AIC values between a the null and non-null model
diffAIC <- AIC(nonE.mv)-AIC(nonE.mv2)
#Big negative values mean that few spp have a linear response to EHR

names(diffAIC) <- colnames(nonE.m)

diffAIC

#Adenalcolen_bicolor                   Caprifoliaceae_sp 
#-0.2350582                           4.0000000 
#Carex_globosa               Clinopodium_douglasii 
#1.4798722                          -5.5842145 
#Cotoneaster_pannosus                   Dryopteris_arguta 
#-5.4368243                          -1.5804368 
#Eurybia_radulina                      Fragaria_vesca 
#1.8032497                          -7.9739120 
#Galium_californicum                    Galium_triflorum 
#-4.5210334                           4.0000000 
#Geranium_molle                    Iris_douglasiana 
#1.8032497                           1.4786535 
#Juncus_patens                  Lonicera_hispidula 
#-2.9692466                          -5.6489893 
#Maianthemum_racemosa                     Melica_subulata 
#-1.0169457                           1.4696444 
#Polystichum_munitum               Psuedotsuga_menziezii 
#-1.9866713                         -11.1744848 
#Pteridium_aquilinum_var_pubescens                   Quercus_agrifolia 
#4.0000000                           3.7276823 
#Rubus_ursinus                     Stachys_bullata 
#-11.3991356                          -1.1179657 
#Symphoricarpos_albus_var_laevigatus               Symphoricarpos_mollis 
#-4.8709285                           4.0000000 
#Symphoricarpos_sp          Toxicodendron_diversilobum 
#-3.6098574                           4.0000000 
#Umbellularia_californica                               Unk_1 
#-9.7530656                           1.2276739 
#Unk_2                               Unk_3 
#1.0190781                          -3.1353253 
#Unk_Aster_1                         Unk_Aster_2 
#-0.1051964                          -0.1051964 
#Viola_ocellata 
#-0.1051964 

##############
#   Ancova   #
##############

#Way of doing a linear model to see how abundance of EHR affects native species

nonE.ancova <- manyglm(nonE.m ~ Site.f + Ehrharta_erecta, data=tot.cover2)

#Note: If you wanted EHR and treatment in the same ANCOVA we would have needed a base line, pretreatment. What we have is incomplete since we didn't survey all sites

AIC(nonE.mv) - AIC(nonE.ancova)
#[1]   1.99883354   2.00000000   1.57807472  -6.38141966  -4.71911741   2.00032255
#[7]  -0.02622818  -9.45987034  -4.60166532   2.00000000   1.99709466   3.70183315
#[13]   1.92687915  -7.59613981  -1.96053714   0.54579615   1.35172361 -12.83892472
#[19]   2.00000000   2.91588729 -10.94480199  -2.97486854  -6.71600397   2.00000000
#[25]  -1.12498617   2.00000000 -11.73803196  -0.08730097   1.94551816  -3.95736814
#[31]   1.99678950  -2.07893716   1.85851528

AIC(nonE.ancova) - AIC(nonE.mv2)
#[1] -2.23389178  2.00000000 -0.09820255  0.79720515 -0.71770688 -3.58075935  1.82947783
#[8]  1.48595836  0.08063191  2.00000000 -0.19384501 -2.22317962 -4.89612575  1.94715054
#[15]  0.94359147  0.92384829 -3.33839492  1.66443992  2.00000000  0.81179501 -0.45433359
#[22]  1.85690280  1.84507542  2.00000000 -2.48487123  2.00000000  1.98496634  1.31497483
#[29] -0.92644007  0.82204280 -2.10198587  1.97374078 -1.96371166

anova(nonE.ancova)

#Time elapsed: 0 hr 3 min 20 sec
#Analysis of Deviance Table

#Model: manyglm(formula = nonE.m ~ Site.f + Ehrharta_erecta, data = tot.cover2)

#Multivariate test:
#Res.Df Df.diff    Dev Pr(>Dev)   
#(Intercept)         32                           
#Site.f              22      10 300.81    0.004 **
#Ehrharta_erecta     21       1  60.93    0.314   
---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  Arguments:
#  Test statistics calculated assuming uncorrelated response (for faster computation)  P-value calculated using 999 resampling iterations via PIT-trap resampling (to account for correlation in testing).

