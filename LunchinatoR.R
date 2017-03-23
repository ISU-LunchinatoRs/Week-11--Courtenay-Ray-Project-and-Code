############

#Libraries
library(MASS)
library(lme4)
library(ggplot2)
library(dplyr)

##############################
# GLMM of just EHR data 2014 #
##############################

#go back one folder to get wd to the folder Ehrharta Project
setwd("../")

setwd("Main plots/Data/Tidy")

#### Analyzing October 2014 collection

oct2014<-read.csv("oct2014_tidy_no_refs.csv", header=T)

#Can use glmer.nb to do GLMM with a negative binomial family

oct2014model<-glmer.nb(Ehrharta_erecta~Treatment + (1|Plot) + (1|Site), data=oct2014)

#glht #lsmeans #or releveling 

summary(oct2014model)

#Results: Herbicide and Pull are both significantly different from the control

#### Analyzing All_years (collections from Dec 2012, Jan/Feb 2013, & OCt 2014)

#go back one folder to get wd to the folder Ehrharta Project
setwd("../")

setwd("Main plots/Data/Tidy")


All_years<-read.csv("All_Years_tidy.csv", header=T)

#####
# Need to add 3 columns for cover native, non-native and unknown
######

#Make new dataset called APlot
APlot<-All_years


#Add column that will be the sum Native species cover

APlot <- mutate(APlot,Natives=c(Rubus_ursinus + Fragaria_vesca + Stachys_bullata + Clinopodium_douglasii + Symphoricarpos_albus_var_laevigatus + Symphoricarpos_sp + Carex_sp + Carex_globosa + Iris_douglasiana + Juncus_patens + Juncus_sp + Sanicula_crassicaulis + Pteridium_aquilinum_var_pubescens + Polystichum_munitum + Dryopteris_arguta + Quercus_agrifolia + Quercus_sp + Umbellularia_californica + Cardamine_oligosperma + Galium_californicum + Unk_conifer + Symphoricarpos_mollis + Lonicera_sp + Claytonia_perfoliata + Psuedotsuga_menziezii + Polygala_californica + Oxalis_oregana + Galium_aparine + Caprifoliaceae_sp + Lonicera_hispidula + Toxicodendron_diversilobum + Adenalcolen_bicolor + Viola_ocellata + Maianthemum_racemosa + Melica_subulata + Eurybia_radulina + Galium_triflorum))

#Add column that will be the sum Non-native species cover                                

APlot <- mutate(APlot,Nonnatives=c(Ehrharta_erecta + Geranium_molle + Geranium_sp + Cotoneaster_pannosus + Cotoneaster_sp + Myosotis_latifolia + Melissa_officinalis + Euphorbia_peplus + Oxalis_pes_capre + Galium_parisiense))

#Add column that will be the sum Unk species cover

APlot <- mutate(APlot,Unknowns=c(Unk_A + Unk_B + Poaceae_sp + Unk_C + Unk_D + Galium_sp + Unk_E + Unk_F + Unk_G + Unk_H + Unk_I + Unk_J + Unk_Aster3 + Unk_K + Unk_L + Lemon.Cots + Lemon.leafed.cot + Unknown.Rd.Cot.See.Pic + X.Lima.bean.shaped.Cot. + Unk_M + Unk_N + Unk_O + Unk_Aster_1 + Unk_1 + Unk_2 + Unk_Aster_2 + Unk_3))


##Graphing using ggplot

############################################
#box plots of Native species perecent cover#
############################################

APlot$Collection_Date <- factor(APlot$Collection_Date,labels = c("Dec 2012", "Jan/Feb 2013", "Oct 2014"))

#Make Treatment catagorical
as.character(APlot$Treatment)
#boxplot                
#Here, the y-axis label is split with the \n 


ggplot(APlot, aes(x = Collection_Date, y = Natives, fill=Treatment)) +
  geom_boxplot() + scale_x_discrete(name = "") +
  scale_y_continuous(name = "Percent cover\nNative species") + theme_bw()+geom_vline(xintercept = 1.5, linetype = "dashed") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

#save graphic
setwd("../../")
setwd("Graphics")

ggsave("Main_Percent_Native.png")

#################################
#Box plots of EHR perecent cover#
#################################

###convert Month/Year into a labelled factor in order to use it as a grouping variable.

All_years$Collection_Date <- factor(All_years$Collection_Date,labels = c("Dec 2012", "Jan/Feb 2013", "Oct 2014"))


#Make Treatment catagorical
as.character(All_years$Treatment)
#boxplot                
#Here, the y-axis label is split with the \n 


ggplot(All_years, aes(x = Collection_Date, y = Ehrharta_erecta, fill=Treatment)) +
  geom_boxplot() + scale_x_discrete(name = "") +
  scale_y_continuous(name = "Percent cover\nEhrharta erecta") + theme_bw()+geom_vline(xintercept = 1.5, linetype = "dashed") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 


#save graphic
setwd("../../")
setwd("Graphics")

ggsave("Main_Percent_EHR.png")

############################
# Calculating SPP Richness #
############################

#Steps
#1 Convert to long using dplyr.
#2 Convert to presence absence data
#3 Sum columns

#Want to use the data set All_Years_tidy

SR_all<-read.csv("All_Years_tidy.csv", header=T)

#Want to Trim so that my only columns are treatment, Collection Date, Plot, and Species

SR_all <- SR_all[c(2, 6:81)]

#gather(cases, "year", "n", 2:4) and average by plot

SR_gather<-gather(SR_all, key=Species, value=Percent_cover, 4:77)

SR_gather<-SR_gather%>%group_by(Collection_Date, Treatment, Plot, Species)%>% summarise(Avg_per_cov=mean(Percent_cover))

#Convert Species data with a count > 0 to 1
SR_gather$Avg_per_cov[SR_gather$Avg_per_cov >0] =1


#Sum "Avg_per_cov" which is really now presence/absence data by species

SR_gather<-SR_gather%>%group_by(Collection_Date, Treatment, Plot)%>% summarise(Spp_richness=sum(Avg_per_cov))

#View(SR_gather)

############################
#Box plots of Spp richness #
############################

###convert Month/Year into a labelled factor in order to use it as a grouping variable.

SR_gather$Collection_Date <- factor(SR_gather$Collection_Date,labels = c("Dec 2012", "Jan/Feb 2013", "Oct 2014"))

#Make Treatment catagorical
as.character(SR_gather$Treatment)
#boxplot                
#Here, the y-axis label is split with the \n 


ggplot(SR_gather, aes(x = Collection_Date, y = Spp_richness, fill=Treatment)) +
  geom_boxplot() + scale_x_discrete(name = "") +
  scale_y_continuous(name = "Species Richness") + theme_bw()+geom_vline(xintercept = 1.5, linetype = "dashed") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 


#save graphic
setwd("../../")
setwd("Graphics")

ggsave("Species_richness_main.png")



