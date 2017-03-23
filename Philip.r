# Philip's code for 'non-Ehrharta' spp comp analysis

library(dplyr)
library(tidyr)
temp <- gather(oct2014, 'species', 'cover', 8:41)
tot.cover <- temp %>% group_by(Treatment, Site, species) %>% summarize(cover=sum(cover))
tot.cover2 <- spread(tot.cover, species, cover)
  
nonE.m <- as.matrix(tot.cover2[,3:36])
  
library(mvabund)

nonE.mv <- manyglm(nonE.m ~ Site + Treatment, data=tot.cover2)
anova(nonE.mv)

rownames(nonE.m) <- paste(tot.cover2$Treatment, tot.cover2$Site)

library(vegan)
nonE.nmds <- metaMDS(nonE.m, k=2, autotransform=F, expand=F, dist='manh')
plot(nonE.nmds, display='sites', type='t')
