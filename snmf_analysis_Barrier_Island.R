#Florida Barrier Island frogs 
#snmf (structure) analysis R script
#Last modification 05/01/2024
#Luis Amador

#####Load packages
library("adegenet")
library("dartR") 
library("LEA") 
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
library("vcfR")

####LEA population structure analysis#####

### Anaxyrus terrestris ###
LEA::struct2geno("Anaxyrus_terrestris_structure_noextrarows_new.str", FORMAT = 2, ploidy = 2 , extra.row = 0, extra.col = 2)
project_0 = NULL
project_0 = snmf("Anaxyrus_terrestris_structure_noextrarows_new.str.geno",
                 K = 1:5,
                 entropy = TRUE,
                 repetitions = 10,
                 project = "new")

plot(project_0, col = "blue", pch = 19, cex = 1.2, main = "Anaxyrus terrestris")
best3At = which.min(cross.entropy(project_0, K = 3)) 
best2At = which.min(cross.entropy(project_0, K = 2))

my.colors <- c("tomato","lightblue")
colores <- c("tomato", "lightblue", "springgreen3")

barchart(project_0, K = 3, run = best3At,
         border = NA, space = 0,
         col = colores,
         xlab = "Individuals",
         ylab = "Ancestry proportions", main = "Anaxyrus terrestris") -> b3At
axis(1, at = 1:length(b3At$order), labels = b3At$order, las=1, cex.axis = .3)

barchart(project_0, K = 2, run = best2At,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions", main = "Anaxyrus terrestris") -> b2At
axis(1, at = 1:length(b2At$order), labels = b2At$order, las=1, cex.axis = .3)


### Hyla cinerea ###
LEA::struct2geno("Hyla_cinerea_structure_noextrarows_new.str", FORMAT = 2, ploidy = 2 , extra.row = 0, extra.col = 2)
project_2 = NULL
project_2 = snmf("Hyla_cinerea_structure_noextrarows_new.str.geno",
                 K = 1:5,
                 entropy = TRUE,
                 repetitions = 10,
                 project = "new")

plot(project_2, col = "blue", pch = 19, cex = 1.2, main = "Hyla cinerea")
best3 = which.min(cross.entropy(project_2, K = 3)) 
best2 = which.min(cross.entropy(project_2, K = 2))

barchart(project_2, K = 3, run = best3,
         border = NA, space = 0,
         col = colores,
         xlab = "Individuals",
         ylab = "Ancestry proportions", main = "Hyla cinerea") -> b3Hc
axis(1, at = 1:length(b3Hc$order), labels = b3Hc$order, las=1, cex.axis = .3)

barchart(project_2, K = 2, run = best2,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions", main = "Hyla cinerea") -> b2Hc
axis(1, at = 1:length(b2Hc$order), labels = b2Hc$order, las=1, cex.axis = .3)


#############Hyla squirella##############
LEA::struct2geno("Hyla_squirella_structure_noextrarows_new.str", FORMAT = 2, ploidy = 2 , extra.row = 0, extra.col = 2)
project_3 = NULL
project_3 = snmf("Hyla_squirella_structure_noextrarows_new.str.geno",
                 K = 1:5,
                 entropy = TRUE,
                 repetitions = 10,
                 project = "new")

plot(project_3, col = "blue", pch = 19, cex = 1.2, main = "Hyla squirella")
best3_Hs = which.min(cross.entropy(project_3, K = 3)) 
best2_Hs = which.min(cross.entropy(project_3, K = 2))

barchart(project_3, K = 3, run = best3_Hs,
         border = NA, space = 0,
         col = colores,
         xlab = "Individuals",
         ylab = "Ancestry proportions", main = "Hyla squirella") -> bp_3.1
axis(1, at = 1:length(bp_3.1$order), labels = bp_3.1$order, las=1, cex.axis = .3)

barchart(project_3, K = 2, run = best2_Hs,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions", main = "Hyla squirella") -> bp_3
axis(1, at = 1:length(bp_3$order), labels = bp_3$order, las=1, cex.axis = .3)


#############Rana sphenocephala#################################
LEA::struct2geno("Rana_sphenocephala_structure_noextrarows_new.str", FORMAT = 2, ploidy = 2 , extra.row = 0, extra.col = 2)
project_1 = NULL
project_1 = snmf("Rana_sphenocephala_structure_noextrarows_new.str.geno",
                 K = 1:5,
                 entropy = TRUE,
                 repetitions = 10,
                 project = "new")

plot(project_1, col = "blue", pch = 19, cex = 1.2, main = "Rana sphenocephala")
best3Rs = which.min(cross.entropy(project_1, K = 3)) 
best2Rs = which.min(cross.entropy(project_1, K = 2))

barchart(project_1, K = 3, run = best3Rs,
         border = NA, space = 0,
         col = colores,
         xlab = "Individuals",
         ylab = "Ancestry proportions", main = "Rana sphenocephala") -> best3Rs
axis(1, at = 1:length(best3Rs$order), labels = best3Rs$order, las=1, cex.axis = .3)

barchart(project_2, K = 2, run = best2,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions", main = "Hyla cinerea") -> b2Rs
axis(1, at = 1:length(b2Rs$order), labels = b2Rs$order, las=1, cex.axis = .3)
