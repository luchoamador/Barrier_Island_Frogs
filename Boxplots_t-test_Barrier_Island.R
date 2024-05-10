### Florida Barrier Island frogs ###
##PGenetic (Nucleotide) diversity box plots - R script
##Last modification 02/14/2024
##Luis Amador##

#setwd("Dropbox/Gen_diversity_amphians_NSF/US-amphibians/")

library(tidyverse)

#Hyla squirella
squi <- read.csv("Hyla_squirella/Hyla_squirella_phylo_summarystats.csv")
t.test(squi$Pi ~ squi$Island_Mainland)
boxplot(squi$Pi ~ squi$Island_Mainland, xlab = "Populations", ylab = "Nucleotide diversity", main="Hyla squirella", sub= "t(4.84) = -2.8676, p = 0.03641")

H_squibx <- squi %>%
  ggplot() +
  aes(x=Island_Mainland, y=Pi)+
  geom_boxplot(notch = FALSE, fill="gold")+labs(x="Populations", y="Nucleotide diversity")+
  scale_shape_manual(values = c(21,21))+
  theme_minimal()
H_squibx

#Rana sphenocephala
sphe <- read.csv("Rana_sphenocephala/Rana_sphenocephala_phylo_summarystats.csv")
t.test(sphe$Pi ~ sphe$Island_Main)
boxplot(sphe$Pi ~ sphe$Island_Main, xlab = "Populations", ylab = "Nucleotide diversity", main="Rana sphenocephala", sub= "t(4.94) = -3.1485, p = 0.02586")

R_sphenbx <- sphe %>%
  ggplot() +
  aes(x = Island_Main, y = Pi) + geom_boxplot(notch = FALSE, fill="gold") + labs(x = "Populations", y="Nucleotide diversity") +
  scale_shape_manual(values = c(21,21)) +
  theme_minimal()
R_sphenbx

#Anaxyrus terrestris
terr <- read.csv("Anaxyrus_terrestris/Anaxyrus_terrestris_phylo_summarystats.csv")
terr
boxplot(terr$Pi ~ terr$Island_main, xlab = "Populations", ylab = "Nucleotide diversity", main="Anaxyrus terrestris", sub="t(6.38) = -2.8787, p = 0.02625")
t.test(terr$Pi ~ terr$Island_main)

A_terrbx <- terr %>%
  ggplot() +
  aes(x = Island_main, y = Pi) + geom_boxplot(notch = FALSE, fill="gold") + labs(x = "Populations", y="Nucleotide diversity") +
  scale_shape_manual(values = c(21,21)) +
  theme_minimal()
A_terrbx


#Hyla cinerea
cine <- read.csv("Hyla_cinerea/Hyla_cinerea_phylo_summarystats.csv")
boxplot(cine$Pi ~ cine$Mainland_Island, xlab = "Populations", ylab = "Nucleotide diversity", main="Hyla cinerea", sub="t(4.86) = -4.0613, p = 0.01031")
t.test(cine$Pi ~ cine$Mainland_Island)

H_cinebx <- cine %>%
  ggplot() +
  aes(x = Mainland_Island, y = Pi) + geom_boxplot(notch = FALSE, fill="gold") + labs(x = "Populations", y="Nucleotide diversity") +
  scale_shape_manual(values = c(21,21)) +
  theme_minimal()
H_cinebx

library(ggpubr)
ggarrange(A_terrbx, R_sphenbx, H_cinebx, H_squibx, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)

