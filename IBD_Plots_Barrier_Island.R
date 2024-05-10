###########Florida Barrier Island frogs########### 
###Isolation by Distance (IBD) analyses R script
###Last modification 04/30/2024
###Luis Amador###

#Load packages
library(ape)
library(tidyverse)
library(sf)
library(stringr)
library(geosphere)
library(dartR)
library(geodist)
library(vcfR)
library(MASS)


#################################
### IBD - Anaxyrus terrestris ###
#################################

#setwd("~/Dropbox/Gen_diversity_amphians_NSF/US-amphibians/Anaxyrus_terrestris")
#setwd("/Users/Luis Amador/Dropbox/Gen_diversity_amphians_NSF/US-amphibians/Anaxyrus_terrestris/")
At_coord <- read.csv("At_coords.csv", sep = "")
coords_At_ll <- data.frame(Lon=At_coord$Longitude, Lat=At_coord$Latitude)

#Geographic distance
At_geodist <- distm(coords_At_ll)
max(At_geodist)
class(At_geodist)
A_terrestris <- At_geodist[, 1:ncol(At_geodist)]/1000
class(A_terrestris)

#Genetic distance
at_vcf <- read.vcfR(file = "Anaxyrus_terrestris_snps.vcf")
at_gl <- vcfR2genlight(at_vcf)
at_divergence <- gl.dist.ind(at_gl, method = "Euclidean", scale = TRUE)
class(at_divergence)
At_gendist <- as.matrix(at_divergence)

class(At_geodist)
class(At_gendist)
plot(A_terrestris, At_gendist)

#Replace values upper diagonal wth NAs
Ante_geod <- replace(A_terrestris, upper.tri(A_terrestris)[ ], NA)
class(Ante_geod)
Ante_geod <- as.dist(Ante_geod)
Ante_gend <- replace(At_gendist, upper.tri(At_gendist)[ ], NA)
class(Ante_gend)
Ante_gend <- as.dist(Ante_gend)
Ante_gend <- unname(Ante_gend)
plot(Ante_geod, Ante_gend)

#IBD 
At_IBD <- mantel.randtest(Ante_gend, Ante_geod, nrepet = 1000)
At_IBD
plot(At_IBD)

#Two-Dimensional Kernel Density Estimation IBD
dens <- kde2d(Ante_geod, Ante_gend, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
At_kd <- plot(Ante_geod, Ante_gend, pch=20,cex=1.2, ylab="Genetic distance", xlab="Geographic distance (km)", main=expression(italic("Anaxyrus terrestris"))) +
image(dens, col=transp(myPal(300),.7), add=TRUE)
#abline(a=coef(Ante_geod~Ante_gend), b=0, col="blue")
#text(20, 0.38, "p-value: 0.2298")
#text(20, 0.37, "r: 0.1436")


#########################
#### Hyla cinerea IBD ###
#########################

#Set working directory
#setwd("/home/luis/Dropbox/Gen_diversity_amphians_NSF/US-amphibians/Hyla_cinerea/")
#setwd("/Users/Luis Amador/Dropbox/Gen_diversity_amphians_NSF/US-amphibians/Hyla_cinerea/")

#Read coordinates
Hc_coord <- read.table("coords_Hc.txt", sep = "", header = TRUE)
#Create a small dataframe with only latitude and longitude
coords_Hc_ll <- data.frame(Lon=Hc_coord$Longitude, Lat=Hc_coord$Latitude)

##Geographic distances##
#Euclidean geographic distance
dgeo_hc <- distm(coords_Hc_ll)
class(dgeo_hc)
max(dgeo_hc)
H_cinerea <- dgeo_hc[, 1:ncol(dgeo_hc)]/1000
class(H_cinerea)

Hc_geod <- replace(dgeo_hc, upper.tri(At_geodist)[ ], NA)
class(Hc_geod)
Hc_geod <- as.dist(H_cinerea)

##Genetic distance##
#read vcf file
hc_vcf <- read.vcfR(file = "Hyla_cinerea_snps.vcf")
hc_gl <- vcfR2genlight(hc_vcf) #convert format to genlight from vcf
hc_divergence <- gl.dist.ind(hc_gl, method = "Euclidean", scale = TRUE) #Genetic distance
class(hc_divergence)
hc_divergence <- unname(hc_divergence)

#Check distances
plot(Hc_geod, hc_divergence)

#IBD
Hc_IBD <- mantel.randtest(hc_divergence, Hc_geod, nrepet = 1000)
Hc_IBD
plot(Hc_IBD)

#Two-Dimensional Kernel Density Estimation IBD
dens <- kde2d(Hc_geod, hc_divergence, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
Hc_kd <- plot(Hc_geod, hc_divergence, pch=20,cex=1.2, ylab="Genetic distance", xlab="Geographic distance (km)", main=expression(italic("Hyla cinerea"))) +
   image(dens, col=transp(myPal(300),.7), add=TRUE)
#text(20, 0.385, "p-value: 0.0010")
#text(20, 0.375, "r: 0.6473")


############################
#### Hyla squirella IBD ####
############################

#Set working directory
#setwd("/home/luis/Dropbox/Gen_diversity_amphians_NSF/US-amphibians/Hyla_squirella/")
#setwd("/home/luis/Dropbox/Gen_diversity_amphians_NSF/US-amphibians/Hyla_squirella/")
#setwd("/Users/Luis Amador/Dropbox/Gen_diversity_amphians_NSF/US-amphibians/Hyla_squirella/")

#Read coordinates
Hs_coord <- read.table("coords_Hs.txt", sep = "", header = TRUE)
#Create a small dataframe with only latitude and longitude
coords_Hs_ll <- data.frame(Lon=Hs_coord$Longitude, Lat=Hs_coord$Latitude)

##Geographic distances##
#Euclidean geographic distance
dgeo_hs <- distm(coords_Hs_ll)
class(dgeo_hs)
max(dgeo_hs)
#meters to km
H_squirella <- dgeo_hs[, 1:ncol(dgeo_hs)]/1000
class(H_squirella)
dgeo_hs <- as.dist(H_squirella)

##Genetic distance##
#read vcf file
hs_vcf <- read.vcfR(file = "Hyla_squirella_snps.vcf")
hs_gl <- vcfR2genlight(hs_vcf) #convert format to genlight from vcf
hs_divergence <- gl.dist.ind(hs_gl, method = "Euclidean", scale = TRUE) #Genetic distance
class(hs_divergence)
hs_divergence <- unname(hs_divergence)
#Check distances
plot(dgeo_hs, hs_divergence)

#IBD
Hs_IBD <- mantel.randtest(hs_divergence, dgeo_hs, nrepet = 1000)
Hs_IBD
plot(Hs_IBD)

#Two-Dimensional Kernel Density Estimation IBD
dens <- kde2d(dgeo_hs, hs_divergence, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(dgeo_hs, hs_divergence, pch=20,cex=1.2, ylab="Genetic distance", xlab="Geographic distance (km)", main=expression(italic("Hyla squirella"))) +
image(dens, col=transp(myPal(300),.7), add=TRUE)
#text(100, 0.38, "r: 0.7134")
#text(100, 0.37, "p-value: 0.0010")


################################
#### Rana sphenocephala IBD ####
################################

#setwd("/home/luis/Dropbox/Gen_diversity_amphians_NSF/US-amphibians/Rana_sphenocephala/")
#setwd("/Users/Luis Amador/Dropbox/Gen_diversity_amphians_NSF/US-amphibians/Rana_sphenocephala/")

#Read coordinates
Rs_coord <- read.table("coords_Rs.txt", sep = "", header = TRUE)
#Create a small dataframe with only latitude and longitude
coords_Rs_ll <- data.frame(Lon=Rs_coord$Longitude, Lat=Rs_coord$Latitude)

##Geographic distances##
dgeo_Rs <- distm(coords_Rs_ll)
class(dgeo_Rs)
max(dgeo_Rs)
#meters to km
R_spheno <- dgeo_Rs[, 1:ncol(dgeo_Rs)]/1000
dgeo_Rs <- as.dist(R_spheno)

##Genetic distance##
#read vcf file
Rs_vcf <- read.vcfR(file = "Rana_sphenocephala_snps.vcf")
Rs_gl <- vcfR2genlight(Rs_vcf) #convert format to genlight from vcf
Rs_divergence <- gl.dist.ind(Rs_gl,  method = "Euclidean", scale = TRUE) #Genetic distance
class(Rs_divergence)
Rs_divergence <- unname(Rs_divergence)
#Check distances
plot(dgeo_Rs, Rs_divergence)

#IBD
Rs_IBD <- mantel.randtest(Rs_divergence, dgeo_Rs, nrepet = 1000)
Rs_IBD
plot(Rs_IBD)

#Two-Dimensional Kernel Density Estimation IBD

par(mfrow = c(2,2), oma=c(1,3,1,1))

line = 1
cex = 2
side = 3
adj=-0.05

#Plotting IBD kernel density 
dens1 <- kde2d(Ante_geod, Ante_gend, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Ante_geod, Ante_gend, pch=20,cex=1.2, ylab="Genetic distance", xlab="Geographic distance (km)", main=expression(italic("Anaxyrus terrestris", scientific=FALSE))) +
image(dens1, col=transp(myPal(300),.7), add=TRUE) +
 # abline(lm(Ante_gend~Ante_geod))
text(20, 0.37, "p-value: 0.2298")
text(20, 0.38, "r: 0.1436")
mtext("A", side=side, line=line, cex=cex, adj=adj)

dens2 <- kde2d(dgeo_Rs, Rs_divergence, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(dgeo_Rs, Rs_divergence, pch=20,cex=1.2, ylab="Genetic distance", xlab="Geographic distance (km)", main=expression(italic("Rana sphenocephala"))) +
image(dens2, col=transp(myPal(300),.7), add=TRUE) +
text(95, 0.383, "r: 0.7011 ")
text(95, 0.373, "p-value: 0.0010")
mtext("B", side=side, line=line, cex=cex, adj=adj)

dens3 <- kde2d(Hc_geod, hc_divergence, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Hc_geod, hc_divergence, pch=20,cex=1.2, ylab="Genetic distance", xlab="Geographic distance (km)", main=expression(italic("Hyla cinerea"))) +
image(dens3, col=transp(myPal(300),.7), add=TRUE) +
text(20, 0.375, "p-value: 0.0010")
text(20, 0.385, "r: 0.6473")
mtext("C", side=side, line=line, cex=cex, adj=adj)

dens4 <- kde2d(dgeo_hs, hs_divergence, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(dgeo_hs, hs_divergence, pch=20,cex=1.2, ylab="Genetic distance", xlab="Geographic distance (km)", main=expression(italic("Hyla squirella"))) +
image(dens4, col=transp(myPal(300),.7), add=TRUE) +
text(100, 0.38, "r: 0.7134")
text(100, 0.37, "p-value: 0.0010")
mtext("D", side=side, line=line, cex=cex, adj=adj)

#Plotting IBD histograms
plot(At_IBD)
mtext("A", side=side, line=line, cex=cex, adj=adj)
plot(Rs_IBD)
mtext("B", side=side, line=line, cex=cex, adj=adj)
plot(Hc_IBD)
mtext("C", side=side, line=line, cex=cex, adj=adj)
plot(Hs_IBD)
mtext("D", side=side, line=line, cex=cex, adj=adj)