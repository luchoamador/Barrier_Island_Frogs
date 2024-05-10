#Florida Barrier Island frogs 
#PCA plots R script
#Last modification 05/01/2024
#Luis Amador

#Load packages
library(dartR)
library(vcfR)

#Set working directory
#setwd("~/Dropbox/Gen_diversity_amphians_NSF/US-amphibians")

#Anaxyrus terresris PCA
At_vcf <- read.vcfR("Anaxyrus_terrestris/Anaxyrus_terrestris_snps.vcf")
At_gl <- vcfR2genlight(At_vcf)
#define populations
at_gl1 <- gl.define.pop(At_gl, ind.list = c("Aterr_LNB00214", "Aterr_LNB00270", "Aterr_LNB00274", "Aterr_LNB01461"), new = "Island")
at_gl2 <- gl.define.pop(at_gl1, ind.list = c("Aterr_LNB01390", "Aterr_LNB01194", "Aterr_LNB01197", "Aterr_LNB01649", "Aterr_LNB01720"), new = "Mainland")
#check the population names
popNames(at_gl2)
at_pca <- gl.pcoa(at_gl2)
gl.pcoa.plot(at_pca, at_gl2, xaxis = 1, yaxis = 2, ellipse = TRUE, scale = TRUE)
#gl.pcoa.plot(at_pca, at_gl2, xaxis = 1, yaxis = 2, ellipse = TRUE, scale = TRUE, interactive = TRUE)
PCA1_3_At <- gl.pcoa.plot(at_pca, at_gl2, xaxis = 1, yaxis = 3, ellipse = TRUE, pop.labels = "pop")

#Hyla cinerea PCA
hc_vcf <- read.vcfR(file = "Hyla_cinerea/Hyla_cinerea_snps.vcf")
hc_gl <- vcfR2genlight(hc_vcf) #convert format to genlight from vcf
#define populations
hc_gl1 <- gl.define.pop(hc_gl, ind.list = c("Hcinerea_LNB00114", "Hcinerea_LNB00296",	"Hcinerea_LNB00113",	"Hcinerea_LNB00297"), new = "Island")
hc_gl2 <- gl.define.pop(hc_gl1, ind.list = c("Hcinerea_LNB00371", "Hcinerea_LNB00540", "Hcinerea_LNB01209",	"Hcinerea_LNB01845", "Hcinerea_LNB01846"), new = "Mainland")
popNames(hc_gl2)
hc_pca <- gl.pcoa(hc_gl2)
gl.pcoa.plot(hc_pca, hc_gl2, xaxis = 1, yaxis = 2, interactive = T, ellipse = T)
PCA1_3_Hc <- gl.pcoa.plot(hc_pca, hc_gl2, xaxis = 1, yaxis = 3, ellipse = TRUE, pop.labels = "pop")

# Rana sphenocephala PCA
rs_vcf <- read.vcfR(file = "Rana_sphenocephala/Rana_sphenocephala_snps.vcf")
rs_gl <- vcfR2genlight(rs_vcf) #convert format to genlight from vcf
#define populations
rs_gl$ind.names
rs_gl1 <- gl.define.pop(rs_gl, ind.list = c("Rspheno_LNB00233", "Rspheno_LNB00248",	"Rspheno_LNB00263",	"Rspheno_LNB00338"), new = "Island")
rs_gl2 <- gl.define.pop(rs_gl1, ind.list = c("Rspheno_LNB01247", "Rspheno_LNB00146", "Rspheno_LNB00147",	"Rspheno_LNB01836", "Rspheno_LNB01837"), new = "Mainland")
popNames(rs_gl2)
rs_pca <- gl.pcoa(rs_gl2)
gl.pcoa.plot(rs_pca, rs_gl2, xaxis = 1, yaxis = 2, interactive = TRUE, ellipse = TRUE)
PCA1_3_Rs <- gl.pcoa.plot(rs_pca, rs_gl2, xaxis = 1, yaxis = 3, ellipse = TRUE, pop.labels = "pop")
#gl.pcoa.plot(rs_pca, rs_gl2, xaxis = 1, yaxis = 2, zaxis = 3)

#Hyla squirella PCA
hs_vcf <- read.vcfR("Hyla_squirella/Hyla_squirella_snps.vcf")
hs_gl <- vcfR2genlight(hs_vcf)
#define populations
hs_gl$ind.names
Hs_gl1 <- gl.define.pop(hs_gl, ind.list = c("Hsquirella_LNB00119", "Hsquirella_LNB00126", "Hsquirella_LNB00198", "Hsquirella_LNB00199"), new = "Island")
Hs_gl2 <- gl.define.pop(Hs_gl1, ind.list = c("Hsquirella_LNB00353", "Hsquirella_LNB00549",  "Hsquirella_LNB00550", "Hsquirella_LNB01693",
"Hsquirella_LNB01694"), new = "Mainland")
hs_pca <- gl.pcoa(hs_gl)
gl.pcoa.plot(hs_pca, Hs_gl2, xaxis = 1, yaxis = 2, ellipse = TRUE, interactive = TRUE)
PCA1_3_Hs <- gl.pcoa.plot(hs_pca, Hs_gl2, xaxis = 1, yaxis = 3, ellipse = TRUE, pop.labels = "pop")

#Create a four PCAs plot
line = 1
cex = 2
side = 3
adj=-0.05

gl.pcoa.plot(at_pca, at_gl2, xaxis = 1, yaxis = 3, ellipse = TRUE, scale = TRUE, pop.labels = "pop")
mtext("A", side=side, line=line, cex=cex, adj=adj)
PCA1_3_Rs
class(PCA1_3_At)
mtext("B", side=side, line=line, cex=cex, adj=adj)
PCA1_3_Hc
mtext("C", side=side, line=line, cex=cex, adj=adj)
PCA1_3_Hs
mtext("D", side=side, line=line, cex=cex, adj=adj)

