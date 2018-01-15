#### ConvergeAnt | PGLS ####
### Mark Wright, 13 December 2017

setwd("E:/Mark/_project/Analysis/")
rm(list=ls())
dev.off()

# Libraries used
#install.packages("ape")
#install.packages("geiger")
#install.packages("nlme")
#install.packages("phytools")
#install.packages("car")
#install.packages("AICcmodavg")
library("ape")
library("geiger")
library("nlme")
library("phytools")
library("car")
library("AICcmodavg")



#### Loading and formatting data ####
raw_data <- read.csv(file = "Data_SurfaceArea.csv", header = TRUE, sep = ";")
raw_data <- raw_data[raw_data$Ad.Juv != "Juv",]
raw_data <- raw_data[order(raw_data$Species_21),]

tree <- read.tree(file = "tree_Fred_29mitos_mafft_gapex_gblocks_reltime.tre.txt")
tree <- drop.tip(phy = tree, tip = c("Otocyon_megalotis", "Chaetophractus_vellerosus", "Zaedyus_pichiy", "Chaetophractus_villosus", "Tolypeutes_matacus", "Tolypeutes_tricinctus"))
plot(tree)

body.size <- lm(log(raw_data$Nas.FM) ~ raw_data$LogMass)
summary(body.size)
plot(log(raw_data$Nas.FM) ~ raw_data$LogMass, col = raw_data$Diet, pch = 16)
abline(a = body.size$coefficients[1], b = body.size$coefficients[2])
Anova(body.size)
plot(residuals(body.size), col = raw_data$Diet, pch = 16)
Anova(lm(residuals(body.size) ~ raw_data$Diet))


# Organize dataframe by species rather than by individuals
tdata <- aggregate(Diet~Species_21, raw_data, mean)
row.names(tdata) <- tdata$Species_21
colnames(tdata)[1] <- "Order"
tdata$Order <- raw_data$Order[duplicated(raw_data$Species_21)==FALSE]
tdata$Diet <- raw_data$Diet[duplicated(raw_data$Species_21)==FALSE]
tdata[,3] <- aggregate(Nas.FM~Species_21, raw_data, mean)[2]
#tdata[,4] <- aggregate(RSA.1~Species_21, raw_data, mean)[2]
#tdata[,5] <- aggregate(OSA.1~Species_21, raw_data, mean)[2]
tdata[,4] <- aggregate(Endos~Species_21, raw_data, mean)[2]
tdata[,5] <- aggregate(Ectos~Species_21, raw_data, mean)[2]
tdata[,6] <- aggregate(Ethmos~Species_21, raw_data, mean)[2]
tdata[,7] <- aggregate(Maxillo~Species_21, raw_data, mean)[2]
#tdata[,10] <- aggregate(cNT~Species_21, raw_data, mean)[2]
#tdata[,11] <- aggregate(rNT~Species_21, raw_data, mean)[2]
tdata[,8] <- aggregate(Endo1~Species_21, raw_data, mean)[2]
tdata[,9] <- aggregate(C_Max~Species_21, raw_data, mean)[2]
tdata[,10] <- aggregate(C_Endo1r~Species_21, raw_data, mean)[2]
tdata[,11] <- aggregate(C_Endos~Species_21, raw_data, mean)[2]
tdata[,12] <- aggregate(C_Ectos~Species_21, raw_data, mean)[2]
tdata[,13] <- aggregate(TreeTip~Species_21, raw_data, mean)[2]
tdata <- tdata[tdata$TreeTip != 16,]             # Need mass for Manis culionensis
tdata[,14] <- aggregate(LogMass~Species_21, raw_data, mean)[2]

#colnames(tdata)[c(3:5,17)] <- c("Body", "RSA", "OSA", "Tip")
colnames(tdata)[c(3,13:14)] <- c("Body", "Tip", "Mass")

raw_data <- raw_data[order(raw_data$TreeTip),]
tdata <- tdata[order(tdata$Tip),]



#### Surface Area Linear Models ####
mod_OSA <- lm(tdata$OSA ~ tdata$Body * tdata$Diet)
mod_RSA <- lm(tdata$RSA ~ tdata$Body * tdata$Diet)
mod_SA_ratio <- lm(tdata$OSA/tdata$RSA ~ tdata$Body * tdata$Diet)

mod_bro_OSA <- gls(OSA ~ Body * Diet, correlation = corBrownian(phy = tree), data = tdata)
mod_bro_RSA <- gls(RSA ~ Body * Diet, correlation = corBrownian(phy = tree), data = tdata)
mod_bro_SA_ratio <- gls(OSA/RSA ~ Body * Diet, correlation = corBrownian(phy = tree), data = tdata)

# Martins models not working
mod_mar_OSA <- gls(OSA ~ Body * Diet, correlation = corMartins(1, phy = tree), data = tdata)
mod_mar_RSA <- gls(RSA ~ Body * Diet, correlation = corMartins(1, phy = tree), data = tdata)
mod_mar_SA_ratio <- gls(OSA/RSA ~ Body * Diet, correlation = corMartins(1, phy = tree), data = tdata)

mod_pag_OSA <- gls(OSA ~ Body * Diet, correlation = corPagel(1, phy = tree), data = tdata)
mod_pag_RSA <- gls(OSA ~ Body * Diet, correlation = corPagel(1, phy = tree), data = tdata)
mod_pag_SA_ratio <- gls(OSA/RSA ~ Body * Diet, correlation = corPagel(1, phy = tree), data = tdata)


# Model selection
extractAIC(mod_OSA)[2]
summary(mod_bro_OSA)$AIC
summary(mod_pag_OSA)$AIC

extractAIC(mod_RSA)[2]
summary(mod_bro_RSA)$AIC
summary(mod_pag_RSA)$AIC

extractAIC(mod_SA_ratio)[2]
summary(mod_bro_SA_ratio)$AIC
summary(mod_pag_SA_ratio)$AIC


# Results & Plots
summary(mod_bro_OSA)
summary(mod_RSA)
summary(mod_SA_ratio)

plot(log(tdata$OSA) ~ log(tdata$Body), col = tdata$Order, pch = 16, xlab = "Body Size", ylab = "Olfactory SA")
plot(log(tdata$RSA) ~ log(tdata$Body), col = tdata$Order, pch = 16, xlab = "Body Size", ylab = "Respiratory SA")
plot(tdata$OSA/tdata$RSA ~ log(tdata$Body), col = tdata$Order, pch = 16, xlab = "Body Size", ylab = "SA Ratio (Olf/Res)")



#### Complexity Linear Models ####
mod_Max <- lm(tdata$C_Max ~ tdata$Mass * tdata$Diet)
#mod_En1 <- lm(tdata$C_Endo1r ~ tdata$Body * tdata$Diet)
mod_ENs <- lm(tdata$C_Endos ~ tdata$Mass * tdata$Diet)
mod_ECs <- lm(tdata$C_Ectos ~ tdata$Mass * tdata$Diet)

mod_bro_C_Max <- gls(C_Max ~ Mass * Diet, correlation = corBrownian(phy = tree), data = tdata)
#mod_bro_C_EN1 <- gls(C_Endo1r ~ Body * Diet, correlation = corBrownian(phy = tree), data = tdata)
mod_bro_C_ENs <- gls(C_Endos ~ Mass * Diet, correlation = corBrownian(phy = tree), data = tdata)
mod_bro_C_ECs <- gls(C_Ectos ~ Mass * Diet, correlation = corBrownian(phy = tree), data = tdata)

mod_pag_C_Max <- gls(C_Max ~ Mass * Diet, correlation = corPagel(1, phy = tree), data = tdata)
#mod_pag_C_EN1 <- gls(C_Endo1r ~ Body * Diet, correlation = corPagel(1, phy = tree), data = tdata)
mod_pag_C_ENs <- gls(C_Endos ~ Mass * Diet, correlation = corPagel(1, phy = tree), data = tdata)
mod_pag_C_ECs <- gls(C_Ectos ~ Mass * Diet, correlation = corPagel(1, phy = tree), data = tdata)


# Model selection
extractAIC(mod_Max)[2]
summary(mod_bro_C_Max)$AIC
summary(mod_pag_C_Max)$AIC

extractAIC(mod_ENs)[2]
summary(mod_bro_C_ENs)$AIC
summary(mod_pag_C_ENs)$AIC

extractAIC(mod_ECs)[2]
summary(mod_bro_C_ECs)$AIC
summary(mod_pag_C_ECs)$AIC


# Results & Plots
summary(mod_Max)
summary(mod_ENs) # *****
summary(mod_ECs)

plot(tdata$C_Max ~ log(tdata$Mass), col = tdata$Order, pch = 16, xlab = "Body Size", ylab = "Maxilloturbinal Complexity")
plot(tdata$C_Endos ~ log(tdata$Mass), col = tdata$Order, pch = 16, xlab = "Body Size", ylab = "Endoturbinals Complexity")
plot(tdata$C_Ectos ~ log(tdata$Mass), col = tdata$Order, pch = 16, xlab = "Body Size", ylab = "Ectoturbinals Complexity")

plot(tdata$C_Max ~ tdata$Diet)
plot(tdata$C_Endos ~ tdata$Diet) # *****
plot(tdata$C_Ectos ~ tdata$Diet)



#### Complexity Linear Models without Carnivora ####
tdata_SC <- tdata[tdata$Order != "Carnivora",]
tree_SC <- drop.tip(phy = tree, tip = c("Proteles_cristata", "Vulpes_vulpes", "Alopex_lagopus", "Felis_silvestris", "Hyaena_hyaena"))
plot(tree_SC)

mod_SC_ENs <- lm(tdata_SC$C_Endos ~ tdata_SC$Body * tdata_SC$Diet)
mod_SC_ECs <- lm(tdata_SC$C_Ectos ~ tdata_SC$Body * tdata_SC$Diet)

mod_SC_bro_C_ENs <- gls(C_Endos ~ Body * Diet, correlation = corBrownian(phy = tree_SC), data = tdata_SC)
mod_SC_bro_C_ECs <- gls(C_Ectos ~ Body * Diet, correlation = corBrownian(phy = tree_SC), data = tdata_SC)

mod_SC_pag_C_ENs <- gls(C_Endos ~ Body * Diet, correlation = corPagel(1, phy = tree), data = tdata_SC)
mod_SC_pag_C_ECs <- gls(C_Ectos ~ Body * Diet, correlation = corPagel(1, phy = tree), data = tdata_SC)


# Model selection
extractAIC(mod_SC_ENs)[2]
summary(mod_SC_bro_C_ENs)$AIC
summary(mod_SC_pag_C_ENs)$AIC

extractAIC(mod_SC_ECs)[2]
summary(mod_SC_bro_C_ECs)$AIC
summary(mod_SC_pag_C_ECs)$AIC


# Results & Plots
summary(mod_SC_ENs)
summary(mod_SC_ECs)

plot(tdata_SC$C_Endos ~ log(tdata_SC$Body), col = tdata_SC$Order, pch = 16, xlab = "Body Size", ylab = "Endoturbinals Complexity")
plot(tdata_SC$C_Ectos ~ log(tdata_SC$Body), col = tdata_SC$Order, pch = 16, xlab = "Body Size", ylab = "Ectoturbinals Complexity")

plot(tdata_SC$C_Endos ~ tdata_SC$Diet)
plot(tdata_SC$C_Ectos ~ tdata_SC$Diet)




#### Notes/Tests ####


















# Loading and formatting data for phylogenetic models

OSA <- tdata$OSA
RSA <- tdata$RSA
SA_ratio <- tdata$OSA/tdata$RSA
Order <- tdata$Order
Diet <-  tdata$Diet
Body <- tdata$Body
Body_adj <- tdata$Body[tdata$C_Endo1r > 0]
SA_Endos <- tdata$Endos
SA_Ectos <- tdata$Ectos
SA_Ethmos <- tdata$Ethmos
SA_Max <- tdata$Maxillo
SA_Nas <- tdata$rNT + tdata$cNT
SA_Endo1 <- tdata$Endo1
SA_Endo1_adj <- tdata$Endo1[tdata$C_Endo1r > 0]
Cmplx_2D_Max <- tdata$C_Max
Cmplx_2D_Endo1r <- tdata$C_Endo1r
Cmplx_2D_Endo1r_adj <- tdata$C_Endo1r[tdata$C_Endo1r > 0]
Cmplx_2D_Endos <- tdata$C_Endos
Cmplx_2D_Ectos <- tdata$C_Ectos



# Exploring the data visually
plot(tdata$RSA/tdata$Body ~ tdata$Order, main = "Respiratory SA adjusted for body size", ylim = c(0,140))
plot(tdata$OSA/tdata$Body ~ tdata$Order, main = "Olfactory SA adjusted for body size", ylim = c(0,140))
plot(log(SA_ratio) ~ log(Body))
plot(log(Cmplx_2D_Ectos) ~ log(Body))
plot(log(Cmplx_2D_Endos) ~ log(Body))
plot(log(Cmplx_2D_Endo1r) ~ log(Body))
plot(log(Cmplx_2D_Max) ~ log(Body))



#### Linear models ####
# Does complexity explain the variance in surface area for each turbinal?
model1 <- lm(log(SA_Ectos) ~ log(Cmplx_2D_Ectos))
summary(model1)

model2 <- lm(log(SA_Endos) ~ log(Cmplx_2D_Endos))
summary(model2)

model3 <- lm(log(SA_Max) ~ log(Cmplx_2D_Max))
summary(model3)

model4 <- lm(log(SA_Endo1_adj) ~ log(Cmplx_2D_Endo1r_adj))
summary(model4)
# Answer: yes; in all four cases, but to varying degrees



# Does adding body size to the model change the results?
model5 <- lm(log(SA_Ectos) ~ log(Cmplx_2D_Ectos)*log(Body))
summary(model5)

model6 <- lm(log(SA_Endos) ~ log(Cmplx_2D_Endos)*log(Body))
summary(model6)

model7 <- lm(log(SA_Max) ~ log(Cmplx_2D_Max)*log(Body))
summary(model7)

model8 <- lm(log(SA_Endo1_adj) ~ log(Cmplx_2D_Endo1r_adj)*log(Body_adj))
summary(model8)
# Answer: Yes, body size explains variation in Max and Endo1, but not in ethmos


Cmplx_2D_Ethmos <- (Cmplx_2D_Endos + Cmplx_2D_Ectos)/2
Cmplx_ratio <- Cmplx_2D_Ethmos/Cmplx_2D_Max
# Now with Surface Area ratios
model9 <- lm(log(SA_ratio) ~ log(Body)*Cmplx_ratio*Diet)
summary(model9)

Anova(lm(SA_ratio ~ Cmplx_ratio))


# Now take into account tree
GLS.model <- gls(SA_ratio ~ Body*Diet, correlation = corBrownian(phy = tree), data = tdata)
summary(GLS.model)

GLS.model.martins <- gls(SA_ratio ~ Body + Diet, correlation = corMartins(1, phy = tree), data = tdata)
summary(GLS.model.martins)

GLS.model.pagel <- gls(SA_ratio ~ Body*Diet, correlation = corPagel(1, phy = tree), data = tdata)
summary(GLS.model.pagel)


# Drop carnivores completely
tdata_sansCarn <- tdata[tdata$Order != "Carnivora",]

model10 <- lm(tdata_sansCarn$OSA/tdata_sansCarn$RSA ~ tdata_sansCarn$Body*tdata_sansCarn$Diet)
summary(model10)






#### Surface Area correlation to Complexity ####
SAdata.comp <- SAdata[is.na(SAdata$Endo1c) == FALSE,]
lm.OSA1.E1C <- lm(log(SAdata.comp$OSA.1) ~ SAdata.comp$Endo1c)
lm.OSA1.E1R <- lm(log(SAdata.comp$OSA.1) ~ SAdata.comp$Endo1r)
lm.OSA1.ETH <- lm(log(SAdata.comp$OSA.1) ~ SAdata.comp$Ethmo)
lm.OSA1.MAX <- lm(log(SAdata.comp$OSA.1) ~ SAdata.comp$Max)
lm.OSA1.NAS <- lm(log(SAdata.comp$OSA.1) ~ SAdata.comp$Nas)
lm.RSA1.E1C <- lm(log(SAdata.comp$RSA.1) ~ SAdata.comp$Endo1c)
lm.RSA1.E1R <- lm(log(SAdata.comp$RSA.1) ~ SAdata.comp$Endo1r)
lm.RSA1.ETH <- lm(log(SAdata.comp$RSA.1) ~ SAdata.comp$Ethmo)
lm.RSA1.MAX <- lm(log(SAdata.comp$RSA.1) ~ SAdata.comp$Max)
lm.RSA1.NAS <- lm(log(SAdata.comp$RSA.1) ~ SAdata.comp$Nas)

SAdata.pgls.comp <- SAdata.pgls
SAdata.pgls.comp[,5] <- c(1.6480, 1.2861, 1.5405, 0, 0, 0, 1.4962, 1.5878, 0, 1.4911, 1.4521, 1.5888, 0, 0, 0, 1.4141, 0, 1.6327)
SAdata.pgls.comp[,6] <- c(1.4857, 1.4117, 1.5457, 0, 0, 0, 1.5742, 1.4735, 0, 1.5016, 1.3674, 1.5687, 0, 0, 0, 1.4314, 0, 1.5742)
SAdata.pgls.comp[,7] <- c(1.6766, 1.4072, 1.3374, 0, 0, 0, 1.3484, 1.5276, 0, 1.5946, 1.2706, 1.5792, 0, 0, 0, 1.2852, 0, 1.6437)
SAdata.pgls.comp <- SAdata.pgls.comp[SAdata.pgls.comp$V5 > 0.01,]

model <- lm(log(SAdata.pgls.comp$V7) ~ log(SAdata.pgls.comp$RSA)*log(SAdata.pgls.comp$NasFM)*SAdata.pgls.comp$Diet)
summary(model)

model2 <- lm(log(SAdata.pgls.comp$V6) ~ log(SAdata.pgls.comp$OSA))
summary(model2)

model3 <- lm(log(SAdata.pgls.comp$V6) ~ log(SAdata.pgls.comp$NasFM))
summary(model3)

boxplot(residuals ~ SAdata.pgls.comp$Diet)
plot(log(SAdata.pgls.comp$V5) ~ log(SAdata.pgls.comp$NasFM), xlab = "Body Size Proxy", ylab = "Complexity")
boxplot(log(SAdata.pgls.comp$V5) ~ SAdata.pgls.comp$Diet, xlab = "Diet", ylab = "Complexity")

OSA2 <- SAdata.pgls.comp$OSA
Diet2 <-  SAdata.pgls.comp$Diet
NasFM2 <- SAdata.pgls.comp$NasFM
Comp2 <- SAdata.pgls.comp$V5
tree <- compute.brlen(tree, 1)

GLS.model <- gls(OSA2 ~ NasFM2 + Diet2 + Comp2, correlation = corBrownian(phy = tree), data = SAdata.pgls)
summary(GLS.model)

GLS.model.martins <- gls(OSA ~ NasFM + Diet, correlation = corMartins(1, phy = tree), data = SAdata.pgls)
summary(GLS.model.martins)

GLS.model.pagel <- gls(OSA ~ NasFM + Diet, correlation = corPagel(1, phy = tree), data = SAdata.pgls)
summary(GLS.model.pagel)

Anova(lm.OSA1.NasFM)
Anova(GLS.model)
Anova(GLS.model.martins)



#### Non-phylogenetic models ####
lm.OSA1.NasFM <- lm(log(SAdata.pgls$OSA) ~ log(SAdata.pgls$NasFM) + SAdata.pgls$Diet)
lm.RSA1.NasFM <- lm(SAdata$RSA.1 ~ SAdata$Nas.FM + SAdata$Diet)
summary(lm.OSA1.NasFM)
summary(lm.RSA1.NasFM)

lm.OSA1.NasFM.resid <- lm(resid(lm.OSA1.NasFM) ~ SAdata$Diet)
summary(lm.OSA1.NasFM.resid)
plot(resid(lm.OSA1.NasFM) ~ SAdata$Diet)

plot(log(SAdata.pgls$OSA) ~ log(SAdata.pgls$NasFM), xlab = "Body Size")
abline(a = lm.OSA1.NasFM$coefficients[1], b = lm.OSA1.NasFM$coefficients[2])

model <- lm(log(SAdata.pgls$OSA) ~ log(SAdata.pgls$NasFM))
residuals <- resid(model)
model2 <- lm(residuals ~ SAdata.pgls$Diet)
summary(model2)
Anova(model2)
boxplot(residuals ~ SAdata.pgls$Diet)



#### Phylogenetic models ####
GLS.model <- gls(OSA ~ NasFM + Diet, correlation = corBrownian(phy = tree), data = SAdata.pgls)
summary(GLS.model)

GLS.model.martins <- gls(OSA ~ NasFM + Diet, correlation = corMartins(1, phy = tree), data = SAdata.pgls)
summary(GLS.model.martins)

GLS.model.pagel <- gls(OSA ~ NasFM + Diet, correlation = corPagel(1, phy = tree), data = SAdata.pgls)
summary(GLS.model.pagel)

Anova(lm.OSA1.NasFM)
Anova(GLS.model)
Anova(GLS.model.martins)

extractAIC(lm.OSA1.NasFM)
extractAIC(GLS.model)
extractAIC(GLS.model.martins)

