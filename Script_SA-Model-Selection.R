#### ConvergeAnt | SA - Model Selection ####
### Mark Wright, 6 December 2017

setwd("E:/Mark/_project/Analysis/")
rm(list=ls())
dev.off()

# Libraries used
#install.packages("nlme")
#install.packages("AICcmodavg")
library("nlme")
library("AICcmodavg")

# Load Data
SAdata <- read.csv(file = "Data_SurfaceArea.csv", header = TRUE, sep = ";")



#### Compare Proxies for predicting Surface Area ####
# Olfactory SA - Formula 1
lm.OSA1.Nasal <- lm(SAdata$OSA.1 ~ SAdata$Nasal)
lm.OSA1.NasOrb <- lm(SAdata$OSA.1 ~ SAdata$Nas.Orb)
lm.OSA1.NasFM <- lm(SAdata$OSA.1 ~ SAdata$Nas.FM)
lm.OSA1.PreFM <- lm(SAdata$OSA.1 ~ SAdata$PM.FM)

AIC.F1.OSA <- aictab(list(lm.OSA1.Nasal, lm.OSA1.NasOrb, lm.OSA1.NasFM, lm.OSA1.PreFM), modnames = c("Nasal", "NasOrb", "NasFM", "PreFM"))

#### Respiratory SA - Formula 1 ####
lm.RSA1.Nasal <- lm(SAdata$RSA.1 ~ SAdata$Nasal)
lm.RSA1.NasOrb <- lm(SAdata$RSA.1 ~ SAdata$Nas.Orb)
lm.RSA1.NasFM <- lm(SAdata$RSA.1 ~ SAdata$Nas.FM)
lm.RSA1.PreFM <- lm(SAdata$RSA.1 ~ SAdata$PM.FM)

AIC.F1.RSA <- aictab(list(lm.RSA1.Nasal, lm.RSA1.NasOrb, lm.RSA1.NasFM, lm.RSA1.PreFM), modnames = c("Nasal", "NasOrb", "NasFM", "PreFM"))



#### Olfactory SA - Formula 2 ####
lm.OSA2.Nasal <- lm(SAdata$OSA.2 ~ SAdata$Nasal)
lm.OSA2.NasOrb <- lm(SAdata$OSA.2 ~ SAdata$Nas.Orb)
lm.OSA2.NasFM <- lm(SAdata$OSA.2 ~ SAdata$Nas.FM)
lm.OSA2.PreFM <- lm(SAdata$OSA.2 ~ SAdata$PM.FM)

AIC.F2.OSA <- aictab(list(lm.OSA2.Nasal, lm.OSA2.NasOrb, lm.OSA2.NasFM, lm.OSA2.PreFM), modnames = c("Nasal", "NasOrb", "NasFM", "PreFM"))

#### Respiratory SA - Formula 2 ####
lm.RSA2.Nasal <- lm(SAdata$RSA.2 ~ SAdata$Nasal)
lm.RSA2.NasOrb <- lm(SAdata$RSA.2 ~ SAdata$Nas.Orb)
lm.RSA2.NasFM <- lm(SAdata$RSA.2 ~ SAdata$Nas.FM)
lm.RSA2.PreFM <- lm(SAdata$RSA.2 ~ SAdata$PM.FM)

AIC.F2.RSA <- aictab(list(lm.RSA2.Nasal, lm.RSA2.NasOrb, lm.RSA2.NasFM, lm.RSA2.PreFM), modnames = c("Nasal", "NasOrb", "NasFM", "PreFM"))



#### Olfactory SA - Formula 3 ####
lm.OSA3.Nasal <- lm(SAdata$OSA.3 ~ SAdata$Nasal)
lm.OSA3.NasOrb <- lm(SAdata$OSA.3 ~ SAdata$Nas.Orb)
lm.OSA3.NasFM <- lm(SAdata$OSA.3 ~ SAdata$Nas.FM)
lm.OSA3.PreFM <- lm(SAdata$OSA.3 ~ SAdata$PM.FM)

AIC.F3.OSA <- aictab(list(lm.OSA3.Nasal, lm.OSA3.NasOrb, lm.OSA3.NasFM, lm.OSA3.PreFM), modnames = c("Nasal", "NasOrb", "NasFM", "PreFM"))

#### Respiratory SA - Formula 3 ####
lm.RSA3.Nasal <- lm(SAdata$RSA.3 ~ SAdata$Nasal)
lm.RSA3.NasOrb <- lm(SAdata$RSA.3 ~ SAdata$Nas.Orb)
lm.RSA3.NasFM <- lm(SAdata$RSA.3 ~ SAdata$Nas.FM)
lm.RSA3.PreFM <- lm(SAdata$RSA.3 ~ SAdata$PM.FM)

AIC.F3.RSA <- aictab(list(lm.RSA3.Nasal, lm.RSA3.NasOrb, lm.RSA3.NasFM, lm.RSA3.PreFM), modnames = c("Nasal", "NasOrb", "NasFM", "PreFM"))





#### Summary Table ####
summ.surf <- data.frame(matrix(nrow = 16, ncol = 6))
colnames(summ.surf) <- c("F1.OSA", "F1.RSA", "F2.OSA", "F2.RSA", "F3.OSA", "F3.RSA")
rownames(summ.surf) <- c("Nasal.p", "Nasal.R2", "Nasal.AICc", "Nasal.wt", "NasOrb.p", "NasOrb.R2", "NasOrb.AICc", "NasOrb.wt", "NasFM.p", "NasFM.R2", "NasFM.AICc", "NasFM.wt", "PreFM.p", "PreFM.R2", "PreFM.AICc", "PreFM.wt")

summary(lm.OSA1.Nasal)$r.squared
summary(lm.OSA1.Nasal)$coefficients[2,4]
AIC.F1.OSA$AICc[AIC.F1.OSA$Modnames == "Nasal"]
AIC.F1.OSA$AICcWt[AIC.F1.OSA$Modnames == "Nasal"]

options(scipen = 999)
spec_dec <- function(x, k) trimws(format(round(x, k), nsmall=k))

### Summary Table - Formula 1
summ.surf[,1] <- c(spec_dec(summary(lm.OSA1.Nasal)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA1.Nasal)$r.squared, 4),
                   spec_dec(AIC.F1.OSA$AICc[AIC.F1.OSA$Modnames == "Nasal"], 2),
                   spec_dec(AIC.F1.OSA$AICcWt[AIC.F1.OSA$Modnames == "Nasal"], 2),
                   
                   spec_dec(summary(lm.OSA1.NasOrb)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA1.NasOrb)$r.squared, 4),
                   spec_dec(AIC.F1.OSA$AICc[AIC.F1.OSA$Modnames == "NasOrb"], 2),
                   spec_dec(AIC.F1.OSA$AICcWt[AIC.F1.OSA$Modnames == "NasOrb"], 2),
                   
                   spec_dec(summary(lm.OSA1.NasFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA1.NasFM)$r.squared, 4),
                   spec_dec(AIC.F1.OSA$AICc[AIC.F1.OSA$Modnames == "NasFM"], 2),
                   spec_dec(AIC.F1.OSA$AICcWt[AIC.F1.OSA$Modnames == "NasFM"], 2),
                   
                   spec_dec(summary(lm.OSA1.PreFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA1.PreFM)$r.squared, 4),
                   spec_dec(AIC.F1.OSA$AICc[AIC.F1.OSA$Modnames == "PreFM"], 2),
                   spec_dec(AIC.F1.OSA$AICcWt[AIC.F1.OSA$Modnames == "PreFM"], 2)
                   )

summ.surf[,2] <- c(spec_dec(summary(lm.RSA1.Nasal)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA1.Nasal)$r.squared, 4),
                   spec_dec(AIC.F1.RSA$AICc[AIC.F1.RSA$Modnames == "Nasal"], 2),
                   spec_dec(AIC.F1.RSA$AICcWt[AIC.F1.RSA$Modnames == "Nasal"], 2),
                   
                   spec_dec(summary(lm.RSA1.NasOrb)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA1.NasOrb)$r.squared, 4),
                   spec_dec(AIC.F1.RSA$AICc[AIC.F1.RSA$Modnames == "NasOrb"], 2),
                   spec_dec(AIC.F1.RSA$AICcWt[AIC.F1.RSA$Modnames == "NasOrb"], 2),
                   
                   spec_dec(summary(lm.RSA1.NasFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA1.NasFM)$r.squared, 4),
                   spec_dec(AIC.F1.RSA$AICc[AIC.F1.RSA$Modnames == "NasFM"], 2),
                   spec_dec(AIC.F1.RSA$AICcWt[AIC.F1.RSA$Modnames == "NasFM"], 2),
                   
                   spec_dec(summary(lm.RSA1.PreFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA1.PreFM)$r.squared, 4),
                   spec_dec(AIC.F1.RSA$AICc[AIC.F1.RSA$Modnames == "PreFM"], 2),
                   spec_dec(AIC.F1.RSA$AICcWt[AIC.F1.RSA$Modnames == "PreFM"], 2)
)



### Summary Table - Formula 2
summ.surf[,3] <- c(spec_dec(summary(lm.OSA2.Nasal)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA2.Nasal)$r.squared, 4),
                   spec_dec(AIC.F2.OSA$AICc[AIC.F2.OSA$Modnames == "Nasal"], 2),
                   spec_dec(AIC.F2.OSA$AICcWt[AIC.F2.OSA$Modnames == "Nasal"], 2),
                   
                   spec_dec(summary(lm.OSA2.NasOrb)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA2.NasOrb)$r.squared, 4),
                   spec_dec(AIC.F2.OSA$AICc[AIC.F2.OSA$Modnames == "NasOrb"], 2),
                   spec_dec(AIC.F2.OSA$AICcWt[AIC.F2.OSA$Modnames == "NasOrb"], 2),
                   
                   spec_dec(summary(lm.OSA2.NasFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA2.NasFM)$r.squared, 4),
                   spec_dec(AIC.F2.OSA$AICc[AIC.F2.OSA$Modnames == "NasFM"], 2),
                   spec_dec(AIC.F2.OSA$AICcWt[AIC.F2.OSA$Modnames == "NasFM"], 2),
                   
                   spec_dec(summary(lm.OSA2.PreFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA2.PreFM)$r.squared, 4),
                   spec_dec(AIC.F2.OSA$AICc[AIC.F2.OSA$Modnames == "PreFM"], 2),
                   spec_dec(AIC.F2.OSA$AICcWt[AIC.F2.OSA$Modnames == "PreFM"], 2)
)

summ.surf[,4] <- c(spec_dec(summary(lm.RSA2.Nasal)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA2.Nasal)$r.squared, 4),
                   spec_dec(AIC.F2.RSA$AICc[AIC.F2.RSA$Modnames == "Nasal"], 2),
                   spec_dec(AIC.F2.RSA$AICcWt[AIC.F2.RSA$Modnames == "Nasal"], 2),
                   
                   spec_dec(summary(lm.RSA2.NasOrb)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA2.NasOrb)$r.squared, 4),
                   spec_dec(AIC.F2.RSA$AICc[AIC.F2.RSA$Modnames == "NasOrb"], 2),
                   spec_dec(AIC.F2.RSA$AICcWt[AIC.F2.RSA$Modnames == "NasOrb"], 2),
                   
                   spec_dec(summary(lm.RSA2.NasFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA2.NasFM)$r.squared, 4),
                   spec_dec(AIC.F2.RSA$AICc[AIC.F2.RSA$Modnames == "NasFM"], 2),
                   spec_dec(AIC.F2.RSA$AICcWt[AIC.F2.RSA$Modnames == "NasFM"], 2),
                   
                   spec_dec(summary(lm.RSA2.PreFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA2.PreFM)$r.squared, 4),
                   spec_dec(AIC.F2.RSA$AICc[AIC.F2.RSA$Modnames == "PreFM"], 2),
                   spec_dec(AIC.F2.RSA$AICcWt[AIC.F2.RSA$Modnames == "PreFM"], 2)
)



### Summary Table - Formula 3
summ.surf[,5] <- c(spec_dec(summary(lm.OSA3.Nasal)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA3.Nasal)$r.squared, 4),
                   spec_dec(AIC.F3.OSA$AICc[AIC.F3.OSA$Modnames == "Nasal"], 2),
                   spec_dec(AIC.F3.OSA$AICcWt[AIC.F3.OSA$Modnames == "Nasal"], 2),
                   
                   spec_dec(summary(lm.OSA3.NasOrb)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA3.NasOrb)$r.squared, 4),
                   spec_dec(AIC.F3.OSA$AICc[AIC.F3.OSA$Modnames == "NasOrb"], 2),
                   spec_dec(AIC.F3.OSA$AICcWt[AIC.F3.OSA$Modnames == "NasOrb"], 2),
                   
                   spec_dec(summary(lm.OSA3.NasFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA3.NasFM)$r.squared, 4),
                   spec_dec(AIC.F3.OSA$AICc[AIC.F3.OSA$Modnames == "NasFM"], 2),
                   spec_dec(AIC.F3.OSA$AICcWt[AIC.F3.OSA$Modnames == "NasFM"], 2),
                   
                   spec_dec(summary(lm.OSA3.PreFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.OSA3.PreFM)$r.squared, 4),
                   spec_dec(AIC.F3.OSA$AICc[AIC.F3.OSA$Modnames == "PreFM"], 2),
                   spec_dec(AIC.F3.OSA$AICcWt[AIC.F3.OSA$Modnames == "PreFM"], 2)
)

summ.surf[,6] <- c(spec_dec(summary(lm.RSA3.Nasal)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA3.Nasal)$r.squared, 4),
                   spec_dec(AIC.F3.RSA$AICc[AIC.F3.RSA$Modnames == "Nasal"], 2),
                   spec_dec(AIC.F3.RSA$AICcWt[AIC.F3.RSA$Modnames == "Nasal"], 2),
                   
                   spec_dec(summary(lm.RSA3.NasOrb)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA3.NasOrb)$r.squared, 4),
                   spec_dec(AIC.F3.RSA$AICc[AIC.F3.RSA$Modnames == "NasOrb"], 2),
                   spec_dec(AIC.F3.RSA$AICcWt[AIC.F3.RSA$Modnames == "NasOrb"], 2),
                   
                   spec_dec(summary(lm.RSA3.NasFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA3.NasFM)$r.squared, 4),
                   spec_dec(AIC.F3.RSA$AICc[AIC.F3.RSA$Modnames == "NasFM"], 2),
                   spec_dec(AIC.F3.RSA$AICcWt[AIC.F3.RSA$Modnames == "NasFM"], 2),
                   
                   spec_dec(summary(lm.RSA3.PreFM)$coefficients[2,4], 3),
                   spec_dec(summary(lm.RSA3.PreFM)$r.squared, 4),
                   spec_dec(AIC.F3.RSA$AICc[AIC.F3.RSA$Modnames == "PreFM"], 2),
                   spec_dec(AIC.F3.RSA$AICcWt[AIC.F3.RSA$Modnames == "PreFM"], 2)
)

write.table(summ.surf, file = "Data_SA_Model-Selection.csv", sep = ";", col.names = NA)


