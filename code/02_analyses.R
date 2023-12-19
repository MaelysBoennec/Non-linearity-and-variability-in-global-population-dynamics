################################################################################
### Non-linearity and temporal variability are overlooked components of 
### global population dynamics
### 12/2023

#### Part 2 - Analyses of population trajectories and variability

################################################################################
# Libraries ---------------------------------------------------------------
library(ggplot2)
library(ggpattern)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(tidyr)
library(treemapify)
library(sjstats) #cv function
library(forcats)
library(ggExtra)
library(moments)
library(lme4)
library(multcomp)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

################################################################################
# Load the data -----------------------------------------------------------

final.nl <- read.csv("data/non_linear_models.csv")
nl.modelsrd <- read.csv("data/random_models.csv")
nl.modsrd.sum <- read.csv("data/random_models_summary.csv")

################################################################################
########################### NON-LINEAR TRAJECTORIES ############################
################################################################################

################################################################################
# Non linearity proportion ------------------------------------------------
final.nl <- final.nl %>% mutate(LoNL = ifelse(acceleration=="constant","linear","non linear"))
table(final.nl$LoNL)/6437
table(final.nl$shape_class)

################################################################################
# Taxonomic and biogeographic patterns ------------------------------------
### Pearson's khi 2 tests

### HABITAT TYPE
TA <- xtabs(~ (acceleration=="constant") +  System, data=final.nl)
chisq.test(TA, correct=FALSE) # Significant

# For each habitat type, we test if the proportion of non-linearity within the
# habitat is significantly different from the overall proportion (44.8%, p=0.448)
all.habitat.tests <- final.nl %>% 
  dplyr::group_by(System) %>% 
  dplyr::do(data.frame(broom::tidy(prop.test(x=sum(.$LoNL=="non linear"),n=nrow(.),p=0.448,alternative="two.sided",correct = TRUE)))) %>% 
  dplyr::ungroup()

### TAXONOMIC GROUP
TA <- xtabs(~ (acceleration=="constant") + taxa, data=final.nl)
chisq.test(TA[,-4], correct=FALSE) # Signif (without column 4 = invertebrates bc n<5)

# For each taxonomic group, we test if the proportion of non-linearity within the
# group is significantly different from the overall proportion (44.8%, p=0.448)
all.taxa.tests <- final.nl %>% 
  filter(taxa!="Invertebrates") %>% # same : invertebrates bc n<5)
  dplyr::group_by(taxa) %>% 
  dplyr::do(data.frame(broom::tidy(prop.test(x=sum(.$LoNL=="non linear"),n=nrow(.),p=0.448,alternative="two.sided",correct = TRUE)))) %>% 
  dplyr::ungroup()

#### REGIONS
TA <- xtabs(~ (acceleration=="constant") + Region, data=final.nl)
chisq.test(TA, correct=FALSE)

# For each region, we test if the proportion of non-linearity within the
# region is significantly different from the overall proportion (44.8%, p=0.448)
all.region.tests <- final.nl %>% 
  dplyr::group_by(Region) %>% 
  dplyr::do(data.frame(broom::tidy(prop.test(x=sum(.$LoNL=="non linear"),n=nrow(.),p=0.448,alternative="two.sided",correct = TRUE)))) %>% 
  dplyr::ungroup()

### REALMS
TA <- xtabs(~ (acceleration=="constant") + full_realm, data=final.nl)
chisq.test(TA, correct=FALSE) # Significant

# For each realm, we test if the proportion of non-linearity within the
# realm is significantly different from the overall proportion (44.8%, p=0.448)
all.realm.tests <- final.nl %>% 
  dplyr::group_by(full_realm) %>% 
  dplyr::do(data.frame(broom::tidy(prop.test(x=sum(.$LoNL=="non linear"),n=nrow(.),p=0.448,alternative="two.sided",correct = TRUE)))) %>% 
  dplyr::ungroup()

### RLC
final.nl$redlistCategory<-as.factor(final.nl$redlistCategory)
levels(final.nl$redlistCategory)<-c("Critically Endangered",
                                    "Data Deficient",
                                    "Endangered",
                                    "Extinct in the Wild", 
                                    "Least Concern",
                                    "Least Concern",
                                    "Near Threatened",
                                    "Vulnerable")

redlist_data <- final.nl %>% filter(!(is.na(redlistCategory)))

TA <- xtabs(~ (acceleration=="constant") + redlistCategory, data=redlist_data)
chisq.test(TA, correct=FALSE) # Significant

# For each RL category, we test if the proportion of non-linearity within the
# category is significantly different from the overall proportion (44.8%, p=0.448)
all.rlc.tests <- redlist_data %>% 
  dplyr::group_by(redlistCategory) %>% 
  dplyr::do(data.frame(broom::tidy(prop.test(x=sum(.$LoNL=="non linear"),n=nrow(.),p=0.448,alternative="two.sided",correct = TRUE)))) %>% 
  dplyr::ungroup()

################################################################################
########################### TEMPORAL VARIABILITY ###############################
################################################################################

### Some preliminary transformations
final.nl$System<-as.factor(final.nl$System)
final.nl$taxa<-as.factor(final.nl$taxa)
final.nl$Region<-as.factor(final.nl$Region)
final.nl$full_realm<-as.factor(final.nl$full_realm)
final.nl$dir2<-as.factor(final.nl$dir2)

no.trend.all<-final.nl %>% filter(dir2=="No trend\nlinear")
inc.nl.all<-final.nl %>% filter(dir2=="Increase\nnon linear")

### Null models
null.D <- lmer(D ~ 1|Binomial, data = final.nl) 
null.CV <- lmer(CV ~ 1|Binomial, data = final.nl) 
null.MSE <- lmer(MSE ~ 1|Binomial, data = final.nl) 

null.NTD <- lmer(D ~ 1|Binomial, data = no.trend.all) 
null.NTCV <- lmer(CV ~ 1|Binomial, data = no.trend.all) 
null.NTMSE <- lmer(MSE ~ 1|Binomial, data = no.trend.all) 

################################################################################
# Temporal variability among habitat types --------------------------------

### Variability calculated with D
hab.D <- lmerTest::lmer(D ~ System + (1|Binomial), data = final.nl)
anova(hab.D,null.D) # Signif

D.hab.summary <- data.frame(summary(hab.D)$coefficients) %>% 
  mutate(system=levels(final.nl$System),
         variability=rep("D",nrow(.)))

D.hab.summary$mean <- c(D.hab.summary$Estimate[1],
                        D.hab.summary$Estimate[2]+D.hab.summary$Estimate[1],
                        D.hab.summary$Estimate[3]+D.hab.summary$Estimate[1])

D.hab.summary <- D.hab.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

### Variability calculated with CV
hab.CV <- lmerTest::lmer(CV ~ System + (1|Binomial), data = final.nl)
anova(hab.CV,null.CV) # Signif

CV.hab.summary <- data.frame(summary(hab.CV)$coefficients) %>% 
  mutate(system=levels(final.nl$System),variability=rep("CV",nrow(.)))

CV.hab.summary$mean <- c(CV.hab.summary$Estimate[1],
                         CV.hab.summary$Estimate[2]+CV.hab.summary$Estimate[1],
                         CV.hab.summary$Estimate[3]+CV.hab.summary$Estimate[1])

CV.hab.summary <- CV.hab.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)


### Variability calculated with MSE
hab.MSE <- lmerTest::lmer(MSE ~ System + (1|Binomial), data = final.nl)
anova(hab.MSE,null.MSE) # Signif

MSE.hab.summary <- data.frame(summary(hab.MSE)$coefficients) %>% 
  mutate(system=levels(final.nl$System),variability=rep("MSE",nrow(.))) 

MSE.hab.summary$mean <- c(MSE.hab.summary$Estimate[1],
                              MSE.hab.summary$Estimate[2]+MSE.hab.summary$Estimate[1],
                              MSE.hab.summary$Estimate[3]+MSE.hab.summary$Estimate[1])

MSE.hab.summary <- MSE.hab.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

full.hab.summary<-rbind(D.hab.summary,CV.hab.summary,MSE.hab.summary)

### Pairwise comparisons
mchabD<-glht(hab.D, linfct = mcp(System = "Tukey"))
mchabCV<-glht(hab.CV, linfct = mcp(System = "Tukey"))
mchabMSE<-glht(hab.MSE, linfct = mcp(System = "Tukey"))

### letters corresponding to post hoc tests
tuk.cld <- cld(mchabD) 
hab.letD <- tuk.cld$mcletters$Letters
hab.letD_df <- data.frame(System=levels(final.nl$System),letters=hab.letD)

tuk.cld <- cld(mchabCV) 
hab.letCV <- tuk.cld$mcletters$Letters
hab.letCV_df <- data.frame(System=levels(final.nl$System),letters=hab.letCV)

tuk.cld <- cld(mchabMSE) 
hab.letMSE <- tuk.cld$mcletters$Letters
hab.letMSE_df <- data.frame(System=levels(final.nl$System),letters=hab.letMSE)

################################################################################
# Temporal variability among habitat types - No trend focus ---------------

### Variability calculated with D
hab.NTD <- lmerTest::lmer(D ~ System + (1|Binomial), data = no.trend.all)
anova(hab.NTD,null.NTD) # Signif

NTD.hab.summary <- data.frame(summary(hab.NTD)$coefficients) %>% 
  mutate(system=levels(final.nl$System),
         variability=rep("D",nrow(.)))

NTD.hab.summary$mean <- c(NTD.hab.summary$Estimate[1],
                        NTD.hab.summary$Estimate[2]+NTD.hab.summary$Estimate[1],
                        NTD.hab.summary$Estimate[3]+NTD.hab.summary$Estimate[1])

NTD.hab.summary <- NTD.hab.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

### Variability calculated with CV
hab.NTCV <- lmerTest::lmer(CV ~ System + (1|Binomial), data = no.trend.all)
anova(hab.NTCV,null.NTCV) # Signif

NTCV.hab.summary <- data.frame(summary(hab.NTCV)$coefficients) %>% 
  mutate(system=levels(final.nl$System),variability=rep("CV",nrow(.)))

NTCV.hab.summary$mean <- c(NTCV.hab.summary$Estimate[1],
                         NTCV.hab.summary$Estimate[2]+NTCV.hab.summary$Estimate[1],
                         NTCV.hab.summary$Estimate[3]+NTCV.hab.summary$Estimate[1])

NTCV.hab.summary <- NTCV.hab.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

### Variability calculated with MSE
hab.NTMSE <- lmerTest::lmer(MSE ~ System + (1|Binomial), data = no.trend.all)
anova(hab.NTMSE,null.NTMSE) # Signif

NTMSE.hab.summary <- data.frame(summary(hab.NTMSE)$coefficients) %>% 
  mutate(system=levels(final.nl$System),variability=rep("MSE",nrow(.)))

NTMSE.hab.summary$mean <- c(NTMSE.hab.summary$Estimate[1],
                           NTMSE.hab.summary$Estimate[2]+NTMSE.hab.summary$Estimate[1],
                           NTMSE.hab.summary$Estimate[3]+NTMSE.hab.summary$Estimate[1])

NTMSE.hab.summary <- NTMSE.hab.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

full.NThab.summary<-rbind(NTD.hab.summary,NTCV.hab.summary,NTMSE.hab.summary)

## Pairwise comparisons
mchabNTD<-glht(hab.NTD, linfct = mcp(System = "Tukey"))
mchabNTCV<-glht(hab.NTCV, linfct = mcp(System = "Tukey"))
mchabNTMSE<-glht(hab.NTMSE, linfct = mcp(System = "Tukey"))

# letters corresponding to post hoc tests
tuk.cld <- cld(mchabNTD) 
hab.letNTD <- tuk.cld$mcletters$Letters
hab.letNTD_df <- data.frame(System=levels(final.nl$System),letters=hab.letNTD)

tuk.cld <- cld(mchabNTCV) 
hab.letNTCV <- tuk.cld$mcletters$Letters
hab.letNTCV_df <- data.frame(System=levels(final.nl$System),letters=hab.letNTCV)

tuk.cld <- cld(mchabNTMSE) 
hab.letNTMSE <- tuk.cld$mcletters$Letters
hab.letNTMSE_df <- data.frame(System=levels(final.nl$System),letters=hab.letNTMSE)

################################################################################
# Temporal variability among taxa -----------------------------------------

# Variability calculated with D
tax.Dbis <- lmerTest::lmer(D ~ taxa + (1|Binomial), data = final.nl %>% filter(taxa!="Invertebrates"))
null.Dbis <- lmer(D ~ 1|Binomial, data = final.nl %>% filter(taxa!="Invertebrates"))
anova(tax.Dbis,null.Dbis)

# tax.D <- lmerTest::lmer(D ~ taxa + (1|Binomial), data = final.nl)
# anova(tax.D,null.D) # Signif

D.tax.summary <- data.frame(summary(tax.Dbis)$coefficients) %>% 
  mutate(taxa=levels(final.nl$taxa)[-4],variability=rep("D",nrow(.))) 

D.tax.summary$mean <- c(D.tax.summary$Estimate[1],
                         D.tax.summary$Estimate[2]+D.tax.summary$Estimate[1],
                         D.tax.summary$Estimate[3]+D.tax.summary$Estimate[1],
                         D.tax.summary$Estimate[4]+D.tax.summary$Estimate[1],
                         D.tax.summary$Estimate[5]+D.tax.summary$Estimate[1],
                         D.tax.summary$Estimate[6]+D.tax.summary$Estimate[1])

D.tax.summary <- D.tax.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

# Variability calculated with CV
tax.CVbis <- lmerTest::lmer(CV ~ taxa + (1|Binomial), data = final.nl %>% filter(taxa!="Invertebrates"))
null.CVbis <- lmer(CV ~ 1|Binomial, data = final.nl %>% filter(taxa!="Invertebrates"))
anova(tax.CVbis,null.CVbis)

# tax.CV <- lmerTest::lmer(CV ~ taxa + (1|Binomial), data = final.nl)
# anova(tax.CV,null.CV) 

CV.tax.summary <- data.frame(summary(tax.CVbis)$coefficients) %>% 
  mutate(taxa=levels(final.nl$taxa)[-4],variability=rep("CV",nrow(.))) 

CV.tax.summary$mean <- c(CV.tax.summary$Estimate[1],
                          CV.tax.summary$Estimate[2]+CV.tax.summary$Estimate[1],
                          CV.tax.summary$Estimate[3]+CV.tax.summary$Estimate[1],
                          CV.tax.summary$Estimate[4]+CV.tax.summary$Estimate[1],
                          CV.tax.summary$Estimate[5]+CV.tax.summary$Estimate[1],
                          CV.tax.summary$Estimate[6]+CV.tax.summary$Estimate[1])

CV.tax.summary <- CV.tax.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

# Variability calculated with MSE
tax.MSEbis <- lmerTest::lmer(MSE ~ taxa + (1|Binomial), data = final.nl %>% filter(taxa!="Invertebrates"))
null.MSEbis <- lmer(MSE ~ 1|Binomial, data = final.nl %>% filter(taxa!="Invertebrates"))
anova(tax.MSEbis,null.MSEbis)

# tax.MSE <- lmerTest::lmer(MSE ~ taxa + (1|Binomial), data = final.nl)
# anova(tax.MSE,null.MSE) # NOT signif

MSE.tax.summary <- data.frame(summary(tax.MSEbis)$coefficients) %>% 
  mutate(taxa=levels(final.nl$taxa)[-4],variability=rep("MSE",nrow(.))) 

MSE.tax.summary$mean <- c(MSE.tax.summary$Estimate[1],
                          MSE.tax.summary$Estimate[2]+MSE.tax.summary$Estimate[1],
                          MSE.tax.summary$Estimate[3]+MSE.tax.summary$Estimate[1],
                          MSE.tax.summary$Estimate[4]+MSE.tax.summary$Estimate[1],
                          MSE.tax.summary$Estimate[5]+MSE.tax.summary$Estimate[1],
                          MSE.tax.summary$Estimate[6]+MSE.tax.summary$Estimate[1])

MSE.tax.summary <- MSE.tax.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

full.tax.summary<-rbind(D.tax.summary,CV.tax.summary,MSE.tax.summary)

## Pairwise comparisons --> without invertebrates (too low)
## Pairwise comparisons
mctaxD<-glht(tax.Dbis, linfct = mcp(taxa = "Tukey"))
mctaxCV<-glht(tax.CVbis, linfct = mcp(taxa = "Tukey"))
mctaxMSE<-glht(tax.MSEbis, linfct = mcp(taxa = "Tukey"))

# letters corresponding to post hoc tests
tuk.cld <- cld(mctaxD) 
tax.letD <- tuk.cld$mcletters$Letters
tax.letD_df <- data.frame(taxa=levels(final.nl$taxa)[-4],variability="D",letters=tax.letD)

tuk.cld <- cld(mctaxCV) 
tax.letCV <- tuk.cld$mcletters$Letters
tax.letCV_df <- data.frame(taxa=levels(final.nl$taxa)[-4],variability="CV",letters=tax.letCV)

tuk.cld <- cld(mctaxMSE) 
tax.letMSE <- tuk.cld$mcletters$Letters
tax.letMSE_df <- data.frame(taxa=levels(final.nl$taxa)[-4],variability="MSE",letters=tax.letMSE)

full.tax.letters<-rbind(tax.letD_df,tax.letCV_df,tax.letMSE_df)
full.tax.models<-left_join(full.tax.summary,full.tax.letters,by=c("taxa","variability"))

################################################################################
# Temporal variability among regions --------------------------------------

# Variability calculated with D
reg.D <- lmerTest::lmer(D ~ Region + (1|Binomial), data = final.nl)
anova(reg.D,null.D) # Signif

D.reg.summary <- data.frame(summary(reg.D)$coefficients) %>% 
  mutate(region=levels(final.nl$Region),variability=rep("D",nrow(.))) 

D.reg.summary$mean <- c(D.reg.summary$Estimate[1],
                        D.reg.summary$Estimate[2]+D.reg.summary$Estimate[1],
                        D.reg.summary$Estimate[3]+D.reg.summary$Estimate[1],
                        D.reg.summary$Estimate[4]+D.reg.summary$Estimate[1],
                        D.reg.summary$Estimate[5]+D.reg.summary$Estimate[1],
                        D.reg.summary$Estimate[6]+D.reg.summary$Estimate[1],
                        D.reg.summary$Estimate[7]+D.reg.summary$Estimate[1],
                        D.reg.summary$Estimate[8]+D.reg.summary$Estimate[1])

D.reg.summary <- D.reg.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

# standardised effect size
D.reg.sd <- final.nl %>% group_by(Region) %>%
  summarise(sd = sd(D))

D.reg.global <- cbind(D.reg.summary,D.reg.sd)
D.reg.global$std.mean <- D.reg.global$mean/D.reg.global$sd
D.reg.global$std.lower <- D.reg.global$CIinf/D.reg.global$sd
D.reg.global$std.upper <- D.reg.global$CIsup/D.reg.global$sd

# Variability calculated with CV
reg.CV <- lmerTest::lmer(CV ~ Region + (1|Binomial), data = final.nl)
anova(reg.CV,null.CV) # Signif

CV.reg.summary <- data.frame(summary(reg.CV)$coefficients) %>% 
  mutate(region=levels(final.nl$Region),variability=rep("CV",nrow(.))) 

CV.reg.summary$mean <- c(CV.reg.summary$Estimate[1],
                        CV.reg.summary$Estimate[2]+CV.reg.summary$Estimate[1],
                        CV.reg.summary$Estimate[3]+CV.reg.summary$Estimate[1],
                        CV.reg.summary$Estimate[4]+CV.reg.summary$Estimate[1],
                        CV.reg.summary$Estimate[5]+CV.reg.summary$Estimate[1],
                        CV.reg.summary$Estimate[6]+CV.reg.summary$Estimate[1],
                        CV.reg.summary$Estimate[7]+CV.reg.summary$Estimate[1],
                        CV.reg.summary$Estimate[8]+CV.reg.summary$Estimate[1])

CV.reg.summary <- CV.reg.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

# standardised effect size
CV.reg.sd <- final.nl %>% group_by(Region) %>%
  summarise(sd = sd(CV))

CV.reg.global <- cbind(CV.reg.summary,CV.reg.sd)
CV.reg.global$std.mean <- CV.reg.global$mean/CV.reg.global$sd
CV.reg.global$std.lower <- CV.reg.global$CIinf/CV.reg.global$sd
CV.reg.global$std.upper <- CV.reg.global$CIsup/CV.reg.global$sd

# Variability calculated with MSE
reg.MSE <- lmerTest::lmer(MSE ~ Region + (1|Binomial), data = final.nl)
anova(reg.MSE,null.MSE) # Signif

MSE.reg.summary <- data.frame(summary(reg.MSE)$coefficients) %>% 
  mutate(region=levels(final.nl$Region),variability=rep("MSE",nrow(.))) 

MSE.reg.summary$mean <- c(MSE.reg.summary$Estimate[1],
                         MSE.reg.summary$Estimate[2]+MSE.reg.summary$Estimate[1],
                         MSE.reg.summary$Estimate[3]+MSE.reg.summary$Estimate[1],
                         MSE.reg.summary$Estimate[4]+MSE.reg.summary$Estimate[1],
                         MSE.reg.summary$Estimate[5]+MSE.reg.summary$Estimate[1],
                         MSE.reg.summary$Estimate[6]+MSE.reg.summary$Estimate[1],
                         MSE.reg.summary$Estimate[7]+MSE.reg.summary$Estimate[1],
                         MSE.reg.summary$Estimate[8]+MSE.reg.summary$Estimate[1])

MSE.reg.summary <- MSE.reg.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

# Standardised effect size
MSE.reg.sd <- final.nl %>% group_by(Region) %>%
  summarise(sd = sd(MSE))

MSE.reg.global <- cbind(MSE.reg.summary,MSE.reg.sd)
MSE.reg.global$std.mean <- MSE.reg.global$mean/MSE.reg.global$sd
MSE.reg.global$std.lower <- MSE.reg.global$CIinf/MSE.reg.global$sd
MSE.reg.global$std.upper <- MSE.reg.global$CIsup/MSE.reg.global$sd

full.reg.global<-rbind(D.reg.global,CV.reg.global,MSE.reg.global)

## Pairwise comparisons
mcregD<-glht(reg.D, linfct = mcp(Region = "Tukey"))
mcregCV<-glht(reg.CV, linfct = mcp(Region = "Tukey"))
mcregMSE<-glht(reg.MSE, linfct = mcp(Region = "Tukey"))

# letters corresponding to post hoc tests
tuk.cld <- cld(mcregD) 
reg.letD <- tuk.cld$mcletters$Letters
reg.letD_df <- data.frame(Region=levels(final.nl$Region),letters=reg.letD)

tuk.cld <- cld(mcregCV) 
reg.letCV <- tuk.cld$mcletters$Letters
reg.letCV_df <- data.frame(Region=levels(final.nl$Region),letters=reg.letCV)

tuk.cld <- cld(mcregMSE) 
reg.letMSE <- tuk.cld$mcletters$Letters
reg.letMSE_df <- data.frame(Region=levels(final.nl$Region),letters=reg.letMSE)

################################################################################
# Temporal variability among realms ---------------------------------------

# Variability calculated with D
realm.D <- lmerTest::lmer(D ~ full_realm + (1|Binomial), data = final.nl)
anova(realm.D,null.D) # Signif

D.realm.summary <- data.frame(summary(realm.D)$coefficients) %>% 
  mutate(realm=levels(final.nl$full_realm),variability=rep("D",nrow(.))) 

D.realm.summary$mean <- c(D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[2]+D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[3]+D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[4]+D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[5]+D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[6]+D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[7]+D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[8]+D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[9]+D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[10]+D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[11]+D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[12]+D.realm.summary$Estimate[1],
                        D.realm.summary$Estimate[13]+D.realm.summary$Estimate[1])

D.realm.summary <- D.realm.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

# Variability calculated with CV
realm.CV <- lmerTest::lmer(CV ~ full_realm + (1|Binomial), data = final.nl)
anova(realm.CV,null.CV) # Signif

CV.realm.summary <- data.frame(summary(realm.CV)$coefficients) %>% 
  mutate(realm=levels(final.nl$full_realm),variability=rep("CV",nrow(.))) 

CV.realm.summary$mean <- c(CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[2]+CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[3]+CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[4]+CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[5]+CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[6]+CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[7]+CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[8]+CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[9]+CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[10]+CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[11]+CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[12]+CV.realm.summary$Estimate[1],
                          CV.realm.summary$Estimate[13]+CV.realm.summary$Estimate[1])

CV.realm.summary <- CV.realm.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

# Variability calculated with MSE
realm.MSE <- lmerTest::lmer(MSE ~ full_realm + (1|Binomial), data = final.nl)
anova(realm.MSE,null.MSE) # Signif

MSE.realm.summary <- data.frame(summary(realm.MSE)$coefficients) %>% 
  mutate(realm=levels(final.nl$full_realm),variability=rep("MSE",nrow(.))) 

MSE.realm.summary$mean <- c(MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[2]+MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[3]+MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[4]+MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[5]+MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[6]+MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[7]+MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[8]+MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[9]+MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[10]+MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[11]+MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[12]+MSE.realm.summary$Estimate[1],
                           MSE.realm.summary$Estimate[13]+MSE.realm.summary$Estimate[1])

MSE.realm.summary <- MSE.realm.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

full.realm.summary<-rbind(D.realm.summary,CV.realm.summary,MSE.realm.summary)

## Pairwise comparisons
mcrealmD<-glht(realm.D, linfct = mcp(full_realm = "Tukey"))
mcrealmCV<-glht(realm.CV, linfct = mcp(full_realm = "Tukey"))
mcrealmMSE<-glht(realm.MSE, linfct = mcp(full_realm = "Tukey"))

# letters corresponding to post hoc tests
tuk.cld <- cld(mcrealmD) 
realm.letD <- tuk.cld$mcletters$Letters
realm.letD_df <- data.frame(realm=levels(final.nl$full_realm),letters=realm.letD)

tuk.cld <- cld(mcrealmCV) 
realm.letCV <- tuk.cld$mcletters$Letters
realm.letCV_df <- data.frame(realm=levels(final.nl$full_realm),letters=realm.letCV)

tuk.cld <- cld(mcrealmMSE) 
realm.letMSE <- tuk.cld$mcletters$Letters
realm.letMSE_df <- data.frame(realm=levels(final.nl$full_realm),letters=realm.letMSE)

################################################################################
# Temporal variability among RLC ------------------------------------------

data.rlc <- final.nl %>% 
  filter(!is.na(redlistCategory)) %>% 
  mutate(redlistCategory=fct_relevel(redlistCategory,c("Extinct in the Wild","Critically Endangered","Endangered",
                                 "Vulnerable","Near Threatened","Least Concern","Data Deficient"))) %>% 
  filter(redlistCategory!="Extinct in the Wild")

null.rlc.D <- lmerTest::lmer(D ~ 1|Binomial, data = data.rlc)
null.rlc.CV <- lmerTest::lmer(CV ~ 1|Binomial, data = data.rlc)
null.rlc.MSE <- lmerTest::lmer(MSE ~ 1|Binomial, data = data.rlc)

# Variability calculated with D
rlc.D <- lmerTest::lmer(D ~ redlistCategory + (1|Binomial), data = data.rlc)
anova(rlc.D,null.rlc.D) # Signif

D.rlc.summary <- data.frame(summary(rlc.D)$coefficients) %>% 
  mutate(redlistCategory=c("Critically Endangered","Endangered",
                           "Vulnerable","Near Threatened","Least Concern","Data Deficient"),
         variability=rep("D",nrow(.))) 

D.rlc.summary$mean <- c(D.rlc.summary$Estimate[1],
                        D.rlc.summary$Estimate[2]+D.rlc.summary$Estimate[1],
                        D.rlc.summary$Estimate[3]+D.rlc.summary$Estimate[1],
                        D.rlc.summary$Estimate[4]+D.rlc.summary$Estimate[1],
                        D.rlc.summary$Estimate[5]+D.rlc.summary$Estimate[1],
                        D.rlc.summary$Estimate[6]+D.rlc.summary$Estimate[1])

D.rlc.summary <- D.rlc.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

# Variability calculated with CV
rlc.CV <- lmerTest::lmer(CV ~ redlistCategory + (1|Binomial), data = data.rlc)
anova(rlc.CV,null.rlc.CV) # Signif

CV.rlc.summary <- data.frame(summary(rlc.CV)$coefficients) %>% 
  mutate(redlistCategory=c("Critically Endangered","Endangered",
                           "Vulnerable","Near Threatened","Least Concern","Data Deficient"),
         variability=rep("CV",nrow(.))) 

CV.rlc.summary$mean <- c(CV.rlc.summary$Estimate[1],
                        CV.rlc.summary$Estimate[2]+CV.rlc.summary$Estimate[1],
                        CV.rlc.summary$Estimate[3]+CV.rlc.summary$Estimate[1],
                        CV.rlc.summary$Estimate[4]+CV.rlc.summary$Estimate[1],
                        CV.rlc.summary$Estimate[5]+CV.rlc.summary$Estimate[1],
                        CV.rlc.summary$Estimate[6]+CV.rlc.summary$Estimate[1])

CV.rlc.summary <- CV.rlc.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

# Variability calculated with MSE
rlc.MSE <- lmerTest::lmer(MSE ~ redlistCategory + (1|Binomial), data = data.rlc)
anova(rlc.MSE,null.rlc.MSE) # Signif

MSE.rlc.summary <- data.frame(summary(rlc.MSE)$coefficients) %>% 
  mutate(redlistCategory=c("Critically Endangered","Endangered",
                           "Vulnerable","Near Threatened","Least Concern","Data Deficient"),
         variability=rep("MSE",nrow(.))) 

MSE.rlc.summary$mean <- c(MSE.rlc.summary$Estimate[1],
                        MSE.rlc.summary$Estimate[2]+MSE.rlc.summary$Estimate[1],
                        MSE.rlc.summary$Estimate[3]+MSE.rlc.summary$Estimate[1],
                        MSE.rlc.summary$Estimate[4]+MSE.rlc.summary$Estimate[1],
                        MSE.rlc.summary$Estimate[5]+MSE.rlc.summary$Estimate[1],
                        MSE.rlc.summary$Estimate[6]+MSE.rlc.summary$Estimate[1])

MSE.rlc.summary <- MSE.rlc.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

full.rlc.summary<-rbind(D.rlc.summary,CV.rlc.summary,MSE.rlc.summary) %>% 
  mutate(redlistCategory=fct_relevel(redlistCategory,c("Critically Endangered","Endangered",
                                                       "Vulnerable","Near Threatened","Least Concern","Data Deficient")))

## Pairwise comparisons
mcrlcD<-glht(rlc.D, linfct = mcp(redlistCategory = "Tukey"))
mcrlcCV<-glht(rlc.CV, linfct = mcp(redlistCategory = "Tukey"))
mcrlcMSE<-glht(rlc.MSE, linfct = mcp(redlistCategory = "Tukey"))

# letters corresponding to post hoc tests
tuk.cld <- cld(mcrlcD) 
rlc.letD <- tuk.cld$mcletters$Letters
rlc.letD_df <- data.frame(redlistCategory=c("Critically Endangered","Endangered","Vulnerable","Near Threatened","Least Concern","Data Deficient"),
                          variability="D",letters=rlc.letD)

tuk.cld <- cld(mcrlcCV) 
rlc.letCV <- tuk.cld$mcletters$Letters
rlc.letCV_df <- data.frame(redlistCategory=c("Critically Endangered","Endangered","Vulnerable","Near Threatened","Least Concern","Data Deficient"),
                           variability="CV",letters=rlc.letCV)

tuk.cld <- cld(mcrlcMSE) 
rlc.letMSE <- tuk.cld$mcletters$Letters

rlc.letMSE_df <- data.frame(redlistCategory=c("Critically Endangered","Endangered","Vulnerable","Near Threatened","Least Concern","Data Deficient"),
                            variability="MSE",letters=rlc.letMSE)

full.rlc.letters<-rbind(rlc.letD_df,rlc.letCV_df,rlc.letMSE_df)
full.rlc.models<-left_join(full.rlc.summary,full.rlc.letters,by=c("redlistCategory","variability"))


################################################################################
# Temporal variability among trajectories ---------------------------------
library(TeachingDemos)

# Variability calculated with D
shape.D <- lmerTest::lmer(D ~ dir2 + (1|Binomial), data = final.nl)
anova(shape.D,null.D) # Signif

lnl.D <- lmerTest::lmer(D ~ LoNL + (1|Binomial), data = final.nl)
anova(lnl.D,null.D) # Signif
summary(lnl.D)

D.shape.summary <- data.frame(summary(shape.D)$coefficients) %>% 
  mutate(dir2=levels(final.nl$dir2),variability=rep("D",nrow(.))) 

D.shape.summary$mean <- c(D.shape.summary$Estimate[1],
                          D.shape.summary$Estimate[2]+D.shape.summary$Estimate[1],
                          D.shape.summary$Estimate[3]+D.shape.summary$Estimate[1],
                          D.shape.summary$Estimate[4]+D.shape.summary$Estimate[1],
                          D.shape.summary$Estimate[5]+D.shape.summary$Estimate[1],
                          D.shape.summary$Estimate[6]+D.shape.summary$Estimate[1])

D.shape.summary <- D.shape.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

# Variability calculated with CV
shape.CV <- lmerTest::lmer(CV ~ dir2 + (1|Binomial), data = final.nl)
anova(shape.CV,null.CV) # Signif

CV.shape.summary <- data.frame(summary(shape.CV)$coefficients) %>% 
  mutate(dir2=levels(final.nl$dir2),variability=rep("CV",nrow(.))) 

CV.shape.summary$mean <- c(CV.shape.summary$Estimate[1],
                          CV.shape.summary$Estimate[2]+CV.shape.summary$Estimate[1],
                          CV.shape.summary$Estimate[3]+CV.shape.summary$Estimate[1],
                          CV.shape.summary$Estimate[4]+CV.shape.summary$Estimate[1],
                          CV.shape.summary$Estimate[5]+CV.shape.summary$Estimate[1],
                          CV.shape.summary$Estimate[6]+CV.shape.summary$Estimate[1])

CV.shape.summary <- CV.shape.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)


# Variability calculated with MSE
shape.MSE <- lmerTest::lmer(MSE ~ dir2 + (1|Binomial), data = final.nl)
anova(shape.MSE,null.MSE) # Signif

MSE.shape.summary <- data.frame(summary(shape.MSE)$coefficients) %>% 
  mutate(dir2=levels(final.nl$dir2),variability=rep("MSE",nrow(.))) 

MSE.shape.summary$mean <- c(MSE.shape.summary$Estimate[1],
                            MSE.shape.summary$Estimate[2]+MSE.shape.summary$Estimate[1],
                            MSE.shape.summary$Estimate[3]+MSE.shape.summary$Estimate[1],
                            MSE.shape.summary$Estimate[4]+MSE.shape.summary$Estimate[1],
                            MSE.shape.summary$Estimate[5]+MSE.shape.summary$Estimate[1],
                            MSE.shape.summary$Estimate[6]+MSE.shape.summary$Estimate[1])

MSE.shape.summary <- MSE.shape.summary %>% 
  mutate(CIinf=mean-Std..Error,CIsup=mean+Std..Error)

full.shape.summary<-rbind(D.shape.summary,CV.shape.summary,MSE.shape.summary)

## Pairwise comparisons
mcshapeD<-glht(shape.D, linfct = mcp(dir2 = "Tukey"))
mcshapeCV<-glht(shape.CV, linfct = mcp(dir2 = "Tukey"))
mcshapeMSE<-glht(shape.MSE, linfct = mcp(dir2 = "Tukey"))

# letters corresponding to post hoc tests
tuk.cld <- cld(mcshapeD) 
shape.letD <- tuk.cld$mcletters$Letters
shape.letD_df <- data.frame(dir2=levels(final.nl$dir2),letters=shape.letD)

tuk.cld <- cld(mcshapeCV) 
shape.letCV <- tuk.cld$mcletters$Letters
shape.letCV_df <- data.frame(dir2=levels(final.nl$dir2),letters=shape.letCV)

tuk.cld <- cld(mcshapeMSE) 
shape.letMSE <- tuk.cld$mcletters$Letters
shape.letMSE_df <- data.frame(dir2=levels(final.nl$dir2),letters=shape.letMSE)

################################################################################
# shape.ttests <- final.nl %>% 
#   dplyr::group_by(dir2) %>% 
#   dplyr::do(data.frame(broom::tidy(z.test(x=.$D,mu=mean(final.nl$D),stdev=sd(.$D),alternative="two.sided")))) %>% 
#   dplyr::ungroup()

#### Modele complet
# modPop <- glm(MSE ~ abs(Latitude) + System + shape_class, data=final.nl)
# modPop <- lmerTest::lmer(MSE ~ abs(Latitude) + System + dir2 + (1|Binomial), data = final.nl)
# summary(modPop)$coefficients[-1,]
