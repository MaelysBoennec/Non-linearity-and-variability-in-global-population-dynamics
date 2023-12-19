################################################################################
### Non-linearity and temporal variability are overlooked components of 
### global population dynamics
### 12/2023

#### Part 1 - Modeling population trajectories and variability

################################################################################
# Libraries ---------------------------------------------------------------
library(ggplot2)
library(ggpattern)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(tidyr)
library(treemapify)
library(forcats)
library(ggExtra)
library(moments)
library(lme4)
library(multcomp)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
source("https://raw.githubusercontent.com/larmarange/JLutils/master/R/lm_right.R")

################################################################################
# Load the data -----------------------------------------------------------

# Get the raw data (Living Planet Database)
lpd <- readr::read_csv(
  here::here("data","LPD2022_public.csv")
) 

# IUCN Red List Categories
RLC <- readr::read_csv(
  here::here("data","IUCNredlist.csv")
) 

# Combine genus and species name to make comparable species names between 
# the two databases
RLC <- RLC %>% 
  unite(col=Binomial, genusName, speciesName, remove=FALSE) %>% 
  dplyr::select(Binomial,redlistCategory,populationTrend)

################################################################################
# Preliminary functions ---------------------------------------------------

# "fun_taxa" transforms the taxonomic Classes originally implemented in the LPD
# into more general taxonomic groups
fun_taxa <- function(L) {
  newL<-list()
  for (i in (1:length(L))){
    if (L[i]=='Aves'){newL<-append(newL,"Birds")}
    if (L[i]=='Mammalia'){newL<-append(newL,"Mammals")}
    if (L[i]=='Reptilia'){newL<-append(newL,"Reptiles")}
    if (L[i]=='Amphibia'){newL<-append(newL,"Amphibians")}
    if (L[i]=='Actinopteri'|L[i]=='Coelacanthi'|L[i]=='Petromyzonti'|L[i]=='Dipneusti'){newL<-append(newL,"Fish")}
    if (L[i]=='Elasmobranchii'|L[i]=='Holocephali'){newL<-append(newL,"Sharks_Rays")}
    if (L[i]=='Myxini'){newL<-append(newL,"Invertebrates")}
  }
  return(as.character(unlist(newL)))
}

# d calculates the consecutive disparity index
# Code from Fernandez-Martinez et al. (2018)
d <- function(ts){
  k <- 0.01 * mean(ts, na.rm=TRUE)
  n <- length(ts)
  ret <- c()
  for(i in 1:(length(ts)-1)){
    ret[i] <- abs(log((ts[i+1] + k) / (ts[i]+k)))
  }
  sum(ret)/(n-1)
}

# "class.trajectory" classifies time series into one of the nine following classes:
# linear increase, accelerated increase, decelerated increase, linear decrease, 
# accelerated decrease, decelarated decrease, linear no trend, concave, convex trajectories
# Code adapted from Rigal et al., 2020
class.trajectory <- function (Y = NULL, X = NULL, dataset = NULL, interval_size = 0.5){
  if (is.null(Y) == TRUE & is.null(Y) == TRUE & is.null(dataset) == TRUE){
    stop("either 'dataset' or at least 'Y' and 'X' must be specified")
  }
  if (is.null(Y) == TRUE & is.null(Y) == TRUE) {
    Y <- dataset[,1]
    X <- dataset[,2]
  }else{
    if (class(Y) == "character" & class(X) == "character") {
      if (is.null(dataset) == TRUE) {
        stop("if 'Y' and 'X' are character, 'dataset' must exist")
      }else{
        Y <- dataset[, Y]
        X <- dataset[, X]
      }
    }else{
      if (!(class(Y) %in% c("numeric","integer")) == TRUE & !(class(X) %in% c("numeric","integer")) == TRUE) {stop("'Y' and 'X' must be either characters or vector but 'class' must be similar")}
    }
  }
  
  data <- data.frame(cbind(Y, X))
  data <- data[order(data$X),]                                                                      # ordering the X values
  
  if (length(X)<4){
    stop("time series length must be at least 4")
  }
  
  Y <- data$Y
  X <- data$X
  
  linear.model <- lm(Y~X)
  
  orthogonal_polynomial <- lm(Y~poly(X,2, raw=F))                                                   
  # After getting Y = gamma*chi + delta*X' + epsilon with orthogonal polynomial
  # we have to perform a variable change to obtain relevant values in the X interval 
  # for first_order_coefficient, second_order_coefficient and intercept,
  # knowing that X'= alpha*X + beta 
  # and chi = eta*X'^2 + theta
  
  gammab  <-  orthogonal_polynomial$coefficients[3]
  delta  <-  orthogonal_polynomial$coefficients[2]
  epsilon  <-  orthogonal_polynomial$coefficients[1]
  
  alpha  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[2]
  beta  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[1]
  
  eta  <-  1/lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[2]
  theta  <-  (-lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[1])*eta
  
  res <- orthogonal_polynomial$residuals
  linear_residuals <- lm(res^2 ~ X)
  
  Y2<-Y*(max(X)-min(X))/(max(Y)-min(Y))                                                             
  # p2 and p3 are relevant when Y and X amplitudes are equivalent,
  # in particular when studying scaled-to-1 indices, Y and X amplitudes
  # may be very different, so we scaled the amplitudes to calculate p2 and p3 
  polynomial_orthonormal_basis<-lm(Y2~poly(X,2, raw=T))$coefficients
  
  if(summary(orthogonal_polynomial)$coefficients[3, 4] <= 0.05){                                     # non linear case
    classification <- data.frame(first_order_coefficient = (delta+2*beta*gammab*eta)*alpha,
                                 first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
                                 second_order_coefficient = (alpha^2)*gammab*eta,
                                 second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
                                 strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
                                 intercept = epsilon+beta*delta+(beta^2)*gammab*eta+gammab*theta,
                                 x_m = (X[length(X)]-X[1])/2+X[1],
                                 p1 = -(delta+2*beta*gammab*eta)/(2*alpha*gammab*eta),                    # points of interest
                                 p2 = (-polynomial_orthonormal_basis[2]+1)/(2*polynomial_orthonormal_basis[3]),
                                 p3 = (-polynomial_orthonormal_basis[2]-1)/(2*polynomial_orthonormal_basis[3]),
                                 CV = sd(Y)/mean(Y),
                                 MSE = mean(res^2),
                                 D = d(Y)
    )
  }else{                                                                                            # linear case
    classification <- data.frame(first_order_coefficient = delta*alpha,
                                 first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
                                 second_order_coefficient = 0,
                                 second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
                                 strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
                                 intercept = epsilon+delta*beta,
                                 x_m = (X[length(X)]-X[1])/2+X[1],
                                 p1 = NA,
                                 p2 = NA,
                                 p3 = NA,
                                 CV = sd(Y)/mean(Y),
                                 MSE = mean(res^2),
                                 D = d(Y)
    )
  }
  
  classification$r.sq <- summary(orthogonal_polynomial)$adj.r.squared                                # retrieve the adjusted coefficient of determination
  
  # compute the derivaive at xm-delta and at xm + delta with delta being half of the input interval size
  derivative  <-  2*(classification$x_m-(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient+classification$first_order_coefficient
  derivative2  <-  2*(classification$x_m+(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient+classification$first_order_coefficient
  
  
  if(sign(derivative) != sign(derivative2)){                                                        # non consistent direction around x_m
    classification$derivative  <-  NA
    classification$intercept_derivative  <-  NA
  }else{                                                                                            # consistent direction around x_m
    classification$derivative  <-  mean(c(derivative, derivative2))
    classification$intercept_derivative  <-  (classification$second_order_coefficient*classification$x_m^2+classification$first_order_coefficient*classification$x_m+classification$intercept)-classification$x_m*classification$derivative
  }
  
  # compute the derivative of the curvature function
  classification$derivated_curvature  <-  -12*(classification$second_order_coefficient^2)*(2*classification$second_order_coefficient*classification$x_m+classification$first_order_coefficient)*(classification$second_order_coefficient/abs(classification$second_order_coefficient))/
    ((1+(2*classification$second_order_coefficient*classification$x_m+classification$first_order_coefficient)^2)^(2.5))
  
  if(classification$second_order_pvalue>0.05){classification$derivated_curvature <- NA}
  
  classification$direction <- NA                                                                    # classify the direction
  classification$direction[which(classification$derivative > 0)] <- "increase"
  classification$direction[which(classification$derivative < 0)] <- "decrease"
  classification$direction[which(is.na(classification$derivative))] <- "stable"
  classification$direction[which(as.numeric(classification$first_order_pvalue)>0.05 & as.numeric(classification$second_order_pvalue)>0.05)] <- "stable"
  
  classification$acceleration <- NA                                                                 # classify the acceleration
  classification$acceleration[which(classification$derivated_curvature < 0)] <- "accelerated"
  classification$acceleration[which(classification$derivated_curvature > 0)] <- "decelerated"
  classification$acceleration[which(classification$direction == "stable" &
                                      classification$second_order_coefficient < 0)] <- "concave"
  classification$acceleration[which(classification$direction == "stable" &
                                      classification$second_order_coefficient > 0)] <- "convex"
  classification$acceleration[which(is.na(classification$derivated_curvature))] <- "constant"
  
  classification$shape_class <- paste(classification$direction,                                       # give the final classification combining direction and acceleration
                                      classification$acceleration,
                                      sep="_")
  
  linear.model.summary <- summary(linear.model)                                                       # provide the linear approach results for comparison
  
  classification$linear_slope <- linear.model.summary$coefficients[2, 1]
  classification$linear_slope_pvalue <- linear.model.summary$coefficients[2, 4]
  classification$linear_intercept <- linear.model.summary$coefficients[1, 1]
  
  classification$first_X_value <- X[1]
  classification$last_X_value <- X[length(X)]
  
  row.names(classification) <- "Y"
  
  return(classification)
  
}

################################################################################
# Reshape data into long form ---------------------------------------------

lpd.long <- lpd %>%
  tidyr::pivot_longer(cols=c(33:103), names_to = "year", values_to = "abundance") %>% 
  dplyr::mutate(year = as.numeric(year),
                abundance = as.numeric(abundance)) %>% 
  drop_na(abundance) # remove NAs in the abundance column

lpd.long <- lpd.long  %>%
  dplyr::distinct() %>%     # remove duplicate rows
  dplyr::group_by(ID) %>%   # group rows so that each group is one population
  dplyr::mutate(minyear = min(year),
                maxyear = max(year),
                duration = (maxyear - minyear)+1,
                scaleab = (abundance - min(abundance))/(max(abundance) - min(abundance)),
                logab = log(abundance+1),
                taxa = fun_taxa(Class)) %>%
  dplyr::filter(is.finite(scaleab)) %>% 
  ungroup() %>% 
  dplyr::group_by(ID,Binomial,Common_name) %>%
  dplyr::mutate(points=length(year))

# save the data into a csv file
write.csv(lpd.long,file="data/lpd_long.csv")

################################################################################
# Calculate population change for each population -------------------------

########## LINEAR MODELS ##########
linear.models <- lpd.long %>%
  # Group by the key variables that we want to iterate over
  dplyr::group_by(ID,Binomial,Common_name,
                  Location,Country,Region,Latitude,Longitude,
                  Class,Order,System,
                  T_realm,T_biome,FW_realm,FW_biome,M_realm,M_biome,
                  minyear,maxyear,duration,points,
                  taxa) %>%
  # Create a linear model for each group
  # Extract model coefficients using tidy() from the broom package
  dplyr::do(broom::tidy(lm(logab ~ year, .))) %>%
  # Filter out slopes and remove intercept values
  dplyr::filter(term == "year") %>%
  # Get rid of the column term as we don't need it any more
  dplyr::select(-term) %>%
  # Remove groupings
  dplyr::ungroup() 


# Add a column to classify population trends according to the slopes and their
# significance
linear.models <- linear.models %>%
  # sorting the dataframe according the mean trend of taxonomic groups
  dplyr::group_by(taxa) %>%
  dplyr::mutate(mean_trend = mean(estimate)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(taxa = forcats::fct_reorder(taxa, -mean_trend)) %>% 
  dplyr::mutate(classif=
                  dplyr::case_when(
                    estimate>0 & p.value<0.05 ~ "positive",
                    estimate<0 & p.value<0.05 ~ "negative",
                    p.value>=0.05 ~ "no trend"))


########## NON-LINEAR MODELS ##########

# According to Rigal et al. (2020) methodology, the non-linear classification of 
# the time series is possible when time series consist of at least 4 points
nl.sub <- lpd.long[which(lpd.long$points>4),]

# Apply the non linear classification on the logged abundances time series
nl.models <- nl.sub %>% 
  dplyr::group_by(ID,Binomial,Common_name,
                  Location,Country,Region,Latitude,Longitude,
                  Class,Order,System,
                  T_realm,T_biome,FW_realm,FW_biome,M_realm,M_biome,
                  minyear,maxyear,duration,points,
                  taxa) %>%
  dplyr::do(class.trajectory(.$logab,.$year)) %>%
  dplyr::ungroup() 

# Calculate the magnitude of the curvature (using the radius of curvature)
# Attribute the usual linear classification (lin_class)
# Attribute "good news" or "bad news" depending on the non-linear trajectories
#     considering declines and concave trajectories are "bad news" and increases and convex trajectories are "good news"
# Re-classify the 9 types of trajectories into 6 more simple classes (dir2)
nl.models <- nl.models %>% 
  dplyr::mutate(magnitude=
                  dplyr::case_when(
                    direction=="increase" & second_order_pvalue<=0.05 ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p3+first_order_coefficient)^2)^(3/2))),
                    direction=="decrease" & second_order_pvalue<=0.05 ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p2+first_order_coefficient)^2)^(3/2))),
                    acceleration=="concave"|acceleration=="convex" ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p1+first_order_coefficient)^2)^(3/2))),
                    second_order_pvalue>0.05 ~ 0),
                lin_class=
                  dplyr::case_when(
                    linear_slope>0 & linear_slope_pvalue<0.05 ~ "positive",
                    linear_slope<0 & linear_slope_pvalue<0.05 ~ "negative",
                    linear_slope_pvalue>=0.05 ~ "no trend"),
                gnbn=dplyr::case_when((direction=="decrease"|acceleration=="concave") ~ "bad",
                                      (direction=="increase"|acceleration=="convex") ~ "good",
                                      shape_class=="stable_constant" ~ "no"),
                dir2=dplyr::case_when((direction=="decrease" & acceleration=="constant")~"decrease linear",
                                      (direction=="decrease" & acceleration!="constant")~"decrease non linear",
                                      (direction=="increase" & acceleration=="constant")~"increase linear",
                                      (direction=="increase" & acceleration!="constant")~"increase non linear",
                                      (direction=="stable" & acceleration=="constant")~"no trend linear",
                                      (direction=="stable" & acceleration!="constant")~"no trend non linear")
  )


# Same methodology for scaled abundance
nl.models.scaled <- nl.sub %>% 
  dplyr::group_by(ID,Binomial,Common_name,
                  Location,Country,Region,Latitude,Longitude,
                  Class,Order,System,
                  T_realm,T_biome,FW_realm,FW_biome,M_realm,M_biome,
                  minyear,maxyear,duration,points,
                  taxa) %>%
  dplyr::do(class.trajectory(.$scaleab,.$year)) %>%
  dplyr::ungroup() 

nl.models.scaled <- nl.models.scaled %>% 
  dplyr::mutate(magnitude=
                  dplyr::case_when(
                    direction=="increase" & second_order_pvalue<=0.05 ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p3+first_order_coefficient)^2)^(3/2))),
                    direction=="decrease" & second_order_pvalue<=0.05 ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p2+first_order_coefficient)^2)^(3/2))),
                    acceleration=="concave"|acceleration=="convex" ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p1+first_order_coefficient)^2)^(3/2))),
                    second_order_pvalue>0.05 ~ 0),
                lin_class=
                  dplyr::case_when(
                    linear_slope>0 & linear_slope_pvalue<0.05 ~ "positive",
                    linear_slope<0 & linear_slope_pvalue<0.05 ~ "negative",
                    linear_slope_pvalue>=0.05 ~ "no trend"),
                gnbn=dplyr::case_when((direction=="decrease"|acceleration=="concave") ~ "bad",
                                      (direction=="increase"|acceleration=="convex") ~ "good",
                                      shape_class=="stable_constant" ~ "no"),
                dir2=dplyr::case_when((direction=="decrease" & acceleration=="constant")~"decrease linear",
                                      (direction=="decrease" & acceleration!="constant")~"decrease non linear",
                                      (direction=="increase" & acceleration=="constant")~"increase linear",
                                      (direction=="increase" & acceleration!="constant")~"increase non linear",
                                      (direction=="stable" & acceleration=="constant")~"no trend linear",
                                      (direction=="stable" & acceleration!="constant")~"no trend non linear")
  )


# Same for raw abundance
nl.models.raw <- nl.sub %>% 
  dplyr::group_by(ID,Binomial,Common_name,
                  Location,Country,Region,Latitude,Longitude,
                  Class,Order,System,
                  T_realm,T_biome,FW_realm,FW_biome,M_realm,M_biome,
                  minyear,maxyear,duration,points,
                  taxa) %>%
  dplyr::do(class.trajectory(.$abundance,.$year)) %>%
  dplyr::ungroup() 

nl.models.raw <- nl.models.raw %>% 
  dplyr::mutate(magnitude=
                  dplyr::case_when(
                    direction=="increase" & second_order_pvalue<=0.05 ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p3+first_order_coefficient)^2)^(3/2))),
                    direction=="decrease" & second_order_pvalue<=0.05 ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p2+first_order_coefficient)^2)^(3/2))),
                    acceleration=="concave"|acceleration=="convex" ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p1+first_order_coefficient)^2)^(3/2))),
                    second_order_pvalue>0.05 ~ 0),
                lin_class=
                  dplyr::case_when(
                    linear_slope>0 & linear_slope_pvalue<0.05 ~ "positive",
                    linear_slope<0 & linear_slope_pvalue<0.05 ~ "negative",
                    linear_slope_pvalue>=0.05 ~ "no trend"),
                gnbn=dplyr::case_when((direction=="decrease"|acceleration=="concave") ~ "bad",
                                      (direction=="increase"|acceleration=="convex") ~ "good",
                                      shape_class=="stable_constant" ~ "no"),
                dir2=dplyr::case_when((direction=="decrease" & acceleration=="constant")~"decrease linear",
                                      (direction=="decrease" & acceleration!="constant")~"decrease non linear",
                                      (direction=="increase" & acceleration=="constant")~"increase linear",
                                      (direction=="increase" & acceleration!="constant")~"increase non linear",
                                      (direction=="stable" & acceleration=="constant")~"no trend linear",
                                      (direction=="stable" & acceleration!="constant")~"no trend non linear")
  )

########## RANDOM MODELS ##########
nl.modelsrd <- data.frame()

for (i in (1:100)){
  nl.modsi <- nl.sub20 %>%
    dplyr::group_by(ID,Binomial,Common_name,
                    Location,Country,Region,Latitude,Longitude,
                    Class,Order,System,
                    T_realm,T_biome,FW_realm,FW_biome,M_realm,M_biome,
                    minyear,maxyear,duration,points,
                    taxa) %>%
    # for each population time series, abundances are re-sampled and the trajectory
    # and temporal variability are extracted 
    # this process is repeated 100 times (loop for)
    dplyr::do(class.trajectory(sample(.$logab,length(.$logab)),.$year)) %>%
    dplyr::ungroup()

  nl.modelsrd <- bind_rows(nl.modelsrd,nl.modsi)
}

# For each population time series, we calculated the average randomized D, CV and MSE values
nl.modsrd.sum <- nl.modelsrd %>%
  group_by(ID,Binomial) %>%
  summarize(meanD=mean(D),
            meanCV=mean(CV),
            meanMSE=mean(MSE))

write.csv(nl.modelsrd, "data/random_models.csv")
write.csv(nl.modsrd.sum, "data/random_models_summary.csv")

########## SUB-SELECTION MODELS ##########

# We selected population time series (TS) with 20 years of data at least
nl.sub20 <- lpd.long[which(lpd.long$points>=20),]

# linear models # 32,211 TS
length(unique(linear.models$Binomial)) # 4883 unique species

# non linear models # 21,841 TS
# populations w/ 10 years data
final.nl10 <- nl.models %>% filter(points>=10) # 14,695 TS

# populations w/ 20 years data - log abundance
final.nl <- nl.models %>% filter(points>=20)   # 6,437 TS

# populations w/ 20 years data - scaled abundance
final.nl.sc <- nl.models.scaled %>% filter(points>=20)

# populations w/ 20 years data - raw abundance
final.nl.raw <- nl.models.raw %>% filter(points>=20)

# Some transformation to make the dataset manageable
final.nl$dir2<-as.factor(final.nl$dir2)
final.nl$System<-as.factor(final.nl$System)
final.nl$taxa<-as.factor(final.nl$taxa)
final.nl$Latitude<-as.numeric(final.nl$Latitude)
final.nl <- final.nl %>% 
  mutate(full_realm=case_when(T_realm!="NULL"~T_realm,
                              FW_realm!="NULL"~FW_realm,
                              M_realm!="NULL"~M_realm),
         LoNL = ifelse(acceleration=="constant","linear","non linear")) 

final.nl<-left_join(final.nl,RLC,by="Binomial")
levels(final.nl$dir2)<-c("Decrease\nlinear","Decrease\nnon linear","Increase\nlinear",
                         "Increase\nnon linear","No trend\nlinear","No trend\nnon linear")

final.nl$redlistCategory<-as.factor(final.nl$redlistCategory) 

length(unique(final.nl$Binomial))  # 1,257 unique species in the final dataset

################################################################################
### First small exploration
table(final.nl$shape_class) # Number of TS within each trajectory type
table(final.nl$acceleration=="constant")    # 2,887 non-linear trajectories : 44.8 % (log abundance)
table(final.nl.sc$acceleration=="constant") # 45.8 % (scaled abundance)
table(final.nl.raw$acceleration=="constant")# 45.8 % (raw abundance)

table(final.nl10$acceleration=="constant") # With 10 years data : 30 %

# Add Z scores : we calculated a D (resp. CV, MSE) corrected for temporal variability 
# that would be expected by chance. We calculated the D statistics (resp. CV, MSE) 
# for the randomized abundance time series. This process was repeated 100 times 
# for each time series, and the average randomized D (resp. CV, MSE) value was calculated. 
# The corrected D (resp. CV, MSE) was obtained by subtracting the value expected by
# chance from the value calculated on the true series. Thus, a positive corrected D 
# (resp. CV, MSE) means that the population exhibits greater temporal variability 
# than would be expected by chance.
nl.modstestD <- left_join(final.nl,nl.modsrd.sum,by=c("ID","Binomial")) %>% 
  mutate(ZscoreD=D-meanD,
         ZscoreCV=CV-meanCV,
         ZscoreMSE=MSE-meanMSE) 

write.csv(nl.modstestD, "data/non_linear_models.csv")
