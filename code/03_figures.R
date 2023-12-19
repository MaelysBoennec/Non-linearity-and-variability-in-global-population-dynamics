################################################################################
### Non-linearity and temporal variability are overlooked components of 
### global population dynamics
### 12/2023
### Maelys Boennec (maelys.boennec@umontpellier.fr)

#### Part 3 - Figures

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
library(ggridges)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
#source("https://raw.githubusercontent.com/larmarange/JLutils/master/R/lm_right.R")

################################################################################
# Load the data -----------------------------------------------------------

final.nl <- read.csv("data/non_linear_models.csv")
nl.modelsrd <- read.csv("data/random_models.csv")
nl.modsrd.sum <- read.csv("data/random_models_summary.csv")

################################################################################
# Functions ---------------------------------------------------------------

plot_classif_tot <- function(seed,data.long,data,col,xsp,ysp,lin=T){
  set.seed(seed)
  data<-data.long[data.long$ID==sample(data$ID,1),]
  sp<-stringr::str_replace(data$Binomial[1],"_"," ")
  abtype<-data.long$Units
  p<-ggplot(data=data,aes(x=year,y=logab))+
    geom_line(linewidth=1.5,color="grey50")
  if (lin){p<-p+geom_smooth(method="lm",se=FALSE,color=col,size=2)}
  else{p<-p+geom_smooth(method="lm",formula=y~poly(x,2),se=FALSE,color=col,size=2)}
  p<-p+theme_LPI()+ylab(paste("log(",abtype,")"))+
    annotate(geom="text", x=xsp, y=ysp, label=sp, size=5, family="Montserrat",fontface="italic")
}

theme_LPI <- function(){
  theme_classic() +
    theme(plot.title = element_text(size=14,family="Montserrat",face="bold",hjust=0.5),
          axis.text = element_text(size = 14,family="Montserrat"), 
          axis.title = element_text(size=14,family="Montserrat"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          # panel.grid.major.x = element_blank(),                                          
          # panel.grid.minor.x = element_blank(),
          # panel.grid.minor.y = element_blank(),
          # panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(0.7, 0.7, 0.7, 0.7), units = , "cm"))
}

drawWorld<-function(lats) {
  world_map<-map_data("world")
  g1<-ggplot()+coord_fixed()+xlab("")+ylab("")
  g1<-g1+geom_polygon(data=world_map, aes(x=long, y=lat, group=group), colour="gray80", fill="gray80")
  g1<-g1+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
               panel.background=element_rect(fill="white", colour="white"), axis.line=element_line(colour="white"),
               axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
  return(g1)
}

plot_ex <- function(seed,data.long,data,lin=T){
  set.seed(seed)
  data<-data.long[data.long$ID==sample(data$ID,1),]
  p<-ggplot(data=data,aes(x=year,y=logab))
  if (lin){p<-p+geom_smooth(method="lm",se=FALSE,color="black",size=1.5)}
  else{p<-p+geom_smooth(method="lm",formula=y~poly(x,2),se=FALSE,color="black",size=1.5)}
  p<-p+theme_minimal()+theme(panel.grid.major = element_blank())
}

plot_prop <- function(p,col,line=F){
  p <- p + geom_bar(position = "fill") +
    geom_text(aes(label = after_stat(count)), stat = "count", position = position_fill(.5),family="Montserrat",fontface="bold",color="white") +
    xlab("") + ylab("") +
    scale_fill_manual(values=col)
  if(line==T){p <- p+ geom_hline(yintercept = 0.448, linetype="dashed",color="grey30")}
  p <- p +  coord_flip() +
    theme_classic() +
    theme(text = element_text(family="Montserrat",face="bold"),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank())
  p
}

plot_prop_traj <- function(p,col){
  p <- p + geom_bar(position = "fill") +
    geom_text(aes(label = after_stat(count)), stat = "count", position = position_fill(.5),family="Montserrat",fontface="bold",color="white") +
    xlab("") + ylab("") +
    scale_fill_manual(values=col)
  #if(line==T){p <- p+ geom_hline(yintercept = sum(dframe_taxa$LoNL=="non linear")/6437, linetype="dashed",color="grey30")}
  p <- p +  coord_flip() +
    theme_classic() +
    theme(text = element_text(family="Montserrat",face="bold"),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank())
  p
}

plot_prop2 <- function(p,col,line=F){
  p <- p + geom_bar(position = "fill") +
    geom_text(aes(label = after_stat(count)), stat = "count", position = position_fill(.5),family="Montserrat",fontface="bold",color="white") +
    xlab("") + ylab("") +
    scale_fill_manual(values=col)
  if(line==T){p <- p+ geom_hline(yintercept = sum(dframe_taxa$LoNL=="non linear")/6437, linetype="dashed",color="grey30")}
  p <- p +  
    theme_classic() +
    theme(text = element_text(family="Montserrat",face="bold"),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank())
  p
}

plot_number <- function(df,ctgr,col){
  p <- ggplot(df, aes(x = fct_rev(fct_infreq(ctgr)), fill = ctgr)) +
     geom_bar() +
     geom_text(aes(label = after_stat(count)), hjust=1,stat = "count",family="Montserrat",fontface="bold",color="black") +
     xlab("") + ylab("Number of population") +
     scale_fill_manual(values=col)+
     coord_flip()+ theme_classic() +
     theme(text = element_text(family="Montserrat",face="bold"),
           legend.position = "none",
           axis.ticks.y = element_blank(),
           axis.line.y = element_blank())
  p
}

plot_rain <- function(p){
  p <- p + 
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), fill="lightslategray", color="lightslategray",scale="width") +
    geom_point(aes(y = D), color="lightslategray",position = position_jitter(width = .15), size = .5, alpha = 0.5) +
    geom_boxplot(width = .2, outlier.shape = NA, alpha = 0.8,color="lightslategray") +
    xlab(element_blank()) + #coord_flip() + 
    #scale_y_continuous(trans="log10")+
    #stat_summary(fun=mean, colour="black", geom="point", shape=16, size=4) +
    guides(fill=FALSE, colour=FALSE) +
    theme_classic() + ylab("D") +
    theme(legend.position = "none",text = element_text(family="Monserrat",size=15),
          axis.line.x = element_blank(),axis.ticks.x = element_blank(),
          axis.text.x=element_text(face="bold"))
  p
}

################################################################################
# Colors ------------------------------------------------------------------
cols <- c("#8F0606", "#D1AB41", "#426962")
cols22 <- c("#8F0606", "#426962", "#D1AB41")
cols2 <- c("#6E2E2E", "#AD8080", "#4D806A", "#739989", "#8F7232", "#E0BD72")
colstraj <- c("#A84646","#C4993D","#4D806A","#739989","#E0BD72","#AD8080")
colprop<-c("#4c8cac", "#cf994c")
cols3 <- c("#156085", "#FF9721","#156085", "#FF9721","#156085", "#FF9721")
cols4 <- c("#A84646","#6E2E2E","#AD8080",
           "#4D806A","#38574A", "#739989",
           "#C4993D","#E0BD72","#8F7232")
cols5 <- c("#6E2E2E","#A84646","#38574A",
           "#4D806A","#C4993D","#8F7232",
           "#AD8080", "#739989","#E0BD72")
colssyst <- c("#4D806A","#273045","#8F7232")
colssyst <- c("#0b775e", "#273046", "#a2a475")
colstaxa <- c("#4D806A","#8F7232","#156085","#AD8080","#FF9721","#2FA1A1")
colsregion <- c("#AD8080"s,"#72B5DB","#4D806A","#8F7232","#156085","#FF9721","#AA9AB8","#FFD16E")
colsrealms=c("#AD8080",'#94b6d2','#704404','#34405e',
             '#d8b25c','#f7b615','#968c8c','#dd8047',
             '#ebddc3','#7ba79d','#a5ab81','#38784d',
             '#775f55')
colstaxawhole <- c("#4D806A","#2FA1A1","#8F7232","#FF9721","#156085","#AD8080",
                   "#4D806A","#2FA1A1","#8F7232","#FF9721","#156085","#AD8080")
colstaxa2 <- c("#4D806A","#8F7232","#156085","#AD8080","#FF9721","#2FA1A1",
               "#8F7232","#156085","#AD8080","#FF9721","#2FA1A1")
colstaxa3 <- c("#4D806A","#8F7232","#156085","#AD8080","#FF9721","#2FA1A1",
               "#4D806A","#8F7232","#156085","#AD8080","#FF9721","#2FA1A1")
colstaxa3 <- c("#4D806A","#8F7232","#156085","#AD8080","#FF9721","#2FA1A1",
               "darkseagreen","#C9AF72","skyblue3","rosybrown2","#FFBC70","#88D4D4")
colstaxa4 <- c("darkseagreen","#C9AF72","skyblue3","rosybrown2","#FFBC70","#88D4D4",
               "#C9AF72","skyblue3","rosybrown2","#FFBC70","#88D4D4")

################################################################################
############################### MAIN FIGURES ###################################
################################################################################

################################################################################
# Figure 2 | All classifications ------------------------------------------

# A) Linear decrease
lin_dec <- final.nl  %>% filter(shape_class=="decrease_constant" & lin_class=="negative")
lin_dec.long <- lpd.long %>% dplyr::filter(ID %in% lin_dec$ID)
(A<-plot_classif_tot(8,lin_dec.long,lin_dec,col="#A84646",xsp=1995,ysp=0.55)+
    annotate(geom="text", x=2011, y=0.78, label="N=1170\n18.2%", family="Montserrat",fontface="bold",size=5,col="#A84646")+
    ggtitle("Linear decrease"))

# B) No trend (linear)
no_trend <- final.nl  %>% filter(shape_class=="stable_constant" & lin_class=="no trend")
no_trend.long <- lpd.long %>% dplyr::filter(ID%in% no_trend$ID)
(B<-plot_classif_tot(33,no_trend.long,no_trend,col="#C4993D",xsp=1985,ysp=1)+
    annotate(geom="text", x=2004.5, y=1.68, label="N=1337\n20.8%", family="Montserrat",fontface="bold",size=5,col="#C4993D")+
    ggtitle("No trend"))

# C) Linear increase
lin_inc <- final.nl  %>% filter(shape_class=="increase_constant" & lin_class=="positive")
lin_inc.long <- lpd.long %>% dplyr::filter(ID%in% lin_inc$ID)
(C<-plot_classif_tot(13,lin_inc.long,lin_inc,col="#4D806A",xsp=2005,ysp=4.45)+
    annotate(geom="text", x=1980, y=4.93, label="N=1043\n16.2%", family="Montserrat",fontface="bold",size=5,col="#4D806A")+
    ggtitle("Linear increase"))

# D) Accelerated decrease
accel_dec <- final.nl %>% filter(shape_class=="decrease_accelerated" & lin_class=="negative")
accel_dec.long <- lpd.long %>% dplyr::filter(ID%in% accel_dec$ID)
(D<-plot_classif_tot(42,accel_dec.long,accel_dec,col="#6E2E2E",xsp=1982,ysp=2.2,lin=F)+
    annotate(geom="text", x=2010, y=2.76, label="N=215\n3.3%", family="Montserrat",fontface="bold",size=5,col="#6E2E2E")+
    ggtitle("Decrease accelerated"))

# E) Concave
concave <- final.nl  %>% filter(shape_class=="stable_concave" & lin_class=="no trend")
concave.long <- lpd.long %>% dplyr::filter(ID%in% concave$ID)
(E<-plot_classif_tot(4,concave.long,concave,col="#E0BD72",xsp=1990,ysp=4.5,lin=F)+
    annotate(geom="text", x=1993, y=8.5, label="N=757\n11.8%", family="Montserrat",fontface="bold",size=5,col="#E0BD72")+
    ggtitle("Concave"))

# F) Decelerated increase
decel_inc <- final.nl  %>% filter(shape_class=="increase_decelerated" & lin_class=="positive")
decel_inc.long <- lpd.long %>% dplyr::filter(ID%in% decel_inc$ID)
(Fi<-plot_classif_tot(6,decel_inc.long,decel_inc,col="#739989",xsp=1985,ysp=6.7,lin=F)+
    annotate(geom="text", x=1955, y=8.75, label="N=276\n4.3%", family="Montserrat",fontface="bold",size=5,col="#739989")+
    ggtitle("Increase decelerated"))

# G) Decelerated decrease
decel_dec <- final.nl %>% filter(shape_class=="decrease_decelerated" & lin_class=="negative")
decel_dec.long <- lpd.long %>% dplyr::filter(ID%in% decel_dec$ID)
(G<-plot_classif_tot(57,decel_dec.long,decel_dec,col="#AD8080",xsp=1965,ysp=6.4,lin=F)+
    annotate(geom="text", x=2006, y=9, label="N=554\n8.6%", family="Montserrat",fontface="bold",size=5,col="#AD8080")+
    ggtitle("Decrease decelerated"))

# H) Convex
convex <- final.nl  %>% filter(shape_class=="stable_convex" & lin_class=="no trend")
convex.long <- lpd.long %>% dplyr::filter(ID%in% convex$ID)
(H<-plot_classif_tot(15,convex.long,convex,col="#8F7232",xsp=1972,ysp=5.7,lin=F)+
    annotate(geom="text", x=1997.5, y=5.8, label="N=621\n9.6%", family="Montserrat",fontface="bold",size=5,col="#8F7232")+
    ggtitle("Convex"))

# I) Accelerated increase
accel_inc <- final.nl  %>% filter(shape_class=="increase_accelerated" & lin_class=="positive")
accel_inc.long <- lpd.long %>% dplyr::filter(ID%in% accel_inc$ID)
(I<-plot_classif_tot(61,accel_inc.long,accel_inc,col="#38574A",xsp=1989,ysp=0,lin=F)+
    annotate(geom="text", x=1972, y=1.15, label="N=464\n7.2%", family="Montserrat",fontface="bold",size=5,col="#38574A")+
    ggtitle("Increase accelerated"))

(all<-cowplot::plot_grid(A,B,C,D,E,Fi,G,H,I,labels=c("A","B","C","D","E","F","G","H","I")))
ggsave("outputs/all_classif.png",all,type='cairo',width=15,height=12,dpi=200)

################################################################################
# Figure 3 | Cross-classification -----------------------------------------
final.nl$shape_class <- as.factor(final.nl$shape_class)

double_classif <- final.nl %>% 
  mutate(shape_class=fct_relevel(shape_class,
                                 "decrease_constant","decrease_accelerated","decrease_decelerated",
                                 "increase_constant","increase_accelerated","increase_decelerated",  
                                 "stable_constant","stable_concave","stable_convex"))

(cross_classif <- ggplot(double_classif, aes(x = lin_class, fill = shape_class)) +
    geom_bar(position = "fill",width = 0.25) +
    geom_text(aes(label = ifelse(after_stat(count)>10,after_stat(count),"")), stat = "count", position = position_fill(.5),family="Montserrat",fontface="bold",color="white") +
    xlab("") + ylab("Proportion") +
    scale_fill_manual(values=cols4)+ theme_void() + labs(fill = "Linearity") + 
    theme(text = element_text(family="Montserrat",face="bold"),legend.position = "none"))

ggsave("outputs/cross_classif.png",cross_classif,width = 7, height = 6)


################################################################################
# Figure 4 | Temporal variability across trajectory types -----------------

# Comparison across shape class
(rain_D_shape <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(dir2, D, mean), y = D, fill = dir2)) +
   geom_flat_violin(aes(color=dir2), position = position_nudge(x = 0.2, y = 0),scale="width") +
   geom_point(aes(y = D, color = dir2),position = position_jitter(width = .15), size = .5, alpha = 0.2) +
   geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.7) +
   #geom_pointrange()+
   xlab(element_blank()) + coord_flip() +
   scale_colour_manual(values = cols2) +
   scale_fill_manual(values = cols2) +
   #stat_summary(aes(fill=dir2), fun=mean, geom="point", shape=16, size=3) +
   geom_text(data = shape.letD_df, aes(label = letters, x=dir2, y = 3.5 ), size=4, family="Montserrat",fontface="italic",colour="black", size=5)  +
   guides(fill=FALSE, colour=FALSE) +
   geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
   theme_classic() + ylab("D") +
   theme(legend.position = "none",text = element_text(family="Monserrat"),
         axis.line.y = element_blank(),axis.ticks.y = element_blank(),
         axis.text.y=element_text(face="bold")))


full.shape.summary$dir2<-as.factor(full.shape.summary$dir2)
full.shape.summary$dir2<-fct_relevel(full.shape.summary$dir2,
                                     levels(forcats::fct_reorder(final.nl$dir2, final.nl$MSE, mean)))
(var_plot_shape <- ggplot(full.shape.summary, aes(x = dir2, y = mean)) +
    geom_pointrange(aes(ymin = CIinf, ymax = CIsup, colour = variability,
                        shape = variability),
                    position = position_dodge(0.7), size = 0.5) + 
    scale_colour_manual(values = c("grey70", "grey40", "black")) +
    theme_classic() + coord_flip() +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    theme(axis.title.x = element_text(family="Montserrat"),
          axis.text.y = element_text(family="Monserrat",hjust=1,face="bold"),
          axis.text.x = element_text(family="Monserrat"),
          legend.text = element_text(family="Monserrat", face = "italic"),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.9, 0.2),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5)) +
    labs(y = "Effect size"))

(co_var <- cowplot::plot_grid(rain_D_shape,var_plot_shape,rel_widths=c(1,0.8),labels=c("A","B"),align="hv"))
ggsave("outputs/variability_full.png",co_var,height = 6,width=12)

################################################################################
# Figure 5 | Non linearity and temporal variability among habitat  --------

### Linear vs Non linear Habitats
dframe_reg <- final.nl %>% 
  mutate(Region = forcats::fct_reorder(.f = Region,
                                       .x = LoNL,
                                       .fun = function(.x) mean(.x == "non linear"),
                                       .desc = FALSE),
         System = forcats::fct_reorder(.f = System,
                                       .x = LoNL,
                                       .fun = function(.x) mean(.x == "non linear"),
                                       .desc = FALSE))

prop_hab <- ggplot(dframe_reg, aes(x = System, fill = LoNL)) 
(prop_hab <- plot_prop(prop_hab,colprop,line=T)+theme(legend.position = "none")+
    #ggtitle("Relative proportions of non-linear and linear trajectories")+
    theme(plot.title = element_text(hjust=0.5)))

### Linear vs Non linear Regions
prop_reg <- ggplot(dframe_reg, aes(x = Region, fill = LoNL))
(prop_reg <- plot_prop(prop_reg,colprop,line=T)+theme(legend.position = "bottom")+
    #ggtitle("Region")+
    theme(plot.title = element_text(hjust=0.5)))

### Temporal variability habitats
full.hab.summary$system<-as.factor(full.hab.summary$system)
full.hab.summary$system<-fct_relevel(full.hab.summary$system,rev(levels(dframe_reg$System)))
(var_plot_hab <- ggplot(full.hab.summary, aes(x = system, y = mean)) +
    geom_pointrange(aes(ymin = CIinf, ymax = CIsup, colour = variability,
                        shape = variability),
                    position = position_dodge(0.7), size = 0.5) +
    scale_colour_manual(values = c("grey70", "grey40", "black")) +
    scale_y_continuous(breaks = c(0.5, 0),
                       limits = c(-0.1, 0.6),
                       labels = c("0.5", "0"))  +
    theme_classic() +
    theme(axis.title.y = element_text(family="Montserrat",face="bold"),
          axis.text.x = element_text(family="Montserrat",angle=40,hjust=1,face="bold"),
          axis.text.y = element_text(family="Monserrat"),
          legend.text = element_text(family="Monserrat", face = "italic"),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.2, 0.8),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Effect size"))

### Temporal variability regions
full.reg.global$region<-as.factor(full.reg.global$region)
full.reg.global$region<-fct_relevel(full.reg.global$region,rev(levels(dframe_reg$Region)))
(var_plot_reg <- ggplot(full.reg.global, aes(x = region, y = mean)) +
    geom_pointrange(aes(ymin = CIinf, ymax = CIsup, colour = variability,
                        shape = variability),
                    position = position_dodge(0.7), size = 0.5) +
    scale_colour_manual(values = c("grey70", "grey40", "black")) +
    scale_y_continuous(breaks = c(1, 0.5, 0),
                       limits = c(-0.1, 1.1),
                       labels = c("1.0", "0.5", "0")) +
    theme_classic() +
    theme(axis.title.y = element_text(family="Montserrat",face="bold"),
          axis.text.x = element_text(family="Montserrat",angle=40,hjust=1,face="bold"),
          axis.text.y = element_text(family="Montserrat"),
          legend.text = element_text(family="Montserrat", face = "italic"),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.05, 0.8),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Effect\nsize"))

### Final plot 
(ppp1<-cowplot::plot_grid(prop_hab,prop_reg,nrow=2,labels = c("A","B"),align = "hv",rel_heights = c(0.35,1)))
(ppp2<-cowplot::plot_grid(var_plot_hab,var_plot_reg,labels = c("C","D"),nrow=1,align="hv",rel_widths = c(0.3,0.8)))
(coploreg<-cowplot::plot_grid(ppp1,ppp2,ncol=1))
ggsave("outputs/nl_var_reg.png",coploreg,height = 9,width=11)


################################################################################
# Figure 6 | Non-linearity and temporal variability among taxonomi --------

### Linear vs Non linear Taxa
dframe_taxa <- final.nl %>% 
  filter(taxa != "Invertebrates") %>% 
  mutate(taxa = forcats::fct_reorder(.f = taxa,
                                     .x = LoNL,
                                     .fun = function(.x) mean(.x == "linear"),
                                     .desc = TRUE))

dframe_taxa$shape_class<-fct_relevel(dframe_taxa$shape_class,
                                     levels(c("decrease_accelerated","decrease_decelerated","increase_accelerated","increase_decelerated",
                                              "stable_concave","stable_convex","decrease_constant","increase_constant","stable_constant")))

prop_taxa <- ggplot(dframe_taxa, aes(x = taxa, fill = LoNL)) 
(prop_taxa <- plot_prop(prop_taxa,colprop,line=T)+theme(legend.position = "none")+
    #ggtitle("Taxonomic group")+
    theme(plot.title = element_text(hjust=0.5)))

### Linear vs Non linear Red List Category
dframe_rlc <- final.nl %>% 
  filter(!(is.na(redlistCategory)))

dframe_rlc$redlistCategory<-factor(dframe_rlc$redlistCategory,
                                   levels=c("Extinct in the Wild",
                                            "Critically Endangered",
                                            "Endangered",
                                            "Vulnerable",
                                            "Near Threatened",
                                            "Least Concern",
                                            "Data Deficient"))

prop_rlc <- ggplot(dframe_rlc %>% filter(redlistCategory!="NA"&redlistCategory!="Extinct in the Wild"), aes(x = redlistCategory, fill = LoNL)) 
(prop_rlc <- plot_prop(prop_rlc,colprop,line=T)+#ggtitle("IUCN Red List Category")+
    theme(plot.title = element_text(hjust=0.5)))

### Temporal variability taxa
full.tax.summary$taxa<-as.factor(full.tax.summary$taxa)
full.tax.summary$taxa<-fct_relevel(full.tax.summary$taxa,rev(levels(dframe_taxa$taxa)))
(var_plot_tax <- ggplot(full.tax.summary %>% filter(taxa!="Invertebrates"), aes(x = taxa, y = mean)) +
    geom_pointrange(aes(ymin = CIinf, ymax = CIsup, colour = variability,
                        shape = variability),
                    position = position_dodge(0.7), size = 0.5) +
    scale_colour_manual(values = c("grey70", "grey40", "black")) +
    # scale_y_continuous(breaks = c(1, 0.5, 0),
    #                    limits = c(-0.1, 1),
    #                    labels = c("1.0", "0.5", "0")) +
    theme_classic() +
    theme(axis.title.y = element_text(family="Montserrat",face="bold"),
          axis.text.x = element_text(family="Monserrat",angle=30,hjust=1,face="bold"),
          axis.text.y = element_text(family="Monserrat"),
          legend.text = element_text(family="Monserrat", face = "italic"),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.06, 0.8),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Effect\nsize"))

### Temporal variability RLC
full.rlc.summary$redlistCategory<-as.factor(full.rlc.summary$redlistCategory)
full.rlc.summary$redlistCategory<-fct_relevel(full.rlc.summary$redlistCategory,rev(levels(dframe_rlc$redlistCategory)))
(var_plot_rlc <- ggplot(full.rlc.summary,
                        aes(x = redlistCategory, y = mean)) +
    geom_pointrange(aes(ymin = CIinf, ymax = CIsup, colour = variability,
                        shape = variability),
                    position = position_dodge(0.7), size = 0.5) +
    scale_colour_manual(values = c("grey70", "grey40", "black")) +
    theme_classic() +
    theme(axis.title.y = element_text(family="Montserrat",face="bold"),
          axis.text.x = element_text(family="Monserrat",angle=30,hjust=1,face="bold"),
          axis.text.y = element_text(family="Monserrat"),
          legend.text = element_text(family="Monserrat", face = "italic"),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.06, 0.8),
          #legend.position="none",
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Effect\nsize"))

### Final plot
(coplotaxorlc<-cowplot::plot_grid(prop_taxa,prop_rlc,var_plot_tax,var_plot_rlc,
                                 nrow=4,labels = c("A","B","C","D"),align = "hv",
                                 rel_heights = c(0.6,0.7,0.8,0.8)))
ggsave("outputs/nl_var_tax.png",coplotaxorlc,height = 10,width=11)


################################################################################
########################## SUPPLEMENTARY FIGURES ###############################
################################################################################

################################################################################
# SM1 | Temporal, geographical, and taxonomic extent ----------------------

### FS1.1 | DURATION 
(duration.all <- ggplot(final.nl, aes(points)) +
   geom_histogram(binwidth = 1, alpha = 0.2, position = "identity") +
   geom_line(stat = "density", aes(y = ..count..), size = 1.5) +
   theme_classic() +
   scale_y_continuous(limits = c(0, 450)) +
   geom_vline(xintercept = mean(final.nl$points), size = 1) +
   scale_x_continuous(breaks = c(20, 25, 30, 35, 40, 45, 50, 55,60,65,70)) +
   theme(axis.line.x = element_line(color = "black"),
         axis.line.y = element_line(color = "black"),
         text = element_text(family="Montserrat",size=15,face="bold")) +
   labs(x = "Number of monitoring years", y = "Number of time series") +
   geom_vline(xintercept = 19.5, colour = "grey30", size = 1, linetype = "dashed") +
   guides(colour = F, fill = F))
ggsave("outputs/temporal_extent.png",duration.all,width=15,height = 10)

### FS1.2 | DISTRIBUTIONS ACROSS BIOGEOGRAPHIC AND TAXONOMIC GROUPS
barhab<-plot_number(final.nl,final.nl$System,colssyst)
bartax<-plot_number(final.nl %>% filter(taxa!="Invertebrates"),(final.nl %>% filter(taxa!="Invertebrates"))$taxa,colstaxa)
barreg<-plot_number(final.nl,final.nl$Region,colsregion)
barrealm<-plot_number(final.nl,final.nl$full_realm,colsrealms)
bariucn<-plot_number(data.rlc,data.rlc$redlistCategory,c("#ff000e","#ffa15b","#ffef00", "#d5db00","#4bb33c",'grey'))

co11<-cowplot::plot_grid(NULL,barhab,ncol=2,rel_widths = c(0.2,1))
co21<-cowplot::plot_grid(NULL,bartax,ncol=2,rel_widths = c(0.15,1))
co1<-cowplot::plot_grid(co11,barreg,labels=c("A","C"),ncol=1,rel_heights = c(0.3,0.8))
co2<-cowplot::plot_grid(co21,bariucn,labels=c("B","D"),ncol=1,rel_heights = c(0.6,0.7))
co3<-cowplot::plot_grid(co1,co2,ncol=2)
(co4<-cowplot::plot_grid(co3,barrealm,labels=c("","E"),ncol=1))
ggsave("outputs/extents.png",co4,height = 12,width=15)

### FS1.3 | MAPS
(world.map.hab <- drawWorld() +
   geom_point(data = final.nl, aes(x = Longitude, y = Latitude, color=System),size=1,alpha = 0.7)+ 
   scale_colour_manual(values=colssyst) +
   scale_y_continuous(limits = c(-80, 80)) +
   #theme(legend.position = "bottom") +
   guides(color = guide_legend(override.aes = list(size = 4,alpha=1))))

(world.map.tax <- drawWorld() +
    geom_point(data = final.nl %>% filter(taxa!="Invertebrates"), aes(x = Longitude, y = Latitude, color=taxa),size=1,alpha = 0.7)+ 
    scale_colour_manual(values=colstaxa) +
    scale_y_continuous(limits = c(-80, 80)) +
    #theme(legend.position = "bottom",legend.direction = "horizontal") +
    guides(color = guide_legend(ncol=1,override.aes = list(size = 4,alpha=1))))

(world.map.reg <- drawWorld() +
    geom_point(data = final.nl, aes(x = Longitude, y = Latitude, color=Region),size=1,alpha = 0.7)+ 
    scale_colour_manual(values=colsregion) +
    scale_y_continuous(limits = c(-80, 80)) +
    #theme(legend.position = "bottom",legend.direction = "horizontal") +
    guides(color = guide_legend(ncol=1,override.aes = list(size = 4,alpha=1))))

(world.map.realm <- drawWorld() +
    geom_point(data = final.nl, aes(x = Longitude, y = Latitude, color=full_realm),size=1,alpha = 0.7)+ 
    scale_colour_manual(values=colsrealms) +
    scale_y_continuous(limits = c(-80, 80)) +
    #theme(legend.position = "bottom",legend.direction = "horizontal") +
    guides(color = guide_legend(ncol=1,override.aes = list(size = 4,alpha=1))))

(world.map.iucn <- drawWorld() +
    geom_point(data = data.rlc, aes(x = Longitude, y = Latitude, color=redlistCategory),size=1,alpha = 0.7)+ 
    scale_colour_manual(values=c("#ff000e","#ffa15b","#ffef00", "#d5db00","#4bb33c",'grey')) +
    scale_y_continuous(limits = c(-80, 80)) +
    #theme(legend.position = "bottom",legend.direction = "horizontal") +
    guides(color = guide_legend(ncol=1,override.aes = list(size = 4,alpha=1))))

(world.map.nl <- drawWorld() +
    geom_point(data = final.nl, aes(x = Longitude, y = Latitude, color=LoNL),alpha = 0.7)+ 
    scale_colour_manual(values=colprop) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(size = guide_legend(title = "Number of years monitored"),
           color = guide_legend(ncol=1,override.aes = list(size = 4,alpha=1))) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.justification = "top"))

# ggsave(world.map.nl, filename = "world.map.nl.pdf", device = "pdf", dpi = 300,
#        height = 15, width = 15)

(maxiplot<-cowplot::plot_grid(world.map.hab,
                              world.map.tax,
                              world.map.reg,
                              world.map.realm,
                              world.map.iucn,
                              world.map.nl,
                              ncol=2, labels=c("A","B","C","D","E","F")))

ggsave("outputs/maps.png",maxiplot,height = 12,width=20)

### TS1.1 | CROSS-DISTRIBUTION HABITAT TYPES AND REGIONS
table(final.nl$Region,final.nl$System)

################################################################################
# SM2 | Impact of the duration, number of years samples, and start --------

final.nl$points <- as.numeric(final.nl$points)
final.nl$points5cl <- cut(final.nl$points, seq(15,70,by=5))
final.nl$duration5cl <- cut(final.nl$duration, seq(15,70,by=5))

final.nl <- final.nl %>% 
  mutate(LoNL = ifelse(acceleration=="constant","linear","non linear")) 

ggplot(final.nl,aes(x=first_X_value,y=points))+geom_point()

p1supp <- ggplot(test_firstX,aes(x=first_X_value,y=prop))+geom_line()+
  geom_smooth(method="lm",color="black")+ theme_classic()+
  xlab("Starting year of time series") + ylab("Proportion of non-linearity") +
  theme(text = element_text(family="Montserrat",face="bold"),legend.position = "bottom")

p2supp <- ggplot(test_points,aes(x=points,y=prop))+geom_line()+
  geom_smooth(method="lm",color="black")+ theme_classic()+
  xlab("Number of points") + ylab("Proportion of non-linearity") +
  theme(text = element_text(family="Montserrat",face="bold"),legend.position = "bottom")

psupp <- cowplot::plot_grid(p1supp,p2supp,ncol=2,labels=c("A","B"))

p3supp <- ggplot(final.nl, aes(x = first_X_value, y = points5cl, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, gradient_lwd = 1, bandwidth=2) +
  scale_fill_viridis(name = "first_X_value", option = "C") + theme_classic() +
  theme(legend.position = "none", text = element_text(family="Monserrat",face="bold",size=12))+
  xlab("Starting year of the time series")+ylab("Number of points within the time series")

(supp<-cowplot::plot_grid(psupp,p3supp,ncol=1,labels=c("","C"),rel_heights = c(1,1.5)))
ggsave("outputs/side_analysis.png",supp,height = 10,width=10)

################################################################################
# SM3 | Detailed analysis of non-linearity among biogeographic ------------

### Donuts Taxonomic groups
taxa_df <- as.data.frame(table(final.nl$taxa,final.nl$dir2)) %>% 
  filter(Freq>0) %>% 
  group_by(Var1) %>% 
  mutate(ntot=sum(Freq)) %>% 
  ungroup() %>% 
  mutate(perc=round((Freq/ntot)*100,2),
         LNL=ifelse((Var2 %in% c("Decrease\nnon linear","Increase\nnon linear","No trend\nnon linear")),"Non linear","Linear"),
         Var1=fct_relevel(Var1,c("Amphibians","Invertebrates","Sharks_Rays","Birds","Reptiles","Fish","Mammals")),
         Var2=fct_relevel(Var2,c("Decrease\nlinear","No trend\nlinear","Increase\nlinear",
                                 "Increase\nnon linear","No trend\nnon linear","Decrease\nnon linear")),
         hab=as.factor(paste(Var1,"\nN =",ntot)))

(taxa_donuts <- ggplot(data=taxa_df %>% filter(Var1!="Invertebrates"), 
                       aes(x=2, y=perc, group=Var2, colour=LNL, fill=Var2)) +
    geom_bar(width = 1, stat = "identity") + xlim(0.5, 2.5) +
    coord_polar("y", start=0) + 
    facet_wrap(.~ hab,ncol=2) +theme_void()+ 
    theme(legend.position = "none", 
          # legend.box = "vertical",
          text = element_text(family="Monserrat",face="bold",size=12)) +
    geom_text(aes(label = paste(perc,"%"),color=LNL),fontface="bold",
              position = position_stack(vjust = 0.5),
              show.legend = FALSE) +
    guides(fill=guide_legend(title="Type of trajectory",nrow=1)) +
    scale_color_manual(values=c("black","white")) +
    scale_fill_manual(values=colstraj))

### Donuts Habitat types
hab_df <- as.data.frame(table(final.nl$System,final.nl$dir2)) %>% 
  filter(Freq>0) %>% 
  group_by(Var1) %>% 
  mutate(ntot=sum(Freq)) %>% 
  ungroup() %>% 
  mutate(perc=round((Freq/ntot)*100,2),
         LNL=ifelse((Var2 %in% c("Decrease\nnon linear","Increase\nnon linear","No trend\nnon linear")),"Non linear","Linear"),
         Var2=fct_relevel(Var2,c("Decrease\nlinear","No trend\nlinear","Increase\nlinear",
                                 "Increase\nnon linear","No trend\nnon linear","Decrease\nnon linear")),
         hab=as.factor(paste(Var1,"\nN =",ntot)))

(hab_donuts <- ggplot(data=hab_df, aes(x=2, y=perc, group=Var2, colour=LNL, fill=Var2)) +
    geom_bar(width = 1, stat = "identity") + xlim(0.5, 2.5) +
    coord_polar("y", start=0) + 
    facet_wrap(.~ hab,nrow=1) +theme_void()+ 
    theme(legend.position = "none", text = element_text(family="Monserrat",face="bold",size=12)) +
    geom_text(aes(label = paste(perc,"%"),color=LNL),fontface="bold",
              position = position_stack(vjust = 0.5),
              show.legend = FALSE) +
    scale_color_manual(values=c("black","white")) +
    scale_fill_manual(values=colstraj))

### Donuts Regions
reg_df <- as.data.frame(table(final.nl$Region,final.nl$dir2)) %>% 
  filter(Freq>0) %>% 
  group_by(Var1) %>% 
  mutate(ntot=sum(Freq)) %>% 
  ungroup() %>% 
  mutate(perc=round((Freq/ntot)*100,2),
         LNL=ifelse((Var2 %in% c("Decrease\nnon linear","Increase\nnon linear","No trend\nnon linear")),"Non linear","Linear"),
         Var2=fct_relevel(Var2,c("Decrease\nlinear","No trend\nlinear","Increase\nlinear",
                                 "Increase\nnon linear","No trend\nnon linear","Decrease\nnon linear")),
         hab=as.factor(paste(Var1,"\nN =",ntot)))

(regions_donuts <- ggplot(data=reg_df, aes(x=2, y=perc, group=Var2, colour=LNL, fill=Var2)) +
    geom_bar(width = 1, stat = "identity") + xlim(0.5, 2.5) +
    coord_polar("y", start=0) + 
    facet_wrap(.~ hab,nrow=2) +theme_void()+ 
    theme(legend.position = "bottom", legend.box = "vertical",
          text = element_text(family="Monserrat",face="bold",size=12)) +
    geom_text(aes(label = paste(perc,"%"),color=LNL),fontface="bold",
              position = position_stack(vjust = 0.5),
              show.legend = FALSE) +
    
    guides(fill=guide_legend(title="",nrow=1),color="none") +
    scale_color_manual(values=c("black","white")) +
    scale_fill_manual(values=colstraj))

### Donuts RLC
rlc_df <- as.data.frame(table(final.nl$redlistCategory,final.nl$dir2)) %>% 
  filter(Freq>0) %>% 
  group_by(Var1) %>% 
  mutate(ntot=sum(Freq)) %>% 
  ungroup() %>% 
  mutate(perc=round((Freq/ntot)*100,2),
         LNL=ifelse((Var2 %in% c("Decrease\nnon linear","Increase\nnon linear","No trend\nnon linear")),"Non linear","Linear"),
         Var1=fct_relevel(Var1,c("Extinct in the Wild","Critically Endangered","Endangered",
                                 "Vulnerable","Near Threatened","Least Concern","Data Deficient")),
         Var2=fct_relevel(Var2,c("Decrease\nlinear","No trend\nlinear","Increase\nlinear",
                                 "Increase\nnon linear","No trend\nnon linear","Decrease\nnon linear")),
         hab=as.factor(paste(Var1,"\nN =",ntot)))

(rlc_donuts <- ggplot(data=rlc_df %>% filter(Var1!="Extinct in the Wild"), 
                      aes(x=2, y=perc, group=Var2, colour=LNL, fill=Var2)) +
    geom_bar(width = 1, stat = "identity") + xlim(0.5, 2.5) +
    coord_polar("y", start=0) + 
    facet_wrap(.~ hab,ncol=2) +theme_void()+ 
    theme(legend.position = "none", legend.box = "vertical",
          text = element_text(family="Monserrat",face="bold",size=12)) +
    geom_text(aes(label = paste(perc,"%"),color=LNL),fontface="bold",
              position = position_stack(vjust = 0.5),
              show.legend = FALSE) +
    guides(fill=guide_legend(title="",nrow=1),color="none") +
    scale_color_manual(values=c("black","white")) +
    scale_fill_manual(values=colstraj))

(codoreg <- cowplot::plot_grid(hab_donuts,regions_donuts,labels=c("A","B"),ncol=1,rel_heights = c(0.75,2)))
(codotax <- cowplot::plot_grid(taxa_donuts,rlc_donuts,labels=c("A","B"),ncol=2))
ggsave("outputs/donuts_reg.png",codoreg,height=10,width=10)
ggsave("outputs/donuts_tax.png",codotax,height=10,width=12)

################################################################################
# SM4 | Detailed analysis of populations' temporal variability ------------

### Temporal variability in population change according to trajectory types
(rain_MSE_shape <- ggplot(data = final.nl %>% filter(MSE<15), aes(x = forcats::fct_reorder(dir2, MSE, mean), y = MSE, fill = dir2)) +
   geom_flat_violin(aes(color=dir2), position = position_nudge(x = 0.2, y = 0),scale="width") +
   geom_point(aes(y = MSE, color = dir2),position = position_jitter(width = .15), size = .5, alpha = 0.2) +
   geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.7) +
   #geom_pointrange()+
   xlab(element_blank()) + coord_flip() +
   scale_colour_manual(values = cols2) +
   scale_fill_manual(values = cols2) +
   #stat_summary(aes(fill=dir2), fun=mean, geom="point", shape=16, size=3) +
   geom_text(data = shape.letMSE_df, aes(label = letters, x=dir2, y = 12.5 ), size=4, family="Montserrat",fontface="italic",colour="black", size=5)  +
   guides(fill=FALSE, colour=FALSE) +
   geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
   theme_classic() + ylab("MSE") +
   theme(legend.position = "none",text = element_text(family="Monserrat"),
         axis.line.y = element_blank(),axis.ticks.y = element_blank(),
         axis.text.y=element_text(face="bold")))

(rain_CV_shape <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(dir2, CV, mean), y = CV, fill = dir2)) +
    geom_flat_violin(aes(color=dir2), position = position_nudge(x = 0.2, y = 0),scale="width") +
    geom_point(aes(y = CV, color = dir2),position = position_jitter(width = .15), size = .5, alpha = 0.2) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.7) +
    #geom_pointrange()+
    xlab(element_blank()) + coord_flip() +
    scale_colour_manual(values = cols2) +
    scale_fill_manual(values = cols2) +
    #stat_summary(aes(fill=dir2), fun=mean, geom="point", shape=16, size=3) +
    geom_text(data = shape.letCV_df, aes(label = letters, x=dir2, y = 7), size=4, family="Montserrat",fontface="italic",colour="black", size=5)  +
    guides(fill=FALSE, colour=FALSE) +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    theme_classic() + ylab("CV") +
    theme(legend.position = "none",text = element_text(family="Monserrat"),
          axis.line.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text.y=element_text(face="bold")))

(covar<-cowplot::plot_grid(rain_MSE_shape,rain_CV_shape,labels=c("A","B")))
ggsave("outputs/rain_variability_supp.png",covar,height=6,width=10)

### Temporal variability in population change according to trajectory types
(var_syst_shape_D <- ggplot(final.nl,aes(x=System,y=D,color=System,fill=System))+
  geom_jitter(size=0.5,width=0.15,alpha=0.3)+
  #geom_boxplot(width=0.1,outlier.shape = NA)+
  geom_flat_violin(position = position_nudge(x=0.2,y=0),scale="width")+
  scale_color_manual(values=colssyst)+
  scale_fill_manual(values=colssyst)+
  theme_classic()+theme(axis.line.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")+coord_flip()+
  facet_wrap(~dir2))

(var_reg_shape_D <- ggplot(final.nl,aes(x=Region,y=D,color=Region,fill=Region))+
  geom_jitter(size=0.5,width=0.15,alpha=0.3)+
  #geom_boxplot(width=0.1,outlier.shape = NA)+
  geom_flat_violin(position = position_nudge(x=0.2,y=0),scale="width")+
  scale_color_manual(values=colsregion)+
  scale_fill_manual(values=colsregion)+
  theme_classic()+theme(axis.line.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")+coord_flip()+
  facet_wrap(~dir2))

ggsave("outputs/var_syst_shape.png",var_syst_shape_D,height = 6,width = 10)
ggsave("outputs/var_reg_shape.png",var_reg_shape_D,height = 6,width = 10)


################################################################################
#### Non linear proportion and effect sizes of temporal variability among realms

# Linear vs Non linear Realms
dframe_realm <- final.nl %>% 
  mutate(LoNL = ifelse(acceleration=="constant","linear","non linear")) %>% 
  mutate(full_realm = forcats::fct_reorder(.f = full_realm,
                                           .x = LoNL,
                                           .fun = function(.x) mean(.x == "linear"),
                                           .desc = TRUE))

prop_realm <- ggplot(dframe_realm, aes(x = full_realm, fill = LoNL))
(prop_realm <- plot_prop(prop_realm,colprop,line=T)+ggtitle("Realm")+theme(plot.title = element_text(hjust=0.5)))

full.realm.summary$realm<-as.factor(full.realm.summary$realm)
full.realm.summary$realm<-fct_relevel(full.realm.summary$realm,rev(levels(dframe_realm$full_realm)))
(var_plot_realm <- ggplot(full.realm.summary, aes(x = realm, y = mean)) +
    geom_pointrange(aes(ymin = CIinf, ymax = CIsup, colour = variability,
                        shape = variability),
                    position = position_dodge(0.7), size = 0.5) + 
    scale_colour_manual(values = c("grey70", "grey40", "black")) +
    theme_classic() +
    theme(axis.title.y = element_text(family="Montserrat",face="bold"),
          axis.text.x = element_text(family="Monserrat",angle=40,hjust=1,face="bold"),
          axis.text.y = element_text(family="Monserrat"),
          legend.text = element_text(family="Monserrat", face = "italic"),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.06, 0.8),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Effect size"))

(coplo_realm_nlvar <- cowplot::plot_grid(prop_realm,var_plot_realm,ncol=1,labels=c("A","B")))
ggsave("outputs/nl_var_realm.png",coplo_realm_nlvar,height=10,width=12)

################################################################################
### OLDER PLOTS

# # Temporal variability ----------------------------------------------------
# 
# # Habitat type
# # rain_hab <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(System, D, mean), y = D))
# # (rain_D_hab <- plot_rain(rain_hab) +
# #     geom_text(data = hab.letters_df, aes(label = letters, x=System, y = 3.5 ), colour="black", size=5, family="Monserrat", fontface="italic"))
# # 
# # rain_habCV <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(System, CV, mean), y = CV))
# # (rain_CV_hab <- plot_rain(rain_habCV)+ylab("CV")+
# #     geom_text(data = hab.letCV_df, aes(label = letters, x=System, y = 6 ), colour="black", size=5, family="Monserrat", fontface="italic"))
# 
# rain_habMSE <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(System, MSE, mean), y = MSE))
# (rain_MSE_hab <- plot_rain(rain_habMSE)+ylab("MSE")+
#     geom_text(data = hab.letMSE_df, aes(label = letters, x=System, y = 15), 
#               colour="black", size=5, family="Monserrat", fontface="italic"))
# 
# 
# rain_habD <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(System, D, mean), y = D))
# (rain_D_hab <- plot_rain(rain_habD)+ylab("D"))
# rain_habD <- ggplot(data = no.trend.all, aes(x = forcats::fct_reorder(System, D, mean), y = D))
# (rain_D_hab <- plot_rain(rain_habD)+ylab("D"))
# rain_habD <- ggplot(data = inc.nl.all, aes(x = forcats::fct_reorder(System, D, mean), y = D))
# (rain_D_hab <- plot_rain(rain_habD)+ylab("D"))
# 
# rain_regD <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(Region, D, mean), y = D))
# (rain_D_reg <- plot_rain(rain_regD)+ylab("D"))
# rain_regD <- ggplot(data = no.trend.all, aes(x = forcats::fct_reorder(Region, D, mean), y = D))
# (rain_D_reg <- plot_rain(rain_regD)+ylab("D"))
# 
# # Taxa
# # rain_tax <- ggplot(data = final.nl %>% filter(taxa!="Invertebrates"), aes(x = forcats::fct_reorder(taxa, D, mean), y = D)) 
# # (rain_D_tax <- plot_rain(rain_tax) +
# #     geom_text(data = tax.letters_df, aes(label = letters, x=taxa, y = 3.5 ), colour="black", size=5, family="Monserrat", fontface="italic"))
# # 
# # rain_taxCV <- ggplot(data = final.nl  %>% filter(taxa!="Invertebrates"), aes(x = forcats::fct_reorder(taxa, CV, mean), y = CV))
# # (rain_CV_tax <- plot_rain(rain_taxCV)+ylab("CV"))
# 
# rain_taxMSE <- ggplot(data = final.nl %>% filter(taxa!="Invertebrates"), aes(x = forcats::fct_reorder(taxa, MSE, mean), y = MSE))
# (rain_MSE_tax <- plot_rain(rain_taxMSE)+ylab("MSE")+
#     geom_text(data = tax.letMSE_df %>% filter(taxa!="Invertebrates"), aes(label = letters, x=taxa, y = 15), 
#               colour="black", size=5, family="Monserrat", fontface="italic"))
# 
# # Regions
# rain_reg <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(Region, D, mean), y = D)) 
# (rain_D_reg <- plot_rain(rain_reg) +
#     geom_text(data = reg.letters_df, aes(label = letters, x=Region, y = 3.5 ), colour="black", size=5, family="Monserrat", fontface="italic"))
# 
# rain_regCV <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(Region, CV, mean), y = CV))
# (rain_CV_reg <- plot_rain(rain_regCV)+ylab("CV"))
# 
# rain_regMSE <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(Region, CV, mean), y = MSE))
# (rain_MSE_reg <- plot_rain(rain_regMSE)+ylab("MSE")+theme(axis.text.y = element_blank()))
# 
# # Realms
# rain_realm <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(full_realm, D, mean), y = D)) 
# (rain_D_realm <- plot_rain(rain_realm) +
#     geom_text(data = realm.letters_df, aes(label = letters, x=full_realm, y = 3.5 ), colour="black", size=5, family="Monserrat", fontface="italic"))
# 
# rain_realmCV <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(full_realm, CV, mean), y = CV))
# (rain_CV_realm <- plot_rain(rain_realmCV)+ylab("CV"))
# 
# rain_realmMSE <- ggplot(data = final.nl, aes(x = forcats::fct_reorder(full_realm, CV, mean), y = MSE))
# (rain_MSE_realm <- plot_rain(rain_realmMSE)+ylab("MSE")+theme(axis.text.y = element_blank()))
# 
# # IUCN Red List Categories
# data.rlc$redlistCategory<-factor(data.rlc$redlistCategory,
#                                    levels=c("Extinct in the Wild",
#                                             "Critically Endangered",
#                                             "Endangered",
#                                             "Vulnerable",
#                                             "Near Threatened",
#                                             "Least Concern",
#                                             "Data Deficient"))
# 
# rain_rlc <- ggplot(data = data.rlc, aes(x = redlistCategory, y = D))
# (rain_D_rlc <- plot_rain(rain_rlc) +
#   geom_text(data = rlc.letters_df, aes(label = letters, x=redlistCategory, y = 3.5 ), colour="black", size=5, family="Monserrat", fontface="italic"))
# 
# 
# rain_rlcCV <- ggplot(data = data.rlc, aes(x = redlistCategory, y = CV))
# (rain_CV_rlc <- plot_rain(rain_rlcCV)+ylab("CV"))
# 
# rain_rlcMSE <- ggplot(data = data.rlc, aes(x = redlistCategory, y = MSE))
# (rain_MSE_rlc <- plot_rain(rain_rlcMSE)+ylab("MSE")+theme(axis.text.y = element_blank()))
# 
# # Final
# rain_D1 <- cowplot::plot_grid(NULL,rain_D_hab,rel_widths = c(0.08,1))
# rain_D2 <- cowplot::plot_grid(NULL,rain_D_tax,rel_widths = c(0.08,1))
# (rain_D <- cowplot::plot_grid(rain_D1,rain_D2,rain_D_rlc,labels=c("A","B","C"),ncol=1, rel_heights = c(0.5,1,1)))
# ggsave("temporal_variability.png", rain_D, height = 10, width = 12)
# 
# rain_Dsupp1 <- cowplot::plot_grid(NULL,rain_D_reg,rel_widths = c(0.08,1))
# (rain_D_supp<-cowplot::plot_grid(rain_Dsupp1 ,rain_D_realm,ncol=1,labels=c("A","B"),rel_heights = c(0.8,1)))
# ggsave("supp_variability.png",rain_D_supp,height = 10, width = 12)
# 
# rain_comp1 <- cowplot::plot_grid(rain_CV_hab,rain_MSE_hab,rel_widths = c(1,0.8))
# rain_comp2 <- cowplot::plot_grid(rain_CV_tax,rain_MSE_tax,rel_widths = c(1,0.8))
# rain_comp3 <- cowplot::plot_grid(rain_CV_reg,rain_MSE_reg,rel_widths = c(1,0.8))
# rain_comp4 <- cowplot::plot_grid(rain_CV_realm,rain_MSE_realm,rel_widths = c(1,0.8))
# rain_comp5 <- cowplot::plot_grid(rain_CV_rlc,rain_MSE_rlc,rel_widths = c(1,0.8))
# 
# ggsave("comp_hab.png", rain_comp1, height = 4, width = 15)
# ggsave("comp_tax.png", rain_comp2, height = 6, width = 15)
# ggsave("comp_reg.png", rain_comp3, height = 8, width = 15)
# ggsave("comp_realm.png", rain_comp4, height = 10, width = 15)
# ggsave("comp_rlc.png", rain_comp5, height = 6, width = 15)
# 
# 
# 
# #ggsave("outputs/comparisonDDZ.png",coDDZ,height = 6,width=12)



################################################################################


################################################################################
# # Comparisons D, CV and MSE
# res<-cor(final.nl[,c("D","CV","MSE","bres")])
# col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# corrplot::corrplot(res, method="color", col=col(200),  
#                    type="upper", order="hclust", 
#                    addCoef.col = "black", # Add coefficient of correlation
#                    tl.col="black", tl.srt=45, #Text label color and rotation
#                    # Combine with significance
#                    p.mat = p.mat, sig.level = 0.01, insig = "blank", 
#                    # hide correlation coefficient on the principal diagonal
#                    diag=FALSE 
# )
# corrplot::corrplot(res, method="color", col=col(200),
#                    type = "upper", order = "hclust", 
#                    addCoef.col = "black",
#                    p.mat = p.mat, sig.level = 0.01, insig = "blank", 
#                    tl.col = "black",diag=FALSE)
# 
# pairs(final.nl[,c("D","CV","MSE")],
#       upper.panel = panel.cor,    # Correlation panel
#       lower.panel = panel.smooth) # Smoothed regression lines
# 
# 
# ggplot(final.nl,aes(x=D,y=CV,color=LoNL))+
#   scale_color_manual(values=cols3) + 
#   geom_point()+geom_abline(slope=1) + theme_minimal()+geom_smooth(method="lm",linetype="dashed")+
#   facet_wrap(~direction)
# 
# ggplot(final.nl,aes(x=D,y=MSE,color=LoNL))+
#   scale_color_manual(values=cols3) + 
#   geom_point()+geom_abline(slope=1) + theme_minimal()+geom_smooth(method="lm",linetype="dashed")+
#   facet_wrap(~direction)
# 
# ggplot(final.nl,aes(x=CV,y=MSE,color=LoNL))+
#   scale_color_manual(values=cols3) + 
#   geom_point()+geom_abline(slope=1) + theme_minimal()+geom_smooth(method="lm",linetype="dashed")+
#   facet_wrap(~direction)
# 
# dframe_taxa[is.na(dframe_taxa$derivative),]$derivative<-0
# dframe_rlc[is.na(dframe_rlc$derivative),]$derivative<-0
# final.nl.sc[is.na(final.nl.sc$derivative),]$derivative<-0
# final.nl[is.na(final.nl$derivative),]$derivative<-0
# 
# 
# MSE_taxa <- ggplot(dframe_taxa)+
#   geom_point(aes(x=MSE,y=derivative,color=taxa))+theme_bw()+
#   scale_color_manual(values=colstaxa)+
#   geom_abline(slope=0,intercept = 0)+ theme_bw()+
#   # geom_abline(slope=0,intercept = 0.1, linetype="dashed") + 
#   # geom_abline(slope=0,intercept = -0.1, linetype="dashed") +
#   geom_smooth(aes(x=MSE,y=derivative),formula=y~poly(x,2),color="grey30")+
#   facet_wrap(~dir2,scales="free")
# 
# D_taxa<-ggplot(dframe_taxa)+
#   geom_point(aes(x=D,y=derivative,color=taxa))+theme_bw()+
#   scale_color_manual(values=colstaxa)+
#   geom_abline(slope=0,intercept = 0)+ theme_bw()+
#   # geom_abline(slope=0,intercept = 0.1, linetype="dashed") + 
#   # geom_abline(slope=0,intercept = -0.1, linetype="dashed") +
#   geom_smooth(aes(x=D,y=derivative),formula=y~poly(x,2),color="grey30")+
#   facet_wrap(~dir2,scales="free")
# 
# cowplot::plot_grid(MSE_taxa,D_taxa,nrow=2)
# 
# ggplot(dframe_rlc %>% filter(redlistCategory!="Extinct in the Wild" & redlistCategory!="NA"))+
#   geom_point(aes(x=MSE,y=derivative,color=dir2))+theme_bw()+
#   scale_color_manual(values=cols2)+
#   geom_abline(slope=0,intercept = 0)+ theme_bw()+
#   # geom_abline(slope=0,intercept = 0.1, linetype="dashed") + 
#   # geom_abline(slope=0,intercept = -0.1, linetype="dashed") +
#   geom_smooth(aes(x=MSE,y=derivative),method="lm",color="grey30")+
#   facet_wrap(~redlistCategory,scales="free")
# 
# hist(final.nl$MSE,breaks=60)
# range(final.nl$MSE)
# 
# essai <- final.nl  %>% 
#   group_by(direction,LoNL) %>% 
#   summarize(mean_D=mean(D),mean_MSE=mean(MSE),mean_CV=mean(CV))
# 
# A<-ggplot(essai, aes(direction, LoNL, fill= mean_D)) + 
#   geom_tile() + theme_minimal() +
#   scale_fill_viridis(option="magma")
# 
# B<-ggplot(essai, aes(direction, LoNL, fill= mean_MSE)) + 
#   geom_tile() + theme_minimal() +
#   scale_fill_viridis(option="magma")
# 
# C<-ggplot(essai, aes(direction, LoNL, fill= mean_CV)) + 
#   geom_tile() + theme_minimal() +
#   scale_fill_viridis(option="magma")
# 
# cowplot::plot_grid(A,B,C,ncol=1)