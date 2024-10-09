################################################################################
library(dplyr)
library(haven)
library(SIMR)
library(ggplot2)
library(reshape2)
library(lme4)
library(misty)
library(labelled)
library(see)
library(gtsummary)
library(ggstatsplot)
library(tidyr)
library(corrplot)
library(ggalluvial)
library(labelled)
library(interactions)
library(margins)
library(prediction)
library(jtools)
library(fastDummies)
library(emmeans)
library(sjPlot)
library(pwr)
library(WebPower)
library(labelled)
#New function for violin plot
#GEOM VIOLIN SPLIT
geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  GeomSplitViolin <- ggplot2::ggproto(
    "GeomSplitViolin", 
    GeomViolin, 
    draw_group = function(self, data, ..., draw_quantiles = NULL) {
      data <- transform(data, 
                        xminv = x - violinwidth * (x - xmin), 
                        xmaxv = x + violinwidth * (xmax - x))
      grp <- data[1,'group']
      newdata <- plyr::arrange(
        transform(data, x = if(grp%%2==1) xminv else xmaxv), 
        if(grp%%2==1) y else -y
      )
      newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
      newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
      if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
        quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
        aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
        ggplot2:::ggname("geom_split_violin", 
                         grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...), quantile_grob))
      } else {
        ggplot2:::ggname("geom_split_violin", ggplot2::GeomPolygon$draw_panel(newdata, ...))
      }
    }
  )
  
  
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  # define flat violin geom
  GeomFlatViolin <- ggplot2::ggproto(
    "Violinist", 
    Geom,
    setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||% (ggplot2::resolution(data$x, FALSE) * 0.9)
      
      # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
      data %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(ymin = min(y),
                      ymax = max(y),
                      xmin = x,
                      xmax = x + width / 2) %>%
        dplyr::ungroup()
      
    },
    draw_group = function(data, panel_scales, coord) {
      # Find the points for the line to go all the way around
      data <- transform(data, xminv = x,
                        xmaxv = x + violinwidth * (xmax - x))
      
      # Make sure it's sorted properly to draw the outline
      newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                       plyr::arrange(transform(data, x = xmaxv), -y))
      
      # Close the polygon: set first and last point the same
      # Needed for coord_polar and such
      newdata <- rbind(newdata, newdata[1,])
      
      ggplot2:::ggname("geom_flat_violin", 
                       ggplot2::GeomPolygon$draw_panel(newdata, panel_scales, coord))
    },
    draw_key = draw_key_polygon,
    default_aes = aes(weight = 1, colour = "grey20", 
                      fill = "white", size = 0.5,
                      alpha = NA, linetype = "solid"),
    required_aes = c("x", "y")
  )
  
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

################################################################################
# 1)                                                                           
# The survey was administered to 3652 undergraduate biology students.                                                                              
# Students with no consent and not check honest questions were removed.    
# Provided data shows these reductions and has an n =2625   
################################################################################
# Load up the data to be used in the analyses
data <- read.csv("C:/Users/rqa2a/OneDrive - Middle Tennessee State University/Research - ReCCEE Projects/ReCCEE randomize controlled trial/Analysis (Fall 2022 & Spring 2023)/24 04 18 randomize study data.csv")

data$religion2<-to_factor(data$religion2)
data$religion4 <- as.factor(data$religion4)
data$religion4 <- recode(data$religion4, '1' = "Christian", '2'= "Other", '3'= "Agnostic", '4'="Atheist")
data$religion4 <- relevel (data$religion4, ref =  "Christian")

data$cond_cc<-as.factor(data$cond_cc)
data$cond_cc<- relevel (data$cond_cc, ref =  "Control")
data$cond<-as.factor(data$cond)
data$cond<- relevel (data$cond, ref =  "Christian Instructor")
data$relcat<-as.factor(data$relcat)
data$class<-as.factor(data$class)

################################################################################
#Variable naming information 
#cond: conditions between control,  Christian Instructor, and Non-Religious Instructor 
#cond_cc: conditions control and cultural competence
#religion4: religious affiliation broken down by Christian (n=1325), Agnostic (n=598), Atheist (n=235), and Other (n=467)
#religion2: religious affiliation broken down by Christian, and Non-Christian
#class: Course informaion
#relcat: Religiosity level broken down by High, Moderate, and Low
#rel: Average scores religiosity level (1-5 point)
#pconf: Average PRE-score of perceived conflict (1-5 point)
#ppconf: Average POST-score of perceived conflict (1-5 point)
#pcomp: Average PRE-score of perceived compatibility (1-5 point)
#ppcomp: Average POST-score of perceived compabicompatibilitylity (1-5 point)
################################################################################
# Power analyses
#Power calculation for two proportions (different sample sizes) 
#How much power is to detect small effect sizes? #all sample
#effect sizes of 0.2 (small), 0.5 (medium), or 0.8 (large). 
pwr.2p2n.test(h = 0.2 , n1 =1325 , n2 =235 , sig.level = 0.05, power = ) #christ #atheist #0.80
pwr.2p2n.test(h = 0.2, n1 =1325 , n2 =465 , sig.level = 0.05, power = ) #christ #other #0.95
pwr.2p2n.test(h = 0.2, n1 =1325 , n2 =598 , sig.level = 0.05, power = ) #christ agnostic#0.98

#How much power is to detect small effect sizes? #only cultural competence condition
data%>%filter(!cond=="Control") %>% select(religion4) %>% table()
pwr.2p2n.test(h = 0.2 , n1 =865 , n2 =157 , sig.level = 0.05, power = ) #christ #atheist #0.63
pwr.2p2n.test(h = 0.2 , n1 =405 , n2 =157 , sig.level = 0.05, power = ) #agnostic #atheist #0.56
pwr.2p2n.test(h = 0.2, n1 =865 , n2 =308 , sig.level = 0.05, power = ) #christ #other #0.85
pwr.2p2n.test(h = 0.2, n1 =865 , n2 =405 , sig.level = 0.05, power = ) #christ agnostic
#0.91

#check the sample in each class
class<-as.data.frame(table(data$class)) 
summary(class$Freq)
table(data$class) #we have observation with missing clas

#delete sample with missing class
data<-data%>%filter(!class=="") #n=2625 #number of class/cluster=17
table(data$class) #now 20 is the smallest sample


#chrombha alpha
cronbach.alpha(data[c("pcomp1", "pcomp2", "pcomp3")], na.rm = TRUE) #alpha: 0.73
cronbach.alpha(data[c("pconf1", "pconf3", "pconf3.0")], na.rm = TRUE) #alpha: 0.66

cronbach.alpha(data[c("mate1", "mate2", "mate3", "mate4", "mate5", "mate6", "mate7", "mate8","mate9")], na.rm = TRUE) #alpha: 0.901
cronbach.alpha(data[c("human1", "human2", "human3", "human4", "human5", "human6", "human7", "human8")], na.rm = TRUE) #alpha: 0.9
cronbach.alpha(data[c("evound1", "evound2", "evound3", "evound4", "evound5", "evound6", "evound7", "evound8", "evound9", "evound10", "evound11", "evound12", "evound13")], na.rm = TRUE) #alpha: 0.583



#To calculate the statistical power given sample size and effect size
# how many sample sizes within each cluster to detect large effect size
# n = sample size within each cluster
#f = effect size
# J = number of clusters
wp.crt2arm(f = 0.40, n = NULL, J =17, icc = 0.02, alpha = 0.05, power = 0.8)
#n=17



################################################################################
# 2)  STRUCTURE VALIDITY 
pre<-data[c("pcomp1", "pcomp2", "pcomp3", "pconf1", "pconf3", "pconf3.0")]
post<-data[c("ppcomp1", "ppcomp2", "ppcomp3", "ppconf1", "ppconf2", "ppconf3")]
colnames(pre)<-colnames(post)
prepost<-rbind(pre,post)
mod.2 <- ' compatibility =~ ppcomp1 + ppcomp2 + ppcomp3
conflict =~ ppconf1 + ppconf2 + ppconf3'
mod2.fit <- cfa(mod.2, data=prepost, std.lv=TRUE, estimator = "WLS")
summary(mod2.fit, standardized=TRUE, fit.measures=TRUE)
fitMeasures(mod2.fit, c("cfi", "rmsea", "srmr"))
#cfi rmsea 
#0.868 0.095 

################################################################################
#3) 
#need to determine whether or not including a random effect structure is justified
#by comparing the deff (design effect) value 
#The design effect quantifies the extent to which the expected sampling error 
#in a survey departs from the sampling error that can be expected under simple 
#random sampling.
#The design effect can be equivalent defined as the the actual sample size divided 
#by the effective sample size. Thus, where the true sampling variance is twice 
#that computed under the assumption of simple random sampling the design effect is 2.0.
#source: https://docs.displayr.com/wiki/Design_Effects_and_Effective_Sample_Size#The_design_effect_.28deff.29

library(Hmisc)
round(deff(data$ppcomp, data$class), 2)
#       n clusters      rho     deff 
# 2620.00    17.00     0.01     2.51 
#the deff is 2.51, which means that the variance of our estimate is 2.51 times 
#larger than what it would be if I had used a simple random sample.
#This suggests that there is a considerable amount of clustering in our data, 
#and that the variable class may indeed be a significant source of this clustering.
round(deff(data.rq2$ppcomp, data.rq2$class), 2)
#n clusters      rho     deff 
#1732.00    17.00     0.01     2.05 

round(deff(data$ppconf, data$class), 2)
#n clusters      rho     deff 
#2620.00    17.00     0.00     1.74 
round(deff(data.rq2$ppconf, data.rq2$class), 2)
#n clusters      rho     deff 
#1732.0     17.0      0.0      0.9  
#QUESTION for Dr. Yi, should I not include the class as random effect for this variable then?

################################################################################
# ANALYSES 
data.rq2<-data%>%filter(!cond=="Control") #n=1733 (without control, for RQ2)

#RQ 1. Does religious cultural competence in evolution instruction lead to better 
#student outcomes than evolution instruction without religious cultural competence?

#Compatibility
pcomp1 <- lmer(ppcomp ~ pcomp + cond_cc + (1 | class), data = data)
plot_model(pcomp1, show.values = TRUE, value.offset = .5, title = "RQ1: Does cultural competence work?")  + theme_sjplot()
summ(pcomp1) #TABLE 1

library(dplyr)
ppcomp.g <- data[c("ppcomp", "pcomp", "cond_cc")] %>%
  pivot_longer(cols = c(ppcomp, pcomp), names_to = "Pre.post", values_to = "Score") %>%
  mutate(
    PrePost = if_else(Pre.post == "pcomp", "Pre", "Post"),
    Variable = Pre.post,
    PrePost = factor(PrePost, levels = c("Pre", "Post"))
  )
#FIGURE 2 B
ggplot(ppcomp.g, aes(cond_cc, Score, fill = PrePost))+ geom_split_violin(alpha=0.7) +
  geom_boxplot(aes(fill=PrePost),fatten=NULL, position = position_dodge(width = 0.1),width = 0.1, alpha=0.7) +
  stat_summary(fun=mean, geom="crossbar", aes(group = factor(PrePost,level=c("Pre","Post"))),position=position_dodge(.1), width=.1) + theme_bw() +
  labs(y= "Perceived Compatibility", x="Conditions") +
  scale_fill_manual(values=c('#F3DC76', '#F4AE4C'), name = "Conditions")

#Conflict
pconf1 <- lmer(ppconf ~ pconf + cond_cc + (1 | class), data = data)
plot_model(pconf1, show.values = TRUE, value.offset = .5, title = "RQ1: Does cultural competence work?")  + theme_sjplot()
summ(pconf1) #TABLE 1

ppconf.g <- data %>% select(ppconf, pconf, cond_cc) %>%
  pivot_longer(cols = c(ppconf, pconf), names_to = "Pre.post", values_to = "Score") %>%
  mutate(
    PrePost = if_else(Pre.post == "pconf", "Pre", "Post"),
    Variable = Pre.post,
    PrePost = factor(PrePost, levels = c("Pre", "Post"))
  )
#FIGURE 2 B
ggplot(ppconf.g, aes(cond_cc, Score, fill = PrePost))+ geom_split_violin(alpha=0.7) +
  geom_boxplot(aes(fill=PrePost),fatten=NULL, position = position_dodge(width = 0.1),width = 0.1, alpha=0.7) +
  stat_summary(fun=mean, geom="crossbar", aes(group = factor(PrePost,level=c("Pre","Post"))),position=position_dodge(.1), width=.1) + theme_bw() +
  labs(y= "Perceived conflict", x="Conditions") +
  scale_fill_manual(values=c('#F3DC76', '#F4AE4C'), name = "Conditions")


#mate
mate1 <- lmer(pmate ~ mate + cond_cc + (1 | class), data = data)
plot_model(mate1, show.values = TRUE, value.offset = .5, title = "RQ1: Does cultural competence work?")  + theme_sjplot()
summ(mate1) #TABLE 1

pmate.g <- data %>% select(pmate, mate, cond_cc) %>%
  pivot_longer(cols = c(pmate, mate), names_to = "Pre.post", values_to = "Score") %>%
  mutate(
    PrePost = if_else(Pre.post == "mate", "Pre", "Post"),
    Variable = Pre.post,
    PrePost = factor(PrePost, levels = c("Pre", "Post"))
  )
#FIGURE 2 B
ggplot(pmate.g, aes(cond_cc, Score, fill = PrePost))+ geom_split_violin(alpha=0.7) +
  geom_boxplot(aes(fill=PrePost),fatten=NULL, position = position_dodge(width = 0.1),width = 0.1, alpha=0.7) +
  stat_summary(fun=mean, geom="crossbar", aes(group = factor(PrePost,level=c("Pre","Post"))),position=position_dodge(.1), width=.1) + theme_bw() +
  labs(y= "General Evolution Acceptance", x="Conditions") +
  scale_fill_manual(values=c('#F3DC76', '#F4AE4C'), name = "Conditions")

#human
human1 <- lmer(phuman ~ human + cond_cc + (1 | class), data = data)
plot_model(human1, show.values = TRUE, value.offset = .5, title = "RQ1: Does cultural competence work?")  + theme_sjplot()
summ(human1) #TABLE 1

phuman.g <- data %>% select(phuman, human, cond_cc) %>%
  pivot_longer(cols = c(phuman, human), names_to = "Pre.post", values_to = "Score") %>%
  mutate(
    PrePost = if_else(Pre.post == "human", "Pre", "Post"),
    Variable = Pre.post,
    PrePost = factor(PrePost, levels = c("Pre", "Post"))
  )
#FIGURE 2 B
ggplot(phuman.g, aes(cond_cc, Score, fill = PrePost))+ geom_split_violin(alpha=0.7) +
  geom_boxplot(aes(fill=PrePost),fatten=NULL, position = position_dodge(width = 0.1),width = 0.1, alpha=0.7) +
  stat_summary(fun=mean, geom="crossbar", aes(group = factor(PrePost,level=c("Pre","Post"))),position=position_dodge(.1), width=.1) + theme_bw() +
  labs(y= "Human Evolution Acceptance", x="Conditions") +
  scale_fill_manual(values=c('#F3DC76', '#F4AE4C'), name = "Conditions")

#human
human1 <- lmer(phuman ~ human + cond_cc + (1 | class), data = data)
plot_model(human1, show.values = TRUE, value.offset = .5, title = "RQ1: Does cultural competence work?")  + theme_sjplot()
summ(human1) #TABLE 1

phuman.g <- data %>% select(phuman, human, cond_cc) %>%
  pivot_longer(cols = c(phuman, human), names_to = "Pre.post", values_to = "Score") %>%
  mutate(
    PrePost = if_else(Pre.post == "human", "Pre", "Post"),
    Variable = Pre.post,
    PrePost = factor(PrePost, levels = c("Pre", "Post"))
  )
#FIGURE 2 B
ggplot(phuman.g, aes(cond_cc, Score, fill = PrePost))+ geom_split_violin(alpha=0.7) +
  geom_boxplot(aes(fill=PrePost),fatten=NULL, position = position_dodge(width = 0.1),width = 0.1, alpha=0.7) +
  stat_summary(fun=mean, geom="crossbar", aes(group = factor(PrePost,level=c("Pre","Post"))),position=position_dodge(.1), width=.1) + theme_bw() +
  labs(y= "Human Evolution Acceptance", x="Conditions") +
  scale_fill_manual(values=c('#F3DC76', '#F4AE4C'), name = "Conditions")
################################################################################
# ANALYSES 
#RQ 2.Does an instructor revealing a secular or religious identity 
#impact the effectiveness of religious cultural competence in evolution
#education for students from different religious and non-religious identities?


#Compatibility
pcomp.rq2 <- data.rq2[c("pcomp", "ppcomp", "cond", "class", "religion4", "rel")]
pcomp.rq2<-na.omit(pcomp.rq2)

#MODEL 2
pcomp2.rq2 <- lmer(ppcomp ~ pcomp  + cond + (1 | class), data = pcomp.rq2)
plot_model(pcomp2.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(pcomp2.rq2)#Table 2

#MODEL 3
pcomp3.rq2 <- lmer(ppcomp ~ pcomp  + cond + religion4 + (1 | class), data = pcomp.rq2)
plot_model(pcomp3.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(pcomp3.rq2)#Table 2

#MODEL 4
pcomp4.rq2 <- lmer(ppcomp ~ pcomp  + cond * religion4 + (1 | class), data = pcomp.rq2)
plot_model(pcomp4.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(pcomp4.rq2)#Table 2

# Calculate estimated marginal 
em.pcomp<- emmeans(pcomp4.rq2, ~ cond*religion4)
con.pcomp<-pairs(em.pconf, by="religion4", combine=TRUE)
sum.pcomp<-summary(em.pconf)
con.pcomp #TABLE3
sum.pcomp #FOR FIGURE 3
sum.pcomp$religion4<-factor(sum.pcomp$religion4, levels=c("Atheist", "Agnostic", "Other", "Christian"))
sum.pcomp$cond<-factor(sum.pcomp$cond, levels=c("Christian Instructor", "Non-Religious Instructor"))


# FIGURE 3
ggplot(sum.pcomp, aes(x = religion4,  y = emmean, color = cond, group=cond)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean  +SE, linetype = factor(religion4)), width = 0.2, linewidth = 1, linetype = "solid") +
  #geom_line(linetype="dashed", size=1)+ 
  geom_point(size = 2, stroke = 1) +
  scale_color_manual(values = c("black", "gray"), name = "Cultural Competence with") +
  labs(x = "Students' religious identity", y = "Estimated Marginal Means in \n Perceived Compatibility") +
  theme_classic() + annotate("text",
                             label="*",
                             x=1,
                             y=4.7,
                             size=7) +
  scale_y_continuous(breaks=seq(4, 4.75, 0.25), limits =c(4, 4.75))


#Conflict 
#not included in the manuscript yet
pconf.rq2 <- data.rq2[c("pconf", "ppconf", "cond", "class", "religion4", "rel")]
pconf.rq2<-na.omit(pconf.rq2)

#MODEL 2
pconf2.rq2 <- lmer(ppconf ~ pconf  + cond + (1 | class), data = pconf.rq2)
plot_model(pconf2.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(pconf2.rq2)#Table not showed in manuscript yet

pconf3.rq2 <- lmer(ppconf ~ pconf  + cond + religion4 + (1 | class), data = pconf.rq2)
plot_model(pconf3.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(pconf4.rq2)#Table not showed in manuscript

pconf4.rq2 <- lmer(ppconf ~ pconf  + cond * religion4 + (1 | class), data = pconf.rq2)
plot_model(pconf4.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(pconf4.rq2)#Table not showed in manuscript

# Calculate estimated marginal 
em.pconf <- emmeans(pconf4.rq2, ~ cond*religion4)
con.pconf<-pairs(em.pconf, by="religion4", combine=TRUE)
sum.pconf<-summary(em.pconf)
con.pconf #TABLE3
sum.pconf #FOR FIGURE 3
sum.pconf$religion4<-factor(sum.pconf$religion4, levels=c("Atheist", "Agnostic", "Other", "Christian"))
sum.pconf$cond<-factor(sum.pconf$cond, levels=c("Christian Instructor", "Non-Religious Instructor"))


# FIGURE 3
ggplot(sum.pconf, aes(x = religion4,  y = emmean, color = cond, group=cond)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean  +SE, linetype = factor(religion4)), width = 0.2, linewidth = 1, linetype = "solid") +
  #geom_line(linetype="dashed", size=1)+ 
  geom_point(size = 2, stroke = 1) +
  scale_color_manual(values = c("black", "gray"), name = "Cultural Competence with") +
  labs(x = "Students' religious identity", y = "Estimated Marginal Means in \n Perceived Conflict") +
  theme_classic() +
  scale_y_continuous(breaks=seq(2.25, 3.50, 0.25), limits =c(2.25, 3))

##MATE
mate.rq2 <- data.rq2[c("mate", "pmate", "cond", "class", "religion4", "rel")]
mate.rq2<-na.omit(mate.rq2)

#MODEL 2
mate2.rq2 <- lmer(pmate ~ mate  + cond + (1 | class), data = mate.rq2)
plot_model(mate2.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(mate2.rq2)#Table 2

#MODEL 3
mate3.rq2 <- lmer(pmate ~ mate  + cond + religion4 + (1 | class), data = mate.rq2)
plot_model(mate3.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(mate3.rq2)#Table 2

#MODEL 4
mate4.rq2 <- lmer(pmate ~ mate  + cond * religion4 + (1 | class), data = mate.rq2)
plot_model(mate4.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(mate4.rq2)#Table 2

# Calculate estimated marginal 
em.mate <- emmeans(mate4.rq2, ~ cond*religion4)
con.mate<-pairs(em.mate, by="religion4", combine=TRUE)
sum.mate<-summary(em.mate)
con.mate #TABLE3
sum.mate #FOR FIGURE 3
sum.mate$religion4<-factor(sum.mate$religion4, levels=c("Atheist", "Agnostic", "Other", "Christian"))
sum.mate$cond<-factor(sum.mate$cond, levels=c("Christian Instructor", "Non-Religious Instructor"))


ggplot(sum.mate, aes(x = religion4,  y = emmean, color = cond, group=cond)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean  + SE, linetype = factor(religion4)), width = 0.2, linewidth = 1, linetype = "solid") +
  geom_point(size = 2, stroke = 1) +
  scale_color_manual(values = c("black", "gray"), name = "Cultural Competence with") +
  labs(x = "Students' religious identity", y = "Estimated Marginal Means in \n General Evolution Acceptance") +
  theme_classic() +  
  scale_y_continuous(breaks=seq(4, 4.5, 0.25), limits =c(4, 4.55))
  

#HUMAN
human.rq2 <- data.rq2[c("human", "phuman", "cond", "class", "religion4", "rel")]
human.rq2<-na.omit(human.rq2)

#MODEL 2
human2.rq2 <- lmer(phuman ~ human  + cond + (1 | class), data = human.rq2)
plot_model(human2.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(human2.rq2)#Table 2

#MODEL 3
human3.rq2 <- lmer(phuman ~ human  + cond + religion4 + (1 | class), data = human.rq2)
plot_model(human3.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(human3.rq2)#Table 2

#MODEL 4
human4.rq2 <- lmer(phuman ~ human  + cond * religion4 + (1 | class), data = human.rq2)
plot_model(human4.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(human4.rq2)#Table 2

# Calculate estihumand marginal 
em.human <- emmeans(human4.rq2, ~ cond*religion4)
con.human<-pairs(em.human, by="religion4", combine=TRUE)
sum.human<-summary(em.human)
con.human #TABLE3
sum.human #FOR FIGURE 3
sum.human$religion4<-factor(sum.human$religion4, levels=c("Atheist", "Agnostic", "Other", "Christian"))
sum.human$cond<-factor(sum.human$cond, levels=c("Christian Instructor", "Non-Religious Instructor"))


ggplot(sum.human, aes(x = religion4,  y = emmean, color = cond, group=cond)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean  + SE, linetype = factor(religion4)), width = 0.2, linewidth = 1, linetype = "solid") +
  geom_point(size = 2, stroke = 1) +
  scale_color_manual(values = c("black", "gray"), name = "Cultural Competence with") +
  labs(x = "Students' religious identity", y = "Estihumand Marginal Means in \n Human Evolution Acceptance") +
  theme_classic() +  
  scale_y_continuous(breaks=seq(3.75, 4.25, 0.25), limits =c(3.75, 4.25))


###evound
evound.rq2 <- data.rq2[c("evound", "pevound", "cond", "class", "religion4", "rel")]
evound.rq2<-na.omit(evound.rq2)

#MODEL 2
evound2.rq2 <- lmer(pevound ~ evound  + cond + (1 | class), data = evound.rq2)
plot_model(evound2.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(evound2.rq2)#Table 2

#MODEL 3
evound3.rq2 <- lmer(pevound ~ evound  + cond + religion4 + (1 | class), data = evound.rq2)
plot_model(evound3.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(evound3.rq2)#Table 2

#MODEL 4
evound4.rq2 <- lmer(pevound ~ evound  + cond * religion4 + (1 | class), data = evound.rq2)
plot_model(evound4.rq2,show.values = TRUE, value.offset = .5, title = "RQ2: Does an instructor revealing a christian or non-religious identity impact the efficacy of ReCCEE?")  + theme_sjplot()
summ(evound4.rq2)#Table 2

# Calculate estievoundd marginal 
em.evound <- emmeans(evound4.rq2, ~ cond*religion4)
con.evound<-pairs(em.evound, by="religion4", combine=TRUE)
sum.evound<-summary(em.evound)
con.evound #TABLE3
sum.evound #FOR FIGURE 3
sum.evound$religion4<-factor(sum.evound$religion4, levels=c("Atheist", "Agnostic", "Other", "Christian"))
sum.evound$cond<-factor(sum.evound$cond, levels=c("Christian Instructor", "Non-Religious Instructor"))


ggplot(sum.evound, aes(x = religion4,  y = emmean, color = cond, group=cond)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean  +SE, linetype = factor(religion4)), width = 0.2, linewidth = 1, linetype = "solid") +
  geom_point(size = 2, stroke = 1) +
  scale_color_manual(values = c("black", "gray"), name = "Cultural Competence with") +
  labs(x = "Students' religious identity", y = "Estimated Marginal Means in \n Evolution Understanding") +
  theme_classic() +
  scale_y_continuous(breaks=seq(0.6, 0.75, 0.02), limits =c(0.6, 0.75))


#### revision

# Select the relevant columns from your data
selected_data <- data %>% dplyr::select(cond, religion4, gender_text, race_text, semester)

# Function to perform Chi-square test
chi_square_test <- function(variable) {
  contingency_table <- table(selected_data$cond, selected_data[[variable]])
  test_result <- chisq.test(contingency_table)
  return(test_result)
}

# Perform Chi-square tests for each categorical variable
variables <- c("religion4", "religion2", "relcat", "gender_text", "race_text", "semester")

results <- lapply(variables, function(var) {
  result <- chi_square_test(var)
  cat("\nChi-square test for", var, "\n")
  print(result)
  cat("------------------------------------------------------\n")
  return(result)
})

names(resusummary_table <- lts) <- variables

# Load necessary libraries
library(dplyr)
library(tidyr)

# Calculate counts for each level of the factor variables grouped by condition

# Summary table for religion4
religion4_summary <- data %>%
  group_by(cond, religion4) %>%
  summarise(count = n(), .groups = 'drop') %>%
  spread(religion4, count, fill = 0)

cat("\nSummary table for religion4\n")
print(religion4_summary)
cat("------------------------------------------------------\n")

# Summary table for gender_text
gender_text_summary <- data %>%
  group_by(cond, gender_text) %>%
  summarise(count = n(), .groups = 'drop') %>%
  spread(gender_text, count, fill = 0)

cat("\nSummary table for gender_text\n")
print(gender_text_summary)
cat("------------------------------------------------------\n")

# Summary table for race_text
race_text_summary <- data %>%
  group_by(cond, race_text) %>%
  summarise(count = n(), .groups = 'drop') %>%
  spread(race_text, count, fill = 0)

cat("\nSummary table for race_text\n")
print(race_text_summary)
cat("------------------------------------------------------\n")

# Summary table for semester
semester_summary <- data %>%
  group_by(cond, semester) %>%
  summarise(count = n(), .groups = 'drop') %>%
  spread(semester, count, fill = 0)

cat("\nSummary table for semester\n")
print(semester_summary)
cat("------------------------------------------------------\n")




###Revision 24 07 31
# only christian students, include religisotu, see ploit religiosity level
christ <- data %>%
  filter(religion4 == "Christian" & cond != "Control")

christ2 <- data %>%
  filter(religion4 == "Christian")

par(mfrow=c(1,2))
hist(christ$rel, 
     main = "Religiosity levels for Christian \n in 2 Cultural Competence Conditions (n=864)")

hist(christ2$rel, 
     main = "Religiosity levels for all Christian students (n=1325)")



library(sjlabelled)
set_label(christ$ppcomp) <- "Perceived Compatibility (post)"
set_label(christ$ppconf) <- "Perceived Conflict (post)"

pcomp.ch <- lmer(ppcomp ~ pcomp + cond * rel + (1 | class), data = christ)
pconf.ch <- lmer(ppconf ~ pconf + cond * rel + (1 | class), data = christ)

summ(pcomp.ch)
summ(pconf.ch)

q1<-plot_model(pcomp.ch, show.values = TRUE, value.offset = .5) + theme_sjplot()
q2<-plot_model(pconf.ch, show.values = TRUE, value.offset = .5) + theme_sjplot()


p1<-plot_model(pcomp.ch, type = "pred", terms = c("rel", "cond")) + theme_sjplot() +
  ylab("Perceived Compatibility") +
  xlab("Religious level") 
p2<-plot_model(pconf.ch, type = "pred", terms = c("rel", "cond")) + theme_sjplot() +
  ylab("Perceived Conflict") +
  xlab("Religious level") 


nested <- ((q1/q2)|(q3/q4)) + plot_annotation(tag_levels = "A", 
                                              title = "Model: Post ~ Pre + Cond*Rel (nested by class) | data : Christian students only (864)
                                              A & B : Regression coefficient plot
                                              C & D : Marginal effects of interaction")
nested #view multi-panel figure

sim1<-sim_slopes(pcomp.ch, pred=rel, modx =cond)
sim2<-sim_slopes(pconf.ch, pred=rel, modx =cond)
t.sim1<-tidy(sim1)
t.sim2<-tidy(sim2)
list.slope<-rbind(t.sim1, t.sim2)
list.slope$p.value<-round(list.slope$p.value, 2)
list.slope$estimate <-round(list.slope$estimate , 2)
list.slope$std.error <-round(list.slope$std.error , 2)
list.slope$statistic  <-round(list.slope$statistic  , 2)
list.slope$conf.low  <-round(list.slope$conf.low  , 2)
list.slope$conf.high  <-round(list.slope$conf.high  , 2)

list.slope<-list.slope %>% dplyr::select(estimate, std.error, statistic, p.value, modx.value, conf.low, conf.high)
list.slope$variable<-c("Compatibility", "Compatibility", "Conflict", "Conflict")
list.slope$variable<-as.factor(list.slope$variable)
list.slope$modx.value<-as.factor(list.slope$modx.value)

list.slope[x]<-"modx.value"
# Create a gt table based on preprocessed
list.slope %>% dplyr::select(variable, modx.value, estimate, p.value, conf.low, conf.high) %>%
  rename(condition = modx.value,
         ??=estimate,
         "Outcome variable"=variable) %>%
  gt() %>%
  tab_header(
    title = " ",
    subtitle = "Results from the simple slope analyses testing \n to what extent religiosity predicts outcome variable \n for students in different conditions")

summ(pconf.ch)
conf<-interactions::interact_plot(pconf.ch, pred =rel, 
                                  modx = cond,
                                  data = christ,
                                  interval = TRUE,
                                  x.label = "Religiosity",
                                  y.label = "Predicted perceived conflict") + ylim(1,5) 
box<-rbind(comp,conf)


####
#DESCRIPTIVE STATISTIC
desc <- data %>% dplyr::select(pcomp, pconf, ppconf, ppcomp, cond, rel, mate, pmate, human, phuman)

#PERCEIVED COMPATIBILITY
#Equality of variances (variance of the differen )
library(DescTools)
leveneTest(pcomp ~ cond, data=data)
#Pr(>F) is p value associated with F statistic. p>0.05 so sig evidence that the variance accross grpup are different. assumption met
aov.pcomp <- aov(pcomp ~ cond, data = data)
summary(aov.pcomp)
# Tukey's HSD test for pairwise comparisons
TukeyHSD(aov.pcomp)
data %>%
  group_by(cond) %>%
  summarise(
    mean = mean(pcomp, na.rm=TRUE),
    sd = sd(pcomp, na.rm=TRUE)
  )
# Calculate effect size (eta squared)
eta_sq_pcomp <- EtaSq(aov.pcomp, type = 2)
round(eta_sq_pcomp[1],3)



#PERCEIVED CONFLICT
leveneTest(pconf ~ cond, data=data)
#Pr(>F) is p value associated with F statistic. p>0.05 so sig evidence that the variance accross grpup are different. assumption met
aov.pconf <- aov(pconf ~ cond, data = data)
summary(aov.pconf)
# Tukey's HSD test for pairwise comparisons
TukeyHSD(aov.pconf)
data %>%
  group_by(cond) %>%
  summarise(
    mean = mean(pconf, na.rm=TRUE),
    sd = sd(pconf, na.rm=TRUE)
  )
# Calculate effect size (eta squared)
eta_sq_pconf <- EtaSq(aov.pconf, type = 2)
round(eta_sq_pconf[1],3)


#MATE
leveneTest(mate ~ cond, data=data)
#Pr(>F) is p value associated with F statistic. p>0.05 so sig evidence that the variance accross grpup are different. assumption met
aov.mate <- aov(mate ~ cond, data = data)
summary(aov.mate)
# Tukey's HSD test for pairwise comparisons
TukeyHSD(aov.mate)
data %>%
  group_by(cond) %>%
  summarise(
    mean = mean(mate, na.rm=TRUE),
    sd = sd(mate, na.rm=TRUE)
  )
# Calculate effect size (eta squared)
eta_sq_mate <- EtaSq(aov.mate, type = 2)
round(eta_sq_mate[1],3)



#human
leveneTest(human ~ cond, data=data)
#Pr(>F) is p value associated with F statistic. p>0.05 so sig evidence that the variance accross grpup are different. assumption met
aov.human <- aov(human ~ cond, data = data)
summary(aov.human)

# Tukey's HSD test for pairwise comparisons
TukeyHSD(aov.human)
data %>%
  group_by(cond) %>%
  summarise(
    mean = mean(human, na.rm=TRUE),
    sd = sd(human, na.rm=TRUE)
  )
# Calculate effect size (eta squared)
eta_sq_human <- EtaSq(aov.human, type = 2)
round(eta_sq_human[1],4)


#HUman mate, and evound

set_label(christ$phuman) <- "Human Evolution Acceptance (post)"
set_label(christ$pmate) <- "General Evolution Acceptance (post)"
set_label(christ$pevound) <- "Evolution Understanding"


human.ch <- lmer(phuman ~ human + cond * rel + (1 | class), data = christ)
mate.ch <- lmer(pmate ~ mate + cond * rel + (1 | class), data = christ)
evound.ch <- lmer(pevound ~ evound + cond * rel + (1 | class), data = christ)

summ(human.ch)
summ(mate.ch)
summ(evound.ch)


q1<-plot_model(human.ch, show.values = TRUE, value.offset = .5) + theme_sjplot()
q2<-plot_model(mate.ch, show.values = TRUE, value.offset = .5) + theme_sjplot()


p3<-plot_model(human.ch, type = "pred", terms = c("rel", "cond")) + theme_sjplot() +
  ylab("Human Evolution Acceptance ") +
  xlab("Religious level") 
p4<-plot_model(mate.ch, type = "pred", terms = c("rel", "cond")) + theme_sjplot() +
  ylab("General Evolution Acceptance") +
  xlab("Religious level") 

p5<-plot_model(evound.ch, type = "pred", terms = c("rel", "cond")) + theme_sjplot() +
  ylab("Evolution Undrstanding") +
  xlab("Religious level") 

nested <- ((q1/q2)|(p3/p4)) + plot_annotation(tag_levels = "A", 
                                              title = "Model: Post ~ Pre + Cond*Rel (nested by class) | data : Christian students only (864)
                                              A & B : Regression coefficient plot
                                              C & D : Marginal effects of interaction")

nested2 <- (p1|p2|p3|p3|p5) + plot_annotation(tag_levels = "A", 
                                              title = "Model: Post ~ Pre + Cond*Rel (nested by class) | data : Christian students only (864)")
nested2

sim1<-sim_slopes(human.ch, pred=rel, modx =cond)
sim2<-sim_slopes(mate.ch, pred=rel, modx =cond)
t.sim1<-tidy(sim1)
t.sim2<-tidy(sim2)
list.slope<-rbind(t.sim1, t.sim2)
list.slope$p.value<-round(list.slope$p.value, 2)
list.slope$estimate <-round(list.slope$estimate , 2)
list.slope$std.error <-round(list.slope$std.error , 2)
list.slope$statistic  <-round(list.slope$statistic  , 2)
list.slope$conf.low  <-round(list.slope$conf.low  , 2)
list.slope$conf.high  <-round(list.slope$conf.high  , 2)

list.slope<-list.slope %>% dplyr::select(estimate, std.error, statistic, p.value, modx.value, conf.low, conf.high)
list.slope$variable<-c("Compatibility", "Compatibility", "Conflict", "Conflict")
list.slope$variable<-as.factor(list.slope$variable)
list.slope$modx.value<-as.factor(list.slope$modx.value)
