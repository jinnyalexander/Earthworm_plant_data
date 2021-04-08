# Species accumulation code and plant community matrix code
library(readxl)
sp_raw <- read_excel("~/Master's research/Coronavirus/Species area plot raw data.xlsx", 
                     sheet = "Sheet1")
sp_raw$Site <- as.factor(sp_raw$Site)
vegdata1 <- read_excel("~/Master's research/Coronavirus/Species area plot raw data.xlsx", 
                       sheet = "Species_area_raw")

# transform nested plot occurrence to ranked abundance values (1=1, 2=2, 4=3, 1024=11)
vegdata <- vegdata1[-c(4)]
vegdata$abun <- vegdata$Area
vegdata$abun[vegdata$abun=="1"]<-11
vegdata$abun[vegdata$abun=="2"]<-10
vegdata$abun[vegdata$abun=="4"]<-9
vegdata$abun[vegdata$abun=="8"]<-8
vegdata$abun[vegdata$abun=="16"]<-7
vegdata$abun[vegdata$abun=="32"]<-6
vegdata$abun[vegdata$abun=="64"]<-5
vegdata$abun[vegdata$abun=="128"]<-4
vegdata$abun[vegdata$abun=="256"]<-3
vegdata$abun[vegdata$abun=="512"]<-2
vegdata$abun[vegdata$abun=="1024"]<-1
vegdata <- vegdata[-c(2)]

# transform the raw data into a site x species community matrix
library(reshape)
spp.matrix <- cast(vegdata, Site ~ Species, value='abun')
spp.matrix <- as.data.frame(spp.matrix)
spp.matrix[is.na(spp.matrix)] <- 0
spp.matrix1 <- spp.matrix[,2:189]
rownames(spp.matrix1) <- spp.matrix[,1]

# make a matrix with site x env variables
worm_env <- read_excel("~/Master's research/Coronavirus/spp.env.xlsx", 
                       sheet = "Sheet3")
spp_env <- read_excel("~/Master's research/Coronavirus/spp.env.xlsx", 
                      sheet = "Sheet1")
library(dplyr)
try1 <- left_join(spp_env, worm_env, by = 'site')
env.matrix11 <- as.data.frame(cbind(site=try1$site, month=try1$month_now, 
                                    utm_e=try1$new_utm_e, overstory=try1$overstory_new,
                                    utm_n=try1$new_utm_n,elev_ft=try1$elev_ft, subsec=try1$subsec, lta=try1$lta, 
                                    topopos=try1$topopos, aspect=try1$aspect, perslope=try1$perslope, drclass=try1$drclass, 
                                    stand_age_2017=try1$stand_age_2017, county=try1$county, IERAT_new=try1$IERAT_new, 
                                    wormed_then_pred=try1$wormed_then_pred, changed_to_wormed=try1$changed_to_wormed, 
                                    wormed_now_actual=try1$wormed_now_actual, new_a=try1$new_a, new_o=try1$new_o, 
                                    epigeic_worms_pres=try1$epigeic_worms_pres, 
                                    epiendogeic_worms_pres=try1$epiendogeic_worms_pres, 
                                    endogeic_worms_pres=try1$endogeic_worms_pres, anecic_worms_pres=try1$anecic_worms_pres, 
                                    epiendogeic_anecic_worms_pres=try1$epiendogeic_anecic_worms_pres,
                                    non_epigeic_worms_pres=try1$non_epigeic_worms_pres, ecogroup_sum=try1$ecogroup_sum, 
                                    worm_category=try1$worm_category))

# get rid of spaces
library(stringr)
env.matrix <- env.matrix11 %>% mutate(across(where(is.character), str_remove_all, 
                                             pattern = fixed(" ")))

# make sure variables are numerical or factors
env.matrix <- data.frame(lapply(env.matrix, as.factor), stringsAsFactors=FALSE)
env.matrix$month <- as.numeric(env.matrix$month)
env.matrix$overstory <- as.numeric(env.matrix$month)
env.matrix$utm_e <- as.numeric(env.matrix$utm_e)
env.matrix$utm_n <- as.numeric(env.matrix$utm_n)
env.matrix$elev_ft <- as.numeric(env.matrix$elev_ft)
env.matrix$perslope <- as.numeric(env.matrix$perslope)
env.matrix$stand_age_2017 <- as.numeric(env.matrix$stand_age_2017)
env.matrix$IERAT_new <- as.numeric(env.matrix$IERAT_new)
env.matrix$new_a <- as.numeric(env.matrix$new_a)
env.matrix$new_o <- as.numeric(env.matrix$new_o)
env.matrix$ecogroup_sum <- as.numeric(env.matrix$ecogroup_sum)

# separate spp.matrix1 into u.u u.w and w.w categories
# here I am using the plant matrix including rare spp
u.u <- subset(spp.matrix1,env.matrix$worm_category=='unwormed')
u.w <- subset(spp.matrix1,env.matrix$worm_category=='short_term_wormed')
w.w <- subset(spp.matrix1,env.matrix$worm_category=='long_term_wormed')

# get rid of spp w/ 0 occurrences
keep.spec1 <- names(u.u)[colSums(u.u>0)>0]
u.u1 <- u.u[,keep.spec1]
keep.spec1 <- names(u.w)[colSums(u.w>0)>0]
u.w1 <- u.w[,keep.spec1]
keep.spec1 <- names(w.w)[colSums(w.w>0)>0]
w.w1 <- w.w[,keep.spec1]

# graph the 3 species area curves


library(ggplot2)
ggplot(sp_raw, aes(x=log(Area), y=log(Total_num), color=worm_category)) + geom_point() + 
  geom_smooth(method=glm, 
              formula = y~(x), show.legend = FALSE) +
  theme(
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.8, 0.1),
    legend.background = element_rect(fill="white",
                                     color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3)) +
  labs(y= "Species accumulation (log)", x = "Area sampled (log)" )+ 
  geom_jitter(width=.1)


# How speciose are my communities?
# shannon diversity
# calculate shannon diversity overall
library(vegan)
Dshan <- diversity(spp.matrix1, 'shannon')
mean.dshan <- mean(Dshan)
# calculate shannon diversity in u.u sites
uu.shan <- diversity(u.u1)
m.uu.shan <- mean(uu.shan)
# calculate shannon diversity in u.w sites
uw.shan <- diversity(u.w1)
m.uw.shan <- mean(uw.shan)
# calculate shannon diversity in w.w sites
ww.shan <- diversity(w.w1)
m.ww.shan <- mean(ww.shan)
# paste 4 shannon indices together
shan.c <- data.frame(mean.dshan, m.uu.shan, m.uw.shan, m.ww.shan)
rshan.c <- round(shan.c, 2)
names(rshan.c) <- c("Mean", "Unwormed","Short-term wormed","Long-term wormed")
rshan.c$Diversity <- "Shannon-Wiener index"
rshan.c1 <- rshan.c[c("Diversity", "Mean", "Unwormed",
                      "Short-term wormed","Long-term wormed")]
rshan.c1

# calculate simpson diversity overall
Dsimp <- diversity(spp.matrix1, 'invsimpson')
mean.dsimp <- mean(Dsimp)
# calculate simpson diversity in u.u sites
uu.simp <- diversity(u.u1, 'invsimpson')
m.uu.simp <- mean(uu.simp)
# calculate simpson diversity in u.w sites
uw.simp <- diversity(u.w1, 'invsimpson')
m.uw.simp <- mean(uw.simp)
# calculate simpson diversity in w.w sites
ww.simp <- diversity(w.w1, 'invsimpson')
m.ww.simp <- mean(ww.simp)
# paste 4 simpson indices together
simp.c <- data.frame(mean.dsimp, m.uu.simp, m.uw.simp, m.ww.simp)
rsimp.c <- round(simp.c, 2)
names(rsimp.c) <- c("Mean", "Unwormed","Short-term wormed","Long-term wormed")
rsimp.c$Diversity <- "Simpson's index (1/D)"
rsimp.c1 <- rsimp.c[c("Diversity", "Mean", "Unwormed",
                      "Short-term wormed","Long-term wormed")]
rsimp.c1

# species richness per plot
spp.matrix.pa <- spp.matrix1 #fist convert matrix to presence absence
spp.matrix.pa[spp.matrix.pa > 0] <- 1
spp.matrix.pa[is.na(spp.matrix.pa)] <- 0
rich <- rowSums(sign(spp.matrix.pa))
# add richness to env.var list
env.matrix$richness <- rich
# mean species richness overall
mean.rich <- mean(rich)
# separate out richness by worm status
u.u2 <- subset(env.matrix$richness,env.matrix$worm_category=='unwormed')
u.w2 <- subset(env.matrix$richness,env.matrix$worm_category=='short_term_wormed')
w.w2 <- subset(env.matrix$richness,env.matrix$worm_category=='long_term_wormed')
# mean species richness by worm status
m.uu.rich <- mean(u.u2)
m.uw.rich <- mean(u.w2)
m.ww.rich <- mean(w.w2)
# paste 4 mean richness values together
rich.c <- data.frame(mean.rich, m.uu.rich, m.uw.rich, m.ww.rich)
rrich.c <- round(rich.c, 2)
names(rrich.c) <- c("Mean", "Unwormed","Short-term wormed","Long-term wormed")
rrich.c$Diversity <- "Richness"
rrich.c1 <- rrich.c[c("Diversity", "Mean", "Unwormed",
                      "Short-term wormed","Long-term wormed")]

# Pielou's evenness
J <- Dshan/log(specnumber(spp.matrix1))
mean.J <- mean(J)
# calculate Pielou's evenness in u.u sites
uu.J <- uu.shan/log(specnumber(u.u1))
m.uu.J <- mean(uu.J)
# calculate Pielou's evenness in u.w sites
uw.J <- uw.shan/log(specnumber(u.w1))
m.uw.J <- mean(uw.J)
# calculate Pielou's evenness in w.w sites
ww.J <- ww.shan/log(specnumber(w.w1))
m.ww.J <- mean(ww.J)
# paste 4 Pielou's evenness together
J.c <- data.frame(mean.J, m.uu.J, m.uw.J, m.ww.J)
rJ.c <- round(J.c, 2)
names(rJ.c) <- c("Mean", "Unwormed","Short-term wormed","Long-term wormed")
rJ.c$Diversity <- "Pielou's evenness"
rJ.c1 <- rJ.c[c("Diversity", "Mean", "Unwormed",
                "Short-term wormed","Long-term wormed")]

# combine diversity indices
div.comb <- rbind(rshan.c1, rsimp.c1, rJ.c1, rrich.c1)
# make a better looking table
library(stargazer)
stargazer(div.comb, type = "html", digits=2, rownames=FALSE,
          summary=FALSE,
          title = "Diversity indices for species area plot data",
          out = "~/Master's research/Coronavirus/Tables/diversity.html")

# Simpson's diversity significance testing
## should worm category be a factor or numeric for ANOVA?
anova(lm(J~env.matrix$worm_category))
pairwise.t.test(J, env.matrix$worm_category, p.adjust.method = "fdr")

# Shannon diversity significance testing
anova(lm(Dshan~env.matrix$worm_category))
pairwise.t.test(Dshan, env.matrix$worm_category, p.adjust.method = "fdr")

# Species richness significance testing
anova(lm(rich~env.matrix$worm_category))
pairwise.t.test(rich, env.matrix$worm_category, p.adjust.method = "fdr")

####beta diversity dissimilarity
# are releves more similar w/in a treatment than between treaments?
# get rid of species with less than 3 occurrences (<5% of releves)
keep.spec <- names(spp.matrix1)[colSums(spp.matrix1>0)>3]
spp.matrix2 <- spp.matrix1[,keep.spec]
# spp.matrix2 <- wisconsin(spp.matrix2) # code for wisconsin standardization
# b-c dissimilarity
vdiss.2017 <- vegdist(spp.matrix1,'bray')
# matrix showing dissimilarity (higher values=less similar)

# mean B-C dissimilarity between and within wormed/unwormed plots
mvdiss.2017 <- meandist(vdiss.2017, env.matrix$worm_category) # mean 2017 dissimilarity
# avg between vs within group dissimilarity
summary(mvdiss.2017)
mvdiss.2017 <- mvdiss.2017[c(3,2,1),c(3,2,1)]
stargazer(mvdiss.2017, type = "html", digits=2,
          summary=FALSE,
          title = "Mean dissimilarity between and within worm categories",
          flip = TRUE,
          out = "~/Master's research/Coronavirus/Tables/meandis.html")

#jaccard's index
vJaccard <- vegdist(spp.matrix1, 'jaccard')
# by treatment
vj.uu <- vegdist(u.u, 'jaccard')
vj.uw <- vegdist(u.w, 'jaccard')
vj.ww <- vegdist(w.w, 'jaccard')
mean(vj.uu)
mean(vj.uw)
mean(vj.ww)

# adonis using bray-curtis dissimilarity
# test if treatment significantly affects spp comp
per.plants <- adonis2(vdiss.2017~env.matrix$worm_category, nperms=100000)
per.plants # no significant differences between wormed categories

# lets try NMS with 2 dimensions
nms.1 <- metaMDS(spp.matrix2)
nms.1
stressplot(nms.1)
ordiplot(nms.1,type="n")
ordispider(nms.1, groups = env.matrix$worm_category, col=1:3)

# lets try NMS with 3 dimensions
nms.3 <- metaMDS(spp.matrix2, k=3)
stressplot(nms.3)
ordiplot(nms.3,type="n")
ordihull(nms.3,groups=env.matrix$worm_category,draw="polygon",col="grey90",label=F)
orditorp(nms.3,display="species",col="red",air=0.01)
orditorp(nms.3,display="sites",cex=1.25,air=0.01)

# try re-creating ordination with ggplot
#Using the scores function from vegan to extract the site scores and convert to a df
data.scores <- as.data.frame(scores(nms.1)) 
# create a column of site names, from the rownames of data.scores
data.scores$site <- rownames(data.scores)  
#  add the grp variable created earlier
data.scores$grp <- env.matrix$worm_category  
#Using the scores function from vegan to extract the species scores and convert to a df
species.scores <- as.data.frame(scores(nms.1, "species")) 
species.scores$species <- rownames(species.scores)  
# create a column of species, from the rownames of species.scores
library(ggplot2)
metadata_nmds <- data.scores
star <- metadata_nmds %>%  group_by(grp) %>% mutate(centroid1 = mean(NMDS1),
                                                    centroid2 = mean(NMDS2)) %>% ungroup
centroid <- metadata_nmds %>% group_by(grp) %>% summarise(axis1 = mean(NMDS1), 
                                                          axis2 = mean(NMDS2), .groups="drop")
ggplot(star, aes(x=NMDS1, xend=centroid1, y=NMDS2, yend=centroid2, 
                 color=grp)) + geom_point() +
  geom_point(data=centroid, mapping=aes(x=axis1, y=axis2, color=grp), 
             shape =16, size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment() +
  # coord_fixed(xlim = c(-0.5, 0.9), ylim = c(-0.5, 0.5)) +
  labs(x="NMDS Axis 1",
       y="NMDS Axis 2") +
  scale_color_manual(name=NULL,
                     breaks=c("unwormed",
                              "short_term_wormed",
                              "long_term_wormed"),
                     values=c("#fde0dd", "#fa9fb5", "#c51b8a"),
                     labels=c("Unwormed (n=6)",
                              "Short-term wormed (n=22)",
                              "Long-term wormed (n=12"))+
  theme_dark() +
  ggtitle("Species area plot data")+
  theme(
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.8, 0.05),
    legend.background = element_rect(fill="white",
                                     color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3))

############## Look at differences in species cover between earthworm categories
# calculate mean cover in w.w vs u.w vs u.u areas
# this is including 0s
mcov.uu <- apply(u.u,2,mean)
mcov.uw <- apply(u.w,2,mean)
mcov.ww <- apply(w.w,2,mean)
mcov.n <- colSums(spp.matrix1>0)
# paste 3 rows together
mcov.c <- data.frame(mcov.n,mcov.uu,mcov.uw,mcov.ww)
rmcov.c <- round(mcov.c, 1)
names(rmcov.c) <- c("n", "Unwormed","Short-term wormed","Long-term wormed")
rmcov.c$Species <- rownames(rmcov.c)
rmcov.c <- rmcov.c[,c(5,1:4)]
test<-rmcov.c[!(rmcov.c$n<=4),] # exclude species with less than 5 occurrences
# make a better looking table
stargazer(test, type = "html", digits=1, rownames=FALSE,
          summary=FALSE,
          title = "Mean understory cover between worm categories",
          flip = FALSE,
          out = "~/Master's research/Coronavirus/Tables/meancover.html")

# Check whether 2017 overstory cover was significantly different among worm categories
test <- aov(overstory ~ worm_category, data = env.matrix)
summary(test)
