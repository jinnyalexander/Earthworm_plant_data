# this is the code for 2017 10 x 10 m releve plant data
library(readxl)
rel_spp <- read_excel("~/Master's research/Coronavirus/releve_raw_data.xlsx", 
              sheet = "rel_spp")
rel_env_raw <- read_excel("~/Master's research/Coronavirus/releve_raw_data.xlsx", 
           sheet = "rel_env", col_types = c("guess", "skip", "text", 
           "text", "numeric", "skip", "numeric","numeric", "numeric", "numeric", "numeric", 
           "text", "text", "numeric", "numeric", "text", "text", "numeric", "text", "numeric", 
           "numeric", "numeric",  "skip", "text", "text", "skip", "numeric", "numeric", "numeric", 
           "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",  
           "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
           "numeric", "numeric", "skip", "skip"))

# add worm_category to rel_env df
worm_cat <- read_excel("~/Master's research/Coronavirus/releve_raw_data.xlsx", 
                       sheet = "worm_cat")
library(tidyverse)
rel_env <- left_join(rel_env_raw,worm_cat, by = 'site')

# omit releves from analysis that had been recently harvested
rel_env <- rel_env[c(-11, -14, -19, -31),]

# get rid of any spaces
rel_env <- rel_env %>% mutate(across(where(is.character), str_remove_all, 
                                     pattern = fixed(" ")))

# Check whether 2017 overstory cover was significantly different among worm categories
test <- aov(overstory_new ~ worm_category, data = rel_env)
summary(test) # nope

# percent cover class plant matrix
test <- as.data.frame(cbind(releve=rel_spp$dnr_releve_nbr, species=rel_spp$taxon_new, 
                            cover=rel_spp$cover_100))
test$releve <- as.factor(test$releve)
test$species <- as.factor(test$species)
test$cover <- as.numeric((test$cover))
relspp.mx <- test %>% pivot_wider(names_from = species, values_from = cover,
                                  values_fill = 0)
relspp.mx <- as.data.frame(relspp.mx)
relspp.mx1 <- relspp.mx[,2:147]
rownames(relspp.mx1) <- relspp.mx[,1]

# for some reason if I don't flip the order of rows in the df subsetting doesn't work right
relspp.mx2 <- relspp.mx1 %>% rownames_to_column('site') %>% map_df(rev) %>% column_to_rownames('site')

# subset data by worm category
u.u <- subset(relspp.mx2,rel_env$worm_category=='unwormed')
u.w <- subset(relspp.mx2,rel_env$worm_category=='short_term_wormed')
w.w <- subset(relspp.mx2,rel_env$worm_category=='long_term_wormed')

# get rid of spp w/ 0 occurrences in worm grouping subset
keep.spec1 <- names(u.u)[colSums(u.u>0)>0]
u.u1 <- u.u[,keep.spec1]
keep.spec1 <- names(u.w)[colSums(u.w>0)>0]
u.w1 <- u.w[,keep.spec1]
keep.spec1 <- names(w.w)[colSums(w.w>0)>0]
w.w1 <- w.w[,keep.spec1]

# diversity indices
library(vegan)
Dshan <- diversity(relspp.mx2, 'shannon')
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

# Pielou's evenness
J <- Dshan/log(specnumber(relspp.mx2))
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

# calculate simpson diversity overall
Dsimp <- diversity(relspp.mx2, 'invsimpson')
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

# species richness per plot
relspp.pa1 <- relspp.mx2 #fist convert matrix to presence absence
relspp.pa1[relspp.pa1 > 0] <- 1
rich <- rowSums(sign(relspp.pa1))
rel_env$richness <- rich # add richness to rel_env list
mean.rich <- mean(rich) # mean species richness overall
# separate out richness by worm status
u.u2 <- subset(rel_env$richness,rel_env$worm_category=='unwormed')
u.w2 <- subset(rel_env$richness,rel_env$worm_category=='short_term_wormed')
w.w2 <- subset(rel_env$richness,rel_env$worm_category=='long_term_wormed')
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

# combine diversity indices
div.comb <- rbind(rshan.c1, rsimp.c1, rJ.c1, rrich.c1)
div.comb

# Pielou's evenness significance testing
anova(lm(J~rel_env$worm_category))
pairwise.t.test(J, rel_env$worm_category, p.adjust.method = "fdr")

# Simpson's diversity significance testing
anova(lm(Dsimp~rel_env$worm_category))
pairwise.t.test(Dsimp, rel_env$worm_category, p.adjust.method = "fdr")

# Shannon diversity significance testing
anova(lm(Dshan~rel_env$worm_category))
pairwise.t.test(Dshan, rel_env$worm_category, p.adjust.method = "fdr")

# Species richness significance testing
anova(lm(rich~rel_env$worm_category))
pairwise.t.test(rich, rel_env$worm_category, p.adjust.method = "fdr")

# is species richness correlated with canopy cover?
cor.test(rel_env$richness, rel_env$overstory_new) # nope

####beta diversity dissimilarity
# are releves more similar w/in a treatment than between treatments?
relspp.mx2 <- relspp.mx1[nrow(relspp.mx1):1,] #flip row order
vdiss.2017 <- vegdist(relspp.mx2,'bray')
# matrix showing dissimilarity (higher values=less similar)

# mean B-C dissimilarity between and within wormed/unwormed plots
mvdiss.2017 <- meandist(vdiss.2017, rel_env$worm_category) # mean 2017 dissimilarity
# avg between vs within group dissimilarity
summary(mvdiss.2017)
mvdiss.2017 <- mvdiss.2017[c(3,2,1),c(3,2,1)]
mvdiss.2017

# adonis using bray-curtis dissimilarity
per.plants <- adonis2(vdiss.2017~rel_env$worm_category, nperms=100000)
per.plants

# post-hoc pairwise perMANOVA
# pairwise adonis function 
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', 
                            p.adjust.m ='bonferroni', perm=100000)
{
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],
                                   metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])],
                permutations = 100000 );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
} 
# apply function to my data
p.wise2017 <- pairwise.adonis(relspp.mx2,rel_env$worm_category, p.adjust.m = "fdr")

# output table
p.wise2017

# NMS with 2 dimensions
nms.1 <- metaMDS(relspp.mx2)
nms.1
stressplot(nms.1)

# graph NMS ordination with ggplot
data.scores <- as.data.frame(scores(nms.1))  
# Using the scores function from vegan to extract the site scores and convert to a df
data.scores$site <- rownames(data.scores)  
# create a column of site names, from the rownames of data.scores
data.scores$grp <- rel_env$worm_category 
#  add the grp variable created earlier
species.scores <- as.data.frame(scores(nms.1, "species"))  
#Using the scores function from vegan to extract the species scores and convert to a df
species.scores$species <- rownames(species.scores)  
# create a column of species, from the rownames of species.scores
library(ggplot2)
metadata_nmds <- data.scores
star <- metadata_nmds %>%  group_by(grp) %>% mutate(centroid1 = mean(NMDS1),
                                                    centroid2 = mean(NMDS2)) %>% ungroup
centroid <- metadata_nmds %>% group_by(grp) %>% summarise(axis1 = mean(NMDS1), 
                                                          axis2 = mean(NMDS2), .groups="drop")
ggplot(star, aes(x=NMDS1, xend=centroid1, y=NMDS2, yend=centroid2, 
                                 color=grp)) +
  geom_point() +
  geom_point(data=centroid, mapping=aes(x=axis1, y=axis2, color=grp), 
             shape =16, size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment() +
  # coord_fixed(xlim = c(-0.5, 1), ylim = c(-0.7, 0.5)) +
  labs(x="NMDS Axis 1",
       y="NMDS Axis 2") +
  scale_color_manual(name=NULL,
                     breaks=c("unwormed",
                              "short_term_wormed",
                              "long_term_wormed"),
                     values=c("#0571b0", "#f4a582", "#ca0020"),
                     labels=c("Unwormed (n=6)",
                              "Short-term wormed (n=23)",
                              "Long-term wormed (n=12)"))+
  theme_bw() +
  theme(
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.2, 0.07),
    legend.background = element_rect(fill="white",
                                     color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3))
# library(gridExtra) # code for side-by-side releve & SAP NMS comparison
# grid.arrange(sap.nms.gg, rel.nms.gg, ncol=2) 

# correlation of env variables with NMS axes
# what environmental variables do I want to test?
# utm_e, utm_n, elev_ft, subsec, lta, topopos,
# remove redundant environmental variables
## removed 'month' as an env var bc we focused on finding unwormed sites earlier in summer
myvars <- names(rel_env) %in% c("new_utm_e", "new_utm_n", "elev_ft", "subsec", "lta",
      "overstory_new", "topopos", "aspect", "perslope", "e_inch", "b_inch", "ctop", 
          "avg_tex", "stand_age_2017, IERAT_new", "wormed_then_pred", "wormed_now_actual",
       "new_a", "new_o", "old_a", "old_o", "epigeic_worms_pres", "epiendogeic_worms_pres",
      "endogeic_worms_pres", "anecic_worms_pres", "epiendogeic_anecic_worms_pres", 
          "non_epigeic_worms_pres", "ecogroup_sum", "worm_category")
newdata <- rel_env[myvars]
ef <- envfit(nms.1, newdata, permutations = 10000)
plot(nms.1)
plot(ef, p.max = 0.01)
capture.output(ef$vectors, file = "~/Master's research/Coronavirus/NMSenvcorr.txt")

############## Look at differences in species cover between earthworm categories
# calculate mean cover in w.w vs u.w vs u.u areas
# this is including 0s
mcov.uu <- apply(u.u,2,mean)
mcov.uw <- apply(u.w,2,mean)
mcov.ww <- apply(w.w,2,mean)
mcov.n <- colSums(relspp.mx1>0)
# paste 3 rows together
mcov.c <- data.frame(mcov.n,mcov.uu,mcov.uw,mcov.ww)
rmcov.c <- round(mcov.c, 1)
names(rmcov.c) <- c("n", "Unwormed","Short-term wormed","Long-term wormed")
rmcov.c$Species <- rownames(rmcov.c)
rmcov.c <- rmcov.c[,c(5,1:4)]
test<-rmcov.c[!(rmcov.c$n<=4),] # exclude species with less than 5 occurrences
test

# total cover by life form
life.form <- as.data.frame(rel_spp[!duplicated(rel_spp$taxon_new), ])
write.csv(life.form, "~/Master's research/Coronavirus/lifeform.csv")
lifeform.ed <- read_excel ("~/Master's research/Coronavirus/lifeform_edited.xlsx")
test <- left_join(rel_spp, lifeform.ed, by="taxon_new")
test.1 <- as.data.frame(test %>% select(dnr_releve_nbr.x, taxon_new, physcode.x, cover_100))
test.1$dnr_releve_nbr.x <- as.factor(test.1$dnr_releve_nbr.x)
test.1$physcode.x <- as.factor(test.1$physcode.x)
test.2 <- test.1 %>% group_by(dnr_releve_nbr.x, physcode.x) %>% summarise(sum(cover_100))
test.3 <- dplyr::rename(test.2, releve = dnr_releve_nbr.x, life_form = physcode.x, 
                        cover = 'sum(cover_100)')
new_lf <- test.3
# Adding in 0's for groupings that were missing in a releve
new_lf1 <- data.frame()
for (r in unique(new_lf$releve)){
  comp <- data.frame("releve"=rep(r,4),"life_form"=c("D","E","G","H"))
  comp <- merge(comp,new_lf,all.y=FALSE,all.x=TRUE)
  comp[is.na(comp)] <- 0
  new_lf1 <- rbind(new_lf1,comp) }
# merging releve number with worm change info
worm.chg.df <- as.data.frame(cbind(rel_env$site, rel_env$worm_category))
worm.chg.df <- worm.chg.df %>% dplyr::rename(releve = V1,worm_category = V2)
new_lf2 <- merge(new_lf1,worm.chg.df,by.x = "releve", by.y = "releve")
# make worm change a factor
new_lf2$worm_category <- as.factor(new_lf2$worm_category)
# create 1 column for each life form
new_lf3 <- spread(new_lf2,key=life_form,value=cover)
new_lf4 <- new_lf3 %>% group_by(worm_category) %>% summarize_all(mean)
new_lf5 <- as.data.frame(cbind(new_lf4$D, new_lf4$E, new_lf4$G, new_lf4$H))
new_lf6 <- as.data.frame(t(new_lf5))
new_lf6 <- new_lf6 %>% dplyr::rename('Long-term wormed'=V1,'Short-term wormed'=V2, 'Unwormed'=V3)
row.names(new_lf6) <- c("Woody deciduous", "Woody evergreen", "Graminoids", "Herbaceous")
new_lf6 <- new_lf6[,3:1]
new_lf6

# try the same thing but with family data
test.1 <- as.data.frame(test %>% select(dnr_releve_nbr.x, taxon_new, family, cover_100))
test.1$dnr_releve_nbr.x <- as.factor(test.1$dnr_releve_nbr.x)
test.1$family <- as.factor(test.1$family)
test.2 <- test.1 %>% group_by(dnr_releve_nbr.x, family) %>% summarise(sum(cover_100)) #issue here
test.3 <- dplyr::rename(test.2, releve = dnr_releve_nbr.x, family = family, 
                        cover = 'sum(cover_100)')
new_lf <- test.3
# Adding in 0's for groupings that were missing in a releve
new_lf1 <- data.frame()
for (r in unique(new_lf$releve)){
  comp <- data.frame("releve"=rep(r,52),"family"=c("Adoxaceae", "Apiaceae","Araceae",  
             "Araliaceae","Aristolochiaceae", "Caprifoliaceae", "Colchicaceae", "Cornaceae",
    "Cyperaceae", "Dryopteridaceae",  "Fagaceae", "Grossulariaceae", "Liliaceae", 
              "Myrsinaceae", "Oleaceae","Orphioglossaceae",  "Pinaceae",  "Ranunculaceae",
             "Rosaceae",  "Ruscaceae",  "Salicaceae","Sapindaceae", "Asteraceae",  "Betulaceae",    
           "Diervillaceae", "Juncaceae", "Lycopodiaceae", "Malvaceae", "Melanthiaceae",  
         "Papaveraceae", "Poaceae",  "Smilacaceae",  "Thymelaeaceae",  "Violaceae",
             "Balsaminaceae",  "Equisetaceae", "Fabaceae", "Osmundaceae",  "Rubiaceae",
          "Saxifragaceae",  "Ulmaceae", "Ericaceae",  "Orchidaceae",  "Berberidaceae",
       "Rhamnaceae", "Onagraceae", "Vitaceae", "Dennstaedtiaceae", "Lamiaceae",
            "Apocynaceae",  "Pteridaceae",  "Anacardiaceae"))
  comp <- merge(comp,new_lf,all.y=FALSE,all.x=TRUE)
  comp[is.na(comp)] <- 0
  new_lf1 <- rbind(new_lf1,comp) }
# merging releve number with worm change info
worm.chg.df <- as.data.frame(cbind(rel_env$site, rel_env$worm_category))
worm.chg.df <- worm.chg.df %>% rename(releve = V1,worm_category = V2)
new_lf2 <- merge(new_lf1,worm.chg.df,by.x = "releve", by.y = "releve")
# make worm change a factor
new_lf2$worm_category <- as.factor(new_lf2$worm_category)
# create 1 column for each life form
new_lf3 <- spread(new_lf2,key=family,value=cover)
new_lf4 <- new_lf3 %>% group_by(worm_category) %>% summarize_all(mean)
new_lf5 <- as.data.frame(new_lf4[,-2])
new_lf6 <- as.data.frame(t(new_lf5))
colnames(new_lf6) <- new_lf6[1,]
new_lf6 <- new_lf6[-1, ] 
familynbr <- as.data.frame(table(new_lf$family)) # n by family
new_lf7 <- rownames_to_column(new_lf6, var = "family") %>% as_tibble()
new_lf8 <- merge(new_lf7,familynbr, by.x = "family", by.y = "Var1")
new_lf8$long_term_wormed <- as.numeric(new_lf8$long_term_wormed)
new_lf8$short_term_wormed <- as.numeric(new_lf8$short_term_wormed)
new_lf8$unwormed <- as.numeric(new_lf8$unwormed)
new_lf8 <- new_lf8[, c(1,5,4,3,2)]
new_lf8