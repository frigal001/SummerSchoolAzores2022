######################################################################
######################################################################
######---------------FD summer school Azores 2022---------------######
######################################################################
######################################################################

setwd("/Users/francoisrigal")

install.packages(c("ade4", "BAT", "vegan", "tidyverse", "gawdis", "devtools", 
                   "ggrepel", "hypervolume", "picante", "GGally", "ggpubr", "nlme", 
                   "emmeans", "ggcorrplot", "devtools", "BiocManager", 
                   "ggnewscale", "ggstar", "ggrepel")) # rest omitted

BiocManager::install("ggtree")
devtools::install_github("xiangpin/ggtreeExtra")

lapply(c("ade4", "BAT", "vegan", "tidyverse", "gawdis", "devtools", 
         "ggrepel", "hypervolume", "picante", "GGally", "ggpubr", "nlme", 
         "emmeans", "ggcorrplot", "devtools", "BiocManager", 
         "ggnewscale", "ggstar", "ggtreeExtra", "ggtree", "ggrepel"), require, character.only = TRUE)

######################################################################
######################################################################
######################################################################

#setwd("/Users/francoisrigal")

species_info <-read.table("traits_matrix_arthropods.txt", h=T)
View(species_info)

sites <- read.table("community_matrix_arthropods.txt", h=T)
View(sites)

dim(species_info)

dim(sites)


Hab <- factor(sites$Habitat, levels = c("NAT", "EXO", "SEM", "INT"))

HabOrd <- sites$Hab_Ordinal

comm <- sites[,-c(1:2)]

dim(comm)

taxo <- species_info[,c(4:1)]

trait <- species_info[,c(6:16)]
View(trait)

origin <- species_info$Origin

comm.ind <- comm[,rownames(species_info[species_info$Origin=="IND",])]

comm.exo <- comm[,rownames(species_info[species_info$Origin=="NON_IND",])]

######################################################################
######################################################################

# overview of the trait matrix

hist(trait$bodysize)

apply(trait[,c(2:5)], 2, sum) %>% barplot

table(trait$GetFood) %>% barplot

table(trait$IngestFood) %>% barplot

apply(trait[,c(8:10)], 2, sum) %>% barplot

table(trait$Dispersal) %>% barplot


######################################################################
######################################################################

# FD gowdis()

# set up the framework of Pavoine for the Gower distance
# we create 3 groups: the multichoice, the continuous and the 
# nominal variables
#colnames(trait[,c(2:5, 8:10)])

w <- prep.binary(trait[,c(2:5, 8:10)], col.blocks=c(4, 3), 
                 label = c("Food", "Period")) # create the matrix for the multichoice

bs <- data.frame(bodysize = log10(trait$bodysize));rownames(bs) <- rownames(trait) 

nom <- data.frame(dispersal = trait$Dispersal, getfood=trait$GetFood, ingfood=trait$IngestFood)
rownames(nom) <- rownames(trait)

ktab1 <- ktab.list.df(list(w, bs, nom)) # global data

distw.tot <- dist.ktab(ktab1, type=c("B", "Q", "N")) # define the type of trait

#ListD <- ldist.ktab(ktab1, type=c("B", "Q", "N")) #list of the distance

kc <- kdist.cor(ktab1, type=c("B", "Q", "N")) # correlations

round(as.dist(kc$paircor), 3) #pairwise correlation between pairs of distances

ggcorrplot(kc$paircor, hc.order = TRUE, type = "lower",
           lab = TRUE)

round(kc$glocor, 3) # contribution

# method gawdis

weights_equal <- gawdis(trait, w.type="analytic")
attr(weights_equal, "weights")

cor(weights_equal, distw.tot)
plot(weights_equal, distw.tot) # high correlation but not max

######################################################################
######################################################################

# visualization of the trait space #

# PCoA

pcoa_trait <- pcoa(distw.tot)

pcoa_trait$values$Relative_eig

p <- ggplot(pcoa_trait$vectors%>%as.data.frame, aes(x=Axis.1, y=Axis.2, color=origin)) +
  geom_point(shape=19, size=2, alpha=0.7) +
  scale_color_manual(name="Species origin",
                     values=c("green4", "grey80"), labels=c("Indigenous", "Exotics")) +
  labs(title="PCoA of Gower Distances",
       x="PCo Axis 1 (30%)",
       y="PCo Axis 2 (25%)") +
  theme_bw()
p


# add the traits 

coord_getfood <- cbind(tapply(pcoa_trait$vectors[,1], trait$GetFood, mean),
             tapply(pcoa_trait$vectors[,2], trait$GetFood, mean))

coord_ingest <- cbind(tapply(pcoa_trait$vectors[,1], trait$IngestFood, mean),
                       tapply(pcoa_trait$vectors[,2], trait$IngestFood, mean))

coord_disp <- cbind(tapply(pcoa_trait$vectors[,1], trait$Dispersal, mean),
                      tapply(pcoa_trait$vectors[,2], trait$Dispersal, mean))

multichoice <- trait[,c(2:5, 8:10)]

multiList <- list()

for (i in 1:ncol(multichoice))
{
  a <- tapply(pcoa_trait$vectors[,1], multichoice[,i], mean)[2]
  b <- tapply(pcoa_trait$vectors[,2], multichoice[,i], mean)[2]
  mat<- t(as.matrix(c(a, b), nrow = 1, ncol = 2));rownames(mat) <- colnames(multichoice)[i]
  multiList[[i]] <- mat
}

coord_bs <- cbind(weighted.mean(pcoa_trait$vectors[,1], log(trait$bodysize)),
                    weighted.mean(pcoa_trait$vectors[,2], log(trait$bodysize)))
rownames(coord_bs) <- "bodysize"


meantraits <- rbind(coord_getfood, coord_ingest, coord_disp, 
                    do.call(rbind, multiList), coord_bs)%>%as.data.frame
colnames(meantraits) <- c("A1", "A2")


p + geom_point(data = meantraits, mapping = aes(A1, A2), size = 2, 
               pch=3, color = "blue") +
  geom_label_repel(data = meantraits, aes(A1, A2), label=rownames(meantraits),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', color = "blue", size = 2) 






# functional dendrogram 


#https://yulab-smu.top/treedata-book/chapter10.html

tree <- as.phylo(hclust(distw.tot)) # dendrogram

df1 <- (sapply(by(comm,Hab,colSums),identity))%>%as.data.frame
df1 <- df1/apply(df1, 1, sum) #relative abundance per species across hab
df1$ID <- rownames(df1)
dat1 <- pivot_longer(df1, cols = c("EXO","INT","NAT","SEM"), names_to="Habitats", 
                     values_to = "Relative_abundance")
dat1$Habitats <- factor(dat1$Habitats, levels = c("NAT","EXO","SEM","INT"))

dat2 <- data.frame(ID=rownames(df1), origin, ab=apply(comm, 2, sum))
dat2$LogAb <- log10(dat2$ab) # log transform abundance




g <- ggtree(tree, layout="fan", size=0.3, open.angle=5, color="grey70") + 
  geom_fruit(data=dat1, geom=geom_tile,
             mapping=aes(y=ID, x=Habitats, alpha=Relative_abundance, fill=Habitats),
             color = "grey50", offset = 0.04,size = 0.02) + 
  scale_fill_manual(values = c("green4", "green", "orange", "red")) + 
  new_scale_fill() +
  geom_fruit(data=dat2, geom=geom_bar,
             mapping=aes(y=ID, x=LogAb, fill=origin),
             pwidth=0.5, 
             orientation="y", 
             stat="identity") +
  scale_fill_manual(name="Origin", values = c("grey", "black"), 
                    labels=c("Indigenous", "Exotics")) 
g

g + layout_rectangular() 



# built a convex hull space #

df.pcoa <- pcoa_trait$vectors[,1:2]%>%as.data.frame

convexhull <- BAT::hull.build(rep(1, ncol(comm)), df.pcoa, axes=0)
plot(convexhull, xlab="PCoA 1", ylab="PCoA 2")


# built a kernel density

kernelgaussian <- BAT::kernel.build(rep(1, ncol(comm)), df.pcoa, axes=0)
plot(kernelgaussian)


######################################################################
######################################################################

# evaluate the quality of the space #

cor(distw.tot, cophenetic(tree)%>%as.dist)

cor(distw.tot, dist(pcoa_trait$vectors[,1:2]))

cor(distw.tot, dist(pcoa_trait$vectors[,1:4]))


tree.quality(distw.tot, tree) # between 0 (bad) and 1 (good)

hyper.quality(distw.tot, pcoa_trait$vectors[,1:2]) # between 0 (bad) and 1 (good)

hyper.quality(distw.tot, pcoa_trait$vectors[,1:4])  # between 0 (bad) and 1 (good)

######################################################################
###################### Framework dissimilarity #######################
######################################################################

## Species based metrics ##

# originality

orig_diss <- BAT::originality(apply(comm, 2, sum), distance = distw.tot, 
                              abund=F)%>%as.vector()

# uniqueness 

uniq_diss <- BAT::uniqueness(apply(comm, 2, sum), 
                             distance = distw.tot)%>%as.vector()

# contribution  X

ggplot(data.frame(orig_diss, origin), aes(origin, orig_diss)) + geom_boxplot() + theme_bw()
wilcox.test(orig_diss ~ origin)

ggplot(data.frame(uniq_diss, origin), aes(origin, uniq_diss)) + geom_boxplot() + theme_bw()
wilcox.test(uniq_diss ~ origin)

## Group-based metrics ##

# Functional richness  X

# function dispersion #

diss.disp.ind <- BAT::dispersion(sqrt(comm.ind), distance=distw.tot)
diss.disp.exo <- BAT::dispersion(sqrt(comm.exo), distance=distw.tot)

ggplot() + geom_boxplot(data = data.frame(diss.disp.ind[,1], Hab), 
                        mapping = aes(x=Hab,y=diss.disp.ind), fill="green4") + 
  theme_bw() + theme(legend.position = "none")

ggplot() + geom_boxplot(data = data.frame(diss.disp.exo[,1], Hab), 
                        mapping = aes(x=Hab,y=diss.disp.exo), fill="grey") + 
  theme_bw() + theme(legend.position = "none")

kruskal.test(diss.disp.ind[,1], Hab)
kruskal.test(diss.disp.exo[,1], Hab)


# functional evenness X

## between groups metrics ##

diss.beta.ind <- picante::comdist(comm.ind, as.matrix(distw.tot))
diss.beta.exo <- picante::comdist(comm.exo, as.matrix(distw.tot))

nmds.beta.ind <- vegan::metaMDS(diss.beta.ind, k=3)
nmds.beta.exo <- vegan::metaMDS(diss.beta.exo, k=2)

datascores.ind <- as.data.frame(vegan::scores(nmds.beta.ind, display=c("sites")))  #extract the site scores
datascores.exo <- as.data.frame(vegan::scores(nmds.beta.exo, display=c("sites")))  #extract the site scores

#Add/calculate spider diagram
scores.ind <- cbind(as.data.frame(datascores.ind), Hab)
centroids.ind <- aggregate(cbind(NMDS1, NMDS2) ~ Hab, data = scores.ind, FUN = mean)
seg.ind <- merge(scores.ind, setNames(centroids.ind, c('Hab','oNMDS1','oNMDS2')),by = 'Hab', sort = FALSE)

scores.exo <- cbind(as.data.frame(datascores.exo), Hab)
centroids.exo <- aggregate(cbind(NMDS1, NMDS2) ~ Hab, data = scores.exo, FUN = mean)
seg.exo <- merge(scores.exo, setNames(centroids.exo, c('Hab','oNMDS1','oNMDS2')),by = 'Hab', sort = FALSE)


mds1 <- ggplot(scores.ind, aes(x = NMDS1, y = NMDS2)) +
  geom_segment(data = seg.ind, mapping = aes(xend = oNMDS1, yend = oNMDS2), color = "green4") + 
  geom_point(data = centroids.ind, size = 4) + geom_point(color = "green4") +                                              
  theme_bw()  +geom_label(data=centroids.ind, aes(x=NMDS1, y=NMDS2, label=Hab), color="green4")


mds2 <- ggplot(scores.exo, aes(x = NMDS1, y = NMDS2)) +
  geom_segment(data = seg.exo, mapping = aes(xend = oNMDS1, yend = oNMDS2), color = "grey50") + 
  geom_point(data = centroids.exo, size = 4) + geom_point(color = "grey50") +                                       
  theme_bw()  +geom_label(data=centroids.exo, aes(x=NMDS1, y=NMDS2, label=Hab), color="grey50")


ggarrange(mds1, mds2, nrow=2)


# analyses #

adonis2(diss.beta.ind ~ Hab)
adonis2(diss.beta.exo ~ Hab)

######################################################################
####################### Framework dendrogram #########################
######################################################################

## Species based metrics ##

# originality 

orig_tree <- BAT::originality(apply(comm, 2, sum), tree=tree, abund=F)%>%as.vector()

# uniqueness 

uniq_tree <- BAT::uniqueness(apply(comm, 2, sum), tree=tree)%>%as.vector()

# contribution

cont_tree <- BAT::contribution(apply(comm, 2, sum), tree=tree, abund = F)%>%as.vector()



ggplot(data.frame(orig_tree, origin), aes(origin, orig_tree)) + geom_boxplot() + theme_bw()
wilcox.test(orig_tree ~ origin)

ggplot(data.frame(uniq_tree, origin), aes(origin, uniq_tree)) + geom_boxplot() + theme_bw()
wilcox.test(uniq_tree ~ origin)

ggplot(data.frame(cont_tree, origin), aes(origin, cont_tree)) + geom_boxplot() + theme_bw()
wilcox.test(cont_tree ~ origin)




## Group-based metrics ##

# richness

tree.richness.ind <- BAT::alpha(comm.ind, tree)
tree.richness.exo <- BAT::alpha(comm.exo, tree)

# dispersion

tree.disp.ind <- BAT::dispersion(comm.ind, tree)
tree.disp.exo <- BAT::dispersion(comm.exo, tree)

# evenness

tree.eve.ind <- BAT::evenness(sqrt(comm.ind), tree)
tree.eve.exo <- BAT::evenness(sqrt(comm.exo), tree)

# example with the three metrics


df.metrics <- data.frame(richness.ind=tree.richness.ind[,1],  richness.exo=tree.richness.exo[,1], 
                         disp.ind=tree.disp.ind[,1], disp.exo=tree.disp.exo[,1], 
                         eve.ind=tree.eve.ind[,1], eve.exo=tree.eve.exo[,1], Hab, plot=factor(1:36))


df.metrics$Hab <- factor(df.metrics$Hab, levels= c("NAT", "EXO", "SEM", "INT"))
df.metrics <- pivot_longer(df.metrics, cols = colnames(df.metrics)[1:6], names_to="metric")
df.metrics$origin <- do.call(rbind, strsplit(df.metrics$metric, "[.]"))[,2]
df.metrics$dimensions <- do.call(rbind, strsplit(df.metrics$metric, "[.]"))[,1]

ggplot() + geom_boxplot(data = df.metrics, 
                        mapping = aes(x=Hab,y=value, fill=origin)) + theme_bw() +
  scale_fill_manual(name='Origin', values = c("grey50", "green4"), label=c("Exotics", "Indigenous")) +
  facet_grid(dimensions ~ origin, scales = "free") + theme(legend.position = "none")





## between-groups metrics ##

tree.beta.ind <- BAT::beta(comm.ind, tree, func = "jaccard", abund=F)$Btotal
tree.beta.exo <- BAT::beta(comm.exo, tree, func = "jaccard", abund=F)$Btotal


nmds.beta.ind <- vegan::metaMDS(tree.beta.ind, k=2)
nmds.beta.exo <- vegan::metaMDS(tree.beta.exo, k=2)

datascores.ind <- as.data.frame(vegan::scores(nmds.beta.ind, display=c("sites")))  #extract the site scores
datascores.exo <- as.data.frame(vegan::scores(nmds.beta.exo, display=c("sites")))  #extract the site scores

#Add/calculate spider diagram
scores.ind <- cbind(as.data.frame(datascores.ind), Hab)
centroids.ind <- aggregate(cbind(NMDS1, NMDS2) ~ Hab, data = scores.ind, FUN = mean)
seg.ind <- merge(scores.ind, setNames(centroids.ind, c('Hab','oNMDS1','oNMDS2')),by = 'Hab', sort = FALSE)

scores.exo <- cbind(as.data.frame(datascores.exo), Hab)
centroids.exo <- aggregate(cbind(NMDS1, NMDS2) ~ Hab, data = scores.exo, FUN = mean)
seg.exo <- merge(scores.exo, setNames(centroids.exo, c('Hab','oNMDS1','oNMDS2')),by = 'Hab', sort = FALSE)


mds1 <- ggplot(scores.ind, aes(x = NMDS1, y = NMDS2)) +
  geom_segment(data = seg.ind, mapping = aes(xend = oNMDS1, yend = oNMDS2), color = "green4") + 
  geom_point(data = centroids.ind, size = 4) + geom_point(color = "green4") +                                              
  theme_bw()  +geom_label(data=centroids.ind, aes(x=NMDS1, y=NMDS2, label=Hab), color="green4")


mds2 <- ggplot(scores.exo, aes(x = NMDS1, y = NMDS2)) +
  geom_segment(data = seg.exo, mapping = aes(xend = oNMDS1, yend = oNMDS2), color = "grey50") + 
  geom_point(data = centroids.exo, size = 4) + geom_point(color = "grey50") +                                       
  theme_bw()  +geom_label(data=centroids.exo, aes(x=NMDS1, y=NMDS2, label=Hab), color="grey50")


ggarrange(mds1, mds2, nrow=2)


# analyses #

adonis2(tree.beta.ind ~ Hab)
adonis2(tree.beta.exo ~ Hab)



######################################################################
####################### Framework convex hull ########################
######################################################################

## Species based metrics ##

# originality X

# uniqueness X

# contribution 

cont_hull <- hull.contribution(convexhull) # convexhull of the global FD space

## Group-based metrics ##

# we first create the hull kernel

hull.comm.ind <- BAT::hull.build(comm.ind, df.pcoa, axes=0)
hull.comm.exo <- BAT::hull.build(comm.exo, df.pcoa, axes=0)

# richness

hull.richness.ind <- BAT::hull.alpha(hull.comm.ind)
hull.richness.exo <- BAT::hull.alpha(hull.comm.exo)

# dispersion X

# evenness X

## Between-groups metrics ##

hull.beta.ind <- BAT::hull.beta(hull.comm.ind, func = "jaccard", comp = FALSE)$Btotal
hull.beta.exo <- BAT::hull.beta(hull.comm.exo, func = "jaccard", comp = FALSE)$Btotal

# mds ...

######################################################################
####################### Framework kernel #############################
######################################################################

## Species-based metrics

# originality

orig_kde <- BAT::kernel.originality(kernelgaussian) # kernel density of the global FD space
 
# uniqueness X

# contribution 

cont_kde <- kernel.contribution(kernelgaussian, func = "neighbor")

## Group-based metrics ##

# we first create the object kernel

kernel.comm.ind <- BAT::kernel.build(comm.ind, df.pcoa, axes=0, abund = F) # with abund=T take times
kernel.comm.exo <- BAT::kernel.build(comm.exo, df.pcoa, axes=0, abund = F) # with abund=T take times

# richness 

kernel.richness.ind <- BAT::kernel.alpha(kernel.comm.ind)
kernel.richness.exo <- BAT::kernel.alpha(kernel.comm.exo)

# dispersion 

kernel.disp.ind <- BAT::kernel.dispersion(kernel.comm.ind)
kernel.disp.exo <- BAT::kernel.dispersion(kernel.comm.exo)

# evenness 

kernel.eve.ind <- kernel.evenness(kernel.comm.ind)
kernel.eve.exo <- kernel.evenness(kernel.comm.exo)

# example with the three metrics


df.metrics <- data.frame(richness.ind=kernel.richness.ind[,1],  richness.exo=kernel.richness.exo[,1], 
                         disp.ind=kernel.disp.ind[,1], disp.exo=kernel.disp.exo[,1], 
                         eve.ind=kernel.eve.ind[,1], eve.exo=kernel.eve.exo[,1], Hab, plot=factor(1:36))


df.metrics$Hab <- factor(df.metrics$Hab, levels= c("NAT", "EXO", "SEM", "INT"))
df.metrics <- pivot_longer(df.metrics, cols = colnames(df.metrics)[1:6], names_to="metric")
df.metrics$origin <- do.call(rbind, strsplit(df.metrics$metric, "[.]"))[,2]
df.metrics$dimensions <- do.call(rbind, strsplit(df.metrics$metric, "[.]"))[,1]

ggplot() + geom_boxplot(data = df.metrics, 
                        mapping = aes(x=Hab,y=value, fill=origin)) + theme_bw() +
  scale_fill_manual(name='Origin', values = c("grey50", "green4"), label=c("Exotics", "Indigenous")) +
  facet_grid(dimensions ~ origin, scales = "free") + theme(legend.position = "none")


## Between groups metrics ##

kernel.beta.ind <- BAT::kernel.beta(kernel.comm.ind, func = "jaccard")$Btotal
kernel.beta.exo <- BAT::kernel.beta(kernel.comm.exo, func = "jaccard")$Btotal

# mds ...

##############################################################################
##############################################################################
###########----------------Metric for a single trait----------------##########
##############################################################################
##############################################################################

cwm.ind <- BAT::cwm(sqrt(comm.ind), trait[colnames(comm.ind),], abund=TRUE)%>%as.data.frame
cwm.exo <- BAT::cwm(sqrt(comm.exo), trait[colnames(comm.exo),], abund=TRUE)%>%as.data.frame

colnames(cwm.ind)
colnames(cwm.exo)

cwm.ind$Hab <- factor(Hab, levels = c("NAT", "EXO", "SEM", "INT"))
cwm.exo$Hab <- factor(Hab, levels = c("NAT", "EXO", "SEM", "INT"))

ggplot() + geom_boxplot(data = cwm.ind, 
                        mapping = aes(x=Hab,y=bodysize), fill="green4") + theme_bw() +
  labs(x="", y="CWM")  + ggtitle("Body size") + ylim(c(0, 25))


ggplot() + geom_boxplot(data = cwm.exo, 
                        mapping = aes(x=Hab,y=bodysize), fill="grey50") + theme_bw() +
  labs(x="", y="CWM") + ggtitle("Body size") + ylim(c(0, 25))


########################################################################
########################################################################
##########------------trait-based community assembly---------###########
########################################################################
########################################################################


# randomize tips #

picante::tipShuffle(tree)

picante::taxaShuffle(distw.tot)%>%as.dist

?randomizeMatrix

picante::randomizeMatrix(comm, null.model = "independentswap")

# better to do the randomization within origin #

data.ind <- picante::match.phylo.comm(tree, comm.ind)
data.exo <- picante::match.phylo.comm(tree, comm.exo)


null.tree.ind.mat <- replicate(1000, BAT::alpha(randomizeMatrix(sqrt(data.ind$comm), null.model = "independentswap"), 
                                                    tree=data.ind$phy)[,1])

null.tree.exo.mat <- replicate(1000, BAT::alpha(randomizeMatrix(sqrt(data.exo$comm), null.model = "independentswap"),
                                                    tree=data.exo$phy)[,1])


# check the null distribution for one community #

hist(null.tree.ind.mat[,30], breaks=10)
segments(null.tree.ind.mat[30], 0, null.tree.ind.mat[30], 10)


## compare to the observed values ##

tree.richness.ind # obs 
apply(null.tree.ind.mat, 1, mean) # mean null
apply(null.tree.ind.mat, 1, sd) # sd null

ses.ind <- (tree.richness.ind - apply(null.tree.ind.mat, 1, mean))/apply(null.tree.ind.mat, 1, sd)

ses.exo <- (tree.richness.exo - apply(null.tree.exo.mat, 1, mean))/apply(null.tree.exo.mat, 1, sd)

df.ses.ind <- data.frame(ses.ind, Hab);df.ses.ind$Hab <- factor(df.ses.ind$Hab, levels = c("NAT", "EXO", "SEM", "INT"))
df.ses.exo <- data.frame(ses.exo, Hab);df.ses.exo$Hab <- factor(df.ses.exo$Hab, levels = c("NAT", "EXO", "SEM", "INT"))


ggplot() + geom_boxplot(data = df.ses.ind, 
                        mapping = aes(x=Hab,y=ses.ind), col="green4") + theme_bw() +
  labs(x="", y="SES Dispersion") + geom_hline(yintercept = 0, lty=3)


ggplot() + geom_boxplot(data = df.ses.exo, 
                        mapping = aes(x=Hab,y=ses.exo), col="grey50") + theme_bw() +
  labs(x="", y="SES Dispersion") + geom_hline(yintercept = 0, lty=3)


tapply(ses.ind, Hab, wilcox.test)

tapply(ses.exo, Hab, wilcox.test)





