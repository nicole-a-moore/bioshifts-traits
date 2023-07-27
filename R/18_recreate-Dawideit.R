# recreating models from Dawideit et al.
library(tidyverse)

## read ecomorph data from the supplement 
ecomorph <- read.csv("data-raw/dispersal/Dawideit_et_al_ecomorph.csv")
ecomorph[ecomorph == "-"] <- NA

## read dispersal data from Paradis
paradis <- read.csv("data-raw/dispersal/Paradis_et_al_2002.csv")

## fix names in Paradis so they match those in ecomorph data
paradis[paradis == "Delichon urbica"] <- "Delichon urbicum"
paradis[paradis == "Parus montanus"] <- "Poecile montanus"
paradis[paradis == "Parus ater"] <- "Periparus ater"
paradis[paradis == "Parus caeruleus"] <- "Cyanistes caeruleus"

## join the databases 
join <- inner_join(ecomorph, paradis)

## ln transform traits to normalize, as done in Dawideit
## note: maybe this is where the difference lies? paper reports variables were either ln(x) or ln(x+1) transformed
join$Kipp.s.distanceLn <- log(join$Kipp.s.distance)
join$Bill.depthLn <- log(join$Bill.depth)

## transform geometric mean of dispersal distances, as done in Dawideit
## ln[ln(x + 1) + 1]
#join$GM_natal_transformed <- log(log(join$GM_natal + 1)+1)
join$GM_natal_transformed <- log1p(log1p(join$GM_natal))

# hist(join$GM_natal)
# hist(join$GM_natal_transformed)

## fit the simple linear two-predictor model
mod <- lm(GM_natal_transformed ~ Kipp.s.distanceLn + Bill.depthLn,
          data = join)
summary(mod)
# Kipp.s.distanceLn  0.30421    0.07894   3.854 0.000374 ***
# Bill.depthLn      -0.35762    0.09300  -3.845 0.000384 ***

## they report:
# Kipp's distance = 0.368
# bill depth = -0.402 


## try phylogenetic comparative model 
## get the phylogenetic tree for British Birds 
library(caper) 
data(BritishBirds)

BritishBirds.tree
join$Species <- str_replace_all(join$Species, " ",  ".")

## make sure all species in ecomorph data are in tree 
length(which(join$Species %in% BritishBirds.tree$tip.label))
## yay - they are

## make comparative object
join_comp <- comparative.data(BritishBirds.tree, join, Species, vcv = TRUE)

## fit the GLS using the phylogenetic comparative method
mod_pgls <- pgls(GM_natal_transformed ~ Kipp.s.distanceLn + Bill.depthLn,
                 lambda = "ML",
                 data = join_comp)

summary(mod_pgls)

## still does not match what they report:
# Kipp's distance = 0.368
# bill depth = -0.402 


## July 20
###################################
#### inspecting data differences 
###################################
#reads in the dataset put together by Britta Dawideit
britta<-read.csv("data-raw/AllyPhillimore/mean.morph&dispersal.txt",sep="\t",header=T)

britta$notch[which(is.na(britta$notch)==T)]<-0
britta$bristle[which(is.na(britta$bristle)==T)]<-0

## check for differences in Kipp's distance / bill depth values between Britta and ecomorph 
which(!ecomorph$Bill.depth %in% britta$billd)
which(!ecomorph$Kipp.s.distance %in% britta$Kipp)

## one difference found:
ecomorph$Kipp.s.distance[36] ## 155.56
britta$Kipp[36] ## 176.375


## what about the tree?
library(ape)
tree<-read.tree("data-raw/AllyPhillimore/Brit_birds_all_MCC_noft.tre")
unwanted.tips<-which(is.na(pmatch(tree$tip.label,as.character(britta[,1])))==T)
bbtree<-drop.tip(tree,unwanted.tips)
vcvtree<-vcv.phylo(bbtree)
reorder.vars<-pmatch(bbtree$tip.label,as.character(britta[,1]))

## do the same thing Ally does to my tree:
unwanted.tips <-which(is.na(pmatch(BritishBirds.tree$tip.label,as.character(britta[,1])))==T)
bbtree_mine <- drop.tip(BritishBirds.tree, unwanted.tips)
vcvtree_mine <-vcv.phylo(bbtree_mine)
reorder.vars<-pmatch(bbtree$tip.label,as.character(britta[,1]))
## they look the same 

## fix difference in Kipp's distance and see if anything changes
join$Kipp.s.distanceLn[36] <- log(britta$Kipp[36])

join_comp <- comparative.data(bbtree_mine, join, Species, vcv = TRUE)

## fit the GLS using the phylogenetic comparative method
mod_pgls <- pgls(GM_natal_transformed ~ Kipp.s.distanceLn + Bill.depthLn,
                 lambda = "ML",
                 data = join_comp)

summary(mod_pgls)
# Kipp.s.distanceLn  0.45227    0.11273  4.0118 0.0003017 ***
# Bill.depthLn      -0.65730    0.15086 -4.3569 0.0001099 ***

## run pglmEstLambda on my data
# reformat data
reorder.vars<-pmatch(bbtree_mine$tip.label,as.character(join[,1]))

#need to appropriately transform varables
join1<-cbind(join$GM_natal_transformed, join$Kipp.s.distanceLn, join$Bill.depthLn)[reorder.vars,]
#note the double logging of dispersal

row.names(join1)<-row.names(vcvtree_mine)
colnames(join1) <- c("GM_natal_transformed", "Kipp.s.distanceLn", "Bill.depthLn")

join1<-data.frame(join1)

model_mine <-pglmEstLambda(join1$GM_natal_transformed ~ join1$Kipp.s.distanceLn + join1$Bill.depthLn,
                           join1,
                     phylomat=vcvtree_mine)
summary(model_mine)

# Term		Estimate		Std Err		T-value	P
# (Intercept) 	0.1949582 	 0.1874538 	 1.040033 	 0.3040058 
# join1$Kipp.s.distanceLn 	0.3762444 	 0.06876041 	 5.471817 	 2.002713e-06 
# join1$Bill.depthLn 	-0.4151167 	 0.08067732 	 -5.145396 	 5.956557e-06 

## YAY! now my numbers are the same as Ally's. 
## but.... still different from the paper
## try leaving the difference in Kipp's distance and seeing if numbers match 
join$Kipp.s.distanceLn[36] <- log(ecomorph$Kipp[36]) ## sets value back to the original value in ecomorph 

## run pglmEstLambda on my data
# reformat data
reorder.vars<-pmatch(bbtree_mine$tip.label,as.character(join[,1]))

#need to appropriately transform varables
join1<-cbind(join$GM_natal_transformed, join$Kipp.s.distanceLn, join$Bill.depthLn)[reorder.vars,]
#note the double logging of dispersal

row.names(join1)<-row.names(vcvtree_mine)
colnames(join1) <- c("GM_natal_transformed", "Kipp.s.distanceLn", "Bill.depthLn")

join1<-data.frame(join1)

model_mine <-pglmEstLambda(join1$GM_natal_transformed ~ join1$Kipp.s.distanceLn + join1$Bill.depthLn,
                           join1,
                           phylomat=vcvtree_mine)
summary(model_mine)

# Term		Estimate		Std Err		T-value	P
# (Intercept) 	0.202499 	 0.1914794 	 1.05755 	 0.2960324 
# join1$Kipp.s.distanceLn 	0.3675882 	 0.07035191 	 5.224992 	 4.570975e-06 
# join1$Bill.depthLn 	-0.4020028 	 0.081652 	 -4.923368 	 1.241377e-05 

### YAY x 2 - these estimates are the same as those in the paper 

