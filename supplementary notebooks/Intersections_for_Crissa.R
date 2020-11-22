library(here)
library(RVenn)
library(dplyr)
library(magrittr)

#genes with GO for stress response
stress <- read.csv(here("supplementary notebooks","outputs","geneID_isDescendant.csv"))

timeStamp <- "DE (2020-07-25_141629)"

###
#For interaction coefficients
###
cagVbat <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_BAT.E.csv")) %>% merge(stress,by="gene_id")
cagVbic <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_BIC.E.csv")) %>% merge(stress,by="gene_id")
bicVbat <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_BAT.E - BIC.E.csv")) %>% merge(stress,by="gene_id")

# 2-way intersections
cagVbat_cagVbic <- select(merge(cagVbat,cagVbic, by = "gene_id"),c(gene_id,fly.x,fly_acc.x,isDescendant.x))
cagVbat_bicVbat <- select(merge(cagVbat,bicVbat, by = "gene_id"),c(gene_id,fly.x,fly_acc.x,isDescendant.x))
cagVbic_bicVbat <- select(merge(cagVbic,bicVbat, by = "gene_id"),c(gene_id,fly.x,fly_acc.x,isDescendant.x))

# 3-way intersections
cagVbat_cagVbic_bicVbat <- select(merge(merge(cagVbat,cagVbic, by = "gene_id"),bicVbat, by = "gene_id"),
                                  c(gene_id,fly.x,fly_acc.x,isDescendant.x)) #could've just used one of the smaller sets, but whatever
#write them out
write.csv(cagVbat_cagVbic,file=here("results","DE",timeStamp,"DE gene intersections","interaction","cagVbat_cagVbic.csv"))
write.csv(cagVbat_bicVbat,file=here("results","DE",timeStamp,"DE gene intersections","interaction","cagVbat_bicVbat.csv"))
write.csv(cagVbic_bicVbat,file=here("results","DE",timeStamp,"DE gene intersections","interaction","cagVbic_bicVbat.csv"))
write.csv(cagVbat_cagVbic_bicVbat,file=here("results","DE",timeStamp,"DE gene intersections","interaction","cagVbat_cagVbic_bicVbat.csv"))


###
#For main coefficients
###
cag <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_E.csv")) %>% merge(stress,by="gene_id")
bat <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_E + BAT.E.csv")) %>% merge(stress,by="gene_id")
bic <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_E + BIC.E.csv")) %>% merge(stress,by="gene_id")

#2-way intersections
cag_bat <- select(merge(cag,bat, by = "gene_id"),c(gene_id,fly.x,fly_acc.x,isDescendant.x))
cag_bic <- select(merge(cag,bic, by = "gene_id"),c(gene_id,fly.x,fly_acc.x,isDescendant.x))
bat_bic <- select(merge(bat,bic, by = "gene_id"),c(gene_id,fly.x,fly_acc.x,isDescendant.x))

#3-way
cag_bat_bic <- select(merge(merge(cag,bat, by = "gene_id"),bic, by = "gene_id"),
                                  c(gene_id,fly.x,fly_acc.x,isDescendant.x)) #could've just used one of the smaller sets, but whatever


#write them out
write.csv(cag_bat,file=here("results","DE",timeStamp,"DE gene intersections","main","cag_bat.csv"))
write.csv(cag_bic,file=here("results","DE",timeStamp,"DE gene intersections","main","cag_bic.csv"))
write.csv(bat_bic,file=here("results","DE",timeStamp,"DE gene intersections","main","bat_bic.csv"))
write.csv(cag_bat_bic,file=here("results","DE",timeStamp,"DE gene intersections","main","cag_bat_bic.csv"))





