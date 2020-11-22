library(here)

timeStamp <- "DE (2020-07-25_141629)"

#genes with GO for stress response
stress <- read.csv(here("supplementary notebooks","outputs","geneID_isDescendant.csv"))


###
#For interaction effect
###
cagVbat <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_BAT.E.csv")) %>% merge(stress,by="gene_id")
cagVbic <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_BIC.E.csv")) %>% merge(stress,by="gene_id")
bicVbat <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_BAT.E - BIC.E.csv")) %>% merge(stress,by="gene_id")

write.csv(cagVbat,file=here("results","DE",timeStamp,"DE gene list with stress info","cagVbat_stress.csv"))
write.csv(cagVbic,file=here("results","DE",timeStamp,"DE gene list with stress info","cagVbic_stress.csv"))
write.csv(bicVbat,file=here("results","DE",timeStamp,"DE gene list with stress info","bicVbat_stress.csv"))


###
#For main effect 
cag <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_E.csv")) %>% merge(stress,by="gene_id")
bat <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_E + BAT.E.csv")) %>% merge(stress,by="gene_id")
bic <- read.csv(here("results","DE",timeStamp,"DE gene list (unique)","DEgenes_E + BIC.E.csv")) %>% merge(stress,by="gene_id")


write.csv(cag,file=here("results","DE",timeStamp,"DE gene list with stress info","cag_stress.csv"))
write.csv(bat,file=here("results","DE",timeStamp,"DE gene list with stress info","bat_stress.csv"))
write.csv(bic,file=here("results","DE",timeStamp,"DE gene list with stress info","bic_stress.csv"))



