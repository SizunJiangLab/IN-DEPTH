wdpath <- "./data/DLBCL_run/"
loadingfunction_wdpath <- c("./src/")
library(enrichR)
library(readxl)
library(openxlsx)
library(dplyr)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
if (websiteLive) head(dbs)
grep("reactome",dbs$libraryName,ignore.case = T,value = T)
grep("GO",dbs$libraryName,ignore.case = T,value = T)
grep("KEGG",dbs$libraryName,ignore.case = T,value = T)
dbs <- c("Reactome_2022", "GO_Biological_Process_2023", "KEGG_2021_Human")

file.list <- list.files(path = file.path(wdpath,"GeoMX_differential_expressed_analysis"),pattern = ".csv")

for (i in file.list){
  ReadinFile <- read.csv(file.path(wdpath,"GeoMX_differential_expressed_analysis",i))
  PositiveIn <- unique(ReadinFile$positive_in)
  ReadinFile$X <- gsub("[.]","-",ReadinFile$X)
  wb <- createWorkbook()
  for (j in PositiveIn){
    print(j)
    tmp.gene <- ReadinFile$X[ReadinFile$positive_in == j]
    enriched <- enrichr(tmp.gene, dbs)
    Reactome_2022 <- enriched$Reactome_2022 %>% filter(P.value < 0.05)
    print(Reactome_2022[1:5,])
    GO_Biological_Process_2023 <- enriched$GO_Biological_Process_2023 %>% filter(P.value < 0.05)
    KEGG_2021_Human <- enriched$KEGG_2021_Human %>% filter(P.value < 0.05)
    addWorksheet(wb, paste0(j,"Reactome_2022"))
    addWorksheet(wb, paste0(j,"GO_Biological_Process_2023"))
    addWorksheet(wb, paste0(j,"KEGG_2021_Human"))
    writeData(wb, sheet = paste0(j,"Reactome_2022"), Reactome_2022)
    writeData(wb, sheet = paste0(j,"GO_Biological_Process_2023"), GO_Biological_Process_2023)
    writeData(wb, sheet = paste0(j,"KEGG_2021_Human"), KEGG_2021_Human)
  }
  file.name.output <- gsub(".csv","Pathway.xlsx",i) 
  saveWorkbook(wb, file.path(wdpath,"GeoMX_differential_expressed_analysis",
                             file.name.output), overwrite = TRUE)
}



