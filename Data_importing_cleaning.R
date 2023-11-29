meta17=read.csv("~/Documents/DE/DE_16S/Oyster_16S/Data/metadata_de17.csv")
data=read.csv("~/Documents/DE/DE_16S/Oyster_16S/Data/DE_DATA_ForGenetics.csv")

##renaming treatment codes
data$Treatment
data$Treatment2=ifelse(data$Treatment=="HH", "HIGH_POLY", 
                       ifelse(data$Treatment=="HL", "HIGH_MONO", 
                              ifelse(data$Treatment=="LL", "LOW_MONO", "LOW_POLY")))

#create unique IDs
2017_NW_HIGH_MONO_B11_CV

data$Colornumber=paste0(data$Color, data$Number)

data$UniqueID=paste("2017", data$Site, data$Treatment2, data$Colornumber, data$Species, sep="_")

meta17_data=merge(meta17, data, by="UniqueID")
write.csv(meta17_data, file="~/Documents/DE/DE_16S/ASVs/cleanedmetadata17.csv")
meta17_data=read.csv("~/Documents/DE/DE_16S/ASVs/cleanedmetadata17.csv")
rownames(meta17_data)=meta17_data$UniqueID