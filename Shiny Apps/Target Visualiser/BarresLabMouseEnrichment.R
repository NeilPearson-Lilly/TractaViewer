library(tidyverse)
library(TissueEnrich)
library(SummarizedExperiment)

# Calculate enrichment (equivalent to method used by the Human Protein Atlas) 
# on Barres lab mouse data.

bldata = "../../Data/barreslab_rnaseq (2).xlsx"
sheets <- openxlsx::getSheetNames(bldata)
xl <- lapply(sheets,openxlsx::read.xlsx,xlsxFile=bldata)
# We only actually want the mouse data here though...
df = as.data.frame(xl[1])

# This is for an example from the vignetter, which definitely works.
# data<-system.file("extdata", "test.expressiondata.txt", package = "TissueEnrich")
# expressionData<-read.table(data,header=TRUE,row.names=1,sep='\t')

# This is for my own data.
expressionData = subset(df, select=3:9)
rownames(expressionData) = df$Gene.symbol
se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),
                         rowData = row.names(expressionData),
                         colData = colnames(expressionData))

output<-teGeneRetrieval(se, maxNumberOfTissues = 6)
enrichments = as.data.frame(assay(output))
enrichments$Group <- gsub('-', ' ', enrichments$Group)
enrichments$Tissue <- gsub('[.]', ' ', enrichments$Tissue)
enrichments = enrichments[enrichments$Group %in% c("Tissue Enriched", "Group Enriched", "Tissue Enhanced"),]
enrichments$Group <- gsub('Enriched', 'enriched', enrichments$Group)
enrichments$Group <- gsub('Enhanced', 'enhanced', enrichments$Group)
# Looks good! Write it to a file.
write.csv(enrichments, "BarresLab_Enrichments.csv", row.names = FALSE)
