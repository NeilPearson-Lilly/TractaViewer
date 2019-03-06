library(tidyverse)
library(DBI)
library(RSQLite)

con = dbConnect(SQLite(), dbname="proteome.db")

datadir = "C/Your/Data/Directory"

datadir_files = list.files(datadir, pattern = "*.xlsx")

for (fl in 1:5) {
  print(paste("READING FILE", fl))
  sheets <- openxlsx::getSheetNames(fl)
  xl <- lapply(sheets,openxlsx::read.xlsx,xlsxFile=fl)
  names(xl) <- sheets
  print("FILE READ")
  
  print("WRITING TO DATABASE")
  for (i in sheets) {
    print(i)
    df = as.data.frame(xl[i], stringsAsFactors = FALSE)
    colnames(df) = lapply(colnames(df), function (x) {sub(paste0("^", make.names(i), '.'), '', x)})
    dbWriteTable(con, i, df, row.names = FALSE, append = TRUE)
  }
  print("DATA WRITTEN")
}

# Couple more could do with adding to the DB too
barreslab_enrichments = read.csv("BarresLab_Enrichments.csv")
dbWriteTable(con, "barreslab_enrichments", barreslab_enrichments, row.names = FALSE, append = TRUE)
# interactions = read.csv("IUPHAR_ligand-gene_interactions.csv")
# dbWriteTable(con, "interactions", interactions, row.names = FALSE, append = TRUE)

