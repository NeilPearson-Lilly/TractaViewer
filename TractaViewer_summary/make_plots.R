# A test thing I'm using to try to make sense of the actual code to read files and produce plots.
# library(readxl)
library(openxlsx)

# An example file...
fl = "target_list v7 short_druggability.xlsx"

# This just reads the first sheet
# xl = read_excel(fl)
# This gets everything!
sheets <- openxlsx::getSheetNames(fl)
xl <- lapply(sheets,openxlsx::read.xlsx,xlsxFile=fl)
names(xl) <- sheets

