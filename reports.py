import pandas as pd
import pweave
import bokeh
from pweave.bokeh import output_pweave, show
import os
from PyQt5.QtCore import QThread, pyqtSignal

# output_pweave()

# pweave.weave("FIR_designp.md", doctype='md2html', informat='markdown', output="FIR_design.html")
# Neat! Now try for my example (with some pandas!)
# pweave.weave("report_test.md", doctype='md2html', informat='markdown', output="rt.html")
# pweave.weave("report_test.ipynb", doctype='md2html', informat='notebook', output="rt.html")


# For the moment, we'll use that as an example. Let's say we want to produce an output for each gene in that list. 
class Reporter(QThread):
    progbar_update = pyqtSignal(int)
    def __init__(self, targets, outpath=None):
        super(QThread, self).__init__()
        self.targets = targets
    
    def make_reports(self):
        for i, row in targets['Basic information'].iterrows():
            # Write some temporary CSV files, invoke the report markdown thing, and then delete the temp files.
            # That's messy, but it's one of the quicker ways of doing this - especially when dealing with a lot of 
            # targets. 
            # I thought about doing some input (for the report thing) verification, but a better approach would be to 
            # make it robust against anything that's going to come out of here. Any issues should have been sorted out
            # long ago by this point.
            print(row['HGNC Name'])
            temp_gtex_file = "temp/gtex.csv"
            targets['GTEX'].to_csv(temp_gtex_file)
            temp_basic_info = "temp/basic.csv"
            targets['Basic information'].to_csv(temp_basic_info)
            temp_pharos = "temp/pharos.csv"
            targets['Pharos'].to_csv(temp_pharos)
            temp_barres_mouse = "temp/barres_mouse.csv"
            targets['Barres mouse'].to_csv(temp_barres_mouse)
            temp_barres_human = "temp/barres_human.csv"
            targets['Barres human'].to_csv(temp_barres_human)
            temp_disease_assocs = "temp/disease_assocs.csv"
            targets['Disease associations'][targets['Disease associations']['series'] == row['series']].to_csv(
                temp_disease_assocs)
            temp_drugs = "temp/drugs.csv"
            targets['Existing drugs'][targets['Existing drugs']['series'] == row['series']].to_csv(temp_drugs)
            temp_druggability = "temp/druggability.csv"
            targets['Druggability'].to_csv(temp_druggability)
            temp_biopharm = "temp/biopharm.csv"
            targets['Biopharmability'].to_csv(temp_biopharm)
            temp_risk = "temp/risk.csv"
            targets['Risk factors'].to_csv(temp_risk)
            
            html_report_file = "Example_output/" + row['HGNC Name'] + ".html"
            # pweave.weave("report_test.md", doctype='md2html', informat='markdown', output=html_report_file)
            pweave.weave("summary_report.md", doctype='md2html', informat='markdown', output=html_report_file)
            
            os.remove(temp_gtex_file)
            os.remove(temp_basic_info)
            os.remove(temp_pharos)
            os.remove(temp_barres_mouse)
            os.remove(temp_barres_human)
            os.remove(temp_disease_assocs)
            os.remove(temp_drugs)
            os.remove(temp_druggability)
            os.remove(temp_biopharm)
            os.remove(temp_risk)
            print("HOLD UP - TESTING MODE")
            input()
    

if __name__ == "__main__":
    targets = pd.read_excel(
        "Example_output/sample_output_druggability.xlsx",
        sheetname=None)
    

