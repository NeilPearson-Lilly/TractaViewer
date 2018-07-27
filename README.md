## TractaViewer

A tool to automatically download and aggregate data relating to the suitability of a list of genes as drug development targets from a large number of databases, and present an overview of that data in a comprehensible manner. 

-----

### Installation

A precompiled binary (Windows 64-bit) is available by unzipping lauch.zip.001. (This executable file should be in the project root directory - the same place the zipped files are now). No further installation is required, though you may wish to make a shortcut pointing to this executable. 

If running from source as a conventional Python script, run launch.py. This will require a Python3 installation with all dependencies installed, including Qt and its PyQt bindings. We recommend using the Anaconda distribution of Python, which comes with most of TractaViewer's dependencies pre-loaded.

### API keys

If you plan to query 100 targets or more per day, we recommend acquiring API keys from OpenTargets and the OpenFDA portal. If invalid API keys are supplied and excessive attempts are made to access OpenTargets and OpenFDA information, TractaViewer will continue to mine other sources, but will leave fields corresponding to those sources blank in its output.

API keys can be obtained for free by registering your email address at the following links:

OpenTargets: https://blog.opentargets.org/api-getting-started-1/

OpenFDA: https://open.fda.gov/api/reference/

### Usage

TractaViewer is capable of accepting lists of targets in the form of HGNC names, Ensembl IDs, UniProt IDs, or any combination of those 3. Choose an Excel or CSV file with this information and TractaViewer will automatically pick out the appropriate columns, provided that their column headers are one of "HGNC Name", "Ensembl ID", or "UniProt ID". TractaViewer can recognise some slight variations on these headers. 

You may also input a list of targets by copy-pasting text as input, but in this case you must specify what kind of identifiers you are providing.

It is possible to supply more than one input file or target list. Subsequent target lists will be appended to the bottom of the existing list shown in the main window.

##### New modality check

Users may, optionally, choose to specify whether modalities involving protein degradation and antagonism/inhibition of protein expression are appropriate for a given target. (Since this is determined largely by therapeutic intent, this is best specified by the user). To specify a target as suitable for these methods, create a column named “New modality check” in the input sheet. Any value in this column indicates suitability for new modalities; leaving it blank indicates unsuitability. TractaViewer will function as normal if this column is not supplied, but the range of New modality bucket outcomes will be more limited. 

##### Disease profiles

We generally recommend that you set up a disease profile. This facilitates a better calculation of the relative safety scores for your targets. To set up a disease profile, click the "Change" button in the main window, and select the pertinent tissues and cell types (if any) for your disease in the dialog that appears. Give your disease profile a memorable name. This profile will be used by default until you create a different one; previous profiles are retained, and remain available after you switch to another.   

Once your list of targets is complete, begin data-mining by clicking the "Download druggability data" button. A dialog will appear, asking where you wish to save output. Progress through your target list is then indicated at the bottom of the main window.

Upon completion of data-mining, a large tabbed data-browser dialog will appear. This can be used to re-order the data according to certain findings, to search for data from certain genes on the list, or to select a subset of genes for the production of visual HTML reports. When this happens, mined data is also saved to the location specified earlier, as an Excel file. The data-browser can be re-opened for older output using the "Open existing druggability data" button. This will also re-populate the input area with the original list of targets.
 

### Example output

An example of mined data for a small set of genes (a subset of a publicly accessible list from the AMP-AD consortium) is available in the example_output directory.

### Provenance

The source databases for each item of data in output files are listed and described in "Data/Data manifest.xlsx".

### Compiling

This project can be compiled into a binary executable using PyInstaller. The path to PyQt .dll libraries must be included, because they aren't included by default. The `--noconsole` flag removes the terminal window that would otherwise appear, causing the application to be viewed purely as a GUI tool to the user.

```PyInstaller --path C:\Path\To\Python\Installation\Lib\site-packages\PyQt5\Qt\bin --clean --noconsole --onefile -i logo2.ico launch.py```

Note that if compiling, data and image assets must be copied into the resulting application's root directory, in the same arrangement as present in this repository. 

