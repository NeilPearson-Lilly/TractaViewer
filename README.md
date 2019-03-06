## TractaViewer

A tool to automatically download and aggregate data relating to the suitability of a list of genes as drug development targets from a large number of databases, and present an overview of that data in a comprehensible manner. 

-----

### Installation

A precompiled binary (Windows 64-bit) is available by unzipping lauch.zip.001. (This executable file should be in the project root directory - the same place the zipped files are now). No further installation is required, though you may wish to make a shortcut pointing to this executable. We strongly recommend running TractaViewer in this manner if possible. 

If running from source as a conventional Python script, run launch.py. This will require a Python3 installation with all dependencies installed, including Qt and its PyQt bindings. We recommend using the Anaconda distribution of Python, which comes with most of TractaViewer's dependencies pre-loaded. However, you will still be required to install at least the Qt5 GUI library yourself. The full list of dependencies can be read [here](Requirements.md).

### Data sources

The majority of sources used by TractaViewer provide APIs, allowing queries without having to store any files locally. 

In some cases, however, it is more practical to maintain local versions of data sources. These cases are listed in `Data/Data_manifest.xlsx`.

Generally, this locally-maintained data is simply a file directly downloaded from the source, but a small number of cases require the construction of new datasets from existing files. These cases, and the scripts used to produce these new datasets, are described [here](Constructed_data_sources.md).

### API keys

If you plan to query 100 targets or more per day, we recommend acquiring API keys from the OpenFDA portal. If invalid API keys are supplied and excessive attempts are made to access OpenFDA information, TractaViewer will continue to mine other sources, but will leave fields corresponding to those sources blank in its output.

API keys can be obtained for free by registering your email address at the following links:

[OpenFDA](https://open.fda.gov/api/reference/)

Previous versions of TractaViewer also recommended using an API key to access OpenTargets, but OpenTargets have since removed the need for API keys when making large numbers of queries.  

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
 
### Output

A full data dictionary for the output file produced by TractaViewer can be read [here](Fields.md).

##### Example output

An example of mined data for a small set of genes (a subset of a publicly accessible list from the AMP-AD consortium) is available in the example_output directory.

##### Shiny app visualisations

Two Shiny apps are provided to facilitate visualisations of data retrieved using TractaViewer. These may be set up as web tools given the appropriate infrastructure.

The first, Bucket Visualiser, provides a means of understanding the drug target landscape of all input genes by plotting various bucket scores in a 3D scatter plot.

The second, Target Visualiser, is a dashboard for a variety of information at the level of a single gene. TractaViewer output (covering what is expected to be a large number of genes) is intended to be parsed into an SQLite database using a provided script, which the app then queries as required.

### Provenance

The source databases for each item of data in output files are listed and described in "Data/Data manifest.xlsx".

### Compiling

This project can be compiled into a binary executable using PyInstaller. The path to PyQt .dll libraries must be included, because they aren't included by default. The `--noconsole` flag removes the terminal window that would otherwise appear, causing the application to be viewed purely as a GUI tool to the user.

```PyInstaller --path C:\Path\To\Python\Installation\Lib\site-packages\PyQt5\Qt\bin --clean --noconsole --onefile -i logo2.ico launch.py```

Note that if compiling, data and image assets must be copied into the resulting application's root directory, in the same arrangement as present in this repository. 

