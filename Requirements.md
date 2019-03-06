# Required modules for TractaViewer operation

We strongly recommend that TractaViewer is run in a Windows environment using the supplied executable, since setting up these dependencies may be difficult in some cases. However, we recognise that this is not always possible.

In order to run TractaViewer from source, you must use Python V3.6 or lower. At the time of writing, V3.7 has an unaddressed bug which prevents intermine from working correctly. 

You will also need the following modules installed in your Python environment:

- numpy
- pandas
- pweave
- requests
- xmltodict
- PyQt5
- biomart
- beautifulsoup4
- intermine
- opentargets
- pypdb
- sklearn

The Qt5 library must also be installed and configured to work with PyQt5. 
