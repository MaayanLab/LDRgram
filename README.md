# LDRgram

This clustergram visualizes data from the LINCS Dataset Registry (LDR), which is a DCIC web application built to allow the uploading and approval of experimental data provided by the LINCS Data Signature Generation Centers. The clustergram shows the number of perturbations applied to cell lines for each assay. Blue tiles represent released data and orange tiles represent unreleased data and some tiles are split to represent partially released data. 

The LDR data is processed in the python script process_LDR.py and produces the json static/networks/LDR_as_cl.json. 