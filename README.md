# count-dicty-cells

R script for processing flow cytometry data of Dictyostelium discoideum cells

## Usage

1. Export your flow cytometry data as .fcs (flow cytometry standard) files. 
2. Make sure you have the Bioconductor flow cytometry packages.  If you don't have them, run these commands: 
```
source("http://bioconductor.org/biocLite.R")
biocLite("flowCore")
biocLite("flowViz")
biocLite("flowStats")
```
3. Choose fluorescence channel. Default choices are red and green. 
4. Name results file (optional).  This will be a tab-delimited text file that can be directly read by R or copy/pasted into Excel.  
5. Execute the script. Type cmd-E in Mac OS X or choose menu item "Edit/Source Document". Make sure your working directory is where your fcs files are.  The script will simultaneously analyze all .fcs files in the working directory. You may want to separately analyze files in batches based on experiment, replicate, strain, et cetera. 
6. Look at the diagnostic plots. They'll show you what the script believes is Dicty, how it's choosing between fluorescent and nonfluorescent cells, whether there is debris inflating nonfluorescent counts, and compare the fluorescence of different samples.  

