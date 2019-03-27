## count-dicty-cells.R v1.1
## Script for analyzing fluorescent Dictyostelium cells via flow cytometry
## For use with R statistical computing environment (r-project.org)
## Tested with R version 3.2.3
## 
## jeff smith 2016
## Dept of Biology, Washington Univ in St Louis

## 
## HOW TO USE THIS SCRIPT: 
## 
## 1. Export your flow cytometry data as .fcs (flow cytometry standard) files. 
##    In the accuri CFlow Plus software, you can choose menu item "File/Export 
##    ALL files as FCS..."
## 
## 2. Make sure you have the Bioconductor flow cytometry packages.  If you don't 
##    have them, run these commands: 
## 
	# source("http://bioconductor.org/biocLite.R")
	# biocLite("flowCore")
	# biocLite("flowViz")
	# biocLite("flowStats")
## 
## 3. Choose fluorescence channel
## 
      # myFluorescenceChannel <- "FL1-H"   # Green
      myFluorescenceChannel <- "FL3-H"   # Red
## 
## 4. Name results file (optional).  This will be a tab-delimited text file that 
##    can be directly read by R or copy/pasted into Excel.  
## 
      myResults <- "flow-cytometry-results.txt"
## 
## 5. Execute this script. Type cmd-E in Mac OS X or choose menu item "Edit/Source 
##    Document".  Make sure your working directory is where your fcs files are.  The 
##    script will simultaneously analyze all .fcs files in the working directory.  
##    You may want to separately analyze files in batches based on experiment, 
##    replicate, strain, et cetera. 
## 
## 6. Look at the diagnostic plots.  They'll show you what the script believes is 
##    Dicty, how it's choosing between fluorescent and nonfluorescent cells, whether 
##    there is debris inflating nonfluorescent counts, and compare the fluorescence 
##    of different samples.  
## 



## ================================================
## YOU SHOULDN'T NEED TO CHANGE ANYTHING BELOW HERE
## ================================================

# Load packages
cat("\nLoading bioconductor.org flow cytometry packages...\n")
library(flowCore)      # Tested with version 1.36.5
library(flowViz)       # Tested with version 1.34.0
library(flowStats)     # Tested with version 3.28.1

# Read in data
cat("Reading data from .fcs files in directory ", getwd(), "...\n", sep = "")
myData <- read.flowSet(
	pattern = ".fcs",
	phenoData = list(
		sampleDescription = "#SAMPLE", 
		dataCell          = "$SMNO", 
		dateCollected     = "$DATE", 
		volumeNanoliters  = "$VOL"
	)
) 
# pData(phenoData(myData))  # summary of metadata

# The Bioconductor flow cytometry packages use "workFlow" objects 
# to manage transformations, gates, and naming of intermediate results 
myWorkFlow <- workFlow(myData, name = "myWorkFlow") 

# Transform data (asinh)
cat("Transforming data (asinh)...\n")
myTransformation <- transformList(
	colnames(myData)[1:12], 
	tfun = asinh, 
	transformationId = "myTransformation"
)
add(myWorkFlow, myTransformation)

# Screen out some of the debris
cat("Screening out debris...\n")
myDebrisScreen <- rectangleGate("FSC-H" = c(11, 18), "SSC-H" = c(10, 20), filterId = "myDebrisScreen") 
add(myWorkFlow, myDebrisScreen, parent = "myTransformation")

# Fit bivariate normal gate to the data in the region where we expect Dicty spores: 
cat("Isolating Dicty cells...\n")
myCellGate <- lymphGate(
	Data(myWorkFlow[["myDebrisScreen+"]]), 
	channels     = c("FSC-H","SSC-H"), 
	preselection = list("FSC-H" = c(13.75, 16.5), "SSC-H" = c(10, 15)), 
	eval         = FALSE, 
	filterId     = "DictyCells"
)
add(myWorkFlow, myCellGate$n2gate, parent="myDebrisScreen+")

# Screen out some of the debris in the Dicty region
cat("Screening out more debris...\n")
myBoundaryFilter <- boundaryFilter(filterId = "myBoundaryFilter", x = myFluorescenceChannel)
add(myWorkFlow, myBoundaryFilter, parent = "DictyCells+")

# Normalize fluorescence peaks if there is more than one fcs file
if ( length(myData) > 1 ) {
	cat("Normalizing fluorescence peaks...\n")
	myNormalization <- normalization(
		normFun         = function(x, parameters, ...) warpSet(x, parameters,...), 
		parameters      = c(myFluorescenceChannel), 
		normalizationId = "myNormalization"
	) 
	add(myWorkFlow, myNormalization, parent = "myBoundaryFilter+") 
}

# Automagically create dividing line between fluorescent and nonfluorescent cells
cat("Measuring number of fluorescent and nonfluorescent cells...\n")
if ( length(myData) > 1 ) {
	myFluorescenceGate <- rangeGate(
		Data(myWorkFlow[["myNormalization"]]), 
		myFluorescenceChannel, 
		plot     = FALSE, 
		filterId = "myFluorescenceGate", 
		sd       = 2.5, 
		refLine  = 4
	) 
	add(myWorkFlow, myFluorescenceGate, parent = "myNormalization")
} else {
	myFluorescenceGate <- rangeGate(
		Data(myWorkFlow[["myBoundaryFilter+"]]), 
		myFluorescenceChannel, 
		plot     = FALSE, 
		filterId = "myFluorescenceGate", 
		sd       = 2.5, 
		refLine  = 4
	) 
	add(myWorkFlow, myFluorescenceGate, parent = "myBoundaryFilter+")
}

# Export results
cat("Exporting results as tab-delimited text file \"", myResults, "\"...\n", sep = "")
results <- cbind(
	pData(phenoData(myData)), 
	"fluorescentCells"    = summary(myWorkFlow[["myFluorescenceGate+"]])$true, 
	"nonfluorescentCells" = summary(myWorkFlow[["myFluorescenceGate+"]])$false, 
	"fractionFluorescent" = summary(myWorkFlow[["myFluorescenceGate+"]])$p
)
results$dateCollected <- as.Date(results$dateCollected, "%d-%b-%Y")
names(results)[names(results) == "name"] <- "fileName"
myResultsFile <- file(myResults)
writeLines(paste("# Flow cytometry results analyzed", Sys.Date()), myResultsFile)
close(myResultsFile)
suppressWarnings(
	write.table(
		results, file = myResults, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE
	)
)

cat("Making diagnostic plots...\n\n")

# Plot Dicty cell/spore gates
if ( length(myData) > 1 ) {
	dev.new(width = 6, height = 6)
} else {
	dev.new(width = 4, height = 4)
}
print(xyplot(
	`SSC-H` ~ `FSC-H`, 
	myWorkFlow[["DictyCells+"]], 
	par.settings = list(gate = list(col = "black", fill = "red", alpha = 0.25)), 
	smooth       = FALSE,
	xbin         = 128,
	xlim         = c(11, 18), 
	ylim         = c(8, 18), 
	xlab         = "Forward scatter (FSC-H)", 
	ylab         = "Side scatter (SSC-H)", 
	main         = "Gate for Dicty cells/spores"
))

# Plot relative fluorescence of samples
if ( length(myData) > 1 ) {
	dev.new(width = 4, height = 6)
	print(densityplot(
		~ ., 
		Data(myWorkFlow[["myBoundaryFilter+"]]), 
		channel = c(myFluorescenceChannel),
		xlim    = c(4, 13), 
		xlab    = paste("Fluorescence (", myFluorescenceChannel, ")", sep = ""),
		ylab    = "Relative abundance (smoothed density)",
		main    = "Relative fluorescence of samples"
	))
}

# Plot fluorescence histogram and gate
myFormula <- substitute(formula(~ a), list(a = as.symbol(myFluorescenceChannel)))
if ( length(myData) > 1 ) {
	dev.new(width = 6, height = 6)
	print(densityplot(
		formula(myFormula),
		Data(myWorkFlow[["myNormalization"]]), 
		plotType = "histogram",
		breaks   = 256,
		stack    = F,
		refline  = myFluorescenceGate@min, 
		xlab     = paste("Normalized fluorescence (", myFluorescenceChannel, ")", sep = ""),
		main     = "Reasonable fluorescence gating?", 
		xlim     = c(4, 13) 
	))
} else {
	dev.new(width = 4, height = 4)
	print(densityplot(
		formula(myFormula),
		Data(myWorkFlow[["myBoundaryFilter+"]]), 
		plotType = "histogram",
		breaks   = 256,
		stack    = F,
		refline  = myFluorescenceGate@min, 
		xlab     = paste("Fluorescence (", myFluorescenceChannel, ")", sep = ""),
		main     = "Reasonable fluorescence gating?", 
		xlim     = c(4, 13) 
	))
}

# Plot fluorescence against forward scatter
if ( length(myData) > 1 ) {
	dev.new(width = 6, height = 6)
} else {
	dev.new(width = 4, height = 4)
}
myFormula <- substitute(formula(a ~ `FSC-H`), list(a = as.symbol(myFluorescenceChannel)))
print(xyplot(
	formula(myFormula),
	myWorkFlow[["DictyCells+"]], 
	smooth = FALSE,
	xbin   = 128,
	xlim   = c(11, 18), 
	xlab   = "Forward scatter (FSC-H)",
	ylab   = paste("Fluorescence (", myFluorescenceChannel, ")", sep = ""),
	main   = "Debris inflating nonfluorescent counts?"
))

