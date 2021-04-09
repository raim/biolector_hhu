# Automating BioLector Data Analysis

The goal is to write an R script that allows to analyze
growth experiments pipetted with the Opentron pipetting robot
and run in the BioLector or Clariostar platereaders 
with **minimal user intervention**.

Required data:

* Data exported from BioLector software (either raw or reduced),
* Microplate layout file,
* Media recipe,
* Physico-chemical constants.

# Preparation

* Your layoutfile should have a structure like: **strain;Glc:amount;Ace:amount;aTc:amount**
* (leave out the things you don`t need)
* for blanks: change strain into "blank"

# How to run

Before you start, you have to adapt the path and experiment IDs at the top of the code:
* adapt "setwd" -> specify the folder, where your data and layout file exist
* adapt "cfile" -> specify the path including your calibration file
* adapt "Parameters" -> specify the path including your Parameters file
* adapt "expID" -> specify the ID of the experiment (name of the data file)
* adapt "calID" -> specify the ID of the calibration from your calibrations.csv file you want to use

	Now run the preparation and the analysis block. A report with all plots will be generated and saved to the specified folder. Done.
	
# Problems

*For some experiments, a pH correction was necessary and for others it wasn`t. This have to be corrected.
*Why there is aTc in some layoutfiles if it is 0 in every well?
