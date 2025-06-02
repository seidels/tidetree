# README

This repository contains scripts to generate NEXUS input data for use with TiDeTree.

## Input data
The starting point is the example_data.csv file, where each row represents a cell and contains its associated barcode information. In this example, we handle editing data across 5 independent sites per cell.

## Converting to NEXUS Format
To convert the example_data.csv into NEXUS format, run the conversion script as shown below. The resulting file, here example_output.tidetree, can be loaded into BEAUti for setting up the analysis.

` Rscript write_nexus.R example_data.csv example_output.tidetree`


