# Analyzer
To compile: make
To run: ./analyzeHGCMuons runlist.txt Plots.root data 2000

# Parameters
runlist.txt: Contains a list of all files to run the analyzer code.
Sample content in runlist.txt:
'''
file1.root
file2.root
'''

Pots.root: Name of output file

data: Placeholder

2000: Parameter used for naming files

# Code to modify
HGCNtupleVariables.h : Update branches according to the input file
AnalyzeHGCMuons.cc : Contains main analysis code in the main function
AnalyzeHGCMuons.h : Header file for above where histograms are booked and finally saved in the destructor ( ~AnalyzeHGCMuons() )
