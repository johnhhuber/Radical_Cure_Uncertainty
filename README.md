# Site-specific biases in phase-III clinical trials underestimate the effect of radical cure against *Plasmodium vivax* hypnozoites
## John H. Huber, Cristian Koepfli, Guido España, Narimane Nekkab, Michael T. White, T. Alex Perkins

A brief description of the code and data used to run the analyses of Huber *et al*. (2021) are provided below.

## Getting Started

The scripts for the R and C++ code are written with relative paths, assuming you have the following folder structure:
```
Radical_Cure_Uncertainty
|   README.md
└─── code 
     └─── R
     └─── condor
     └─── scripts
     └─── src	 	
└─── data
└─── output
```

All C++ code for the transmission model can be found in code/src/. Scripts to run the transmission model under various clinical trial settings can be found in code/condor/. Scripts to process the corresponding outputs can be found in code/scripts/. All R scripts to process results and generate figures are located in code/R/.

### Software and Package Requirements

The *Plasmodium vivax* transmission model requires:

* GCC : GNU Compiler Collection (We used gcc/8.3.0)

A makefile is provided in code/src/ to build the software project. You may need to adjust the 'g++' command to match your system requirements.

The R scripts require the following packages:

* doParallel
* foreach
* ggsci
* icenreg
* RColorBrewer
* seqinr
* shape
* survival

The R scripts should automatically install the necessary packages if they are not already installed on your local machine. 

