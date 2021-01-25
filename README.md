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

All C++ code for the transmission model can be found in code/src/. Scripts to run the transmission model for various clinical trial settings can be found in code/condor/. Scripts to process the corresponding outputs can be found in code/scripts/. All R scripts to process results and generate figures are located in code/R/.

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

## Analysis 

### 1. Generate the input files
* The script 'code/R/gen_analysis.R' can be used to generate the various input files needed to run the transmission model for the different clinical trial scenarios,

### 2. Run the transmission model for the different clinical trial scenarios 

* The scripts to run the transmission model using HTCondor (a framework for high throughput computing) are available at 'code/condor/condor.\*'.

### 3. Process the outputs of the clinicial trial scenarios

* To perform the efficacy calculations for each of the clinical trial scenarios, run the corresponding scripts: 'code/scripts/script_process_\*'.
* To calculate efficacy using a genotyping method, run the scripts: 'code/script/script_genotyping\*'.
* To compute the identify the breakdown of recurrent infections, run the scripts: 'code/scripts/script_recurrent_\*'.

### 4. Generate the figures

* Scripts to generate all figures can be found in 'code/R/'. Scripts named 'run_Fig\*.R' process the data necessary to generate the accompanying figure using 'fig_Fig\*.R' and must be run first. 

## Contact Me
If you encounter difficulties running the code or need help adapting it to meet your own needs, please do not hesitate to contact me. I can be reached at jhuber3 AT nd DOT edu

