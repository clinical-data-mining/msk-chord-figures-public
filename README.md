# msk-chord-figures-public
Code and sample data for the MSK-CHORD manuscript (public version)
NOTE: do NOT upload PHI

## Requirements
Python 3.10
Library dependencies: pandas (1.5.3), numpy (1.24.3), matplotlib (3.7.1), seaborn (0.12.2), statsmodels (0.14.0), lifelines (0.27.7), sksurv (scikit-survival 0.20.0). () = tested on this version

Installation: see specific libraries for details. Expected install time ~15min

## Components:
data: Sample flat files containing characteristics for the MSK-CHORD cohort including follow-up time and mortality status (note: full data will be available on cBioPortal on study release)

code: 

--msk_chord_figures.ipynb: jupyter notebook using the files in data to reproduce the figures in the manuscript. Expected run time ~5min

--run_rsf.py: standalone python file to train and 5-fold cross validate several random survival forests (RSFs) grouped by class of variable on the files in data. Expected run time: up to ~2 hours if all variable groupings are included. Using the example (only demographics and genomics for two cancer types) ~15 min

--NLP folder: contains scripts with regex/rule based methods as well as the clinical-longformer used to predict survival ("radLongformer"). NOTE these are meant to be run on PHI so sample data is not provided. For additional examples using BERT and its variant architectures to train and validate NLP classifiers see https://huggingface.co/
