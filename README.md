# AVENGERS
MAVE approach for classification of BRCA2 SNV

**Processing**

_Scripts to process samples_

- BRCA2_pre_processing_pipeline.sh: bash script that will run the following scripts:
  * run_seqprep.py: run SeqPrep
  * run_remove_n_bases.py: remove N bases
  * run_needle_to_sam.py: align .fq reads using NeedleAll
- get_counts.py: Python script to extract HDR read counts
- add_annotation.py: Python script to add ClinVar annotation
- run_predictions.R: R script to filter and normalize the data, then run the model
- functions.R: functions for run_predictions.R
