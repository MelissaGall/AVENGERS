# AVENGERS
MAVE approach for classification of BRCA2 SNV

**Processing**

_Scripts to process samples_

- BRCA2_pre_processing_pipeline.sh: bash script that will run the following scripts. These scripts were originally written by Shendure's lab team (see credits below). 
  * run_seqprep.py: run SeqPrep
  * run_remove_n_bases.py: remove N bases
  * run_needle_to_sam.py: align .fq reads using NeedleAll
- get_counts.py: Python script to extract HDR read counts
- add_annotation.py: Python script to add ClinVar annotation
- run_predictions.R: R script to filter and normalize the data, then run the model
- functions.R: functions for run_predictions.R

**Link to the original scripts for preprocessing:**

[Shendure lab GitHub page](https://github.com/shendurelab/saturationGenomeEditing_pipeline)

_Reference of the corresponding paper:_

Findlay GM, Daza RM, Martin B, Zhang MD, Leith AP, Gasperini M, Janizek JD, Huang X, Starita LM, Shendure J. 
Accurate classification of BRCA1 variants with saturation genome editing. 
Nature. 2018 Oct;562(7726):217-222. 
doi: 10.1038/s41586-018-0461-z. Epub 2018 Sep 12. PMID: 30209399; PMCID: PMC6181777.
