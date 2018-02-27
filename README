Requirements
============
Python 3 is required to run the code, as well as the packages pyfasta, numpy, intervaltree and sklearn. All of these can be installed via pip.

Files
=====
The data files `hg19.fa`, `all_exons.bed`, `hg19_Alu.bed` and `hsa_hg19_Rybak2015.bed` all need to be present in the `data` folder.The `tmp` folder needs to be writeable.

Usage
========
The pipeline consists of three steps:

1. Generating negative samples: 
  `python generate_negatives.py`
  The results will be stored in `tmp/negatives.bed`.
2. Computing features:
  `python calculate_features.py`
  If necessary, the source of the negative samples can be changed in the code (e.g. to use the exons as negatives). The results will be stored in `tmp/features.npy` and `tmp/labels.npy`.
3. Training and cross-validation as well as feature ranking:
  `python validate.py`