SVAFotate
=========================

SVAFotate is currently undergoing some major updates that are under 
construction. Please check back soon.

Overview
=========================
Annotate a structual variant (SV) VCF with allele frequencies 
(AFs) from large population SV cohorts with a simple command line tool. 
This will add to the INFO field new categories corresponding to the 
maximum AF found for SVs from these SV datasets that overlap 
a given SV in your VCF. SVAFotate enables many additional annotation 
options related to AF metrics.


Installation
========================
0) Installing Miniconda

- If Miniconda is not installed on your system, install it from [miniconda](https://conda.io/en/latest/miniconda.html)


1) Set up new conda environment 

```
$ conda create --name svafotate-env python=3
```

```
$ conda activate svafotate-env
```


2) Install package requirements 

```
$ conda install --file https://raw.githubusercontent.com/fakedrtom/SVAFotate/master/requirements.txt
```


3) Install SVAFotate

```
$ pip install git+https://github.com/fakedrtom/SVAFotate.git
```


4) Check that SVAFotate installed Correctly 

```
$ svafotate --version

svafotate 0.0.1
```

```
$ svafotate -h 


Usage
======================== 
## Options

usage: svafotate [-h] [-v] {annotate,pickle-source,custom-annotation} ...

SVAFotate: Structural Variant VCF Annotation tools
==================================================

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Installed version

[sub-commands]:
  {annotate,pickle-source,custom-annotation}
    annotate            Annotate SV VCF File
    pickle-source       Pickle Source Bed
    custom-annotation   Add custom annotation(s) to source annotation file
```
