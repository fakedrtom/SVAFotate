SVAFotate
=========================

Overview
=========================
Annotate a (lumpy) structual variant (SV) VCF with allele frequencies 
(AFs) from large population SV cohorts (currently CCDG and/or gnomAD) 
with a simple command line tool. This will add to the INFO field new 
categories corresponding to the maximum AF frequency found for SVs from 
these SV datasets that overlap a given SV in your VCF. It will also include 
a field for the number (count) of overlaps between a given SV in your 
VCF and those found in the CCDG and/or gnomAD SV datasets.

Installation
========================
This is just a python script so you should only need to clone this repo
or download the svafotate.py script itself.

## Dependencies

* [cyvcf2](https://github.com/brentp/cyvcf2)
* [pybedtools](http://daler.github.io/pybedtools/#)

Usage
======================== 
## Options

```
usage: svafotate.py [-h] [-i INPUT VCF] [-o OUTPUT VCF] [-f MINIMUM OVERLAP]
                    [-ccdg PATH TO CCDG] [-gnomad PATH TO GNOMAD]
                    [-ci USE CI BOUNDARIES]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT VCF, --in INPUT VCF
                        path to VCF to annotate
  -o OUTPUT VCF, --out OUTPUT VCF
                        output VCF name/path
  -f MINIMUM OVERLAP, --minf MINIMUM OVERLAP
                        minimum reciprocal overlap required between SVs
                        (default 0.01, must be between 0 and 1.0)
  -ccdg PATH TO CCDG, --ccdg PATH TO CCDG
                        path to CCDG SV bed file
  -gnomad PATH TO GNOMAD, --gnomad PATH TO GNOMAD
                        path to gnomAD SV bed file
  -ci USE CI BOUNDARIES, --ci USE CI BOUNDARIES
                        option to use out inner or outer confidence intervals
                        (CIPOS95, CIEND95) for SV boundaries, must answer "in"
                        or "out"
```

Since the idea is to annotate a VCF with AF information from CCDG and/or 
gnomAD, you must include a path to at least one of these. If only one is 
included, only those AFs will be added.

```
python svafotate.py -i your.vcf -o your.annotated.vcf -ccdg ccdg_sv_afs.bed.gz -gnomad gnomad_sv_afs.bed.gz
```

The `-f` option allows for customization in the amount of reciprocal overlap 
between query SVs in the provided VCF and SVs in the SV datasets that are
designated. Any value between 0 and 1 may be entered (default value of 
0.01 if `-f` is not used). SVs with overlaps will have maximum AFs from the datasets 
added to the INFO fields. Query SVs without overlaps will be given AF annotations 
of 0. The higher the value of `-f`, the more precise the SV overlap match must be. 
For more information on how the overlaps are measured and determined please 
refer to [this](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html).

Lumpy SVs provide confidence intervals for the location of the SV breakpoints. The 
optional `-ci` parameter allows these breakpoints to be adjusted using the 95% confidence 
intervals (CIPOS95 and CIEND95). Indicating `in` with the `-ci` option will reduce the 
size of the SV by increasing the start with the right value from CIPOS95 and 
decreasing the end with the left value of CIEND95. The `out` answer will increase 
the size of the SV by subtracting the left value from CIPOS95 from the start and 
adding the right value from CIEND95 to the end. These new boundaries will then be 
used to ascertain overlaps with the indicated SV datasets, but the resulting output 
VCF will maintain the original boundaries in the POS column and END INFO field.

SV Datasets
==========================
Currently this repo includes SV datasets from CCDG and gnomAD that contain AFs 
from large population cohorts. Each of these datasets have been summarized and made
available for download.

## CCDG

Detailed information about the CCDG dataset can be found [here](https://www.biorxiv.org/content/10.1101/508515v1).
While the authors have made a hg38 and hg19 VCF of their dataset available,
the hg38 reflects more samples and hence a larger SV dataset. This was selected,
converted to a BED format and was then converted to hg19 using UCSC's [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) 
tool. With those results an hg19 VCF was generated and then converted to a BED including
CHROM, START, END, SVTYPE, and AF. That bed is included as part of this repo 
and is named `ccdg_sv_afs.bed.gz`.

## gnomAD

Detailed information about the gnomAD dataset can be found [here](https://www.biorxiv.org/content/10.1101/578674v1).
A BED format of their results is available [here](https://gnomad.broadinstitute.org/downloads).
This was downloaded and then summarized with the following command:

```
zcat gnomad_v2_sv.sites.bed.gz | awk '$7=="PASS"' | cut -f 1-3,5,31,39 | bgzip -c > gnomad_sv_afs.bed.gz
```

The resulting BED contains CHROM, START, END, SVTYPE, AF, and POPMAX_AF and is 
included in this repo, named `gnomad_sv_afs.bed.gz`. Please note that only SVs with a 
"PASS" FILTER were included here. 
