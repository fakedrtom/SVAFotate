SVAFotate
=========================

Overview
=========================
Annotate a (lumpy) SV VCF with allele frequencies (AFs) from large 
population SV cohorts (currently CCDG and/or gnomAD) with a simple 
command line tool. This will add to the INFO field new categories 
corresponding to the maximum AF frequency found for SVs from these 
SV datasets that overlap a given SV in your VCF. It will also include 
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
usage: svafotate.py [-h] [-i STRING] [-o STRING] [-f FLOAT] [-ccdg STRING]
                    [-gnomad STRING]

optional arguments:
  -h, --help      show this help message and exit
  -i STRING       path to VCF to annotate
  -o STRING       output VCF name/path
  -f FLOAT        minimum reciprocal overlap required between SVs (default
                  0.1)
  -ccdg STRING    path to CCDG SV bed file
  -gnomad STRING  path to gnomAD SV bed file
```

Since the idea is to annotate a VCF with AF information from CCDG and/or 
gnomAD, you must include a path to at least one of these. If only one is 
included, only those annotations will be added.

```
python svafotate.py -i your.vcf -o your.annotated.vcf -ccdg ccdg_sv_afs.bed.gz -gnomad gnomad_sv_afs.bed.gz
```

The `-f` option allows for customization in the amount of reciprocal overlap 
between query SVs in the provided VCF and SVs in the SV datasets that they are 
being compared to. Any value between 0 and 1 may be entered (default value of 
0.1 if `-f` is not used). SVs with overlaps will have AFs from the datasets 
added to the INFO fields. Query SVs without overlaps will be given AF annotations 
of 0. The higher the value of `-f`, the more precise the SV overlap match must be.
 
SV Datasets
==========================
Currently this repo includes SV datasets from CCDG and gnomAD that contain AFs 
from large population cohorts. Each of these datasets have been summarized and made
available download here.

## CCDG

Detailed information about the CCDG dataset can be found [here](https://www.biorxiv.org/content/10.1101/508515v1).
While the authors have made a hg38 and hg19 VCF of their dataset available,
the hg38 reflects more samples and hence a larger SV dataset. This was selected,
converted to a BED format and was converted to hg19 using UCSC's [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) 
tool. With the results an hg19 VCF was generated and then converted to a BED including
CHROM, START, END, SVTYPE, and AF. That bed is included as part of this repo.

## gnomAD

Detailed information about the gnomAD dataset can be found [here](https://www.biorxiv.org/content/10.1101/578674v1).
A BED format of their results is available [here](https://gnomad.broadinstitute.org/downloads).
This was downloaded and then summarized with the following command:

```
zcat gnomad_v2_sv.sites.bed.gz | awk '$7=="PASS"' | cut -f 1-3,5,31,39 | bgzip -c > gnomad_sv_afs.bed.gz
```

The resulting BED contains CHROM, START, END, SVTYPE, AF, and POPMAX_AF and is 
included in this repo. Please note that only SVs with a "PASS" FILTER were 
included here. 