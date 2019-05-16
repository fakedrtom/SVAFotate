SVAFotate
=========================

Overview
=========================
Annotate a (lumpy) SV VCF with various AFs from CCDG and/or gnomAD with a 
simple command line tool. This will add to the INFO field new categories
corresponding to the maximum AF frequency found for SVs from CCDG and/or gnomAD
that overlap a given SV in your VCF. Also includes a field for the number (count) 
of overlaps between a given SV in your VCF and those found in the CCDG and/or 
gnomAD SV datasets.

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

Since the idea is to annotate a VCF with AF information from CCDG and/or gnomAD, 
you must include a path to at least one of these. If only one is included,
only those annotations will be added.

```
python svafotate.py -i your.vcf -o your.annotated.vcf -ccdg ccdg_sv_afs.bed.gz -gnomad gnomad_sv_afs.bed.gz
```

