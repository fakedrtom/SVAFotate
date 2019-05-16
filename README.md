SVAFotate
=========================

Overview
=========================
Annotate a lumpy SV VCF with various AFs from gnomAD and CCDG with a 
simple command line tool.

## Dependencies

* [cyvcf2](https://github.com/brentp/cyvcf2)
* [pybedtools](http://daler.github.io/pybedtools/#)
 
## Options

```
usage: svafotate.py [-h] [-i STRING] [-f FLOAT] [-o STRING] [-ccdg STRING]
                    [-gnomad STRING]

optional arguments:
  -h, --help      show this help message and exit
  -i STRING       path to VCF to annotate
  -f FLOAT        minimum overlap required between SVs (default 0.1)
  -o STRING       output VCF name/path
  -ccdg STRING    path to CCDG SV bed file
  -gnomad STRING  path to gnomAD SV bed file
```
