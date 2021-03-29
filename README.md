SVAFotate
=========================

SVAFotate is currently undergoing some major updates that are under 
construction. Please check back soon.

Overview
=========================
Being able to distinguish whether a mutation or variant is common or
rare in the general population is crucial in many genomic analyses, including
structural variant (SV) analysis. More and more large population SV cohorts,
with population level allele frequencies (AFs), are becoming available, but
each consist of primarily rare variants that are unique to the individual
dataset. SVAFotate provides the means to aggregate and combine AF and related
data from multiple population datasets into simple annotations that are
then added to a SV VCF. These annotations are then readily available and
can be used for different filtering schemes or other analyses. SVAFotate
is a command-line tool and provides a variety of annotation options related
to AF metrics.


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

Usage
========================
There are currently three SVAFotate subcommands: **annotate**,
**pickle-source**, and **custom-annotation**.

```
$ svafotate -h
usage: svafotate [-h] [-v] {annotate,pickle-source,custom-annotation} ...

SVAFotate: Structural Variant Allele Frequency annotator ==================================================

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Installed version

[sub-commands]:
  {annotate,pickle-source,custom-annotation}
    annotate            Annotate SV VCF File
    pickle-source       Pickle Source Bed
    custom-annotation   Add custom annotation(s) to source annotation file
```

**annotate**

Adds AF related metrics to your SV VCF based on overlapping matches
between SVs in an input VCF and SVs in a provided BED file. This is
the main functionality of SVAFotate. 

```
svafotate annotate -h
usage: svafotate annotate [-h] -v INPUT VCF -o OUTPUT VCF [-b SOURCE BED] [-p Pickled Source Data] [-f [MINIMUM OVERLAP FRACTION [MINIMUM OVERLAP FRACTION ...]]]
                          [-s [SOURCES TO ANNOTATE [SOURCES TO ANNOTATE ...]]] [-a [EXTRA ANNOTATIONS [EXTRA ANNOTATIONS ...]]] [-c OBSERVED SV COVERAGE] [-u UNIQUE SV REGIONS] [-l SV SIZE LIMIT]
                          [-t TARGETS BED FILE] [-ci USE CI BOUNDARIES] [-ci95 USE CI BOUNDARIES] [-e EMBIGGEN THE SV SIZE] [-r REDUCE THE SV SIZE] [--cpu CPU Count]
			  
Annotate your SV VCF files with Population Allele Frequency Info

help:
  -h, --help            show this help message and exit

Required Arguments:
  -v INPUT VCF, --vcf INPUT VCF
                        Path and/or name of the VCF file to annotate.
  -o OUTPUT VCF, --out OUTPUT VCF
                        Path and/or name of the output VCF file.
			sources [SOURCES TO ANNOTATE [SOURCES TO ANNOTATE ...]]

Requires Only One Argument:
  -b SOURCE BED, --bed SOURCE BED
                        Path and/or name of the combined sources AF bed file (--pickled-source can be used instead of --bed) (NOTE: --bed takes higher priority than --pickled-source).
  -p Pickled Source Data, --pickled-source Pickled Source Data
                        Path and/or name of the pickled source data to use (--bed can be used instead of --pickled-source) (NOTE: if --bed is provided, --pickled-source is ignored).

Optional Arguments:
  -f [MINIMUM OVERLAP FRACTION [MINIMUM OVERLAP FRACTION ...]], --minf [MINIMUM OVERLAP FRACTION [MINIMUM OVERLAP FRACTION ...]]
                        A space seperated list of minimum reciprocal overlap fractions required between SVs for each source listed with the `-s` option. If `-s` is not used, only the first minf will be
                        used and will be applied to all sources. minf values must be between 0.0 and 1.0 (Default = 0.001).
  -s [SOURCES TO ANNOTATE [SOURCES TO ANNOTATE ...]], --sources [SOURCES TO ANNOTATE [SOURCES TO ANNOTATE ...]]
                        Space seperated list of data sources to use for annotation. If '-s' is not used, all sources available in the source bed file will be used (Example: ' -s CCDG gnomAD ' ).
  -a [EXTRA ANNOTATIONS [EXTRA ANNOTATIONS ...]], --ann [EXTRA ANNOTATIONS [EXTRA ANNOTATIONS ...]]
                        By default, only the Max_AF, Max_Hets and Max_HomAlt counts, and Max_PopMax_AF are annotated in the output VCF file. `-a` can be used to add additional annotations, with each
                        anntotation seperated by a space (Example ' -a mf best pops ' ). Choices = [all, mf, best, pops, AFR, AMR, EAS, EUR, OTH, SAS, full, mis]
  -c OBSERVED SV COVERAGE, --cov OBSERVED SV COVERAGE
                        Add an annotation reflecting how much of the queried SV genomic space has been previously observed with the same SVTYPE. Uses the data sources listed with -s as the previously
                        observed SVs. Please provide minimum AF to exclude all SVs from data sources with a total AF below that value (must be between 0 and 1.0).
  -u UNIQUE SV REGIONS, --uniq UNIQUE SV REGIONS
                        Generate a file of unique SV regions called 'unique.bed'. These regions reflect genomic space within the queried SV region that have not been previously observed with the same
                        SVTYPE. This will also add an annotation regarding the number of unique regions within a given SV. Please provide minimum AF to exclude all SVs from data sources with a total AF
                        below that value (must be between 0 and 1.0).
  -l SV SIZE LIMIT, --lim SV SIZE LIMIT
                        Only include previously observed SVs from data sources with a size less than or equal to this value (only available when using --cov or --uniq).
  -t TARGETS BED FILE, --target TARGETS BED FILE
                        Path to target regions BED file. Expected format is a tab delimited file listing CHROM START END ID where ID is a genomic region identifier that will be listed as an annotation
                        if an overlap exists between a given SV and the target regions.
  -ci USE CI BOUNDARIES, --ci USE CI BOUNDARIES
                        Expects CIPOS and CIEND to be included in the INFO field of the input VCF (--vcf). If argument is selected, use 'inner' or 'outer' confidence intervals (CIPOS, CIEND) for SV
                        boundaries. Choices = [in, out]
  -ci95 USE CI BOUNDARIES, --ci95 USE CI BOUNDARIES
                        Expects CIPOS95 and CIEND95 to be included in the INFO field of the input VCF (--vcf). If argument is selected, use 'inner' or 'outer' confidence intervals (CIPOS95, CIEND95) for
                        SV boundaries. Choices = [in, out]
  -e EMBIGGEN THE SV SIZE, --emb EMBIGGEN THE SV SIZE
                        Increase the size of the SV coordinates in the input VCF (--vcf) by a single integer; Subtract that value from the start and add it to the end of each set of coordinates.
  -r REDUCE THE SV SIZE, --red REDUCE THE SV SIZE
                        Reduce the size of the SV coordinates in the input VCF (--vcf) by a single integer; Add that value to the start and subtract it from the end of each set of coordinates.
  --cpu CPU Count       The number of cpus to use for multi-threading (Default = 1).
```

As stated, this requires an input VCF and output VCF name. SVAFotate was
developed using SV VCFs derived from the Lumpy SV caller, however it has
been tested on VCFs created by other SV callers. As long as SVTYPE and END
(preferably SVLEN as well) are included in the INFO fields, any SV VCF that
follows expected VCF conventions should be usable with SVAFotate.

SVAFotate also requires a BED file corresponding to population SV data. A
BED file is provided along with SVAFotate, called SVAFotate_core_SV_popAFs.GRCh38.bed.gz,
and was compiled using information gathered from the publicly available
CCDG, gnomAD, and 1000G SV callsets. AF related data (including HOM_REF, HET,
and HOM_ALT counts, where available), were parsed from these into a single BED
file. A different BED file could be used or customized data could be added
to the provided BED file, but please note that SVAFotate expects specific
columns and their order placement to be present in the BED file. If a
different BED file is used or additional data is added to the provided BED,
please ensure that it follows the same ordering and column information. Please
note that all columns do need to be populated with actual data and where data is
unavailable an 'NA' should suffice. As a minimum, it is recommended that any other
BED file used or other data added to this provided BED file include: CHROM, START,
END, SVLEN, SVTYPE, SOURCE, SV_ID, and AF.

**pickle-source**

**custom-annotation***

This is not currently used, but is a placeholder for potential
future developments.

Acknowledgements 
========================
Thank you to Michael Cormier for code review, improvements, and contributions.