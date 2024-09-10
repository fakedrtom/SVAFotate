SVAFotate
=========================

Structural Variant Allele Frequency annotate or SVAFotate is a tool for
annotating structural variant VCFs with population level allele frequency
information and other related metrics.

Please see the [SVAFotate publication](https://doi.org/10.1186/s12859-022-05008-y) 
for more information and for citing purposes.

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

## Quick Links
[**Installation**](https://github.com/fakedrtom/SVAFotate#installation)<br>
[**Usage**](https://github.com/fakedrtom/SVAFotate#usage)<br>
  * [**annotate**](https://github.com/fakedrtom/SVAFotate#annotate)<br>
   [Minimum Overlap Fraction](https://github.com/fakedrtom/SVAFotate#minimum-overlap-fraction)<br>
   [Sources to Annotate](https://github.com/fakedrtom/SVAFotate#sources-to-annotate)<br>
   [Extra Annotations](https://github.com/fakedrtom/SVAFotate#extra-annotations)<br>
   [Observed SV Coverage](https://github.com/fakedrtom/SVAFotate#observed-sv-coverage)<br>
   [Unique SV Regions](https://github.com/fakedrtom/SVAFotate#unique-sv-regions)<br>
   [SV Size Limit](https://github.com/fakedrtom/SVAFotate#sv-size-limit)<br>
   [Targets BED File](https://github.com/fakedrtom/SVAFotate#targets-bed-file)<br>
   [Use CI boundaries](https://github.com/fakedrtom/SVAFotate#use-ci-boundaries)<br>
   [Change SV Size](https://github.com/fakedrtom/SVAFotate#change-sv-size)<br>
   [CPU Count](https://github.com/fakedrtom/SVAFotate#cpu-count)<br>
  * [**pickle-source**](https://github.com/fakedrtom/SVAFotate#pickle-source)<br>  
  * [**custom-annotation**](https://github.com/fakedrtom/SVAFotate#custom-annotation)

[**Acknowledgements**](https://github.com/fakedrtom/SVAFotate#acknowledgements)<br>

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

svafotate 0.1.0
```

Usage
========================
There are currently three SVAFotate subcommands: `annotate`,
`pickle-source`, and `custom-annotation`.

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

annotate
------------------------
This is the main functionality of SVAFotate and adds AF related metrics to
your SV VCF based on overlapping matches between SVs in an input VCF and SVs
in a provided BED file.

```
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
```


As stated, this requires an input VCF (`-v`) and output VCF name (`-o`). SVAFotate was
developed using SV VCFs derived from the Lumpy SV caller, however it has
been tested on VCFs created by other SV callers. As long as SVTYPE and END
(preferably SVLEN as well) are included in the INFO fields, any SV VCF that
follows expected VCF conventions should be usable with SVAFotate. Please note that a
unique SV ID (column 3 or the `ID` column of the VCF) is also expected. 

SVAFotate also requires a BED file corresponding to population SV data (`-b`) that you
wish to compare the SVs in the input VCF against. A BED file with currently available
large population SV data is available, called:

```
SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz
```

This file can be found [here](https://zenodo.org/records/11642574) while 
previous versions of similar data can be found in the [supporting_data](https://github.com/fakedrtom/SVAFotate/tree/master/supporting_data) 
directory. The BED file was compiled using information gathered from the publicly available
CCDG, gnomAD, 1000G, and TOPMed SV callsets. AF related data (including HOM_REF, HET,
and HOM_ALT genotype counts, where available), were parsed from these into a single BED
file. All of these data were prepared from GRCh38 alignments. A different BED file could 
be used or customized, user-specific data could be added to the provided BED file, but 
please note that SVAFotate expects specific columns, their header names, and their order 
placement to be present in the BED file. If a different BED file is used or additional 
data is added to the provided BED, please ensure that it follows the same ordering and 
column information as found in `SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz`. As a minimum, 
it is necessary that any other BED file used or other data added to this provided
BED file include: CHROM, START, END, SVLEN, SVTYPE, SOURCE, SV_ID, and AF. 
Please note that all other columns do *not* need to be populated with actual data and where
data is unavailable an 'NA' should be added.

**CURRENTLY UNAVAILABLE**: 
A serialized pickle object of the BED file can be used in place of the BED file (`-p`).
This may result in faster performance when using the SVAFotate `annotate` subcommand.
Included with SVAFotate is the subcommand [pickle-source](https://github.com/fakedrtom/SVAFotate#pickle-source)
which enables a pickle object of a BED file to be generated and then repeatedly used
with the `annotate` subcommand.

The bare minimum to run the `annotate` subcommand is:

```
svafotate annotate -v in.vcf.gz -o out.vcf -b SVAFotate_core_SV_popAFs.GRCh38.bed.gz
```

The `annotate` subcommand considers each SV in the VCF input individually and identifies
all matching SVs from the BED file. **Matching SVs are defined as SVs from the input VCF
that share overlapping genomic coordinates with SVs provided in the BED file and that are
described as having the same SVTYPE**. Under default settings even a single basepair overlap
can qualify as a matching SV which may result in many matches being rather imprecise (and also
an overabundance of matches). There are options to require more stringent overlap requirements
to help avoid this (please see [Minimum Overlap Fraction](https://github.com/fakedrtom/SVAFotate#minimum-overlap-fraction)
). Once SVAFotate collects all matching SVs, it assesses individual metrics from each matching
SV and creates new population AF related annotations for each SV from the input VCF. By
default, SVAFotate will add the following annotations to the INFO field of the input VCF:

```
Max_AF		 The maximum AF from all matching SVs across all specified data sources found in the provided BED file
Max_Het	      	 The maximum Het count from all matching SVs across all specified data sources found in the provided BED file
Max_HomAlt    	 The maximum HomAlt count from all matching SVs across all specified data sources found in the provided BED file
Max_PopMax_AF	 The maximum PopMax_AF from all matching SVs across all specified data sources found in the provided BED file
```

For each of these annotations the maximum values from all matching SVs is returned.
If there are no matching SVs, each of these annotations will be added with a 0 or 0.0 value.
The following figure provides example input SVs compared to SVs of the same SVTYPE derived
from the CCDG, gnomAD, and 1000G datasets, which are included with the input BED file.
Consider the example SV on the far right which has matches with 4 different SVs of the same
SVTYPE from CCDG, gnomAD, and 1000G in the BED file, and these have reported AFs of
0.01, 0.03, 0.63, and 0.09. In this case the `Max_AF` annotation that is added to the
input SV will be 0.63.

![max_af](https://github.com/fakedrtom/SVAFotate/blob/master/images/max_AF_example_fig.png)

Please note that each of these annotations is determined independently of one another.
That is to say, just because the SV with an AF of 0.63 is returned for the `Max_AF`
annotation, the HET or HOM_ALT genotype counts corresponding to the SV with the AF
of 0.63 are not necessarily going to be returned for the `Max_Het` and `Max_HomAlt`
annotations unless they are also the maximum values from all matching SVs (in this
case they likely would be added given the wide disparity between AFs in all matching SVs).
Also by default an annotation reflecting the number of matches per data source in the
BED file is also added. For example, if CCDG is one of the data sources in the input
BED file, then `CCDG_Count` is added and lists the number of CCDG matching SVs. 

#### Matching more complex SVs
As previously mentioned, matching SVs are defined as an SV from the input VCF that shares
an overlap of genomic coordinates with an SV in the BED file, provided that these SVs share
the same SVTYPE. This is fairly straightforward for many SVTYPES, such as, DELs, DUPs,
and INVs. This matching scheme can be more complicted for other SVTYPEs. Insertions
(INSs), for example, are often reported as a single basepair (or 2bp) genomic coordinate
with an accompanying SVLEN that reflects the size of the insertion. INSs are still matched
based on overlapping coordinates, which generally means that INSs from the input VCF only match
when their coordinates are (nearly) the same as those in the BED file (even if SVLENs differ).
Reported SVLENs for INSs are then used if certain options (such as `-f` or `-a best`) are
requested to better refine the matching INSs. Other even more complex SVTYPES may require more
specialized attention. In some of these cases, it may be helpful to include the `-a mis`
parameter which would add annotations regarding overlapping SVs that have different SVTYPEs.
For more information please see the [Extra Annotations](https://github.com/fakedrtom/SVAFotate#extra-annotations)
section.

There are many additional options beyond the defaults of SVAFotate that may result
in improved or more detailed annotations. Some of these are highly recommended for
most uses of SVAFotate. These options are listed here along with explanations:

#### Minimum Overlap Fraction
```
  -f [MINIMUM OVERLAP FRACTION [MINIMUM OVERLAP FRACTION ...]], --minf [MINIMUM OVERLAP FRACTION [MINIMUM OVERLAP FRACTION ...]]
                        A space seperated list of minimum reciprocal overlap fractions required between SVs for each source listed with the `-s` option. If `-s` is not used, only the first minf will be
                        used and will be applied to all sources. minf values must be between 0.0 and 1.0 (Default = 0.001).
```

By default, SVAFotate will consider any overlap of genomic coordinates between SVs from the
input VCF and SVs in the supplied BED file, even those of a single base pair, as matching SVs.
While this represents all possible overlaps, it could result in imprecise matching of SVs
that do not actually reflect similar SV alleles. **It is highly recommended to use the `-f` option
to increase the required amount of reciprocal overlap between putative matching SVs**. The
higher the value of `-f`, the more exact the overlapping match between SVs must be (and
likely the more precise and limited the number of overlapping matches will be). For example,
a `-f` value of 0.9 would require that 90% of the genomic region belonging to an SV from
the input VCF must overlap with at least 90% of the genomic region belonging to an SV in
the BED file in order to be considered a matching SV (provided they have the same SVTYPE).
If this requirement is not met, the SV from the BED file will not be added as a matching SV.
Given the inherent variance of SV breakpoints, it is difficult to presume what the best value
for `-f` is, but a minimum of 0.5 should be a decent starting value to consider.

#### Sources to Annotate
```
  -s [SOURCES TO ANNOTATE [SOURCES TO ANNOTATE ...]], --sources [SOURCES TO ANNOTATE [SOURCES TO ANNOTATE ...]]
                        Space seperated list of data sources to use for annotation. If '-s' is not used, all sources available in the source bed file will be used (Example: ' -s CCDG gnomAD ' ).
```

The required BED file for running the `annotate` subcommand may contain SV data from multiple
datasets. The `SOURCE` column in this BED file refers to the data source for the given SVs. For
example, the supplied `SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz` file contains SVs from the CCDG,
gnomAD, 1000G, and TOPMed dataset sources. By default SVAFotate will consider all SVs in the BED file for
determining potential matches (and can report metrics specifc to each source). If annotations
based on comparisons with select sources only is desired, the `-s` parameter will reduce the
considered SVs from the BED file to only those belonging to the requested sources. A single
source could be selected or any combination of included sources from the BED file.

#### Extra Annotations
```
  -a [EXTRA ANNOTATIONS [EXTRA ANNOTATIONS ...]], --ann [EXTRA ANNOTATIONS [EXTRA ANNOTATIONS ...]]
                        By default, only the Max_AF, Max_Hets and Max_HomAlt counts, and Max_PopMax_AF are annotated in the output VCF file. `-a` can be used to add additional annotations, with each
                        anntotation seperated by a space (Example ' -a mf best pops ' ). Choices = [all, mf, best, pops, AFR, AMR, EAS, EUR, OTH, SAS, full, mis]
```

A number of optional annotations are available and can be included using the `-a` option. This
option requires one or more of the following choices which each add more annotations as described
here:

***all***

Adds all of the annotations described by each of the other choices available with `-a`.

***mf***

Adds male and female annotations including:

```
Max_Male_AF	    The maximum Male AF from all matching SVs across all specified data sources
Max_Male_Het	    The maximum Male Het count from all matching SVs across all specified data sources
Max_Male_HomAlt	    The maximum Male HomAlt count from all matching SVs across all specified data sources
Max_Female_AF	    The maximum Female AF from all matching SVs across all specified data sources
Max_Female_Het	    The maximum Female Het count from all matching SVs across all specified data sources
Max_Female_HomAlt   The maximum Female HomAlt count from all matching SVs across all specified data sources
```

If used alongside `best`, `pops` or any of the individual population choices, those annotations
will also include male and female specific annotations regarding AF and genotype counts.

***best***

Adds the "best" match for each data source included in the BED file (or as specified)
and creates the following annotations:

```
Best_[data_source]_ID		The SV_ID of the best matching SV for that data source
Best_[data_source]_OFP		The OFP of the best matching SV for that data source
Best_[data_source]_AF	    	The AF from the best SV_ID for that data source
Best_[data_source]_Het    	The Het count from the best SV_ID for that data source
Best_[data_source]_HomAlt 	The HomAlt count from the best SV_ID for that data source
Best_[data_source]_PopMax_AF	The PopMax_AF from the best SV_ID for that data source
```

The "best" match here is defined by the matching SVs with the highest overlap fraction
product (OFP) and considers all overlaps even those that may be filtered out using `-f` 
option. Calculating the OFP is straightforward and relies on the genomic size of
overlapping SVs and the amount of overlap that is shared between them.

![calculating_ofp](https://github.com/fakedrtom/SVAFotate/blob/master/images/calculating_ofp.png)

From this example, an overlap fraction is calculated for each SV by dividing the amount
of overlap by the size of each SV, respectively. Then these fractions are multiplied to
create the OFP which will range between 0.0 and 1.0. As illustrated below, high OFP scores
reflect matching SVs that are more identical in terms of both their genomic sizes and the
amount of overlap they share while low OFP scores suggest a disparity in genomic sizes
between matching SVs or a low amount of shared overlap between them (or both a discrepancy
in sizes and low overlap).

![ofp](https://github.com/fakedrtom/SVAFotate/blob/master/images/ofp_range.png)

If used alongside `mf`, `pops` or any of the individual population choices, "best" annotations
will also be included for those annotations.

***pops***

The suggested BED file `SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz` contains data from gnomAD, 1000G, 
and TOPMed which include population specific metrics belonging to the following populations:
AFR, AMI, AMR, ASJ, EAS, EUR, FIN, MID, NFE, OTH, and SAS. The `pops` choice will add the 
following annotations for all populations:

```
Max_[population]_AF	    The maximum population AF from all matching SVs across all specified data sources
Max_[population]_Het	    The maximum population Het count from all matching SVs across all specified data sources
Max_[population]_HomAlt	    The maximum population HomAlt count from all matching SVs across all specified data sources
```

If preferred, individual population designations can be selected instead to add annotations
pertaining to that population only. If used alongside `mf` or `best`, male and female or "best"
annotations will also be added for the designated populations.

***full***

Adds all information available for all matches. This will add the following data source specific
annotations:

```
[data_source]_Matches	Comma-separated list of each SV match with the data source
```

Please note that this annotation will include all information available in the BED file
for all matching SVs. This is the "kitchen sink" annotation for matches and will potentially
add very lengthy annotations to the INFO field. For example, if an SV from the input
VCF is found to have 4 matches with SVs from the gnomAD dataset, all information from
the BED file for all 4 of those matches will be listed.

***mis***

While matching SVs are defined as having the same SVTYPE it is possible to see if any SVs with a
different SVTYPE also overlap with an SV from the input VCF (a mismatched SV) with the `mis`
option. This will add the following annotations:

```
[data_source]_Mismatches		Comma-separated list of the SV_IDs of overlapping SVs from the data source with different SVTYPEs
[data_source]_Mismatches_Count		The number of overlapping SVs from the data source with different SVTYPEs
[data_source]_Mismatch_SVTYPEs		Comma-separated list of the other overlapping SVTYPEs from the data source
```

The same requirements for traditional matching SVs will be applied here. For example, if the `-f`
option is used, mismatched SVs will be required to meet those same thresholds. Additionally, the 
following annotations will also be added reflecting the "best" mismatch based on OFP scores:

```
Best_[data_source]_Mismatch_ID		The SV_ID of the best overlapping SV with different SVTYPE from the data source
Best_[data_source]_Mismatch_OFP		The OFP of the best overlapping SV with different SVTYPE from the data source
Best_[data_source]_Mismatch_SVTYPE	The SVTYPE of the best mismatch SV_ID from the data source
Best_[data_source]_Mismatch_AF		The AF of the best mismatch SV_ID from the data source 
Best_[data_source]_Mismatch_Het		The Het count of the best mismatch SV_ID from the data source 
Best_[data_source]_Mismatch_HomAlt	The HomAlt count of the best mismatch SV_ID from the data source 
```

When using the `-a` option, any combination of the above (`mf`, `best`, `pops`, `full`, and 
`mis`) choices can be selected, but if `all` is included than all of these annotations will 
be added to the VCF.

#### Observed SV Coverage
```
  -c OBSERVED SV COVERAGE, --cov OBSERVED SV COVERAGE
                        Add an annotation reflecting how much of the queried SV genomic space has been previously observed with the same SVTYPE. Uses the data sources listed with -s as the previously
                        observed SVs. Please provide minimum AF to exclude all SVs from data sources with a total AF below that value (must be between 0 and 1.0).
```

To see how much of a given SV in the input VCF has been observed with the same SVTYPE from SVs
in the BED file, the `-c` option is provided. This will add the `SV_Cov` annotation
which reflects the fraction of the genomic coordinates for an SV in the input VCF that
are found to overlap with SVs of the same SVTYPE in the BED file. The following figure provides example
scenarios with all SVs being the same SVTYPE and the resulting `SV_Cov` annotations that would
be added to the given SVs. 

![max_af](https://github.com/fakedrtom/SVAFotate/blob/master/images/SV_Cov_example_fig.png)

All overlaps from all given sources in the BED file are considered when calculating `SV_Cov`, 
provided that the SVs share the same SVTYPE, but source specific `SV_Cov` annotations are also 
added. For example, if `CCDG` is one of the sources in the input BED file (and has not been excluded 
using the `-s` parameter), `CCDG_SV_Cov` will be added and reflects the `SV_Cov` calculation 
using only SVs from the input BED file that come from `CCDG`. All `SV_Cov` annotations measure 
*all* overlaps and not just those restricted by the `-f` option. This parameter does expect 
a minimum AF value to be included which restricts SVs from the BED file to only those meeting 
that value. In this manner, SVs from the BED file, perhaps very rare ones, can be exlcuded from 
the calculation of `SV_Cov`.

#### Unique SV Regions
```
  -u UNIQUE SV REGIONS, --uniq UNIQUE SV REGIONS
                        Generate a file of unique SV regions called 'unique.bed'. These regions reflect genomic space within the queried SV region that have not been previously observed with the same
                        SVTYPE. This will also add an annotation regarding the number of unique regions within a given SV. Please provide minimum AF to exclude all SVs from data sources with a total AF
                        below that value (must be between 0 and 1.0).
```

The opposite of the `-c` option is availble with the `-u` option. This will identify
"unique" regions within a given SV from the input VCF that are defined as sharing no overlap
with any of the SVs of the same SVTYPE in the BED file. This option will create an additional
output file named `unique.bed` which will list the coordinates of these unique regions, their
corresponding SVTYPE, the SV_ID of the SV from which they are found, as well the sample names that
are genotyped as heterozygous (0/1) and homozygous (1/1) for the alternate allele. The following
figure illustrates how these regions are determined in a variety of matching SV scenarios.

![max_af](https://github.com/fakedrtom/SVAFotate/blob/master/images/SV_Uniq_example_fig.png)

Additionally, the `-u` option will also add the `SV_Uniq` annotation to the VCF which represents
the number of unique regions found within a given SV from the input VCF. Similar to the `-c`
option, unique regions do not consider restrictions from the `-f` option and this parameter
also expects a minimum AF value to be provided that can restict the inclusion of SVs from the BED
file.

#### SV Size Limit
```
  -l SV SIZE LIMIT, --lim SV SIZE LIMIT
                        Only include previously observed SVs from data sources with a size less than or equal to this value (only available when using --cov or --uniq).
```

If any SVs in the BED file are exceedingly large, they may have big effects on the determination
of observed SV coverage and unique SV regions. For example, CCDG reports an exceptionally rare
deletion that is over 61Mb in size. This event is likely to overlap with many putative deletions
that likely represent distinct SV events. Considering this may obscure meaningful coverage and
unique region annotations and interpretations, the SV size limit `-l` is available. Without
including a size limit with the `-l` option all of these example events overlapping with the
very large and rare deletion from CCDG would show observed SV coverage of 1.0 and no
unique SV regions regardless of other more precise and potentially more meaningful comparisons.
By including an SV size limit alongside the `-c` or `-u` options, such scenarios can be avoided.

#### Targets BED File
```
  -t TARGETS BED FILE, --target TARGETS BED FILE
                        Path to target regions BED file. Expected format is a tab delimited file listing CHROM START END ID where ID is a genomic region identifier that will be listed as an annotation
                        if an overlap exists between a given SV and the target regions.
```

There may be particular genomic regions where an overlap with any reported SV event would be of
interest. Using the `-t` option and supplying a simple BED file consisting of CHROM, START, END,
and a region identifier will create a `Target_Overlaps` annotation that lists the supplied
region identifier for all overlaps between the SV and the regions in this targets BED file. Regions
of interest may include genes, specific exons, promoters, enhancers, etc. If the `-t` option is
used when the `-u` option has also been invoked, then the resulting `unique.bed` file from the
`-u` option will also include a column corresponding to any overlaps found between the unique
regions and the regions provided by the `-t` option. This column will be populated with the
region identifier found in the BED file used with the `-t` option.

#### Use CI Boundaries
```
  -ci USE CI BOUNDARIES, --ci USE CI BOUNDARIES
                        Expects CIPOS and CIEND to be included in the INFO field of the input VCF (--vcf). If argument is selected, use 'inner' or 'outer' confidence intervals (CIPOS, CIEND) for SV
                        boundaries. Choices = [in, out]
  -ci95 USE CI BOUNDARIES, --ci95 USE CI BOUNDARIES
                        Expects CIPOS95 and CIEND95 to be included in the INFO field of the input VCF (--vcf). If argument is selected, use 'inner' or 'outer' confidence intervals (CIPOS95, CIEND95) for
                        SV boundaries. Choices = [in, out]
```

Lumpy and other SV callers may generate confidence interval metrics for predicted SV breakpoints. To
use these confidence intervals instead of the given START and END coordinates, the `-ci` or `-ci95`
options may be used. There are two choices for each of these parameters: `in` and `out`. The `in` choice
will use the CIPOS and CIEND (or CIPOS95 and CIEND95) that will result in a smaller SV while the `out`
choice will use the intervals resulting in a larger SV. If CIPOS and CIEND or CIPOS95 and CIEND95 are
not present in the INFO field for a given SV in the input VCF these options will not work. Please be
sure to consider the adjustments to SV genomic coordinates when using the `-f` option for requiring
reciprocal overlaps.

#### Change SV Size
```
  -e EMBIGGEN THE SV SIZE, --emb EMBIGGEN THE SV SIZE
                        Increase the size of the SV coordinates in the input VCF (--vcf) by a single integer; Subtract that value from the start and add it to the end of each set of coordinates.
  -r REDUCE THE SV SIZE, --red REDUCE THE SV SIZE
                        Reduce the size of the SV coordinates in the input VCF (--vcf) by a single integer; Add that value to the start and subtract it from the end of each set of coordinates.
```

It may be helpful in some cases to adjust the reported SV breakpoints for the SVs in the
input VCF. The genomic space corresponding to a reported SV may be enlarged using the
`-e` option or reduced using the `-r` option. Please be sure to consider the adjustments to SV
genomic coordinates when using the `-f` option for requiring reciprocal overlaps.

#### CPU Count
```
  --cpu CPU Count       The number of cpus to use for multi-threading (Default = 1).
```

**Currently experiencing technical difficulties; mileage may vary** 

SVAFotate may use additional cpus for multi-threading purposes which can be 
designated using the `--cpu` option.

pickle-source
------------------------
**Please note that pickle-sources are currently unavaible for use with `annotate`; `pickle-source` 
should still work, but its application with `annotate` may render this less useful for now** 

Since SVAFotate may be used repeatedly on different SV VCFs with the same BED file, it
may be advantageous to create a pickle object of the BED file to improve SVAFotate's
performance when running the `annotate` subcommand. This is optional. If SVAFotate or
the BED file used with the `annotate` subcommand is seldom needed, creating a pickle object
of that BED file is probably not necessary. Nevertheless, creating a pickle object out
of a BED file is straightforward and can be accomplished with the following command:

```
svafotate pickle-source --bed SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz --out SVAFotate_core_SV_popAFs.GRCh38.v4.1.pickle
```

custom-annotation
-----------------------

This is not currently used, but is a placeholder for potential
future developments.

Acknowledgements 
========================
Thank you to Michael Cormier for code review, improvements, and contributions.
