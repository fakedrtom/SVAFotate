An updated core BED file is now available that integrates SVs from 
gnomAD v4.1. The size of this file prevents it from being hosted 
in this repo, but it can be downloaded from [here](https://zenodo.org/records/11642574). 
This BED file contains reported SVs from CCDG, gnomAD, 1000G, and TOPMed. 
Because this version of gnomAD includes new populations, the SV calls 
from CCDG, 1000G, and TOPMed have been reformatted to account for these 
additional populations. Please note that 1000G SVs are identified 
with `ThousG` in this file to better aid downstream use of VCFs 
annotated by SVAFotate. Since gnomAD v4.1 was aligned to GRCh38, no 
liftover was required to create this data. Since this file did not 
require any genomic coordiantes liftover, and because it contains more 
up to date SV calls from population datasets, it is recommended that 
this BED file be used over the older files listed (and described below) 
in this repo. 

The core BED files provided here are made up of SV calls 
from CCDG, gnomAD, and 1000G. CCDG and 1000G were called 
using the GRCh38 reference while the gnomAD calls were 
made using the GRCh37 reference. To create reference-specific 
versions using all three data sets, necessary liftovers 
were performed. This was done using [UCSC's liftOver tool and available chain files](https://genome.ucsc.edu/cgi-bin/hgLiftOver)
. For the SVAFotate_core_SV_popAFs.GRCh37.bed.gz file this 
required CCDG and 1000G data to be converted while the 
SVAFotate_core_SV_popAFs.GRCh38.bed.gz file only required 
the gnomAD data to be converted.

Thanks to Adam English and the Sedlazeck Lab for making SV 
calls from TOPMed available to SVAFotate. These SVs and their 
AFs can be found in the TOPMed.GRCH38.bed.gz file. This can be 
used as the sole input BED for SVAFotate or it can be combined with the 
SVAFotate_core_SV_popAFs.GRCh38.bed.gz core BED file for more 
comprehensive SV annotation.
