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
