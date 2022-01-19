The core BED files provided here are made up of SV calls 
from CCDG, gnomAD, and 1000G. CCDG and 1000G were called 
using the GRCh38 reference while the gnomAD calls were 
made using the GRCh37 reference. To create reference-specific 
versions using all three data sets, necessary liftovers 
were performed. This was done using [UCSC's liftOver tool and available chain files](https://genome.ucsc.edu/cgi-bin/hgLiftOver)
. For the SVAFotate_core_SV_popAFs.GRCh37.bed.gz file this 
required CCDG and 1000G data to be converted while the 
SVAFotate_core_SV_popAFs.GRCh37.bed.gz file only required 
the gnomAD data to be converted.
