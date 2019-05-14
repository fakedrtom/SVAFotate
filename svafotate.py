import os
import sys
import cyvcf2
import gzip
from pybedtools import BedTool

vcf = cyvcf2.VCF(sys.argv[1], gts012=True)
ccdg = gzip.open('/uufs/chpc.utah.edu/common/HIPAA/u0055382/genome_ref/ccdg_v37_sv.sites.bed.gz', 'r')
ccdgbed = BedTool(ccdg)
gnomad = gzip.open('/uufs/chpc.utah.edu/common/HIPAA/u0055382/genome_ref/gnomad_v2_sv.sites.bed.gz', 'r')
gnomadbed = BedTool(gnomad)

def ccdg_overlaps(sv):
    intersect = sv.intersect(ccdgbed, wao = True, f = 0.5, r = True)
    for interval in intersect:
        chrom1,start1,end1 = interval[0],interval[1],interval[2]
        svtype1,svtype2 = interval[3],interval[7]
        sv = str(chrom1) + ':' + str(start1) + ':' +  str(end1) + ':' + svtype1
        if svtype1 == svtype2:
            af = '%.6f' % float(interval[8])
            ccdg_AFs[sv].append(af)

def gnomad_overlaps(sv):
    intersect = sv.intersect(gnomadbed, wao = True, f = 0.5, r = True)
    for interval in intersect:
        chrom1,start1,end1 = interval[0],interval[1],interval[2]
        svtype1,svtype2 = interval[3],interval[8]
        gnomad_filter = interval[10]
        sv = str(chrom1) + ':' + str(start1) + ':' +  str(end1) + ':' + svtype1
        if svtype1 == svtype2 and gnomad_filter == 'PASS':
            af = '%.6f' % float(interval[34])
            gnomad_AFs[sv].append(af)
            popaf = '%.6f' % float(interval[42])
            gnomad_popAFs[sv].append(popaf)

tmp = []
ccdg_AFs = {}
gnomad_AFs = {}
gnomad_popAFs = {}
for v in vcf:
    chrom = v.CHROM
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    start = int(v.POS)-1
    end = v.INFO.get('END')
    if v.INFO.get('SVTYPE') is not None:
        svtype=v.INFO.get('SVTYPE')
    if svtype == 'BND':
        end = int(v.POS)
    sv = str(chrom) + ':' + str(start) + ':' + str(end) + ':' +svtype
    if sv not in ccdg_AFs:
        ccdg_AFs[sv] = []
    if sv not in gnomad_AFs:
        gnomad_AFs[sv] = []
    if sv not in gnomad_popAFs:
        gnomad_popAFs[sv] = []
    out = [str(chrom), str(start), str(end), svtype]
    tmp.append(out)

tmpbed = BedTool(tmp)
ccdg_overlaps(tmpbed)
gnomad_overlaps(tmpbed)
for interval in tmpbed:
    chrom,start,end,svtype = interval
    sv = str(chrom) + ':' + str(start) + ':' + str(end) + ':' +svtype
    ccdg_maxAF = 0
    if len(ccdg_AFs[sv]) > 0:
        ccdg_maxAF = max(ccdg_AFs[sv])
    gnomad_maxAF = 0
    gnomad_popmaxAF = 0
    if len(gnomad_AFs[sv]) > 0:
        gnomad_maxAF = max(gnomad_AFs[sv])
        gnomad_popmaxAF = max(gnomad_popAFs[sv])
    print chrom,start,end,svtype,ccdg_maxAF,len(ccdg_AFs[sv]),gnomad_maxAF,gnomad_popmaxAF,len(gnomad_AFs[sv])
