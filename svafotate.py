import sys
import cyvcf2
from cyvcf2 import Writer
import gzip
from pybedtools import BedTool
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-v',
                    metavar='STRING',
                    help='path to VCF to annotate')
parser.add_argument('-f', 
                    metavar='FLOAT', 
                    help='minimum overlap required between SVs (default 0.1)')
parser.add_argument('-o',
                    metavar='STRING',
                    help='output VCF name/path')

args = parser.parse_args()

if args.v is None:
    raise NameError('Must include path to VCF with option -v')
else:
    vcf = cyvcf2.VCF(args.v)
if args.f is None:
    minf = float(0.1)
else:
    minf = float(args.f)
if args.o is None:
    raise NameError('Must include name/path to output VCF with option -o')
else:
    output_vcf = args.o

ccdg = gzip.open('/uufs/chpc.utah.edu/common/HIPAA/u0055382/genome_ref/ccdg_v37_sv.sites.bed.gz', 'r')
ccdgbed = BedTool(ccdg)
gnomad = gzip.open('/uufs/chpc.utah.edu/common/HIPAA/u0055382/genome_ref/gnomad_v2_sv.sites.bed.gz', 'r')
gnomadbed = BedTool(gnomad)

def ccdg_overlaps(sv):
    intersect = sv.intersect(ccdgbed, wao = True, f = minf, r = True)
    for interval in intersect:
        chrom1,start1,end1 = interval[0],interval[1],interval[2]
        svtype1,svtype2 = interval[3],interval[7]
        sv = str(chrom1) + ':' + str(start1) + ':' +  str(end1) + ':' + svtype1
        if svtype1 == svtype2:
            af = '%.6f' % float(interval[8])
            ccdg_AFs[sv].append(af)

def gnomad_overlaps(sv):
    intersect = sv.intersect(gnomadbed, wao = True, f = minf, r = True)
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

vcf.close(); vcf = cyvcf2.VCF(args.v)
tmpbed = BedTool(tmp)
ccdg_overlaps(tmpbed)
gnomad_overlaps(tmpbed)

vcf.add_info_to_header({'ID': 'CCDG_maxAF', 'Description': 'The maximum AF from matching SV overlaps with CCDG', 'Type': 'Float', 'Number': '1'})
vcf.add_info_to_header({'ID': 'CCDG_count', 'Description': 'The number of matching SV overlaps with CCDG', 'Type': 'Integer', 'Number': '1'})
vcf.add_info_to_header({'ID': 'gnomAD_maxAF', 'Description': 'The maximum AF from matching SV overlaps with gnomAD', 'Type': 'Float', 'Number': '1'})
vcf.add_info_to_header({'ID': 'gnomAD_popmaxAF', 'Description': 'The maximum PopMax AF from matching SV overlaps with gnomAD', 'Type': 'Float', 'Number': '1'})
vcf.add_info_to_header({'ID': 'gnomAD_count', 'Description': 'The number of matching SV overlaps with gnomAD', 'Type': 'Integer', 'Number': '1'})

new_vcf = Writer(output_vcf, vcf)
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
    ccdg_maxAF = 0
    if len(ccdg_AFs[sv]) > 0:
        ccdg_maxAF = max(ccdg_AFs[sv])
    gnomad_maxAF = 0
    gnomad_popmaxAF = 0
    if len(gnomad_AFs[sv]) > 0:
        gnomad_maxAF = max(gnomad_AFs[sv])
        gnomad_popmaxAF = max(gnomad_popAFs[sv])
    v.INFO['CCDG_maxAF'] = ccdg_maxAF
    v.INFO['CCDG_count'] = len(ccdg_AFs[sv])
    v.INFO['gnomAD_maxAF'] = gnomad_maxAF
    v.INFO['gnomAD_popmaxAF'] = gnomad_popmaxAF
    v.INFO['gnomAD_count'] = len(gnomad_AFs[sv])
    new_vcf.write_record(v)

new_vcf.close(); vcf.close()

#for interval in tmpbed:
#    chrom,start,end,svtype = interval
#    sv = str(chrom) + ':' + str(start) + ':' + str(end) + ':' +svtype
#    ccdg_maxAF = 0
#    if len(ccdg_AFs[sv]) > 0:
#        ccdg_maxAF = max(ccdg_AFs[sv])
#    gnomad_maxAF = 0
#    gnomad_popmaxAF = 0
#    if len(gnomad_AFs[sv]) > 0:
#        gnomad_maxAF = max(gnomad_AFs[sv])
#        gnomad_popmaxAF = max(gnomad_popAFs[sv])
#    print chrom,start,end,svtype,ccdg_maxAF,len(ccdg_AFs[sv]),gnomad_maxAF,gnomad_popmaxAF,len(gnomad_AFs[sv])
