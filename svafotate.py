import sys
import cyvcf2
from cyvcf2 import Writer
import gzip
from pybedtools import BedTool
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-i',
                    metavar='STRING',
                    help='path to VCF to annotate')
parser.add_argument('-f', 
                    metavar='FLOAT', 
                    help='minimum overlap required between SVs (default 0.1)')
parser.add_argument('-o',
                    metavar='STRING',
                    help='output VCF name/path')
parser.add_argument('-ccdg',
                    metavar='STRING',
                    help='path to CCDG SV bed file')
parser.add_argument('-gnomad',
                    metavar='STRING',
                    help='path to gnomAD SV bed file')

args = parser.parse_args()

if args.i is None:
    raise NameError('Must include path to VCF with option -i')
else:
    vcf = cyvcf2.VCF(args.i)
if args.f is None:
    minf = float(0.1)
else:
    minf = float(args.f)
if args.o is None:
    raise NameError('Must include name/path to output VCF with option -o')
else:
    output_vcf = args.o
if args.ccdg is not None:
    ccdg = gzip.open(args.ccdg, 'r')
    ccdgbed = BedTool(ccdg)
if args.gnomad is not None:
    gnomad = gzip.open(args.gnomad, 'r')
    gnomadbed = BedTool(gnomad)
if args.ccdg is None and args.gnomad is None:
    raise NameError('Please include something to annotate with -ccdg or -gnomad')

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
        svtype1,svtype2 = interval[3],interval[7]
        sv = str(chrom1) + ':' + str(start1) + ':' +  str(end1) + ':' + svtype1
        if svtype1 == svtype2:
            af = '%.6f' % float(interval[8])
            gnomad_AFs[sv].append(af)
            popaf = '%.6f' % float(interval[9])
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

vcf.close(); vcf = cyvcf2.VCF(args.i)
tmpbed = BedTool(tmp)
if args.ccdg is not None:
    ccdg_overlaps(tmpbed)
    vcf.add_info_to_header({'ID': 'CCDG_MaxAF', 'Description': 'The maximum AF from matching SV overlaps with CCDG', 'Type': 'Float', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'CCDG_Count', 'Description': 'The number of matching SV overlaps with CCDG', 'Type': 'Integer', 'Number': '1'})
if args.gnomad is not None:
    gnomad_overlaps(tmpbed)
    vcf.add_info_to_header({'ID': 'gnomAD_MaxAF', 'Description': 'The maximum AF from matching SV overlaps with gnomAD', 'Type': 'Float', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'gnomAD_PopMaxAF', 'Description': 'The maximum PopMax AF from matching SV overlaps with gnomAD', 'Type': 'Float', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'gnomAD_Count', 'Description': 'The number of matching SV overlaps with gnomAD', 'Type': 'Integer', 'Number': '1'})

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
    if args.ccdg is not None:
        ccdg_maxAF = 0
        if len(ccdg_AFs[sv]) > 0:
            ccdg_maxAF = max(ccdg_AFs[sv])
        v.INFO['CCDG_MaxAF'] = ccdg_maxAF
        v.INFO['CCDG_Count'] = len(ccdg_AFs[sv])
    if args.gnomad is not None:
        gnomad_maxAF = 0
        gnomad_popmaxAF = 0
        if len(gnomad_AFs[sv]) > 0:
            gnomad_maxAF = max(gnomad_AFs[sv])
            gnomad_popmaxAF = max(gnomad_popAFs[sv])
        v.INFO['gnomAD_MaxAF'] = gnomad_maxAF
        v.INFO['gnomAD_PopMaxAF'] = gnomad_popmaxAF
        v.INFO['gnomAD_Count'] = len(gnomad_AFs[sv])
    new_vcf.write_record(v)

new_vcf.close(); vcf.close()
