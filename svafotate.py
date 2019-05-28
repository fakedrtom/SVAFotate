import sys
import cyvcf2
from cyvcf2 import Writer
import gzip
from pybedtools import BedTool
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-i', '--in',
                    metavar='INPUT VCF',
                    dest="i",
                    help='path to VCF to annotate')
parser.add_argument('-o', '--out',
                    metavar='OUTPUT VCF',
                    dest="o",
                    help='output VCF name/path')
parser.add_argument('-f', '--minf',
                    metavar='MINIMUM OVERLAP',
                    dest="f",
                    type=float,
                    help='minimum reciprocal overlap required between SVs (default 0.01, must be between 0 and 1.0)')
parser.add_argument('-ccdg', '--ccdg',
                    metavar='PATH TO CCDG',
                    dest="ccdg",
                    help='path to CCDG SV bed file')
parser.add_argument('-gnomad', '--gnomad',
                    metavar='PATH TO GNOMAD',
                    dest="gnomad",
                    help='path to gnomAD SV bed file')
#parser.add_argument('-ci', '--ci',
#                    dest="ci",
#                    action='store_true',
#                    help='option to use out confidence intervals (CIPOS95, CIEND95) for SV boundaries')

args = parser.parse_args()

if args.i is None:
    raise NameError('Must include path to VCF with option -i')
else:
    vcf = cyvcf2.VCF(args.i)
if args.f is None:
    minf = float(0.01)
else:
    minf = float(args.f)
if minf > 1 or minf < 0:
    raise NameError('minimum reciprocal overlap must be between 0 and 1.0')
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

def overlaps(bed,data_dict,columns):
    if data_dict == ccdg_AFs:
        intersect = bed.intersect(ccdgbed, wao = True, f = minf, r = True)
    elif data_dict == gnomad_AFs:
        intersect = bed.intersect(gnomadbed, wao = True, f = minf, r = True)
    for interval in intersect:
        chrom1,start1,end1 = interval[0],interval[1],interval[2]
        svtype1,svtype2 = interval[3],interval[7]
        sv = str(chrom1) + ':' + str(start1) + ':' +  str(end1) + ':' + svtype1
        afs = []
        if svtype1 == svtype2:
            for i in columns:
                af = '%.6f' % float(interval[i])
                afs.append(af)
            data_dict[sv].append(afs)

tmp = []
ccdg_AFs = {}
gnomad_AFs = {}
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
#    if args.ci:
#        cipos = v.INFO.get('CIPOS95')
#        ciend = v.INFO.get('CIEND95')
#        start = start + int(cipos[0])
#        end = end + int(ciend[1])
    sv = str(chrom) + ':' + str(start) + ':' + str(end) + ':' +svtype
    if sv not in ccdg_AFs:
        ccdg_AFs[sv] = []
    if sv not in gnomad_AFs:
        gnomad_AFs[sv] = []
    out = [str(chrom), str(start), str(end), svtype]
    tmp.append(out)

vcf.close(); vcf = cyvcf2.VCF(args.i)
tmpbed = BedTool(tmp)
if args.ccdg is not None:
    overlaps(tmpbed,ccdg_AFs,[8])
    vcf.add_info_to_header({'ID': 'CCDG_MaxAF', 'Description': 'The maximum AF from matching SVs with ' + str(minf) + ' overlaps with CCDG', 'Type': 'Float', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'CCDG_Count', 'Description': 'The number of matching SVs with ' + str(minf) + ' overlaps with CCDG', 'Type': 'Integer', 'Number': '1'})
if args.gnomad is not None:
    overlaps(tmpbed,gnomad_AFs,[8,9])
    vcf.add_info_to_header({'ID': 'gnomAD_MaxAF', 'Description': 'The maximum AF from matching SVs with ' + str(minf) + ' overlaps with gnomAD', 'Type': 'Float', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'gnomAD_PopMaxAF', 'Description': 'The maximum PopMax AF from matching SVs with ' + str(minf) + ' overlaps with gnomAD', 'Type': 'Float', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'gnomAD_Count', 'Description': 'The number of matching SVs with ' + str(minf) + ' overlaps with gnomAD', 'Type': 'Integer', 'Number': '1'})

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
#    if args.ci:
#        cipos = v.INFO.get('CIPOS95')
#        ciend = v.INFO.get('CIEND95')
#        start = start + int(cipos[0])
#        end = end + int(ciend[1])
    sv = str(chrom) + ':' + str(start) + ':' + str(end) + ':' +svtype
    if args.ccdg is not None:
        ccdg_maxAF = 0
        if len(ccdg_AFs[sv]) > 0:
            afs = []
            for i in ccdg_AFs[sv]:
                afs.append(i[0])
            ccdg_maxAF = max(afs)
        v.INFO['CCDG_MaxAF'] = ccdg_maxAF
        v.INFO['CCDG_Count'] = len(ccdg_AFs[sv])
    if args.gnomad is not None:
        gnomad_maxAF = 0
        gnomad_popmaxAF = 0
        if len(gnomad_AFs[sv]) > 0:
            afs = []
            pop_afs = []
            for i in gnomad_AFs[sv]:
                afs.append(i[0])
                pop_afs.append(i[1])
            gnomad_maxAF = max(afs)
            gnomad_popmaxAF = max(pop_afs)
        v.INFO['gnomAD_MaxAF'] = gnomad_maxAF
        v.INFO['gnomAD_PopMaxAF'] = gnomad_popmaxAF
        v.INFO['gnomAD_Count'] = len(gnomad_AFs[sv])
    new_vcf.write_record(v)

new_vcf.close(); vcf.close()
