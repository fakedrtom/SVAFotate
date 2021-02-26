import sys
import os
import pickle
import argparse 
import gzip
import io
import copy
from collections import defaultdict
from .utils import process_bed_source

def add_pickle_source(parser):
    
    p = parser.add_parser("pickle-source",
        help="Pickle Source Bed",
        description = "Pickle the Souce Annotation Bed file for improved annotation loading"
    )

    req = p.add_argument_group("Required Arguments")
    p._optionals.title = "help"

    req.add_argument("--bed",
                     metavar = "Souce BED File",
                     required = True, 
                     help = "Path and/or name of the source bed file to pickle"
    )

    req.add_argument("--out",
                     metavar = "Ouput pickle file",
                     required = True,
                     help = "Path and/or name of the output pickle file to create"
    )

    p.set_defaults(func=pickle_source)



def pickle_source(parser,args):

    ## Get data objects from bed file
    s, d, bl, bh = process_bed_source(args.bed)

    ## Cast them to non local function objects
    sources = list(s)
    datas = dict(d)
    bed_lists = dict(bl)
    bed_headers = list(bh)

    ## Pickle annotation objects
    print("\nPickling data")
    out_fh = open(args.out if args.out.endswith(".pickle") else args.out + ".pickle", "wb")
    pickle.dump({"sources":sources, "datas":datas, "bed_lists":bed_lists, "bed_headers":bed_headers}, out_fh)
    out_fh.close()

    print("DONE")

    
