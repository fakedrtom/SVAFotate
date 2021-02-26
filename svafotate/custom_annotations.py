import sys
import os
import argparse



def add_custom_annotation(parser):
    
    p = parser.add_parser("custom-annotation",
        help="Add custom annotation(s) to source annotation file",
        description="Add custom annotation(s) to the source annotation file, which can then be used to annotate an SV VCF file"
    )

    p._optionals.title = "help"
    req = p.add_argument_group("Required Arguments")

    req.add_argument("-f",
                     metavar="Input File",
                     required = True,
                     help="The Input file to use"
    )

    p.set_defaults(func=custom_anno)



def custom_anno(parser,args):
    
    print("\n\n=======================================")
    print("NOT IMPLEMENTED YET")
    print("=======================================\n\n")



