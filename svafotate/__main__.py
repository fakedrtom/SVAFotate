import argparse
import sys

from .__init__ import __version__
from .svafotate_main import add_annotation
from .pickle_source import add_pickle_source 
from .custom_annotations import add_custom_annotation

def main(args=None):
    if args is None:
        args = sys.argv[1:]


    parser = argparse.ArgumentParser(
        prog="svafotate", 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""SVAFotate: Structural Variant VCF Annotation tools
                       =================================================="""
    )


    parser.add_argument(
        "-v",
        "--version",
        help="Installed version",
        action="version",
        version="%(prog)s " + str(__version__),
    )


    sub = parser.add_subparsers(title="[sub-commands]", dest="command")
    sub.required = True

    add_annotation(sub)

    add_pickle_source(sub)

    add_custom_annotation(sub)

    args = parser.parse_args(args)
    args.func(parser, args)
    

if __name__ == "__main__":
    sys.exit(main() or 0)
