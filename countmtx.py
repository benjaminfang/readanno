#! /usr/bin/env python3
"""
Read in annread output, and creat count matrix file.
"""

import argparse


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("annoread", help="annoread output file", nargs="+")
    parser.add_argument("--out", "-O", help="output file name", default="out")
    args = parser.parse_args()
    return args.annoread, args.out


def struc_annoread(annoread_file):
    """
    Structure Annotation File

    Parameters
    ---------------
    annoread_file: The result file generate by annoread.py

    Returns
    ------------
    dt_out: Structured data of annoread_file.
    """
    
    dt_out = {}
    fin = open(annoread_file, "r")

    fin.close()


if __name__ == "__main__":
    anno_files, out = getargs()

    for annofile in anno_files:
        pass