#! /usr/bin/env python3
"""
Read read postion file. count the unique read number belong to a feature.

There are three types of fetures or feature levels.

    Gene_name
    Gene_id
    Transcript_id

A count matrix file is outputed in following format:

feature_id, sample_1, sample_2, ...
id_1, number_1, number_2
id_2, number_3, number_4
.
.
.
"""

import numpy as np
import os
import argparse
from .annodata import ANNOREAD_DATA


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("anno_files", help="annoread output file", nargs="+")
    parser.add_argument("--out", "-O", help="output file name", default="out")
    parser.add_argument("--count_level", choices={"gene_name", "gene_id", "transcript_id"},
                        help="count feature level", default="gene_name")
    args = parser.parse_args()
    return args.anno_files, args.out, args.count_level 


def counting(annoread_dt, count_level="gene_name"):
    count_res = {}
    data = annoread_dt.get_data()
    if count_level == "gene_name":
        for gene_name in data:
            key = gene_name
            count_res[key] = set()
            for gene_id in data[gene_name]:
                for transcript_id in data[gene_name][gene_id]["transcripts"]:
                    for read in data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]:
                        count_res[key].add(read[0])

    elif count_level == "gene_id":
        for gene_name in data:
            for gene_id in data[gene_name]:
                key = gene_name + "|" + gene_id
                count_res[key] = set()
                for transcript_id in data[gene_name][gene_id]["transcripts"]:
                    for read in data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]:
                        count_res[key].add(read[0])
    
    elif count_level == "transcript_id":
        for gene_name in data:
            for gene_id in data[gene_name]:
                for transcript_id in data[gene_name][gene_id]["transcripts"]:
                    key = gene_name + "|" + gene_id + "|" + transcript_id
                    count_res[key] = set()
                    for read in data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]:
                        count_res[key].add(read[0])
    else:
        pass

    for key in count_res:
        count_res[key] = len(count_res[key])
    return count_res


def output_count_mtx(count_data, anno_files, out):
    all_id_s = set()
    for count_data_ele in count_data:
        set_gene = set()
        for key in count_data_ele:
            if count_data_ele[key] > 0:
                set_gene.add(key)
        all_id_s |= set_gene
    all_id_s = list(all_id_s)
    all_id_s.sort()
    
    mtx = []
    for count_data_ele in count_data:
        one_file_dt = []
        for id in all_id_s:
            one_file_dt.append(count_data_ele.get(id, 0))
        mtx.append(one_file_dt)
    mtx = np.asarray(mtx).T
    
    fout = open(out + '.tsv', "w")
    print("\t".join(["feature_id"] + [os.path.basename(ff) for ff in anno_files]), file=fout)
    for idx, row in enumerate(mtx):
        print("\t".join([all_id_s[idx]] + [str(ele) for ele in row]), file=fout)
    fout.close()


def main():
    anno_file_s, out, count_level = getargs()
    count_res_all = []
    
    for anno_file in anno_file_s:
        print("reading data...", anno_file)
        annodata = ANNOREAD_DATA(anno_file)
        print("counting...")
        count_res = counting(annodata, count_level=count_level)
        count_res_all.append(count_res)
    
    print("writing result...")
    output_count_mtx(count_res_all, anno_file_s, out)
    print("Done.")


if __name__ == "__main__":
    main()
