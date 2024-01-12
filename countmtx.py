#! /usr/bin/env python3
"""
Read in annread output, and creat count matrix file.
"""

import argparse
import os

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("annoread", help="annoread output file", nargs="+")
    parser.add_argument("--out", "-O", help="output file name", default="out")
    args = parser.parse_args()
    return args.annoread, args.out


def struc_annoread(annoread_file):
    read_gene_dic = {}
    fin = open(annoread_file, "r")
    for line in fin:
        line = line.rstrip().split("\t")
        readid, gene = line[0], line[-1]
        if readid not in read_gene_dic:
            read_gene_dic[readid] = set([gene])
        else:
            read_gene_dic[readid].add(gene)
    fin.close()
    return read_gene_dic


def filter_read(read_gene_dic):
    data_out = {}
    uniq_gene_read = 0
    mult_gene_read = 0
    for read in read_gene_dic:
        if len(read_gene_dic[read]) == 1:
            data_out[read] = read_gene_dic[read].pop()
            uniq_gene_read += 1
        else:
            mult_gene_read += 1
    print(f"read mapped to one gene: {uniq_gene_read}")
    print(f"read mapped to multiple genes: {mult_gene_read}")
    return data_out


def count_genes_read(read_gene_dic):
    gene_read_count = {}
    for read, gene in read_gene_dic.items():
        if gene not in gene_read_count:
            gene_read_count[gene] = 1
        else:
            gene_read_count[gene] += 1

    return gene_read_count


def main():
    annoread_files, out = getargs()
    count_dics = []
    for anno_file in annoread_files:
        print(f"Processing {anno_file}")
        read_gene_dic = struc_annoread(anno_file)
        read_gene_dic = filter_read(read_gene_dic)
        gene_count = count_genes_read(read_gene_dic)
        count_dics.append(gene_count)
    
    fout = open(out, "w")
    headline = "\t".join(["geneid"] + [".".join(os.path.basename(ff).split(".")[:-1]) for ff in annoread_files])
    print(headline, file=fout)
    gene_union = set()
    for dic in count_dics:
        gene_union |= set(dic.keys())
    gene_union = list(gene_union)
    gene_union.sort()
    for gene in gene_union:
        counts = []
        for dic in count_dics:
            count_gene = dic.get(gene)
            if count_gene:
                counts.append(str(count_gene))
            else:
                counts.append(str(0))
        print("\t".join([gene] + counts), file=fout)
    fout.close()
    return 0
    

if __name__ == "__main__":
    main()