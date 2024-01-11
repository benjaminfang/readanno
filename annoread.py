#! /usr/bin/env python3
"""
Annotate the mapped reads.

The reference genome can be genome or transcriptome. The corresponding
gtf file is need for the annotation process.

"""
import argparse
from biobrary.bioparse import GTF
from biobrary.bioparse import FASTA
import re


def getargs():
    """
    Parse arguments

    Parameters
    ------------
    
    Retures
    -------------
    ref_type:
    gtf:
    out:
    samfile:

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref-type", "-R", required=True,
                        choices=["geno", "tran"],
                        help="The reference type, genome of transcriptom")
    parser.add_argument("--gtf", "-G", default=None,
                        help="annotation file, gtf file.")
    parser.add_argument("--fasta", "-F", default=None,
                        help="fasta file which roled as reference.")
    parser.add_argument("--out", "-O", default="out",
                        help="name of output file.")
    parser.add_argument("samfile", help="sam file.")
    args = parser.parse_args()
    
    if args.ref_type == "geno" and (not args.gtf):
        print("gtf file required for geno ref-type.")
        parser.print_help()
        exit()
    elif args.ref_type == "tran" and (not args.fasta):
        print("transcriptom fasta is required for tran ref-type.")
        parser.print_help()
        exit()
        
    return args.ref_type, args.gtf, args.fasta, args.out, args.samfile


def parse_gtf(gtf_file):
    """
    Parser GTF file using biobrary

    Parameters
    --------------
    gtf_file: file name of gtf file

    Returns
    --------------
    transid_gene_dic: A dictionary
    {transcriptid: geneid, ...}

    """
    trans_gene_dic = {}
    gtf = GTF(gtf_file)
    for gene in gtf:
        geneid = gene.get_geneid()
        for transid in gene.get_transcriptids():
            if transid not in trans_gene_dic:
                trans_gene_dic[transid] = geneid
            else:
                print(f"{transid} has exists in the dictionary.")
                
    return trans_gene_dic


def parse_seqinfo(info_line):
    re_tmp = re.compile(r'\[(.+?=.+?)\]')
    finding = re_tmp.findall(info_line)
    info = {}
    for ele in finding:
        key, value = ele.split("=")
        info[key] = value
    return info


def parse_fasta(fastqfile):
    seqid_geneid_dic = {}
    fasta = FASTA(fastqfile)
    seqinfo = fasta.seqid_info

    for seqid in seqinfo:
        seq_info_dic = parse_seqinfo(seqinfo[seqid])
        seqid_geneid_dic[seqid] = seq_info_dic["gene"]

    return seqid_geneid_dic


def parse_refname(refname):
    """
    parse refname
    """
    refname = refname.split("|")[1].split("_")
    if len(refname) != 6:
        return None, None, None
    genome_contig = "_".join(refname[:2])
    transtype = refname[2]
    transid = "_".join(refname[3: 5])

    return genome_contig, transtype, transid


def annread(samfile, trans_gene_dic, seqid_geneid_dic, out):
    """
    """
    fin = open(samfile, "r")
    fout = open(out, "w")
    for line in fin:
        if line.startswith("@"):
            continue
        line = line.rstrip().split("\t")
        readid, flg, refname, pos, mapq, cigar = line[: 6]
        if flg == "99" or flg == "355":
            genome_contig, transtype, transid = parse_refname(refname)
            if transid:
                geneid = trans_gene_dic[transid]
            else:
                geneid = seqid_geneid_dic[refname]
            print("\t".join([readid, flg, refname, pos, mapq, cigar, geneid]), file=fout)

    fin.close()
    fout.close()
    return


def main():
    """
    Main entry of the program

    Arguments
    -----------

    Returns
    ------------
    status: 0 for nornal; other for abnormal.

    """
    reftype, gtffile, fastafile, out, samfile = getargs()
    if reftype == "tran":
        trans_gene_dic = parse_gtf(gtffile)
        seqid_geneid_dic = parse_fasta(fastafile)
        annread(samfile, trans_gene_dic, seqid_geneid_dic, out)
    elif reftype == "geno":
        print("The geno reference type supporting is one the way.")
    else:
        print("reference type not recognized.")
    return 0


if __name__ == "__main__":
    main()