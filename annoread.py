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
    parser.add_argument("--ref", "-R", required=True,
                        choices=["geno", "tran"],
                        help="The reference type, genome of transcriptom,\
                            transcriptom fasta is from NCBI")
    parser.add_argument("--ref-type", "-T", required=True,
                        choices=["ncbi-genome", "ensamble-genome", "ncbi-rna", "ncbi-genome-rna"],
                        help="reference data origin")
    
    parser.add_argument("--gtf", "-G", required=True,
                        help="annotation file, gtf file.")
    parser.add_argument("--fasta", "-F", required=True,
                        help="fasta file which roled as reference.")
    
    parser.add_argument("--out", "-O", default="out",
                        help="name of output file.")
    
    parser.add_argument("samfile", help="sam file.")
    args = parser.parse_args()

    return args.ref, args.ref_type, args.gtf, args.fasta, args.out, args.samfile


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
    transid_geneid_dic = {}
    gtf = GTF(gtf_file)
    for gene in gtf:
        geneid = gene.get_geneid()
        for transid in gene.get_transcriptids():
            if transid not in transid_geneid_dic:
                transid_geneid_dic[transid] = geneid
            else:
                print(f"{transid} has exists in the dictionary.")
                
    return transid_geneid_dic


def parse_seqinfo(info_line):
    re_tmp = re.compile(r'\[(.+?=.+?)\]')
    finding = re_tmp.findall(info_line)
    info = {}
    for ele in finding:
        key, value = ele.split("=")
        info[key] = value
    return info


def parse_fasta_genome_rna(fastafile):
    seqid_seqinfo_dic = {}
    fasta = FASTA(fastafile)
    seqinfo = fasta.seqid_info

    for seqid in seqinfo:
        seq_info_dic = parse_seqinfo(seqinfo[seqid])
        seqid_seqinfo_dic[seqid] = seq_info_dic

    return seqid_seqinfo_dic


def annread_tran_genome_rna(samfile, seqid_seqinfo_dic, out):
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
            contig_id = "_".join(refname.split("|")[1].split("_")[:2])
            seqinfo = seqid_seqinfo_dic[refname]
            transid = seqinfo.get("transcript_id")
            gbkey = seqinfo.get("gbkey")
            if not transid:
                transid = "*"
            if not gbkey:
                gbkey = "*"
            geneid = seqinfo["gene"]
            
            print("\t".join([readid, flg, refname, pos, mapq, cigar, contig_id,\
                            transid, gbkey, geneid]), file=fout)
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
    ref, reftype, gtffile, fastafile, out, samfile = getargs()
    if ref == "geno":
        print("reference type not recognized.")
    elif ref == "tran":
        if reftype == "ncbi-rna":
            pass
        elif reftype == "ncbi-genome-rna":
            seqid_seqinfo_dic = parse_fasta_genome_rna(fastafile)
            annread_tran_genome_rna(samfile, seqid_seqinfo_dic, out)
        else:
            pass
    else:
        pass

        
    return 0


if __name__ == "__main__":
    main()