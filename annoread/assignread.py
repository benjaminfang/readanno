#! /usr/bin/env python3
"""
Assign the read which mapped to multiple features.

The feature can devide into four level.

    Gene_name: the name of gene.
    Gene_id: one gene_name may have multiple gene_id feature.
    Transcript_id: one gene_id may master multiple transcript_id

The choice or method to assign read when it mapped to multiple features.
    
    keep: the read is kept for all features. one read count multiple time.
    drop: the read is dropped for any features.
    equal: the read is assign to a feautre randomly in equal probability.
    proportion: the read is assign to a feature randomly in probability proportional
        to unique read of each features.
    largest: The read is assign to the feature which have the largest amount unique
        read.
"""
import argparse
from .annodata import ANNOREAD_DATA


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("anno_file", help="annoread file")
    parser.add_argument("--out", "-O", help="output file name", default="out")
    parser.add_argument("--gene_name", 
        choices={"keep", "drop", "equal", "proportion", "largest"},
        help="the method applied to assign read in gene_name level",
        default=None)
    parser.add_argument("--gene_id",
        choices={"keep", "drop", "equal", "proportion", "largest"},
        help="the method applied to assign read in gene_id level",
        default=None)
    parser.add_argument("--transcript_id",
        choices={"keep", "drop", "equal", "proportion", "largest"},
        help="the method applied to assign read in transcript_id level",
        default=None)
    parser.add_argument("--unique_transcript_read",
        help="flag, indicate whether unique read mapped to same transcripts in multiple position",
        action="store_true")
    args = parser.parse_args()
    return (args.anno_file, args.out, args.gene_name, args.gene_id,
            args.transcript_id, args.unique_transcript_read)


def output_result(anno_data, out):
    fout = open(out + ".rdpos", "w")
    data = anno_data.get_data()
    for gene_name in data:
        print(">" + gene_name, file=fout)
        for gene_id in data[gene_name]:
            ref_id = data[gene_name][gene_id]["ref_id"]
            gene_range = data[gene_name][gene_id]["gene_range"]
            gene_range = ",".join([str(ele) for ele in gene_range])
            ori = data[gene_name][gene_id]["ori"]
            print("\t".join(["&" + gene_id, ref_id, gene_range, ori]), file=fout)
            transcripts = data[gene_name][gene_id]["transcripts"]
            for transcript_id in transcripts:
                gbkey = transcripts[transcript_id]["gbkey"]
                transcript_range = transcripts[transcript_id]["transcript_range"]
                transcript_range = ";".join([",".join([str(ele[0]), str(ele[1])]) for ele in transcript_range])
                start_codon = transcripts[transcript_id]["start_codon"]
                if start_codon:
                    start_codon = ";".join([",".join([str(ele[0]), str(ele[1])]) for ele in start_codon])
                else:
                    start_codon = "*"
                stop_codon = transcripts[transcript_id]["stop_codon"]
                if stop_codon:
                    stop_codon = ";".join([",".join([str(ele[0]), str(ele[1])]) for ele in stop_codon])
                else:
                    stop_codon = "*"
                print("\t".join(["$" + transcript_id, gbkey, transcript_range,
                                start_codon, stop_codon]),
                                file=fout)
                reads = transcripts[transcript_id]["reads"]
                for read in reads:
                    read_id, read_pos, read_len, align_pos = read
                    print("\t".join([read_id, str(read_pos), str(read_len),
                        ";".join([",".join([str(ele[0]), str(ele[1])]) for ele in align_pos])]),
                        file=fout)
    fout.close()


def main():
    anno_file, out, gene_name_lev, gene_id_lev, transcript_id_lev, uniq_transcript \
        = getargs()
    print("Start...", anno_file)

    print("Reading annoread data...")
    anno_data = ANNOREAD_DATA(anno_file)

    if gene_name_lev:
        print("Assigning reads across gene_name...")
        anno_data.assign_reads_across_gene_name(method=gene_name_lev)

    if gene_id_lev:
        print("Assigning reads across gene_id...")
        anno_data.assign_reads_across_gene_id(method=gene_id_lev)

    if transcript_id_lev:
        print("Assigning reads across transcript_id...")
        anno_data.assign_reads_across_transcript_id(method=transcript_id_lev)

    if uniq_transcript:
        print("Uniquing reads mapped to multiple postion of a transcript...")
        anno_data.unique_reads_one_transcript()
    
    print("outputing result...")
    output_result(anno_data, out)

    print("Done.")

if __name__ == "__main__":
    main()
