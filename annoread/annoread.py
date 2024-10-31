#!/usr/bin/env python3
"""
Annotated read associated features based on mapping data.

Output read position file, ".rdpos", with format:
'
>gene_name
&gene_id seqid gene_range orientation
$transcript_id gbkey transcript_blocks_relative_to_gene start_codon/* stop_codon/*
read_id read_pos_in_trans read_length read_pos_align_with_the_gene
'
"""

import argparse
import biobrary.bioparse as bp
import re
import os

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("samfile", help="sam files")
    parser.add_argument("--out", default="out", help="out file")
    parser.add_argument("--fasta_file", help="NCBI rna from genome data.", 
                        required=True)
    parser.add_argument("--gtf_file", help="NCBI GTF file", required=True)
    parser.add_argument("--read_type", choices={"se", "pe"}, required=True, 
                        help="sequence type, single end or pair end")

    args = parser.parse_args()
    return args.samfile, args.out, args.fasta_file, args.gtf_file, args.read_type


def get_gtf_trans_cds_info(gtf_file):
    """
    dt_out: {gene_name: {gene_id: [refname, gene_range, ori, 
    {trans_id: [exon_range, gbkey, cds_range, start_codon_range, stop_codon_range],
    ...}], ...}, ...}
    """
    dt_out = {}
    transcript_id_gene_id = {}
    transcript_id_gene = {} #maybe useful in future
    gtf = bp.parse_gtf(gtf_file)
    for gene in gtf:
        gene_id = gene.get_gene_id()
        gene_name = gene.get_attr("gene")
        ori = gene.get_ori()
        gene_range = gene.get_range()
        refname = gene.get_seq_name()
        if gene_name not in dt_out:
            dt_out[gene_name] = {gene_id: [refname, gene_range, ori, {}]}
        else:
            dt_out[gene_name][gene_id] = [refname, gene_range, ori, {}]
        trans_ids = gene.get_transcript_id_s()
        for trans_id in trans_ids:
            trans = gene.get_transcript(trans_id)
            exon = trans.get_exon()
            cds = trans.get_CDS()
            start_codon = trans.get_start_codon()
            stop_codon = trans.get_stop_codon()

            transcript_id_gene_id[trans_id] = gene_id
            transcript_id_gene[trans_id] = gene_name
            if exon:
                dt_out[gene_name][gene_id][-1][trans_id] = [exon.get_range(), exon.get_attr("gbkey")]
            else:
                dt_out[gene_name][gene_id][-1][trans_id] = [None, None]

            if cds:
                dt_out[gene_name][gene_id][-1][trans_id].append(cds.get_range())
            else:
                dt_out[gene_name][gene_id][-1][trans_id].append(None)

            if start_codon:
                dt_out[gene_name][gene_id][-1][trans_id].append(start_codon.get_range())
            else:
                dt_out[gene_name][gene_id][-1][trans_id].append(None)
            if stop_codon:
                dt_out[gene_name][gene_id][-1][trans_id].append(stop_codon.get_range())
            else:
                dt_out[gene_name][gene_id][-1][trans_id].append(None)

    return dt_out, transcript_id_gene_id, transcript_id_gene


def get_fasta_seq_info(fasta_file):
    """
    dt_out: {seq_id: [refseq, gene, transcript_id, gbkey, left, right, ori], ...}
    """
    fasta = bp.parse_fasta(fasta_file)
    seq_id_s = fasta.get_seq_id_s()
    dt_out = {}
    retmp = re.compile(r'\[(.+?)=(.+?)\]')
    for seq_id in seq_id_s:
        seq_entry = fasta.get_seq_entry(seq_id)
        info = seq_entry.get_seq_info()
        info = {ele[0]: ele[1] for ele in retmp.findall(info)}
        refseq = "_".join(seq_id.split("|")[1].split("_")[:2])
        gene = info.get("gene")
        transcript_id = info.get("transcript_id")
        if not transcript_id:
            transcript_id = gene
        gbkey = info.get("gbkey")
        ori = "+"
        location = info.get("location")
        if location.startswith("complement"):
            location = location[11: -1]
            ori = "-"
        if location.startswith("join"):
            location = location[5: -1]
        location = location.split("..")
        left = location[0]
        right = location[-1]
        if "," in left:
            left = left.split(",")[0]
        if "," in right:
            right = right.split(",")[0]
        if left[0].isdigit():
            left = int(left)
        else:
            left = int(left[1:])
        if right[0].isdigit():
            right = int(right)
        else:
            right = int(right[1:])
        dt_out[seq_id] = [refseq, gene, transcript_id, gbkey, left, right, ori]
    return dt_out


def correct_fasta_info(fasta_info, gtf_info, transcript_id_gene_id, transcript_id_gene, out):
    dt_out = {}
    flog = open(out + "-fasta-info-correctiong.log", "w")
    for seq_id in fasta_info:
        refseq_fasta, gene_name_fasta, transcript_id_fasta, gbkey_fasta, \
            left_fasta, right_fasta, ori_fasta = fasta_info[seq_id]
        gene_name_gtf = transcript_id_gene.get(transcript_id_fasta)
        gene_id_gtf = transcript_id_gene_id.get(transcript_id_fasta)
        #assert gene_name_fasta == gene_name_gtf #note, the name may not same.
        if gene_name_gtf and gene_id_gtf:
            if gene_name_gtf in gtf_info and gene_id_gtf in gtf_info[gene_name_gtf] \
                and transcript_id_fasta in gtf_info[gene_name_gtf][gene_id_gtf][-1]:
                refseq_gtf, range_gtf, ori_gtf  = gtf_info[gene_name_gtf][gene_id_gtf][:-1]
                if (refseq_fasta == refseq_gtf and
                    ori_fasta == ori_gtf and left_fasta >= range_gtf[0][0] and
                     right_fasta <= range_gtf[0][1]):
                    
                    dt_out[seq_id] = [gene_name_gtf, gene_id_gtf, transcript_id_fasta, gbkey_fasta]
                else:
                    print(f"unknown seq (transcript_id): {refseq_fasta} {gene_name_gtf} {gene_id_gtf} {transcript_id_fasta} {gbkey_fasta}", file=flog)
            else:
                print(f"unknown seq (gene_id): {refseq_fasta} {gene_name_gtf} {gene_id_gtf} {transcript_id_fasta} {gbkey_fasta}", file=flog)
        else:
            print(f"{transcript_id_fasta} not found in transcript_id_gene_id dic", file=flog)
    print(len(dt_out), "FASTA sequence is found in GTF file", file=flog)
    flog.close()
    return dt_out


def struc_sam_file(samfile, read_type):
    fin = open(samfile, "r")
    for line in fin:
        if line[0] != "@":
            break

    dt = []
    line = line.rstrip().split('\t')
    readid, flg, refname, pos, read = line[0], line[1], line[2], line[3], line[9]
    if read_type == "pe" and (flg == "99" or flg == "355"):
        dt.append([readid, refname, int(pos), len(read)])
    elif read_type == "se" and (flg == "0" or flg == "256"):
        dt.append([readid, refname, int(pos), len(read)])
    
    if read_type == "pe":
        for line in fin:
            line = line.rstrip().split('\t')
            readid, flg, refname, pos, read = line[0], line[1], line[2], line[3], line[9]
            if flg == "99" or flg == "355":
                dt.append([readid, refname, int(pos), len(read)])
    elif read_type == "se":
        for line in fin:
            line = line.rstrip().split('\t')
            readid, flg, refname, pos, read = line[0], line[1], line[2], line[3], line[9]
            if flg == "0" or flg == "256":
                dt.append([readid, refname, int(pos), len(read)])
    
    fin.close()
    return dt


def annotate_read(samdata, fasta_info, out):
    dt_anno = {}
    flog = open(out + "-annotation.log",  "w")
    anno_count = 0
    for line in samdata:
        readid, seq_id, pos, readlen = line
        fasta_info_ext = fasta_info.get(seq_id)
        if fasta_info_ext:
            #print(">", refname)
            anno_count += 1
            gene_name_gtf, gene_id_gtf, transcript_id_fasta, gbkey_fasta = fasta_info_ext
            if gene_name_gtf in dt_anno:
                if gene_id_gtf in dt_anno[gene_name_gtf]:
                    if transcript_id_fasta in dt_anno[gene_name_gtf][gene_id_gtf]:
                        dt_anno[gene_name_gtf][gene_id_gtf][transcript_id_fasta].append([readid, pos, readlen])
                    else:
                        dt_anno[gene_name_gtf][gene_id_gtf][transcript_id_fasta] = [[readid, pos, readlen]]
                else:
                    dt_anno[gene_name_gtf][gene_id_gtf] = {transcript_id_fasta: [[readid, pos, readlen]]}
            else:
                dt_anno[gene_name_gtf] = {gene_id_gtf: {transcript_id_fasta: [[readid, pos, readlen]]}}
    print(f"{anno_count} reads were annotated.", file=flog)
    flog.close()
    return dt_anno


def change_coordinate_relative_to_gene_range(gene_range, trans_range, ori):
    range_rel = []
    gene_left = gene_range[0]
    gene_right = gene_range[1]
    trans_range.sort(key=lambda x:x[0])
    if ori == "+":
        for block in trans_range:
            range_rel.append((block[0] - gene_left + 1, block[1] - gene_left + 1))
    elif ori == "-":
        for block in trans_range[::-1]:
            range_rel.append((gene_right - block[1] + 1, gene_right - block[0] + 1))
    return range_rel


def align_read_to_transcipt(trans_rel, read_pos, read_len, flog):
    trans_block_len = [ele[1] - ele[0] + 1 for ele in trans_rel]
    trans_block_len_cum = [0]
    read_tail_len = int(read_len * .1)
    for idx, ele in enumerate(trans_block_len):
        trans_block_len_cum.append(trans_block_len_cum[idx] + ele)
    trans_block_len_cum = trans_block_len_cum[1:]
    read_start = read_pos
    read_stop = read_pos + read_len - 1
    read_start_in = 0
    for idx_start, cum_len in enumerate(trans_block_len_cum):
        if read_start <= cum_len:
            read_start_in = 1
            break
    read_stop_in = 0
    for idx_stop, cum_len in enumerate(trans_block_len_cum):
        if read_stop <= cum_len:
            read_stop_in = 1
            break
    
    if not read_start_in:
        print("read left out", file=flog)
        return None
    if not read_stop_in:
        if read_stop - trans_block_len_cum[-1] < read_tail_len:
            read_stop = trans_block_len_cum[-1]
        else:
            print("read right out", trans_rel, read_pos, read_len, file=flog)
            return None
    assert idx_start <= idx_stop

    dt_out = []
    if idx_start == idx_stop:
        dt_out.append([trans_rel[idx_start][1] - (trans_block_len_cum[idx_start] - read_start), \
                      trans_rel[idx_stop][1] - (trans_block_len_cum[idx_stop] - read_stop)])
    else:
        dt_out.append([trans_rel[idx_start][1] - (trans_block_len_cum[idx_start] - read_start), trans_rel[idx_start][1]])
        for idx in range(idx_start + 1, idx_stop):
            dt_out.append(trans_rel[idx])
        dt_out.append([trans_rel[idx_stop][0], trans_rel[idx_stop][1] - (trans_block_len_cum[idx_stop] - read_stop)])

    return dt_out


def output_result(annoed_read, gtf_info, out):
    fout = open(out + ".rda", "w")
    flog = open(out + "-reads-align.log", "w")
    for gene_name in annoed_read:
        print(f">{gene_name}", file=fout)
        gene_name_info = gtf_info[gene_name]
        for gene_id in annoed_read[gene_name]:
            gene_id_info = gene_name_info[gene_id]
            gene_id_ref, gene_id_range, gene_id_ori, gene_id_trans_id_info = gene_id_info
            gene_id_range_str = ",".join([str(ele) for ele in gene_id_range[0]])
            print('\t'.join(["&" + gene_id, gene_id_ref, gene_id_range_str, gene_id_ori]), file=fout)
            for transcript_id in annoed_read[gene_name][gene_id]:
                transcript_id_info = gene_id_trans_id_info[transcript_id]
                trans_range, trans_gbkey, protein_range, start_codon, stop_codon = transcript_id_info
                trans_range_rel = change_coordinate_relative_to_gene_range(gene_id_range[0], trans_range, gene_id_ori)
                trans_range_str = ";".join([",".join([str(e) for e in ele]) for ele in trans_range_rel])
                if start_codon:
                    start_codon_rel = change_coordinate_relative_to_gene_range(gene_id_range[0], start_codon, gene_id_ori)
                    start_codon_str = ";".join([",".join([str(e) for e in ele]) for ele in start_codon_rel])
                else:
                    start_codon_str = "*"
                if stop_codon:
                    stop_codon_rel = change_coordinate_relative_to_gene_range(gene_id_range[0], stop_codon, gene_id_ori)
                    stop_codon_str = ";".join([",".join([str(e) for e in ele]) for ele in stop_codon_rel])
                else:
                    stop_codon_str = "*"
                print('\t'.join(["$" + transcript_id, trans_gbkey, trans_range_str, start_codon_str, stop_codon_str]), file=fout)
                for read in annoed_read[gene_name][gene_id][transcript_id]:
                    read_id, read_pos, read_len = read
                    read_align = align_read_to_transcipt(trans_range_rel, read_pos, read_len, flog)
                    if read_align:
                        read_align_str = ";".join([",".join([str(e) for e in ele]) for ele in read_align])
                        print("\t".join([read_id, str(read_pos), str(read_len), read_align_str]), file=fout)
    fout.close()
    flog.close()

def main():
    samfile, out, fasta_file, gtf_file, read_type = getargs()
    print("\nStart...", os.path.basename(samfile))

    print("Reading annotation information...")
    gtf_info, transcript_id_gene_id, transcript_id_gene = \
        get_gtf_trans_cds_info(gtf_file)
    fasta_info = get_fasta_seq_info(fasta_file)
    fasta_info = correct_fasta_info(fasta_info, gtf_info, transcript_id_gene_id,
                                    transcript_id_gene, out)
    
    print("Reading sam file...")
    sam_dt = struc_sam_file(samfile, read_type)
    
    print("Annotating...")
    annoed_read = annotate_read(sam_dt, fasta_info, out)
    
    print("Outputing result...")
    output_result(annoed_read, gtf_info, out)

    print("Done.")


if __name__ == "__main__":
    main()
