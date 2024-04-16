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

import random
import numpy as np
import os


class ANNOREAD_DATA:
    """
    The class to handle reads annotation data generate by annoread.py.
    """
    def __init__(self, anno_file):
        """
        Structure Annotation File

        Parameters
        ---------------
        anno_file: The result file generate by annoread.py

        Returns
        ------------
        dt_out: Structured data of annoread_file.
            {gene_name: 
                {gene_id: 
                    {
                        'ref_id': string,
                        'gene_range': [left, right],
                        'ori': string,
                        'transcipts': {
                            transcipt_id: {
                                'gbkey': string,
                                'transcript_range': list->[[p1, p2], ...],
                                'start_codon': list/None,
                                'stop_codon': list/None,
                                'reads': [
                                    [read_id, read_pos, read_len, alig_pos],
                                    ...
                                ]
                            }
                        }  

                    }
                }
            }
        """
        self.__data = None
        data = {}
        gene_name = None
        gene_id = None
        transcipt_id = None
        fin = open(anno_file, "r")
        for line in fin:
            line = line.rstrip()
            if line[0] == ">":
                gene_name = line[1:]
                data[gene_name] = {}
            elif line[0] == "&":
                line = line[1:].split('\t')
                gene_id, ref_id, gene_range, ori = line
                gene_range = [int(ele) for ele in gene_range.split(",")]
                data[gene_name][gene_id] = {}
                data[gene_name][gene_id]["ref_id"] = ref_id
                data[gene_name][gene_id]["gene_range"] = gene_range
                data[gene_name][gene_id]["ori"] = ori
                data[gene_name][gene_id]["transcripts"] = {}
            elif line[0] == "$":
                line = line[1:].split('\t')
                transcipt_id, gbkey, trans_range, start_codon, stop_codon = line
                trans_range = [[int(nn) for nn in ele.split(",")] for ele in trans_range.split(";")]
                if start_codon != "*":
                    start_codon = [[int(nn) for nn in ele.split(",")] for ele in start_codon.split(";")]
                else:
                    start_codon = None
                if stop_codon != "*":
                    stop_codon = [[int(nn) for nn in ele.split(",")] for ele in stop_codon.split(";")]
                else:
                    stop_codon = None
                data[gene_name][gene_id]["transcripts"][transcipt_id] = {}
                data[gene_name][gene_id]["transcripts"][transcipt_id]["gbkey"] = gbkey
                data[gene_name][gene_id]["transcripts"][transcipt_id]["transcript_range"] = trans_range
                data[gene_name][gene_id]["transcripts"][transcipt_id]["start_codon"] = start_codon
                data[gene_name][gene_id]["transcripts"][transcipt_id]["stop_codon"] = stop_codon
                data[gene_name][gene_id]["transcripts"][transcipt_id]["reads"] = []
            else:
                line = line.split('\t')
                read_id, read_pos, read_len, align_pos = line
                read_pos = int(read_pos)
                read_len = int(read_len)
                align_pos = [[int(nn) for nn in ele.split(",")] for ele in align_pos.split(";")]
                data[gene_name][gene_id]["transcripts"][transcipt_id]["reads"].append([read_id, read_pos, read_len, align_pos])
        fin.close()
        self.__data = data


    def get_data(self):
        return self.__data


    def __get_read_id_affiliation(self, set_data):
        dt_out = {}
        for key in set_data:
            for read_id in set_data[key]:
                if read_id not in dt_out:
                    dt_out[read_id] = [key]
                else:
                    dt_out[read_id].append(key)

        return dt_out


    def __set_unique_ele_number(self, owner_reads_set):
        """
        Some times most reads is common in all owers. Calculate the unique_element
        firstly, to avoid repease calculation.
        """
        owner_s = list(owner_reads_set.keys())
        owner_s_len = len(owner_s)
        owner_idx = list(range(owner_s_len))
        uniqu_ele_num = []
        for ii in range(owner_s_len):
            idx_pop = owner_idx.pop(0)
            set_pop = owner_reads_set[owner_s[idx_pop]]
            left_set = set()
            for jj in owner_idx:
                left_set |= owner_reads_set[owner_s[jj]]
            uniqu_ele_num.append([len(set_pop - left_set), idx_pop])
            owner_idx.append(idx_pop)

        dt_out = {}
        for ele_num in uniqu_ele_num:
            dt_out[owner_s[ele_num[1]]] = ele_num[0]
        return dt_out


    def __assign_reads_drop(self, read_affiliation_dic):
        """
        Read has multiple owers would dropped.
        """
        read_not_keep_dic = {}
        for read_id in read_affiliation_dic:
            if len(read_affiliation_dic[read_id]) > 1:
                for owner in read_affiliation_dic[read_id]:
                    if owner in read_not_keep_dic:
                        read_not_keep_dic[owner].append(read_id)
                    else:
                        read_not_keep_dic[owner] = [read_id]
        
        return read_not_keep_dic


    def __assign_reads_equal(self, read_affiliation_dic):
        """
        Read has multiple owers will randomly assign to one owers. 
        """
        read_not_keep_dic = {}
        for read_id in read_affiliation_dic:
            if len(read_affiliation_dic[read_id]) > 1:
                owner_s = read_affiliation_dic[read_id]
                owner_s.remove(random.choice(owner_s))
                for owner in owner_s:
                    if owner in read_not_keep_dic:
                        read_not_keep_dic[owner].append(read_id)
                    else:
                        read_not_keep_dic[owner] = [read_id]
        return read_not_keep_dic


    def __assign_reads_largest(self, read_affiliation_dic, owner_reads_set):
        read_not_keep_dic = {}
        set_unique_num = self.__set_unique_ele_number(owner_reads_set)
        all_owner_num = len(owner_reads_set)
        larget_owner = [[set_unique_num[ele], ele] for ele in set_unique_num]
        larget_owner.sort()
        larget_owner = larget_owner[-1][1]
        owner_left = list(set_unique_num.keys())
        owner_left.remove(larget_owner)

        for read_id in read_affiliation_dic:
            if len(read_affiliation_dic[read_id]) > 1:
                if len(read_affiliation_dic[read_id]) != all_owner_num:
                    owner_s = read_affiliation_dic[read_id]
                    owner_s_len = len(owner_s)
                    idx = list(range(owner_s_len))
                    unique_read_number = []
                    for count in range(owner_s_len):
                        pop_idx = idx.pop(0)
                        set_a = owner_reads_set[owner_s[pop_idx]]
                        set_left = set()
                        for ii in idx:
                            set_left |= owner_reads_set[owner_s[ii]]
                        unique_read_number.append([len(set_a - set_left), pop_idx])
                        idx.append(pop_idx)
                    unique_read_number.sort()
                    largest_set_idx = unique_read_number[-1][-1]
                    owner_s.remove(owner_s[largest_set_idx])
                    for owner in owner_s:
                        if owner in read_not_keep_dic:
                            read_not_keep_dic[owner].append(read_id)
                        else:
                            read_not_keep_dic[owner] = [read_id]
                else:
                    for owner in owner_left:
                        if owner in read_not_keep_dic:
                            read_not_keep_dic[owner].append(read_id)
                        else:
                            read_not_keep_dic[owner] = [read_id]                
        
        return read_not_keep_dic


    def __assign_reads_proportion(self, read_affiliation_dic, owner_reads_set):
        read_not_keep_dic = {}
        set_unique_num = self.__set_unique_ele_number(owner_reads_set)
        all_owner_num = len(owner_reads_set)
        all_owner = list(set_unique_num.keys())
        weights_all = [set_unique_num[ele] for ele in all_owner]
        if 0 in weights_all:
            weights_all = [ele + 1 for ele in weights_all]
        
        for read_id in read_affiliation_dic:
            if len(read_affiliation_dic[read_id]) > 1:
                owner_s = read_affiliation_dic[read_id]
                if len(owner_s) != all_owner_num:
                    owner_s_len = len(owner_s)
                    idx = list(range(owner_s_len))
                    unique_read_number = []
                    for count in range(owner_s_len):
                        pop_idx = idx.pop(0)
                        set_a = owner_reads_set[owner_s[pop_idx]]
                        set_left = set()
                        for ii in idx:
                            set_left |= owner_reads_set[owner_s[ii]]
                        unique_read_number.append([len(set_a - set_left), pop_idx])
                        idx.append(pop_idx)
                    #choose a owner with probability be positive to its uniqure_read_number.
                    weights = [0]  * owner_s_len
                    for set_size in unique_read_number:
                        weights[set_size[1]] = set_size[0]
                    if 0 in weights:
                        weights = [ele + 1 for ele in weights]
                    owner_s.remove(random.choices(owner_s, weights=weights, k=1)[0])
                else:
                    owner_s.remove(random.choices(all_owner, weights=weights_all, k=1)[0])

                for owner in owner_s:
                    if owner in read_not_keep_dic:
                        read_not_keep_dic[owner].append(read_id)
                    else:
                        read_not_keep_dic[owner] = [read_id]
        return read_not_keep_dic


    def __assign_reads(self, read_affiliation_dic, owner_reads_set, method):
        read_not_keep_dic = {}

        if method == "assign_drop":
            read_not_keep_dic = self.__assign_reads_drop(read_affiliation_dic)
        elif method == "assign_equal":
            read_not_keep_dic = self.__assign_reads_equal(read_affiliation_dic)
        elif method == "assign_largest":
            read_not_keep_dic = self.__assign_reads_largest(read_affiliation_dic, owner_reads_set)
        elif method == "assign_proportion":
            read_not_keep_dic = self.__assign_reads_proportion(read_affiliation_dic, owner_reads_set)
        else:
            print(f"{method} is not known.")

        return read_not_keep_dic


    def assign_reads_accross_gene_name(self, method="assign_drop"):
        """
        Filter reads which belong to multiple genes.

        Two strategies used to do this,
        1) drop any reads which belong to multiple genes.
        2) keep reads to one gene randomly.

        Parameters
        ---------------
        method: the method applied to filter reads, drop is the default.
            drop, drop the reads.
            equal_assign, keep reads to one gene.
            largest_assign,
            proportion_assign,

        Retures
        ---------------

        """
        data = self.__data
        owner_reads_set = {}
        for gene_name in data:
            owner_reads_set[gene_name] = set()
            for gene_id in data[gene_name]:
                for transcript_id in data[gene_name][gene_id]["transcripts"]:
                    for read in data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]:
                        read_id = read[0]
                        owner_reads_set[gene_name].add(read_id)
        read_affiliation_dic = self.__get_read_id_affiliation(owner_reads_set)
        read_not_keep_dic = self.__assign_reads(read_affiliation_dic, owner_reads_set, method)
        for gene_name in data:
            not_keep_read = read_not_keep_dic.get(gene_name)
            if not_keep_read:
                for gene_id in data[gene_name]:
                    for transcript_id in data[gene_name][gene_id]["transcripts"]:
                        reads = data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]
                        reads_filered = []
                        for read in reads:
                            if read[0] not in not_keep_read:
                                reads_filered.append(read)

                        data[gene_name][gene_id]["transcripts"][transcript_id]["reads"] = reads_filered


    def assign_reads_accross_gene_id(self, method="assign_drop"):
        """
        Parameters
        -------------
        methods: drop, random, largest_assign, proportion_assign

        Retures
        -------------

        """
        data = self.__data
        for gene_name in data:
            owner_reads_set = {}
            for gene_id in data[gene_name]:
                owner_reads_set[gene_id] = set()
                for transcript_id in data[gene_name][gene_id]["transcripts"]:
                    for read in data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]:
                        owner_reads_set[gene_id].add(read[0])
 
            assert len(owner_reads_set) != 0
            if len(owner_reads_set) == 1:
                continue

            read_affiliation_dic = self.__get_read_id_affiliation(owner_reads_set)
            read_not_keep_dic = self.__assign_reads(read_affiliation_dic, owner_reads_set, method)

            for gene_id in data[gene_name]:
                not_keep = read_not_keep_dic.get(gene_id)
                if not_keep:
                    for transcript_id in data[gene_name][gene_id]["transcripts"]:
                        reads = data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]
                        reads_filtered = []
                        for read in reads:
                            if read[0] not in not_keep:
                                reads_filtered.append(read)
                        data[gene_name][gene_id]["transcripts"][transcript_id]["reads"] = reads_filtered


    def assign_reads_accross_transcript_id(self, method="assign_drop"):
        data = self.__data
        for gene_name in data:
            for gene_id in data[gene_name]:
                owner_reads_set = {}
                for transcript_id in data[gene_name][gene_id]["transcripts"]:
                    owner_reads_set[transcript_id] = set()
                    for read in data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]:
                        owner_reads_set[transcript_id].add(read[0])

                assert len(owner_reads_set) != 0
                if len(owner_reads_set) == 1:
                    continue

                read_affiliation_dic = self.__get_read_id_affiliation(owner_reads_set)
                read_not_keep_dic = self.__assign_reads(read_affiliation_dic, owner_reads_set, method)

                for transcript_id in data[gene_name][gene_id]["transcripts"]:
                    not_keep = read_not_keep_dic.get(transcript_id)
                    if not_keep:
                        reads = data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]
                        reads_filtered = []
                        for read in reads:
                            if read[0] not in not_keep:
                                reads_filtered.append(read)
                        data[gene_name][gene_id]["transcripts"][transcript_id]["reads"] = reads_filtered


    def unique_reads_one_transcript(self):
        data = self.__data
        for gene_name in data:
            for gene_id in data[gene_name]:
                for transcript_id in data[gene_name][gene_id]["transcripts"]:
                    read_id_s = []
                    read_filtered = []
                    for read in data[gene_name][gene_id]["transcripts"][transcript_id]["reads"]:
                        if read[0] not in read_id_s:
                            read_id_s.append(read[0])
                            read_filtered.append(read)
                    data[gene_name][gene_id]["transcripts"][transcript_id]["reads"] = read_filtered


def count_reads(annoread_dt):
    gene_count_data = {}
    annoread_dt.assign_reads_accross_gene_name(method="assign_drop")
    print("    reads filtered by assign_drop method")
    data = annoread_dt.get_data()
    for gene_name in data:
        gene_count_data[gene_name] = set()
        for gene_id in data[gene_name]:
            for transcipt_id in data[gene_name][gene_id]["transcripts"]:
                for read in data[gene_name][gene_id]["transcripts"][transcipt_id]["reads"]:
                    gene_count_data[gene_name].add(read[0])

    for gene_name in gene_count_data:
        gene_count_data[gene_name] = len(gene_count_data[gene_name])
    return gene_count_data


def output_count_mtx(count_data, anno_files, out):
    all_gene_names = set()
    for count_data_ele in count_data:
        set_gene = set()
        for key in count_data_ele:
            if count_data_ele[key] > 0:
                set_gene.add(key)
        all_gene_names |= set_gene
    all_gene_names = list(all_gene_names)
    all_gene_names.sort()
    mtx = []
    for count_data_ele in count_data:
        one_file_dt = []
        for gene in all_gene_names:
            one_file_dt.append(count_data_ele.get(gene, 0))
        mtx.append(one_file_dt)
    mtx = np.asarray(mtx).T
    fout = open(out + '.tsv', "w")
    print("\t".join(["geneid"] + [os.path.basename(ff) for ff in anno_files]), file=fout)
    for idx, row in enumerate(mtx):
        print("\t".join([all_gene_names[idx]] + [str(ele) for ele in row]), file=fout)
    fout.close()


if __name__ == "__main__":
    import argparse
    def getargs():
        parser = argparse.ArgumentParser()
        parser.add_argument("anno_files", help="annoread output file", nargs="+")
        parser.add_argument("--out", "-O", help="output file name", default="out")
        args = parser.parse_args()
        return args.anno_files, args.out

    anno_files, out = getargs()
    count_data = []
    for annofile in anno_files:
        print("reading...", annofile)
        annodata = ANNOREAD_DATA(annofile)
        print("    counting...")
        count_data.append(count_reads(annodata))
    print("writing result...")
    output_count_mtx(count_data, anno_files, out)
