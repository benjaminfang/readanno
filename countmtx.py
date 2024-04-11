#! /usr/bin/env python3
"""
Read in annread output, and creat count matrix file.
"""

import argparse

class ANNOTREAD_DATE:
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
                        'gene_range': (left, right),
                        'ori': string,
                        'transcipts': {
                            transcipt_id: {
                                'gbkey': string,
                                'trans_range': list,
                                'start_codon': list,
                                'stop_codon': list,
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
                trans_range = [ele.split(",") for ele in trans_range.split(";")]
                trans_range = [(int(ele[0]), int(ele[1])) for ele in trans_range]
                if start_codon != "*":
                    start_codon = [ele.split(",") for ele in start_codon.split(";")]
                    start_codon = [(int(ele[0]), int(ele[1])) for ele in start_codon]
                else:
                    start_codon = None
                if stop_codon != "*":
                    stop_codon = [ele.split(",") for ele in stop_codon.split(";")]
                    stop_codon = [(int(ele[0]), int(ele[1])) for ele in stop_codon]
                else:
                    stop_codon = None
                data[gene_name][gene_id]["transcripts"][transcipt_id] = {}
                data[gene_name][gene_id]["transcripts"][transcipt_id]["gbkey"] = gbkey
                data[gene_name][gene_id]["transcripts"][transcipt_id]["trans_range"] = trans_range
                data[gene_name][gene_id]["transcripts"][transcipt_id]["start_codon"] = start_codon
                data[gene_name][gene_id]["transcripts"][transcipt_id]["stop_codon"] = stop_codon
                data[gene_name][gene_id]["transcripts"][transcipt_id]["reads"] = []
            else:
                line = line.split('\t')
                read_id, read_pos, read_len, align_pos = line
                read_pos = int(read_pos)
                read_len = int(read_len)
                align_pos = [ele.split(",") for ele in align_pos.split(";")]
                align_pos = [(int(ele[0]), int(ele[1])) for ele in align_pos]
                data[gene_name][gene_id]["transcripts"][transcipt_id]["reads"].append([read_id, read_pos, read_len, align_pos])
        fin.close()
        self.__data = data


    def get_data(self):
        return self.__data




def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("annoread", help="annoread output file", nargs="+")
    parser.add_argument("--out", "-O", help="output file name", default="out")
    args = parser.parse_args()
    return args.annoread, args.out


if __name__ == "__main__":
    anno_files, out = getargs()

    for annofile in anno_files:
        annodata = ANNOTREAD_DATE(annofile)
