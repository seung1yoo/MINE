#!/usr/bin/python

import os
import sys

paths = os.environ["PATH"].split(":")
paths.insert(0, "/BiO/BioPeople/siyoo/MINE/bin")
sys.path.extend(paths)
print(paths)

#py_paths = os.environ["PYTHONPATH"].split(":")
#py_paths.insert(0, "/BiO/BioPeople/siyoo/MINE/bin")
#os.environ["PYTHONPATH"] = ":".join(py_paths)
#print(py_paths)

from lib.MINE_LIB import PARSE_REF_JSON

def gene_list2dic(fn):
    gene_dic = dict()
    for line in open(fn):
        items = line.rstrip('\n').split('\t')
        g_symbol = items[0]
        g_ens = items[1]
        tr_ens = items[2]
        gene_dic.setdefault(tr_ens, [g_ens, g_symbol])
    return gene_dic

def print_gene_dic(gene_dic):
    for tr_ens, infos in gene_dic.iteritems():
        print(tr_ens, infos)


def main(args):
    gene_dic = gene_list2dic(args.gene_list)
    print_gene_dic(gene_dic)
    #
    ref_json_parser = PARSE_REF_JSON(args.ref_json, args.species)

    # Case.1
    region_dic = dict()
    for region in ['UP1K','genebody','DW1K']:
        region_dic.update(ref_json_parser.extract_region(gene_dic, region))
    bed_dic = region_dic2bed(region_dic)
    print region_dic

    # Case.2
    region_dic = dict()
    for region in ['5UTR','CDS','3UTR']: # EXON
        region_dic.update(ref_json_parser.extract_region(gene_dic, region))
        break
    bed_dic = region_dic2bed(region_dic)
    print region_dic

    # Case.3
    region_dic = dict()
    for region in ['Promoter','HCP','ICP','LCP']:
        region_dic.update(ref_json_parser.extract_region(gene_dic, region))
        break
    bed_dic = region_dic2bed(region_dic)
    print region_dic

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('gene_list')
    parser.add_argument('ref_json')
    parser.add_argument('species')
    parser.add_argument('--temp_dir', default='./tmp')
    args = parser.parse_args()
    main(args)
