#!/usr/bin/python

import json

class PARSE_REF_JSON:
    def __init__(self, json_fn, species):
        print(json_fn)
        self.ref_dic = json.load(open(json_fn))[species]

    def cgsite_anno_list(self):
        fn_s = list()
        home = self.ref_dic['HOME']
        for chrom, _f_path in self.ref_dic['CGSITE_ANNO'].iteritems():
            fn_s.append(_f_path.replace('[HOME]', home))
        return fn_s

    def region_extractor(self, anno_fn):
        a_dic = dict()
        for line in open(anno_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['#CHR']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            annos = items[idx_dic['MEMBER.{0}'.format(self.region)]].split(',')
            if '-' in annos:
                continue
            #3_10850348_10892404_genebody_ENSSSCT00000048033
            for anno in annos:
                chrom, start, end, region, annoId = anno.split('_')
                if annoId in self.target_dic:
                    a_dic.setdefault(annoId, {}).\
                          setdefault(region, {}).\
                          setdefault(start, [chrom,start,end])
        return a_dic

    def extract_region(self, target_dic, region):
        self.target_dic = target_dic
        self.region = region
        #
        region_dic = dict()
        for anno_fn in self.cgsite_anno_list():
            print('looking at {0} for {1}'.format(anno_fn, region))
            region_dic.update(self.region_extractor(anno_fn))
        return region_dic


annoTable_haeder='''\
#CHR
START
END
TYPE
NUM.UP1K
NUM.genebody
NUM.DW1K
NUM.5UTR
NUM.CDS
NUM.EXON
NUM.3UTR
NUM.Promoter
NUM.HCP
NUM.ICP
NUM.LCP
NUM.Nshelf
NUM.Nshore
NUM.CGI
NUM.Sshore
NUM.Sshelf
MEMBER.UP1K
MEMBER.genebody
MEMBER.DW1K
MEMBER.5UTR
MEMBER.CDS
MEMBER.EXON
MEMBER.3UTR
MEMBER.Promoter
MEMBER.HCP
MEMBER.ICP
MEMBER.LCP
MEMBER.Nshelf
MEMBER.Nshore
MEMBER.CGI
MEMBER.Sshore
MEMBER.Sshelf
'''
