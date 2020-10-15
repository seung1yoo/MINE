#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE
#
import shutil

class Do_report_sample:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))

        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco

        # variable setting
        #
        self.outdir = os.path.join(self.eco.room['report'], 'samples')
        self.eco.make_dir(self.outdir)
        #
        self.sumup_contents = ['N_READ_RAW',
                               'N_READ_CLEAN','R_READ_CLEAN',
                               'N_READ_MAP','R_READ_MAP',
                               'N_READ_UNIQUE','R_READ_UNIQUE',
                               'NT_C_TOTAL',
                               'NT_CpG_MET','R_CpG_MET',
                               'NT_CpH_MET','R_CpH_MET',
                               'NT_CHH_MET','R_CHH_MET']
        self.sumup_stat_fn = os.path.join(self.outdir,'Sample_Stat.xls')
        #
        self.seq_contents = ['N_BASE','N_BASE_Q30','N_BASE_Q20',
                             'N_READ','N_READ_Q30','N_READ_Q20',
                             'NT_GC','NT_AT','NT_N']
                             #'LEN_MAX','LEN_MIN','LEN_AVG']
        self.seq_stat_fn = os.path.join(self.outdir,'Sequence_Stat.xls')
        #
        self.align_contents = ['N_READ_IN','N_READ_MAP',
                               'N_READ_MAPLQ', 'N_READ_DUP','N_READ_UNIQUE',
                               'NT_C_TOTAL','NT_CpG_MET','NT_CpG_UNMET',
                               'NT_CpH_MET','NT_CpH_UNMET',
                               'NT_CHH_MET','NT_CHH_UNMET']
        self.align_stat_fn = os.path.join(self.outdir,'Alignment_Stat.xls')


        # check_stats & run
        if not self.eco.check_stats(self.name):
            seqstat_inputs = ['do_stat_fastq_raw','do_stat_fastq_clean']
            seq_input_dic = self.select_input(seqstat_inputs)
            self.seqstat_dic = self.seqstat_sum_up(seq_input_dic)
            #
            alignstat_inputs = ['do_bismark_report']
            align_input_dic = self.select_input(alignstat_inputs)
            self.alignstat_dic = self.alignstat_sum_up(align_input_dic)
            #
            self.write_sumup_stat_f()
            self.write_seq_stat_f()
            self.write_align_stat_f()
        else:
            pass

        # finising
        self.update_ecosystem()
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem(self):
        lvs = [self.name, 'all', 'stat.sumup', self.sumup_stat_fn]
        self.eco.add_to_ecosystem(lvs)
        lvs = [self.name, 'all', 'stat.seq', self.seq_stat_fn]
        self.eco.add_to_ecosystem(lvs)
        lvs = [self.name, 'all', 'stat.align', self.align_stat_fn]
        self.eco.add_to_ecosystem(lvs)

    def write_align_stat_f(self):
        out = open(self.align_stat_fn,'w')
        headers = ['SAM_NUM','SAM_NAME']
        headers.extend(self.align_contents)
        out.write('{0}\n'.format('\t'.join(headers)))
        for sam_num, info_dic in sorted(self.alignstat_dic.iteritems()):
            new_items = [sam_num]
            sam_name = self.mine.sample_dic[sam_num]['1']['label']
            new_items.append(sam_name)
            #
            for content in self.align_contents:
                new_items.append(self.comma_into_int(info_dic[content]))
            out.write('{0}\n'.format('\t'.join(new_items)))
        out.close()

    def write_seq_stat_f(self):
        out = open(self.seq_stat_fn,'w')
        headers = ['SAM_NUM','SAM_NAME','STATS']
        headers.extend(self.seq_contents)
        out.write('{0}\n'.format('\t'.join(headers)))
        for sam_num, info_dic in sorted(self.seqstat_dic.iteritems()):
            for _type in ['raw.stat','clean.stat']:
                new_items = [sam_num]
                sam_name = self.mine.sample_dic[sam_num]['1']['label']
                new_items.append(sam_name)
                if _type in ['raw.stat']:
                    new_items.append('RAW')
                elif _type in ['clean.stat']:
                    new_items.append('CLEAN')
                #
                for content in self.seq_contents:
                    if content in ['LEN_AVG']:
                        new_items.append(self.comma_into_float(info_dic[content][_type]))
                    else:
                        new_items.append(self.comma_into_int(info_dic[content][_type]))
                out.write('{0}\n'.format('\t'.join(new_items)))
        out.close()

    def write_sumup_stat_f(self):
        out = open(self.sumup_stat_fn,'w')
        headers = ['SAM_NUM','SAM_NAME']
        headers.extend(self.sumup_contents)
        out.write('{0}\n'.format('\t'.join(headers)))
        for sam_num, info_dic in sorted(self.seqstat_dic.iteritems()):
            new_items = [sam_num]
            sam_name = self.mine.sample_dic[sam_num]['1']['label']
            new_items.append(sam_name)
            #
            for content in self.sumup_contents:
                if content in ['N_READ_RAW']:
                    new_items.append(self.comma_into_int(info_dic['N_READ']['raw.stat']))
                elif content in ['N_READ_CLEAN']:
                    new_items.append(self.comma_into_int(info_dic['N_READ']['clean.stat']))
                elif content in ['R_READ_CLEAN']:
                    numerator = info_dic['N_READ']['clean.stat']
                    denominator = info_dic['N_READ']['raw.stat']
                    new_items.append(self.cal_ratio(numerator, denominator))
                elif content in ['N_READ_MAP']:
                    new_items.append(self.comma_into_int(self.alignstat_dic[sam_num][content]))
                elif content in ['R_READ_MAP']:
                    numerator = self.alignstat_dic[sam_num]['N_READ_MAP']
                    denominator = info_dic['N_READ']['clean.stat']
                    new_items.append(self.cal_ratio(numerator, denominator))
                elif content in ['N_READ_UNIQUE']:
                    new_items.append(self.comma_into_int(self.alignstat_dic[sam_num][content]))
                elif content in ['R_READ_UNIQUE']:
                    numerator = self.alignstat_dic[sam_num]['N_READ_UNIQUE']
                    denominator = info_dic['N_READ']['clean.stat']
                    new_items.append(self.cal_ratio(numerator, denominator))
                elif content in ['NT_C_TOTAL']:
                    new_items.append(self.comma_into_int(self.alignstat_dic[sam_num][content]))
                elif content in ['NT_CpG_MET']:
                    new_items.append(self.comma_into_int(self.alignstat_dic[sam_num][content]))
                elif content in ['R_CpG_MET']:
                    numerator = self.alignstat_dic[sam_num]['NT_CpG_MET']
                    denominator = int(self.alignstat_dic[sam_num]['NT_CpG_MET'])+\
                                  int(self.alignstat_dic[sam_num]['NT_CpG_UNMET'])
                    new_items.append(self.cal_ratio(numerator, denominator))
                elif content in ['NT_CpH_MET']:
                    new_items.append(self.comma_into_int(self.alignstat_dic[sam_num][content]))
                elif content in ['R_CpH_MET']:
                    numerator = self.alignstat_dic[sam_num]['NT_CpH_MET']
                    denominator = int(self.alignstat_dic[sam_num]['NT_CpH_MET'])+\
                                  int(self.alignstat_dic[sam_num]['NT_CpH_UNMET'])
                    new_items.append(self.cal_ratio(numerator, denominator))
                elif content in ['NT_CHH_MET']:
                    new_items.append(self.comma_into_int(self.alignstat_dic[sam_num][content]))
                elif content in ['R_CHH_MET']:
                    numerator = self.alignstat_dic[sam_num]['NT_CHH_MET']
                    denominator = int(self.alignstat_dic[sam_num]['NT_CHH_MET'])+\
                                  int(self.alignstat_dic[sam_num]['NT_CHH_UNMET'])
                    new_items.append(self.cal_ratio(numerator, denominator))
                else:
                    print 'Unknown content {0}'.format(content)
            out.write('{0}\n'.format('\t'.join(new_items)))
        out.close()

    def comma_into_int(self, num):
        return '{0:,}'.format(int(num))

    def comma_into_float(self, num):
        return '{0:,}'.format(round(float(num),2))

    def cal_ratio(self, numerator, denominator):
        ratio = (float(numerator)/float(denominator))*100
        return str(round(ratio,2))

    def alignstat_sum_up(self, align_input_dic):
        alignstat_dic = dict()
        align_fn = align_input_dic['all']['summary.txt']
        for line in open(align_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['File']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            sam_num = self.convert_label2num(items[idx_dic['File']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('N_READ_IN',
                    items[idx_dic['Total Reads']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('N_READ_MAP',
                    items[idx_dic['Aligned Reads']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('N_READ_UNMAP',
                    items[idx_dic['Unaligned Reads']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('N_READ_MAPLQ',
                    items[idx_dic['Ambiguously Aligned Reads']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('N_READ_DUP',
                    items[idx_dic['Duplicate Reads (removed)']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('N_READ_UNIQUE',
                    items[idx_dic['Unique Reads (remaining)']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('NT_C_TOTAL',
                    items[idx_dic['Total Cs']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('NT_CpG_MET',
                    items[idx_dic['Methylated CpGs']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('NT_CpG_UNMET',
                    items[idx_dic['Unmethylated CpGs']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('NT_CpH_MET',
                    items[idx_dic['Methylated CpHs']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('NT_CpH_UNMET',
                    items[idx_dic['Unmethylated CpHs']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('NT_CHH_MET',
                    items[idx_dic['Methylated CHHs']])
            alignstat_dic.setdefault(sam_num, {}).setdefault('NT_CHH_UNMET',
                    items[idx_dic['Unmethylated CHHs']])
        return alignstat_dic

    def convert_label2num(self, _label):
        label = _label.split('_bismark')[0]
        for num, pair_dic in self.mine.sample_dic.iteritems():
            if label in [pair_dic['1']['label']]:
                return num
            else:
                continue
        return False

    def seqstat_sum_up(self, input_dic):
        seqstat_dic = dict()
        for sam_num, type_dic in input_dic.iteritems():
            for _type, seqstat_fn in type_dic.iteritems():
                seqstat_dic.setdefault(sam_num, {})
                seqstat_dic[sam_num] = self.parse_seqstat_fn(seqstat_fn, _type, seqstat_dic[sam_num])
        return seqstat_dic

    def parse_seqstat_fn(self, seqstat_fn, _type, _seqstat_dic):
        n_base = 0
        n_base_q30 = 0
        n_base_q20 = 0
        n_read = 0
        n_read_q30 = 0
        n_read_q20 = 0
        nt_gc = 0
        nt_at = 0
        nt_n = 0
        len_max = 0
        len_min = 0
        len_avg = 0
        #
        for line in open(seqstat_fn):
            items = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            n_base = items[idx_dic['totalLength']]
            n_base_q30 = str(int(items[dix_dic['totalQ30BaseCntR1']]) + \
                             int(items[dix_dic['totalQ30BaseCntR2']]))
            n_base_q20 = str(int(items[dix_dic['totalQ20BaseCntR1']]) + \
                             int(items[dix_dic['totalQ20BaseCntR2']]))
            n_read = items[idx_dic['#totalReadCnt']]
            n_read_q30 = str(int(items[dix_dic['totalQ30ReadCntR1']]) + \
                             int(items[dix_dic['totalQ30ReadCntR2']]))
            n_read_q20 = str(int(items[dix_dic['totalQ20ReadCntR1']]) + \
                             int(items[dix_dic['totalQ20ReadCntR2']]))
            nt_gc = items[idx_dic['totalGCCnt']]
            nt_at = str(int(items[idx_dic['totalLength']]) - \
                        int(items[idx_dic['totalGCCnt']]) - \
                        int(items[idx_dic['totalNCnt']]))
            nt_n = items[idx_dic['totalNCnt']]
        #
        _seqstat_dic.setdefault('N_BASE', {}).setdefault(_type, n_base)
        _seqstat_dic.setdefault('N_BASE_Q30', {}).setdefault(_type, n_base_q30)
        _seqstat_dic.setdefault('N_BASE_Q20', {}).setdefault(_type, n_base_q20)
        _seqstat_dic.setdefault('N_READ', {}).setdefault(_type, n_read)
        _seqstat_dic.setdefault('N_READ_Q30', {}).setdefault(_type, n_read_q30)
        _seqstat_dic.setdefault('N_READ_Q20', {}).setdefault(_type, n_read_q20)
        _seqstat_dic.setdefault('NT_GC', {}).setdefault(_type, nt_gc)
        _seqstat_dic.setdefault('NT_AT', {}).setdefault(_type, nt_at)
        _seqstat_dic.setdefault('NT_N', {}).setdefault(_type, nt_n)
        _seqstat_dic.setdefault('LEN_MAX', {}).setdefault(_type, len_max)
        _seqstat_dic.setdefault('LEN_MIN', {}).setdefault(_type, len_min)
        _seqstat_dic.setdefault('LEN_AVG', {}).setdefault(_type, len_avg)
        #
        return _seqstat_dic

    def iter_seqstat_fh(self, fh):
        lines = list()
        for line in fh:
            if line.startswith('#'):
                yield lines
                lines = list()
                lines.append(line.rstrip('\n'))
            else:
                lines.append(line.rstrip('\n'))
        yield lines

    def _parse_seqstat_fn(self, seqstat_fn, _type, _seqstat_dic):
        n_base = 0
        n_base_q30 = 0
        n_base_q20 = 0
        n_read = 0
        n_read_q30 = 0
        n_read_q20 = 0
        nt_gc = 0
        nt_at = 0
        nt_n = 0
        len_max = 0
        len_min = 0
        len_avg = 0
        #
        for lines in self.iter_seqstat_fh(open(seqstat_fn)):
            if lines[0].startswith('FL'):
                continue
            if lines[0].startswith('# N_BASE'):
                items = lines[1].split('\t')
                n_base = items[1]
                n_base_q30 = items[2]
                n_base_q20 = items[3]
            elif lines[0].startswith('# N_READ'):
                items = lines[1].split('\t')
                n_read = items[1]
                n_read_q30 = items[2]
                n_read_q20 = items[3]
            elif lines[0].startswith('# GC-contents:'):
                items = lines[1].split('\t')
                nt_gc = items[1]
            elif lines[0].startswith('# AT-contents:'):
                items = lines[1].split('\t')
                nt_at = items[1]
            elif lines[0].startswith("# 'N'%"):
                items = lines[1].split('\t')
                nt_n = items[1]
            elif lines[0].startswith('# LEN_TOTAL'):
                items = lines[1].split('\t')
                len_max = items[2]
                len_min = items[3]
                len_avg = items[4]
            else:
                print 'unknown line in {0}'.format(seqstat_fn)
        #
        _seqstat_dic.setdefault('N_BASE', {}).setdefault(_type, n_base)
        _seqstat_dic.setdefault('N_BASE_Q30', {}).setdefault(_type, n_base_q30)
        _seqstat_dic.setdefault('N_BASE_Q20', {}).setdefault(_type, n_base_q20)
        _seqstat_dic.setdefault('N_READ', {}).setdefault(_type, n_read)
        _seqstat_dic.setdefault('N_READ_Q30', {}).setdefault(_type, n_read_q30)
        _seqstat_dic.setdefault('N_READ_Q20', {}).setdefault(_type, n_read_q20)
        _seqstat_dic.setdefault('NT_GC', {}).setdefault(_type, nt_gc)
        _seqstat_dic.setdefault('NT_AT', {}).setdefault(_type, nt_at)
        _seqstat_dic.setdefault('NT_N', {}).setdefault(_type, nt_n)
        _seqstat_dic.setdefault('LEN_MAX', {}).setdefault(_type, len_max)
        _seqstat_dic.setdefault('LEN_MIN', {}).setdefault(_type, len_min)
        _seqstat_dic.setdefault('LEN_AVG', {}).setdefault(_type, len_avg)
        #
        return _seqstat_dic

    def select_input(self,names):
        input_dic = dict()
        for name in names:
            for sam_num, _dic in self.eco.eco_dic[name].iteritems():
                input_dic.setdefault(sam_num, {})
                input_dic[sam_num].update(_dic)
        return input_dic


class Do_report_dmr:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))

        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco

        # variable setting
        #
        self.outdir = os.path.join(self.eco.room['report'], 'DMR')
        self.eco.make_dir(self.outdir)
        #
        self.dmr_stat_contents = ['DMR_NUM','DMR_G1','DMR_G2',
                                  'DMR_CUTOFF',
                                  'N_CANDIDATE','N_DMR','N_DMR_UP','N_DMR_DOWN']
        self.dmr_stat_fn = os.path.join(self.outdir, 'DMR_Stat.xls')
        #
        self.dmr_tags = ['CHR','START','END','q','diff','n.CpG',
                         'p.MWU','p.2D_KS','mean.G1','mean.G2',
                         'DMR.YN.{0}'.format(self.mine.dmr_cut),
                         'DMR.UPDOWN']
        self.anno_tags = ['UP1K','genebody','DW1K',
                          '5UTR','CDS','EXON','3UTR',
                          'Promoter','HCP','ICP','LCP',
                          'Nshelf','Nshore','CGI','Sshore','Sshelf']
        #
        econames = ['do_metilene_anno']
        input_dic = self.select_input(econames)

        # check_stats & run
        if not self.eco.check_stats(self.name):
            dmr_rpt_dic = self.copy_to_report(input_dic)
            self.cmd_stat_dic = self.cal_dmr_stat(input_dic)
            self.write_dmr_stat()
            #
            anno2dmr_f_dic = self.convert_gene_to_dmr(dmr_rpt_dic)
        else:
            dmr_rpt_dic = self.dmr2anno_dic_maker(input_dic)
            anno2dmr_f_dic = self.anno2dmr_dic_maker(dmr_rpt_dic)
            pass

        # finising
        self.update_ecosystem(dmr_rpt_dic, anno2dmr_f_dic)
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem(self, dmr_rpt_dic, anno2dmr_f_dic):
        lvs = [self.name, 'all', 'stat.dmr', self.dmr_stat_fn]
        self.eco.add_to_ecosystem(lvs)
        #
        for dmr_num, dmr_rpt_fn in dmr_rpt_dic.iteritems():
            lvs = [self.name, dmr_num, 'dmr2anno', dmr_rpt_fn]
            self.eco.add_to_ecosystem(lvs)
        for dmr_num, anno2dmr_fn in anno2dmr_f_dic.iteritems():
            lvs = [self.name, dmr_num, 'anno2dmr', anno2dmr_fn]
            self.eco.add_to_ecosystem(lvs)
        #

    def anno2dmr_dic_maker(self, dmr_rpt_dic):
        anno2dmr_f_dic = dict()
        for dmr_num, dmr_rpt_fn in dmr_rpt_dic.iteritems():
            gene2dmr_fn = dmr_rpt_fn.replace('.xls','.anno2dmr.xls')
            anno2dmr_f_dic.setdefault(dmr_num, gene2dmr_fn)
        return anno2dmr_f_dic

    def dmr2anno_dic_maker(self, input_dic):
        dmr_rpt_dic = dict()
        for dmr_num, f_type_dic in input_dic.iteritems():
            dmr_dst_fn = self.make_dmr_report_fn(dmr_num)
            dmr_rpt_dic.setdefault(dmr_num, dmr_dst_fn)
        return dmr_rpt_dic

    def convert_gene_to_dmr(self, dmr_rpt_dic):
        anno2dmr_f_dic = dict()
        for dmr_num, dmr_rpt_fn in dmr_rpt_dic.iteritems():
            gene2dmr_dic, dmrInfo_dic = self.gene2dmr_converter(dmr_rpt_fn)
            #
            gene2dmr_fn = dmr_rpt_fn.replace('.xls','.anno2dmr.xls')
            anno2dmr_f_dic.setdefault(dmr_num, gene2dmr_fn)
            #
            out = open(gene2dmr_fn,'w')
            out_headers = ['#ANNO_TAG','ANNO_ID','ANNO_CHR','ANNO_START','ANNO_END','DMR_ID']
            out_headers.extend(self.dmr_tags)
            out.write('{0}\n'.format('\t'.join(out_headers)))
            for anno_tag, anno_id_dic in gene2dmr_dic.iteritems():
                for anno_id, info_dic in anno_id_dic.iteritems():
                    anno_items = [anno_tag, anno_id]
                    anno_items.append(info_dic['chr'])
                    anno_items.append(info_dic['start'])
                    anno_items.append(info_dic['end'])
                    for dmr_id in info_dic['dmrs']:
                        dmr_items = [dmr_id]
                        for dmr_tag in self.dmr_tags:
                            dmr_items.append(dmrInfo_dic[dmr_id][dmr_tag])
                        out.write('{0}\t{1}\n'.format('\t'.join(anno_items), '\t'.join(dmr_items)))
                        #
                    #
                #
            out.close()
        return anno2dmr_f_dic

    def gene2dmr_converter(self, dmr_rpt_fn):
        gene2dmr_dic = dict()
        dmrInfo_dic = dict()
        for line in open(dmr_rpt_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['#DMR_ID']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            #
            yn = items[idx_dic['DMR.YN.{0}'.format(self.mine.dmr_cut)]]
            if yn in ['N']:
                continue
            #
            dmr_id = items[idx_dic['#DMR_ID']]
            #
            for anno_tag in self.anno_tags:
                members = items[idx_dic['MEMBER.{0}'.format(anno_tag)]].split(',')
                for member in members:
                    if not member:
                        continue
                    (_chr, _start, _end, _anno_tag, _id) = self.member_parser(member, anno_tag)
                    gene2dmr_dic.setdefault(anno_tag, {}).setdefault(_id, {}).setdefault('chr',_chr)
                    gene2dmr_dic.setdefault(anno_tag, {}).setdefault(_id, {}).setdefault('start',_start)
                    gene2dmr_dic.setdefault(anno_tag, {}).setdefault(_id, {}).setdefault('end',_end)
                    gene2dmr_dic.setdefault(anno_tag, {}).setdefault(_id, {}).setdefault('dmrs',{}).setdefault(dmr_id,None)
            #
            for dmr_tag in self.dmr_tags:
                dmrInfo_dic.setdefault(dmr_id, {}).setdefault(dmr_tag, items[idx_dic[dmr_tag]])
            #
        return gene2dmr_dic, dmrInfo_dic

    def member_parser(self, member, mode):
        if mode in ['CGI']:
            _chr, _start, _end, _anno_tag, __id = member.split('_')
            _id = __id.split('|')[0]
        else:
            _chr, _start, _end, _anno_tag, _id = member.split('_')
        return _chr, _start, _end, _anno_tag, _id






    def write_dmr_stat(self):
        out = open(self.dmr_stat_fn,'w')
        out.write('{0}\n'.format('\t'.join(self.dmr_stat_contents)))
        for dmr_num, yn_dic in self.cmd_stat_dic.iteritems():
            new_items = list()
            new_items.append(dmr_num)
            #
            g1_sam_s = self.mine.dmr_dic[dmr_num]['G1']
            g1_name_s = [self.mine.idmatch_dic[x] for x in g1_sam_s]
            new_items.append(','.join(g1_name_s))
            #
            g2_sam_s = self.mine.dmr_dic[dmr_num]['G2']
            g2_name_s = [self.mine.idmatch_dic[x] for x in g2_sam_s]
            new_items.append(','.join(g2_name_s))
            #
            new_items.append(self.mine.dmr_cut)
            #
            n_total = self.cal_n_total(yn_dic)
            new_items.append(n_total)
            #
            n_yes = self.cal_n_yes(yn_dic)
            new_items.append(n_yes)
            #
            n_up = self.cal_n_up(yn_dic)
            new_items.append(n_up)
            #
            n_down = self.cal_n_down(yn_dic)
            new_items.append(n_down)
            #
            out.write('{0}\n'.format('\t'.join([str(x) for x in new_items])))
        out.close()

    def cal_n_total(self, yn_dic):
        n_total = 0
        n_total += yn_dic['Y']['UP']
        n_total += yn_dic['Y']['DOWN']
        if yn_dic['Y'].has_key('FLAT'):
            n_total += yn_dic['Y']['FLAT']
        n_total += yn_dic['N']['UP']
        n_total += yn_dic['N']['DOWN']
        if yn_dic['N'].has_key('FLAT'):
            n_total += yn_dic['N']['FLAT']
        return n_total

    def cal_n_yes(self, yn_dic):
        n_yes = 0
        n_yes += yn_dic['Y']['UP']
        n_yes += yn_dic['Y']['DOWN']
        if yn_dic['Y'].has_key('FLAT'):
            n_yes += yn_dic['Y']['FLAT']
        return n_yes

    def cal_n_up(self, yn_dic):
        n_up = 0
        n_up += yn_dic['Y']['UP']
        return n_up

    def cal_n_down(self, yn_dic):
        n_down = 0
        n_down += yn_dic['Y']['DOWN']
        return n_down

    def cal_dmr_stat(self, input_dic):
        dmr_stat_dic = dict()
        for dmr_num, f_type_dic in input_dic.iteritems():
            dmr_report_fn = self.make_dmr_report_fn(dmr_num)
            dmr_stat_dic.setdefault(dmr_num, {})
            dmr_stat_dic[dmr_num].update(self.parse_dmr_report(dmr_report_fn))
        return dmr_stat_dic

    def parse_dmr_report(self, fn):
        _dic = dict()
        for line in open(fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['#DMR_ID']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            dmr_yn = items[idx_dic['DMR.YN.{0}'.format(self.mine.dmr_cut)]]
            dmr_ud = items[idx_dic['DMR.UPDOWN']]
            #
            _dic.setdefault(dmr_yn,{}).setdefault(dmr_ud,0)
            _dic[dmr_yn][dmr_ud] += 1
        return _dic



    def copy_to_report(self, input_dic):
        dmr_rpt_dic = dict()
        for dmr_num, f_type_dic in input_dic.iteritems():
            dmr_src_fn = f_type_dic['dmr_anno']
            dmr_dst_fn = self.make_dmr_report_fn(dmr_num)
            dmr_rpt_dic.setdefault(dmr_num, dmr_dst_fn)
            if not os.path.exists(dmr_dst_fn):
                shutil.copyfile(dmr_src_fn, dmr_dst_fn)
            else:
                pass
        return dmr_rpt_dic

    def make_dmr_report_fn(self, dmr_num):
        g1_sam_s = self.mine.dmr_dic[dmr_num]['G1']
        g1_name_s = [self.mine.idmatch_dic[x] for x in g1_sam_s]
        #
        g2_sam_s = self.mine.dmr_dic[dmr_num]['G2']
        g2_name_s = [self.mine.idmatch_dic[x] for x in g2_sam_s]
        #
        dmr_report_fn = '{0}-{1}-{2}.xls'.format(dmr_num,
                            '_'.join(g1_name_s),'_'.join(g2_name_s))
        dmr_report_path = os.path.join(self.outdir, dmr_report_fn)
        return dmr_report_path

    def select_input(self,names):
        input_dic = dict()
        for name in names:
            for sam_num, _dic in self.eco.eco_dic[name].iteritems():
                input_dic.setdefault(sam_num, {})
                input_dic[sam_num].update(_dic)
        return input_dic

class Do_report_ref:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))

        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco

        # variable setting
        #
        self.outdir = os.path.join(self.eco.room['report'], 'files')
        self.eco.make_dir(self.outdir)
        #
        self.ref_info_fn = os.path.join(self.outdir, 'Info_Reference.txt')

        # check_stats & run
        if not self.eco.check_stats(self.name):
            self.write_ref_info_f()
            self.copy_promoter_stat()
            self.copy_cpgisland_stat()
            self.copy_cgsite_stat()
        else:
            pass

        # finising
        self.update_ecosystem()
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)
        #

    def update_ecosystem(self):
        lvs = [self.name, 'all', 'ref.info', self.ref_info_fn]
        self.eco.add_to_ecosystem(lvs)

    def copy_cgsite_stat(self):
        src_dic = self.mine.refs['CGSITE_PNG']
        for chr_num, fn in src_dic.iteritems():
            src_fn = fn.replace('[HOME]',self.mine.refs['HOME'])
            outdir = os.path.join(self.outdir, 'CG_site_distribution')
            self.eco.make_dir(outdir)
            dst = os.path.join(outdir, 'CG.site.dist.{0}.png'.format(chr_num))
            if not os.path.exists(dst):
                shutil.copyfile(src_fn, dst)

    def copy_promoter_stat(self):
        src = self.mine.refs['PROMTER_STATS_FN']
        dst = os.path.join(self.outdir, 'Promoter_Annotation_Stat.xls')
        if not os.path.exists(dst):
            shutil.copyfile(src, dst)

    def copy_cpgisland_stat(self):
        src = self.mine.refs['CGI_GC_PNG']
        dst = os.path.join(self.outdir, 'CpGisland_GCpct_Disribution.png')
        if not os.path.exists(dst):
            shutil.copyfile(src, dst)
            shutil.copyfile(src.replace('.png','.xls'),dst.replace('.png','.xls'))
        #
        src = self.mine.refs['CGI_LEN_PNG']
        dst = os.path.join(self.outdir, 'CpGisland_Length_Disribution.png')
        if not os.path.exists(dst):
            shutil.copyfile(src, dst)
            shutil.copyfile(src.replace('.png','.xls'),dst.replace('.png','.xls'))
        #
        src = self.mine.refs['CGI_OBSEXP_PNG']
        dst = os.path.join(self.outdir, 'CpGisland_ObsExp_Disribution.png')
        if not os.path.exists(dst):
            shutil.copyfile(src, dst)
            shutil.copyfile(src.replace('.png','.xls'),dst.replace('.png','.xls'))

    def write_ref_info_f(self):
        out = open(self.ref_info_fn,'w')
        out.write('{0}\n'.format('='.join(['SPECEIS_ID',self.mine.species])))
        out.write('{0}\n'.format('='.join(['SPECIES_NAME',self.mine.refs['SPECIES_NAME']])))
        out.write('{0}\n'.format('='.join(['SPECIES_ALIAS',self.mine.refs['SPECIES_ALIAS']])))
        out.write('{0}\n'.format('='.join(['ASSEMBLY_ACC',self.mine.refs['ASSEMBLY_ACC']])))
        out.write('{0}\n'.format('='.join(['GENESET_VER',self.mine.refs['GENESET_VER']])))
        out.write('\n')
        #
        for line in open(self.mine.refs['CHR_LEN_FN']):
            items = line.rstrip('\n').split('\t')
            _id = items[0]
            _len = int(items[1])
            out.write('{0}\n'.format('='.join(['CHR_LEN:{0}'.format(_id),'{0:,}'.format(_len)])))
        #
        out.close()

    def _write_ref_info_f(self):
        out = open(self.ref_info_fn,'w')
        #
        out.write('\n')
        out.write('{0}\n'.format('#'*50))
        out.write('#{0:^48}#\n'.format('Methylation Analysis Report'))
        out.write('#{0:^48}#\n'.format('Reference Informations'))
        out.write('{0}\n'.format('#'*50))
        out.write('\n')
        #
        #
        out.write('{0}\n'.format('='.join(['SPECEIS_ID',self.mine.species])))
        out.write('{0}\n'.format('='.join(['SPECIES_NAME',self.mine.refs['SPECIES_NAME']])))
        out.write('{0}\n'.format('='.join(['SPECIES_ALIAS',self.mine.refs['SPECIES_ALIAS']])))
        out.write('{0}\n'.format('='.join(['ASSEMBLY_ACC',self.mine.refs['ASSEMBLY_ACC']])))
        out.write('{0}\n'.format('='.join(['GENESET_VER',self.mine.refs['GENESET_VER']])))
        out.write('\n')
        #
        for line in open(self.mine.refs['CHR_LEN_FN']):
            items = line.rstrip('\n').split('\t')
            _id = items[0]
            _len = int(items[1])
            out.write('{0}\n'.format('='.join(['CHR_LEN:{0}'.format(_id),'{0:,}'.format(_len)])))
        #
        out.write('\n')
        out.write('{0}\n'.format('#'*50))
        out.write('\n')
        out.close()


















