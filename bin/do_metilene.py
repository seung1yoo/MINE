#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE

class Do_metilene_anno:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))

        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco

        # variable setting
        self.app_bedtools = App(self.mine.tools, "BEDTOOLS")
        #
        self.input_dic = self.select_input()
        self.outdir = os.path.join(self.eco.room['analysis'], 'do_metilene')
        self.eco.make_dir(self.outdir)
        #
        self.dmr_tags = ['CHR','START','END','q','diff','n.CpG',
                         'p.MWU','p.2D_KS','mean.G1','mean.G2']
        self.anno_tags = ['UP1K','genebody','DW1K',
                          '5UTR','CDS','EXON','3UTR',
                          'Promoter','HCP','ICP','LCP',
                          'Nshelf','Nshore','CGI','Sshore','Sshelf']
        self.anno_fn_dic = dict()
        self.init_anno_fn_dic()

        # check_stats & run
        if not self.eco.check_stats(self.name):
            cmds = self.make_cmds_intersect()
            RunQsub(cmds, self.app_bedtools.que, '5',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    self.name)
        else:
            pass
        dmr_anno_fn_dic = self.annotation()
        self.update_ecosystem(dmr_anno_fn_dic)

        # finising
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem(self, dmr_anno_fn_dic):
        for dmr_num, dmr_anno_fn in dmr_anno_fn_dic.iteritems():
            lvs = [self.name, dmr_num, 'dmr_anno', dmr_anno_fn]
            self.eco.add_to_ecosystem(lvs)

    def annotation(self):
        dmr_anno_fn_dic = dict()
        for dmr_num, f_type_dic in self.input_dic.iteritems():
            #
            out_fn = '{0}.anno'.format(f_type_dic['denovo'])
            dmr_anno_fn_dic.setdefault(dmr_num, out_fn)
            if os.path.exists(out_fn) and os.path.getsize(out_fn) > 1000:
                continue
            #
            dmr_dic = dict()
            for line in open(f_type_dic['denovo']):
                items = line.rstrip('\n').split('\t')
                if items[0] in ['#CHR']:
                    idx_dic = dict()
                    for idx, item in enumerate(items):
                        idx_dic.setdefault(item, idx)
                    continue
                key = self.annotation_key(items)
                # depend on self.dmr_tags
                dmr_dic.setdefault(key,{}).setdefault('CHR',items[idx_dic['#CHR']])
                dmr_dic.setdefault(key,{}).setdefault('START',items[idx_dic['START']])
                dmr_dic.setdefault(key,{}).setdefault('END',items[idx_dic['STOP']])
                dmr_dic.setdefault(key,{}).setdefault('q',items[idx_dic['q-value']])
                dmr_dic.setdefault(key,{}).setdefault('diff',
                        self.p_to_m_to_p(items[idx_dic['mean_methylation_difference']]))
                dmr_dic.setdefault(key,{}).setdefault('n.CpG',items[idx_dic['No.CpGs']])
                dmr_dic.setdefault(key,{}).setdefault('p.MWU',items[idx_dic['p(MWU)']])
                dmr_dic.setdefault(key,{}).setdefault('p.2D_KS',
                        items[idx_dic['p(2D_KS)']])
                dmr_dic.setdefault(key,{}).setdefault('mean.G1',
                        items[idx_dic['mean_{0}:G1'.format(dmr_num)]])
                dmr_dic.setdefault(key,{}).setdefault('mean.G2',
                        items[idx_dic['mean_{0}:G2'.format(dmr_num)]])
                for anno_tag in self.anno_tags:
                    dmr_dic.setdefault(key,{}).setdefault(anno_tag,{})
            #
            for _chr, anno_fn in self.anno_fn_dic[dmr_num].iteritems():
                idx_dic = dict()
                a_headers = self.grep_a_headers(dmr_num)
                b_headers = self.grep_b_headers(_chr)
                headers = a_headers
                headers.extend(b_headers)
                for idx, header in enumerate(headers):
                    idx_dic.setdefault(header, idx)
                #
                for line in open(anno_fn):
                    items = line.rstrip('\n').split('\t')
                    key = self.annotation_key(items)
                    #
                    for anno_tag in self.anno_tags:
                        anno_key = 'MEMBER.{0}'.format(anno_tag)
                        anno_units = items[idx_dic[anno_key]].split(',')
                        for unit in anno_units:
                            if not unit in ['-']:
                                dmr_dic[key][anno_tag].setdefault(unit,None)
                            else:
                                pass
                        #
                    #
                #
            out = open(out_fn, 'w')
            out_headers = list()
            out_headers.append('DMR_ID')
            out_headers.extend(self.dmr_tags)
            out_headers.extend(['DMR.YN.{0}'.format(self.mine.dmr_cut),'DMR.UPDOWN'])
            out_headers.extend(['NUM.{0}'.format(x) for x in self.anno_tags])
            out_headers.extend(['MEMBER.{0}'.format(x) for x in self.anno_tags])
            out.write('#{0}\n'.format('\t'.join(out_headers)))
            dmr_count = 0
            for key, info_dic in sorted(dmr_dic.iteritems()):
                dmr_count += 1
                new_items = list()
                new_items.append('{0}N{1:0>5}'.format(dmr_num, dmr_count))
                for tag in self.dmr_tags:
                    new_items.append(info_dic[tag])
                new_items.append(self.dmr_yes_or_no(info_dic))
                new_items.append(self.dmr_up_or_down(info_dic))
                for tag in self.anno_tags:
                    new_items.append(str(len(info_dic[tag].keys())))
                for tag in self.anno_tags:
                    new_items.append(','.join(info_dic[tag].keys()))
                out.write('{0}\n'.format('\t'.join(new_items)))
            out.close()
        return dmr_anno_fn_dic

    def dmr_yes_or_no(self, info_dic):
        test, cut = self.mine.dmr_cut.split(':')
        if test in ['P']:
            value = info_dic['p.MWU']
        elif test in ['Q']:
            value = info_dic['q']
        #
        if float(value) <= float(cut):
            return 'Y'
        else:
            return 'N'

    def dmr_up_or_down(self, info_dic):
        diff_value = info_dic['diff']
        if float(diff_value) > 0:
            return 'UP'
        elif float(diff_value) < 0:
            return 'DOWN'
        else:
            return 'FLAT'

    def p_to_m_to_p(self, value):
        return str(float(value) * -1)

    def annotation_key(self, items):
        return '_'.join(items[:3])

    def grep_a_headers(self, dmr_num):
        fn = self.input_dic[dmr_num]['denovo']
        fh = open(fn)
        return fh.readline().rstrip('\n').split('\t')

    def grep_b_headers(self, _chr):
        fn = self.mine.refs['CGSITE_ANNO'][_chr]
        fh = open(fn.replace('[HOME]',self.mine.refs['HOME']))
        return fh.readline().rstrip('\n').split('\t')

    def make_cmds_intersect(self):
        cmds = list()
        for dmr_num, f_type_dic in self.input_dic.iteritems():
            opts = list()
            for _chr, _f_path in self.mine.refs['CGSITE_ANNO'].iteritems():
                out_fn = self.intersect_out_fn(f_type_dic, _chr)
                if os.path.exists(out_fn):
                    continue
                opts.append(self.app_bedtools.exe)
                opts.append('intersect')
                opts.append('-a')
                opts.append(f_type_dic['denovo'])
                opts.append('-b')
                opts.append(_f_path.replace('[HOME]', self.mine.refs['HOME']))
                opts.append('-wa')
                opts.append('-wb')
                opts.append('>')
                opts.append(out_fn)
                opts.append('\n')
            cmds.append(' '.join(opts))
        return cmds

    def intersect_out_fn(self, f_type_dic, _chr):
        return os.path.join('{0}.{1}.anno'.format(f_type_dic['denovo'], _chr))

    def init_anno_fn_dic(self):
        for dmr_num, f_type_dic in self.input_dic.iteritems():
            for _chr, _f_path in self.mine.refs['CGSITE_ANNO'].iteritems():
                out_fn = self.intersect_out_fn(f_type_dic, _chr)
                self.anno_fn_dic.setdefault(dmr_num,{}).setdefault(_chr,out_fn)

    def select_input(self):
        input_dic = dict()
        input_dic = self.eco.eco_dic['do_metilene_exe']
        return input_dic

class _Do_metilene_parse:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))

        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco

        # variable setting
        app = App(self.mine.tools, "METILENE")
        #
        self.input_dic = self.select_input()
        self.outdir = os.path.join(self.eco.room['analysis'], 'do_metilene')
        self.eco.make_dir(self.outdir)

        # check_stats & run
        if not self.eco.check_stats(self.name):
            #
            cmds = self.make_cmds_parse(app, self.input_dic, '0.05')
            RunQsub(cmds, app.que, '10',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    self.name)
            #
            self.update_ecosystem_parse(self.input_dic)
        else:
            self.update_ecosystem_parse(self.input_dic)
        # finising
        self.eco.sync_ecosystem()
        #self.eco.add_stats(self.name)

    def make_cmds_parse(self, app, input_dic, q_cut):
        cmds = list()
        for dmr_num, f_type_dic in input_dic.iteritems():
            opts = list()
            #
            opts.append(app.exec_s['filter'])
            opts.append('-q')
            opts.append(f_type_dic['denovo'])
            opts.append('-o')
            opts.append(f_type_dic['denovo'])
            opts.append('-p {0}'.format(q_cut))
            opts.append('-c 10')
            opts.append('-d 0.1')
            opts.append('-l 0')
            opts.append('-a G1')
            opts.append('-b G2')
            opts.append('\n')
            opts.append('\n')
            #
            opts.append("header='#CHR@START@STOP@q-value@mean_methylation_difference@No.CpGs@p(MWU)@p(2D_KS)@mean_{0}:G1@mean_{0}:G2'".format(dmr_num))
            opts.append('\n')
            opts.append('echo $header')
            opts.append("| sed 's/@/\t/g'")
            opts.append('>')
            opts.append('{0}_qval.{1}.xls'.format(f_type_dic['denovo'],q_cut))
            opts.append('\n')
            opts.append('cat')
            opts.append('{0}_qval.{1}.out'.format(f_type_dic['denovo'],q_cut))
            opts.append('| sort -V -k1,1 -k2,2n')
            opts.append('>>')
            opts.append('{0}_qval.{1}.xls'.format(f_type_dic['denovo'],q_cut))
            cmds.append(' '.join(opts))
        return cmds

    def update_ecosystem_parse(self, input_dic):
        for dmr_num, f_type_dic in input_dic.iteritems():
            pass

    def select_input(self):
        input_dic = dict()
        input_dic = self.eco.eco_dic['do_metilene_exe']
        return input_dic


class Do_metilene_exe:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))

        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco

        # variable setting
        app_metilene = App(self.mine.tools, "METILENE")
        #
        self.input_dic = self.select_input()
        self.outdir = os.path.join(self.eco.room['analysis'], 'do_metilene')
        self.eco.make_dir(self.outdir)

        # check_stats & run
        if not self.eco.check_stats(self.name):
            # execute metilene
            #
            cmds = self.make_cmds(app_metilene, self.input_dic)
            RunQsub(cmds, app_metilene.que, '10',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    self.name)
            #
            self.update_ecosystem(self.input_dic)
        else:
            self.update_ecosystem(self.input_dic)

        # finishing
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def make_cmds(self, app, input_dic):
        cmds = list()
        for dmr_num, f_type_dic in input_dic.iteritems():
            outdir = os.path.join(self.outdir, dmr_num)
            self.eco.make_dir(outdir)
            #
            opts = list()
            opts.append("header='#CHR@START@STOP@q-value@mean_methylation_difference@No.CpGs@p(MWU)@p(2D_KS)@mean_{0}:G1@mean_{0}:G2'".format(dmr_num))
            opts.append('\n')
            opts.append('echo $header')
            opts.append("| sed 's/@/\t/g'")
            opts.append('>')
            opts.append(os.path.join(outdir, 'metilene_G1_G2.output.denovo'))
            opts.append('\n')
            #
            opts.append(app.exe)
            opts.append('--maxdist 300')
            opts.append('--mincpgs 10')
            opts.append('--minMethDiff 0.1')
            opts.append('--threads 10')
            opts.append('--mode 1')
            opts.append('--groupA G1')
            opts.append('--groupB G2')
            opts.append('--minNoA {0}'.format(len(self.mine.dmr_dic[dmr_num]['G1'])))
            opts.append('--minNoB {0}'.format(len(self.mine.dmr_dic[dmr_num]['G2'])))
            opts.append('--valley 0.7')
            opts.append(input_dic[dmr_num]['input'])
            opts.append('| sort -V -k1,1 -k2,2n')
            opts.append('>>')
            opts.append(os.path.join(outdir, 'metilene_G1_G2.output.denovo'))
            #
            cmds.append(' '.join(opts))
        return cmds

    def update_ecosystem(self, input_dic):
        for dmr_num, f_type_dic in input_dic.iteritems():
            outdir = os.path.join(self.outdir, dmr_num)
            lvs = [self.name, dmr_num, 'denovo',
                    os.path.join(outdir, 'metilene_G1_G2.output.denovo')]
            self.eco.add_to_ecosystem(lvs)

    def select_input(self):
        input_dic = dict()
        input_dic = self.eco.eco_dic['do_metilene_prepare']
        return input_dic

class Do_metilene_prepare:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))

        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco

        # variable setting
        app_bedtools = App(self.mine.tools, "BEDTOOLS")
        app_metilene = App(self.mine.tools, "METILENE")
        #
        self.input_dic = self.select_input()
        self.outdir = os.path.join(self.eco.room['analysis'], 'do_metilene')
        self.eco.make_dir(self.outdir)

        # check_stats & run
        if not self.eco.check_stats(self.name):
            cmds = self.make_cmds(app_metilene, app_bedtools, self.input_dic)
            RunQsub(cmds, app_bedtools.que, '10',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    self.name)
            self.update_ecosystem_bed2input(self.input_dic)
            #
            self.input_reheader()
            self.update_ecosystem_reheader(self.input_dic)
        else:
            self.update_ecosystem_bed2input(self.input_dic)
            self.update_ecosystem_reheader(self.input_dic)

        # finishing
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def make_cmds(self, app_metilene, app_bedtools, input_dic):
        cmds = list()
        for dmr_num, g_dic in self.mine.dmr_dic.iteritems():
            outdir = os.path.join(self.outdir, dmr_num)
            self.eco.make_dir(outdir)
            group1_ins = [input_dic[sam_num]['call.bedGraph.final'] \
                    for sam_num in g_dic['G1']]
            group2_ins = [input_dic[sam_num]['call.bedGraph.final'] \
                    for sam_num in g_dic['G2']]
            #
            opts = list()
            opts = ['perl']
            opts.append(app_metilene.exec_s['bed2input'])
            opts.append('--in1')
            opts.append(','.join(group1_ins))
            opts.append('--in2')
            opts.append(','.join(group2_ins))
            opts.append('--out')
            opts.append(os.path.join(outdir, 'metilene_G1_G2.input.tmp'))
            opts.append('--h1')
            opts.append('G1')
            opts.append('--h2')
            opts.append('G2')
            opts.append('-b')
            opts.append(app_bedtools.exe)
            cmds.append(' '.join(opts))
        return cmds

    def update_ecosystem_bed2input(self, input_dic):
        for dmr_num, g_dic in self.mine.dmr_dic.iteritems():
            outdir = os.path.join(self.outdir, dmr_num)
            lvs = [self.name, dmr_num, 'input.tmp',
                    os.path.join(outdir, 'metilene_G1_G2.input.tmp')]
            self.eco.add_to_ecosystem(lvs)

    def input_reheader(self):
        for dmr_num, g_dic in self.mine.dmr_dic.iteritems():
            outdir = os.path.join(self.outdir, dmr_num)
            #
            reheaders = list()
            reheaders.append('chrom')
            reheaders.append('pos')
            head_G1s = ['G1_{0}'.format(sam_num) for sam_num in g_dic['G1']]
            reheaders.extend(head_G1s)
            head_G2s = ['G2_{0}'.format(sam_num) for sam_num in g_dic['G2']]
            reheaders.extend(head_G2s)
            #
            out_fh = open(self.eco.eco_dic[self.name][dmr_num]['input.tmp'].replace('.tmp',''), 'w')
            out_fh.write('{0}\n'.format('\t'.join(reheaders)))
            for line in open(self.eco.eco_dic[self.name][dmr_num]['input.tmp']):
                if line.startswith('chrom\tpos'):
                    continue
                out_fh.write(line)
            out_fh.close()

    def update_ecosystem_reheader(self, inputdic):
        for dmr_num, g_dic in self.mine.dmr_dic.iteritems():
            outdir = os.path.join(self.outdir, dmr_num)
            lvs = [self.name, dmr_num, 'input',
                    os.path.join(outdir, 'metilene_G1_G2.input')]
            self.eco.add_to_ecosystem(lvs)

    def select_input(self):
        input_dic = dict()
        input_dic = self.eco.eco_dic['do_targetcheck_grep']
        return input_dic


