#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE


class Do_metilene_parse:
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
            RunQsub(cmds, app.que, '1',
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
            RunQsub(cmds, app_bedtools.que, '1',
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


