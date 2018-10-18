#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE

class Do_metilene:
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
        app_bedtools = App(self.mine.tools, "BEDTOOLS")
        #
        self.outdir = os.path.join(self.eco.room['analysis'], self.name)
        self.eco.make_dir(self.outdir)
        # check_stats & run
        if not self.eco.check_stats(self.name):
            # execute metilene
            input_dic = self.eco.eco_dic['do_metilene_input']
            #
            cmds = self.make_cmds_metilene(app_metilene, input_dic)
            RunQsub(cmds, app_metilene.que, '10',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    self.name)
            #
            self.update_ecosystem_metilene(input_dic)
        else:
            input_dic = self.eco.eco_dic['do_metilene_input']
            self.update_ecosystem_metilene(input_dic)
        # finishing
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def make_cmds_metilene(self, app, input_dic):
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

    def update_ecosystem_metilene(self, input_dic):
        for dmr_num, f_type_dic in input_dic.iteritems():
            outdir = os.path.join(self.outdir, dmr_num)
            lvs = [self.name, dmr_num, 'denovo',
                    os.path.join(outdir, 'metilene_G1_G2.output.denovo')]
            self.eco.add_to_ecosystem(lvs)



class Do_metilene_bismark2input:
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
        self.outdir = os.path.join(self.eco.room['analysis'], self.name)
        self.eco.make_dir(self.outdir)
        # check_stats & run
        if not self.eco.check_stats(self.name):
            input_dic = self.eco.eco_dic['do_bismark']
            cmds = self.make_cmds_bedsort(app_bedtools, input_dic)
            RunQsub(cmds, app_bedtools.que, '1',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    '{0}_bedsort'.format(self.name))
            self.update_ecosystem_bedsort(input_dic)
            #
            cmds = self.make_cmds_bed2input(app_metilene, app_bedtools, input_dic)
            RunQsub(cmds, app_bedtools.que, '1',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    self.name)
            self.update_ecosystem_bed2input(input_dic)
            #
            self.input_reheader()
            self.update_ecosystem_reheader(input_dic)
        else:
            input_dic = self.eco.eco_dic['do_bismark']
            self.update_ecosystem_bedsort(input_dic)
            self.update_ecosystem_bed2input(input_dic)
            self.update_ecosystem_reheader(input_dic)
        # finishing
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

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

    def make_cmds_bed2input(self, app_metilene, app_bedtools, input_dic):
        cmds = list()
        for dmr_num, g_dic in self.mine.dmr_dic.iteritems():
            outdir = os.path.join(self.outdir, dmr_num)
            self.eco.make_dir(outdir)
            group1_ins = [input_dic[sam_num]['call.sort.bedGraph'] \
                    for sam_num in g_dic['G1']]
            group2_ins = [input_dic[sam_num]['call.sort.bedGraph'] \
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

    def make_cmds_bedsort(self, app_bedtools, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            opts = list()
            opts.append(app_bedtools.exe)
            opts.append('sort')
            opts.append('-i')
            opts.append(f_type_dic['call.bedGraph'])
            opts.append('|')
            opts.append('gzip')
            opts.append('-f')
            opts.append('>')
            opts.append(f_type_dic['call.bedGraph'].replace('.bedGraph.gz','.sort.bedGraph.gz'))
            cmds.append(' '.join(opts))
        return cmds

    def update_ecosystem_bedsort(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = ['do_bismark',sam_num,'call.sort.bedGraph',f_type_dic['call.bedGraph'].replace('.bedGraph.gz','.sort.bedGraph.gz')]
            self.eco.add_to_ecosystem(lvs)


