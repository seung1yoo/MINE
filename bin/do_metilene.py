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




