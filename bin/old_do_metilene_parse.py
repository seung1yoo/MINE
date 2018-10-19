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
        self.outdir = os.path.join(self.eco.room['analysis'], 'do_metilene')
        # check_stats & run
        if not self.eco.check_stats(self.name):
            input_dic = self.eco.eco_dic['do_metilene']
            #
            cmds = self.make_cmds_parse(app, input_dic, '0.05')
            RunQsub(cmds, app.que, '1',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    self.name)
            #
            self.update_ecosystem_parse(input_dic)
        else:
            input_dic = self.eco.eco_dic['do_metilene']
            self.update_ecosystem_parse(input_dic)
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



