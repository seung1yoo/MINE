#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE

class Do_filter_fastq:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))
        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco
        # variable setting
        app = App(self.mine.tools, "FILTER_FQ")
        self.input_dic = self.select_input()
        self.outdir = os.path.join(self.eco.room['analysis'], self.name)
        self.eco.make_dir(self.outdir)
        cmds = self.make_cmds(app, self.input_dic)
        # check_stats & run
        if not self.eco.check_stats(self.name):
            RunQsub(cmds, app.que, '1',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    self.name)
        # finishing
        self.update_ecosystem(self.input_dic)
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = [self.name, sam_num, '1.fq.gz',
                    os.path.join(self.outdir,'{0}_1.clean.fq.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, '2.fq.gz',
                    os.path.join(self.outdir,'{0}_2.clean.fq.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)

    def make_cmds(self, app, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            opts = list()
            opts.append(app.exe)
            opts.append('0.1')
            opts.append('20')
            opts.append('0.4')
            opts.append('0')
            opts.append(f_type_dic['1.fq.gz'])
            opts.append(os.path.join(self.outdir,'{0}_1.clean.fq.gz'.format(sam_num)))
            opts.append(f_type_dic['2.fq.gz'])
            opts.append(os.path.join(self.outdir,'{0}_2.clean.fq.gz'.format(sam_num)))
            cmds.append(' '.join(opts))
        return cmds

    def select_input(self):
        input_dic = dict()
        if self.name in ['do_filter_fastq']:
            input_dic = self.eco.eco_dic['do_link_raw']
        else:
            pass
        return input_dic

def main():
    pass

if __name__=='__main__':
    main()
