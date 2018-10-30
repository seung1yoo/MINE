#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE

class Do_stat_fastq:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))
        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco
        # variable setting
        app = App(self.mine.tools, "STAT_FQ")
        self.input_dic = self.select_input()
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
            stat_fn = '{0}.stat'.format(f_type_dic['1.fq.gz'])
            if self.name in ['do_stat_fastq_raw']:
                self.eco.add_to_ecosystem([self.name,sam_num,'raw.stat',stat_fn])
            elif self.name in ['do_stat_fastq_clean']:
                self.eco.add_to_ecosystem([self.name,sam_num,'clean.stat',stat_fn])

    def make_cmds(self, app, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            stat_fn = '{0}.stat'.format(f_type_dic['1.fq.gz'])
            if os.path.exists(stat_fn):
                continue
            opts = list()
            opts.append(app.exe)
            opts.append(f_type_dic['1.fq.gz'])
            opts.append(f_type_dic['2.fq.gz'])
            cmds.append(' '.join(opts))
        return cmds

    def select_input(self):
        input_dic = dict()
        if self.name in ['do_stat_fastq_raw']:
            input_dic = self.eco.eco_dic['do_link_raw']
        elif self.name in ['do_stat_fastq_clean']:
            input_dic = self.eco.eco_dic['do_filter_fastq']
        else:
            pass
        return input_dic

def main():
    print("Not supported")

if __name__=='__main__':
    main()
