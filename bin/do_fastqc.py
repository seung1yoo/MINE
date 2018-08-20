#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE

class Do_fastqc:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))
        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco
        # variable setting
        app = App(self.mine.tools, "FASTQC")
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
            for f_type, f_path in f_type_dic.iteritems():
                base_path = f_path.rstrip('.fq.gz')
                fastqc_dir = '{0}_fastqc'.format(base_path)
                image_dir = '{0}/Images'.format(fastqc_dir)
                seq_content = '{0}/per_base_sequence_content.png'.format(image_dir)
                seq_quality = '{0}/per_base_quality.png'.format(image_dir)
                if f_type.startswith('1'):
                    self.eco.add_to_ecosystem(
                            [self.name,sam_num,'1.qc.content.png',seq_content])
                    self.eco.add_to_ecosystem(
                            [self.name,sam_num,'1.qc.quality.png',seq_quality])
                elif f_type.startswith('2'):
                    self.eco.add_to_ecosystem(
                            [self.name,sam_num,'2.qc.content.png',seq_content])
                    self.eco.add_to_ecosystem(
                            [self.name,sam_num,'2.qc.quality.png',seq_quality])

    def make_cmds(self, app, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            for f_type, f_path in f_type_dic.iteritems():
                opts = list()
                opts.append(app.exe)
                opts.append('-f fastq')
                opts.append('--extract')
                opts.append('-t 7')
                opts.append(f_path)
                cmds.append(' '.join(opts))
        return cmds

    def select_input(self):
        input_dic = dict()
        if self.name in ['do_fastqc_raw']:
            input_dic = self.eco.eco_dic['do_link_raw']
        elif self.name in ['do_fastqc_clean']:
            input_dic = self.eco.eco_dic['do_filter_fastq']
        else:
            pass
        return input_dic

def main():
    pass

if __name__=='__main__':
    main()
