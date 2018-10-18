#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE

class Do_something:
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
        #
        self.input_dic = self.select_input()
        self.outdir = os.path.join(self.eco.room['analysis'], 'methyl_call')
        self.eco.make_dir(self.outdir)
        self.make_samdir()

        # check_stats & run
        if not self.eco.check_stats(self.name):
            cmds = self.make_cmds(app_bedtools, self.input_dic)
            RunQsub(cmds, app_bedtools.que, '1',
                    self.eco.room['logs'],
                    self.eco.room['scripts'], self.name)
            self.update_ecosystem(self.input_dic)
        else:
            self.update_ecosystem(self.input_dic)

        # finising
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = [self.name, sam_num,
                   'call.bedGraph.sort',
                   self.make_out_fn(sam_num)]
            self.eco.add_to_ecosystem(lvs)

    def make_cmds(self, app_bedtools, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            out_fn = self.make_out_fn(sam_num)
            if os.path.exists(out_fn):
                eco.logger.info('Found pre-run : {0}'.format(out_fn))
                continue
            opts = list()
            opts.append(app_bedtools.exe)
            opts.append('sort')
            opts.append('-i')
            opts.append(f_type_dic['call.bedGraph'])
            opts.append('|')
            opts.append('gzip')
            opts.append('-f')
            opts.append('>')
            opts.append(out_fn)
            cmds.append(' '.join(opts))
        return cmds

    def make_out_fn(self, sam_num):
        out_fn = os.path.join(self.outdir, sam_num,
            '{0}.call.bedGraph.sort.gz'.format(sam_num))
        return out_fn

    def make_samdir(self):
        for sam_num, f_type_dic in self.input_dic.iteritems():
            samdir = os.path.join(self.outdir, sam_num)
            self.eco.make_dir(samdir)

    def select_input(self):
        input_dic = dict()
        input_dic = self.eco.eco_dic['do_bismark_call']
        return input_dic


