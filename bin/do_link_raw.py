#!/usr/bin/python

import os
import sys
import time
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE

class Do_link_raw:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))
        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco
        # variable setting
        self.outdir = os.path.join(self.eco.room['analysis'], self.name)
        self.eco.make_dir(self.outdir)
        # check_stats & run
        if not self.eco.check_stats(self.name):
            self.linker()
        # finishing
        self.update_ecosystem()
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem(self):
        for num, pair_dic in self.mine.sample_dic.iteritems():
            for pair, info_dic in pair_dic.iteritems():
                f_ext = self.eco.find_file_extension(info_dic['path'])
                dst = os.path.join(self.outdir, '{0}_{1}.{2}'.format(num, pair, f_ext))
                self.eco.add_to_ecosystem([self.name,num,'{0}.fq.gz'.format(pair),dst])

    def linker(self):
        for num, pair_dic in self.mine.sample_dic.iteritems():
            for pair, info_dic in pair_dic.iteritems():
                src = info_dic['path']
                f_ext = self.eco.find_file_extension(info_dic['path'])
                dst = os.path.join(self.outdir, '{0}_{1}.{2}'.format(num, pair, f_ext))
                if not os.path.exists(dst):
                    os.symlink(src, dst)

