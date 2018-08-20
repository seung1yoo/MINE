#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE

class Do_default:
    def __init__(self, name, mine, eco):
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))
        #
        self.name = name
        self.mine = mine
        self.eco = eco
        #
        if not self.eco.check_stats(self.name):
            # SOME TASK
            pass
        self.eco.add_stats(self.name)
        #
        ## test ################
        #cmds = self.make_cmds()
        #self.run_Popen(cmds)
        ########################

    def _run_Popen(self, cmds):
        procs = list()
        for idx, opts in enumerate(cmds):
            outfh = open(os.path.join(self.eco.room['logs'],'.'.join([self.name,str(idx),'stdout'])), 'a')
            errfh = open(os.path.join(self.eco.room['logs'],'.'.join([self.name,str(idx),'stderr'])), 'a')
            proc = Popen(opts, stdout=outfh, stderr=errfh, shell=True)
            procs.append(proc)
        for proc in procs:
            proc.communicate()
            time.sleep(1)
