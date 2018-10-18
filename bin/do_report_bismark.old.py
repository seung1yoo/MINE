#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE

class Do_report_bismark:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))
        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco
        # variable setting
        app = App(self.mine.tools, "BISMARK")
        #
        self.input_dic = self.select_input()
        self.outdir = os.path.join(self.eco.room['report'], 'bismark')
        self.eco.make_dir(self.outdir)
        #
        if not self.eco.check_stats(self.name):
            report_cmds = self.make_report_cmds(app, self.input_dic)
            RunQsub(report_cmds, app.que, '1',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    '{0}_each'.format(self.name))
            self.update_ecosystem_report(self.input_dic)
            #
            summary_cmds = self.make_summary_cmds(app, self.input_dic)
            RunQsub(summary_cmds, app.que, '1',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    '{0}_summary'.format(self.name))
            self.update_ecosystem_summary()
        else:
            self.update_ecosystem_report(self.input_dic)
            self.update_ecosystem_summary()
            #
        # finishing
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem_summary(self):
        lvs = [self.name, 'all', 'summary.html', os.path.join(self.outdir, 'bismark_summary_report.html')]
        self.eco.add_to_ecosystem(lvs)
        lvs = [self.name, 'all', 'summary.txt', os.path.join(self.outdir, 'bismark_summary_report.txt')]
        self.eco.add_to_ecosystem(lvs)

    def make_summary_cmds(self, app, input_dic):
        cmds = list()
        opts = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            sub_opts = list()
            sub_opts.append('ln -s')
            sub_opts.append(f_type_dic['map.bam'])
            sub_opts.append(os.path.join(self.outdir,'{0}_bismark_bt2_pe.bam'.format(self.mine.idmatch_dic[sam_num])))
            opts.append(' '.join(sub_opts))
            #
            sub_opts = list()
            sub_opts.append('ln -s')
            sub_opts.append(f_type_dic['map.report'])
            sub_opts.append(os.path.join(self.outdir,'{0}_bismark_bt2_PE_report.txt'.format(self.mine.idmatch_dic[sam_num])))
            opts.append(' '.join(sub_opts))
            #
            sub_opts = list()
            sub_opts.append('ln -s')
            sub_opts.append(f_type_dic['dedup.report'])
            sub_opts.append(os.path.join(self.outdir,
                '{0}_bismark_bt2_pe.deduplication_report.txt'.format(self.mine.idmatch_dic[sam_num])))
            opts.append(' '.join(sub_opts))
            #
            sub_opts = list()
            sub_opts.append('ln -s')
            sub_opts.append(f_type_dic['splitting.report'])
            sub_opts.append(os.path.join(self.outdir,
                '{0}_bismark_bt2_pe_splitting_report.txt'.format(self.mine.idmatch_dic[sam_num])))
            opts.append(' '.join(sub_opts))
            #
        opts.append('cd {0}'.format(self.outdir))
        opts.append(app.exec_s['summary'])
        for sam_num, f_type_dic in input_dic.iteritems():
            sub_opts = list()
            sub_opts.append('unlink')
            sub_opts.append(os.path.join(self.outdir,'{0}_bismark_bt2_pe.bam'.format(self.mine.idmatch_dic[sam_num])))
            opts.append(' '.join(sub_opts))
            #
            sub_opts = list()
            sub_opts.append('unlink')
            sub_opts.append(os.path.join(self.outdir,'{0}_bismark_bt2_PE_report.txt'.format(self.mine.idmatch_dic[sam_num])))
            opts.append(' '.join(sub_opts))
            #
            sub_opts = list()
            sub_opts.append('unlink')
            sub_opts.append(os.path.join(self.outdir,
                '{0}_bismark_bt2_pe.deduplication_report.txt'.format(self.mine.idmatch_dic[sam_num])))
            opts.append(' '.join(sub_opts))
            #
            sub_opts = list()
            sub_opts.append('unlink')
            sub_opts.append(os.path.join(self.outdir,
                '{0}_bismark_bt2_pe_splitting_report.txt'.format(self.mine.idmatch_dic[sam_num])))
            opts.append(' '.join(sub_opts))
            #
        cmds.append('\n'.join(opts))
        return cmds

    def update_ecosystem_report(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = [self.name, sam_num, 'report', os.path.join(self.outdir,
                '{0}.bismark_report.html'.format(self.mine.idmatch_dic[sam_num]))]
            self.eco.add_to_ecosystem(lvs)

    def make_report_cmds(self, app, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            opts = list()
            opts.append(app.exec_s['report'])
            opts.append('--dir')
            opts.append(self.outdir)
            opts.append('--output')
            opts.append('{0}.bismark_report.html'.format(self.mine.idmatch_dic[sam_num]))
            opts.append('--alignment_report')
            opts.append(f_type_dic['map.report'])
            opts.append('--dedup_report')
            opts.append(f_type_dic['dedup.report'])
            opts.append('--splitting_report')
            opts.append(f_type_dic['splitting.report'])
            opts.append('--mbias_report')
            opts.append(f_type_dic['M-bias.report'])
            opts.append('--nucleotide_report')
            opts.append(f_type_dic['nucleotide.report'])
            cmds.append(' '.join(opts))
        return cmds



    def select_input(self):
        input_dic = dict()
        if self.name in ['do_report_bismark']:
            input_dic = self.eco.eco_dic['do_bismark']
        else:
            pass
        return input_dic

def main():
    pass

if __name__=='__main__':
    main()

