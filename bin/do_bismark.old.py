#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE

class Do_bismark:
    def __init__(self, name, mine, eco):
        # For debugging
        mine.logger.debug(dir(mine))
        eco.logger.debug(dir(eco))
        # inheritance
        self.name = name
        self.mine = mine
        self.eco = eco
        # variable setting
        app_bsm = App(self.mine.tools, "BISMARK")
        app_sts = App(self.mine.tools, "SAMTOOLS")
        app_bwi = App(self.mine.tools, "BOWTIE2")
        #
        self.input_dic = self.select_input()
        self.outdir = os.path.join(self.eco.room['analysis'], self.name)
        self.eco.make_dir(self.outdir)
        #
        # check_stats & run
        if not self.eco.check_stats(self.name):
            map_cmds = self.make_map_cmds(app_bsm, app_sts, app_bwi, self.input_dic)
            RunQsub(map_cmds, app_bsm.que, '16',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    '{0}_map'.format(self.name))
            self.update_ecosystem_map(self.input_dic)
            #
            dedup_cmds = self.make_dedup_cmds(app_bsm, app_sts, self.input_dic)
            RunQsub(dedup_cmds, app_bsm.que, '1',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    '{0}_dedup'.format(self.name))
            self.update_ecosystem_dedup(self.input_dic)
            #
            call_cmds = self.make_call_cmds(app_bsm, app_sts, self.input_dic)
            RunQsub(call_cmds, app_bsm.que, '8',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    '{0}_call'.format(self.name))
            self.update_ecosystem_call(self.input_dic)
            #
            nucl_cmds = self.make_nucl_cmds(app_bsm, app_sts, self.input_dic)
            RunQsub(nucl_cmds, app_bsm.que, '1',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    '{0}_nucl'.format(self.name))
            self.update_ecosystem_nucl(self.input_dic)
        else:
            self.update_ecosystem_map(self.input_dic)
            self.update_ecosystem_dedup(self.input_dic)
            self.update_ecosystem_call(self.input_dic)
            self.update_ecosystem_nucl(self.input_dic)
        # finishing
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem_nucl(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = [self.name, sam_num, 'nucleotide.report', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated.nucleotide_stats.txt'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)

    def make_nucl_cmds(self, app_bsm, app_sts, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            if os.path.exists(os.path.join(self.outdir,sam_num,'{0}_1.clean_bismark_bt2_pe.deduplicated.nucleotide_stats.txt'.format(sam_num))):
                continue
            opts = list()
            opts.append(app_bsm.exec_s['nucl_cov'])
            opts.append('--dir')
            opts.append(os.path.join(self.outdir, sam_num))
            opts.append('--samtools_path')
            opts.append(app_sts.home)
            opts.append('--genome_folder')
            opts.append(self.mine.refs['BISMARK_DIR'])
            opts.append(self.eco.eco_dic[self.name][sam_num]['dedup.bam'])
            cmds.append(' '.join(opts))
        return cmds

    def update_ecosystem_call(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = [self.name, sam_num, 'call.CHH.OB', os.path.join(self.outdir, sam_num,
                       'CHH_OB_{0}_1.clean_bismark_bt2_pe.deduplicated.txt.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'call.CHH.OT', os.path.join(self.outdir, sam_num,
                       'CHH_OT_{0}_1.clean_bismark_bt2_pe.deduplicated.txt.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'call.CHG.OB', os.path.join(self.outdir, sam_num,
                       'CHG_OB_{0}_1.clean_bismark_bt2_pe.deduplicated.txt.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'call.CHG.OT', os.path.join(self.outdir, sam_num,
                       'CHG_OT_{0}_1.clean_bismark_bt2_pe.deduplicated.txt.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'call.CpG.OB', os.path.join(self.outdir, sam_num,
                       'CpG_OB_{0}_1.clean_bismark_bt2_pe.deduplicated.txt.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'call.CpG.OT', os.path.join(self.outdir, sam_num,
                       'CpG_OT_{0}_1.clean_bismark_bt2_pe.deduplicated.txt.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'call.cov', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated.bismark.cov.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'call.cov.zero', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'call.bedGraph', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated.bedGraph.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'call.cx', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated.CX_report.txt.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'splitting.report', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated_splitting_report.txt'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'M-bias.report', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated.M-bias.txt'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)

    def make_call_cmds(self, app_bsm, app_sts, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            if os.path.exists(os.path.join(self.outdir,sam_num,'{0}_1.clean_bismark_bt2_pe.deduplicated.CX_report.txt.gz'.format(sam_num))):
                continue
            opts = list()
            opts.append(app_bsm.exec_s['call'])
            opts.append('--multicore 8')
            opts.append('--paired-end')
            opts.append('--report')
            opts.append('--gzip')
            opts.append('--output')
            opts.append(os.path.join(self.outdir, sam_num))
            opts.append('--samtools_path')
            opts.append(app_sts.home)
            opts.append('--CX_context')
            opts.append('--zero_based')
            opts.append('--cytosine_report')
            opts.append('--genome_folder')
            opts.append(self.mine.refs['BISMARK_DIR'])
            opts.append(self.eco.eco_dic[self.name][sam_num]['dedup.bam'])
            cmds.append(' '.join(opts))
        return cmds

    def update_ecosystem_dedup(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = [self.name, sam_num, 'dedup.bam',
                   os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated.bam'.format(sam_num)
                       )
                   ]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'dedup.report',
                   os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplication_report.txt'.format(sam_num)
                       )
                   ]
            self.eco.add_to_ecosystem(lvs)

    def make_dedup_cmds(self, app_bsm, app_sts, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            if os.path.exists(os.path.join(self.outdir,sam_num,'{0}_1.clean_bismark_bt2_pe.deduplication_report.txt'.format(sam_num))):
                continue
            opts = list()
            opts.append(app_bsm.exec_s['dedup'])
            opts.append('--paired')
            opts.append('--output_dir')
            opts.append(os.path.join(self.outdir, sam_num))
            opts.append('--bam')
            opts.append('--samtools_path')
            opts.append(app_sts.home)
            opts.append(self.eco.eco_dic[self.name][sam_num]['map.bam'])
            cmds.append(' '.join(opts))
        return cmds

    def update_ecosystem_map(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = [self.name, sam_num, 'map.bam',
                   os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.bam'.format(sam_num)
                       )
                   ]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'map.report',
                   os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_PE_report.txt'.format(sam_num)
                       )
                   ]
            self.eco.add_to_ecosystem(lvs)

    def make_map_cmds(self, app_bsm, app_sts, app_bwi, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            if os.path.exists(os.path.join(self.outdir,sam_num,'{0}_1.clean_bismark_bt2_PE_report.txt'.format(sam_num))):
                continue
            opts = list()
            opts.append(app_bsm.exe)
            opts.append('--parallel 8')
            opts.append('--fastq')
            opts.append('--phred33-quals')
            opts.append('--output_dir')
            opts.append(os.path.join(self.outdir, sam_num))
            opts.append('--bowtie2 -N 0 -L 20 -p 2')
            opts.append('--path_to_bowtie')
            opts.append(app_bwi.home)
            opts.append('--temp_dir')
            opts.append(os.path.join(self.outdir, sam_num, 'tmp'))
            opts.append('--unmapped')
            opts.append('--samtools_path')
            opts.append(app_sts.home)
            opts.append(self.mine.refs['BISMARK_DIR'])
            opts.append('-1')
            opts.append(f_type_dic['1.fq.gz'])
            opts.append('-2')
            opts.append(f_type_dic['2.fq.gz'])
            cmds.append(' '.join(opts))
        return cmds

    def select_input(self):
        input_dic = dict()
        if self.name in ['do_bismark']:
            input_dic = self.eco.eco_dic['do_filter_fastq']
        else:
            pass
        return input_dic

def main():
    pass

if __name__=='__main__':
    main()

