#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE


class Do_bismark_map:
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
        self.outdir = os.path.join(self.eco.room['analysis'], 'do_bismark')
        self.eco.make_dir(self.outdir)

        # check_stats & run
        if not self.eco.check_stats(self.name):
            map_cmds = self.make_cmds(app_bsm, app_sts, app_bwi, self.input_dic)
            RunQsub(map_cmds, app_bsm.que, '40', self.eco.room['logs'], self.eco.room['scripts'], self.name)
            self.update_ecosystem(self.input_dic)
        else:
            self.update_ecosystem(self.input_dic)

        # finishing
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = [self.name, sam_num, 'map.bam',
                   os.path.join(self.outdir, sam_num, '{0}_1.clean_bismark_bt2_pe.bam'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            #
            lvs = [self.name, sam_num, 'map.report',
                   os.path.join(self.outdir, sam_num, '{0}_1.clean_bismark_bt2_PE_report.txt'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)

    def make_cmds(self, app_bsm, app_sts, app_bwi, input_dic):
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
        input_dic = self.eco.eco_dic['do_filter_fastq']
        return input_dic


class Do_bismark_dedup:
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
        self.outdir = os.path.join(self.eco.room['analysis'], 'do_bismark')
        self.eco.make_dir(self.outdir)

        # check_stats & run
        if not self.eco.check_stats(self.name):
            dedup_cmds = self.make_cmds(app_bsm, app_sts, self.input_dic)
            RunQsub(dedup_cmds, app_bsm.que, '40', self.eco.room['logs'], self.eco.room['scripts'], self.name)
            self.update_ecosystem(self.input_dic)
        else:
            self.update_ecosystem(self.input_dic)

        # finishing
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = [self.name, sam_num, 'dedup.bam',
                   os.path.join(self.outdir, sam_num, '{0}_1.clean_bismark_bt2_pe.deduplicated.bam'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            #
            lvs = [self.name, sam_num, 'dedup.report',
                   os.path.join(self.outdir, sam_num, '{0}_1.clean_bismark_bt2_pe.deduplication_report.txt'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)

    def make_cmds(self, app_bsm, app_sts, input_dic):
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
            opts.append(self.eco.eco_dic['do_bismark_map'][sam_num]['map.bam'])
            cmds.append(' '.join(opts))
        return cmds

    def select_input(self):
        input_dic = dict()
        input_dic = self.eco.eco_dic['do_filter_fastq']
        return input_dic

class Do_bismark_call:
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
        self.outdir = os.path.join(self.eco.room['analysis'], 'do_bismark')
        self.eco.make_dir(self.outdir)

        # check_stats & run
        if not self.eco.check_stats(self.name):
            call_cmds = self.make_cmds(app_bsm, app_sts, self.input_dic)
            RunQsub(call_cmds, app_bsm.que, '16', self.eco.room['logs'], self.eco.room['scripts'], self.name)
            self.update_ecosystem(self.input_dic)
        else:
            self.update_ecosystem(self.input_dic)

        # finishing
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem(self, input_dic):
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
            #lvs = [self.name, sam_num, 'call.cov.zero', os.path.join(self.outdir, sam_num,
            #           '{0}_1.clean_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov.gz'.format(sam_num))]
            lvs = [self.name, sam_num, 'call.cov.zero', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'call.bedGraph', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated.bedGraph.gz'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'splitting.report', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated_splitting_report.txt'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)
            lvs = [self.name, sam_num, 'M-bias.report', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated.M-bias.txt'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)

    def make_cmds(self, app_bsm, app_sts, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            call_fn = os.path.join(self.outdir,sam_num,\
                    '{0}_1.clean_bismark_bt2_pe.deduplicated.bismark.cov.gz'.format(sam_num))
            if not os.path.exists(call_fn):
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
                opts.append('--zero_based')
                opts.append('--cytosine_report')
                opts.append('--genome_folder')
                opts.append(self.mine.refs['BISMARK_DIR'])
                opts.append(self.eco.eco_dic['do_bismark_dedup'][sam_num]['dedup.bam'])
            else:
                continue
            cmds.append(' '.join(opts))
        return cmds

    def select_input(self):
        input_dic = dict()
        input_dic = self.eco.eco_dic['do_filter_fastq']
        return input_dic


class Do_bismark_Cytosine:
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
        #
        self.input_dic = self.select_input()
        self.outdir = os.path.join(self.eco.room['analysis'], 'do_bismark')
        self.eco.make_dir(self.outdir)

        # check_stats & run
        if not self.eco.check_stats(self.name):
            nucl_cmds = self.make_cmds(app_bsm, self.input_dic)
            RunQsub(nucl_cmds, app_bsm.que, '16', self.eco.room['logs'], self.eco.room['scripts'], self.name)
            self.update_ecosystem(self.input_dic)
        else:
            self.update_ecosystem(self.input_dic)
        # finishing
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = [self.name, sam_num, 'CX.report', os.path.join(self.outdir, sam_num,'{0}.CX_report.txt'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)

    def make_cmds(self, app_bsm, input_dic):
        cmds = list()
        for sam_num, f_type_dic in input_dic.iteritems():
            if not os.path.exists(os.path.join(self.outdir,sam_num,'{0}.CX_report.txt'.format(sam_num))):
                opts = list()
                opts.append(app_bsm.exec_s['c_cov'])
                opts.append('--dir')
                opts.append(os.path.join(self.outdir, sam_num))
                opts.append('--genome_folder')
                opts.append(self.mine.refs['BISMARK_DIR'])
                opts.append('--CX_context')
                opts.append('-o')
                opts.append('{0}'.format(sam_num))
                opts.append(f_type_dic['call.cov.zero'])
                cmds.append(' '.join(opts))
            else:
                continue
        return cmds

    def select_input(self):
        input_dic = dict()
        input_dic = self.eco.eco_dic['do_bismark_call']
        return input_dic


class Do_bismark_nucl:
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
        self.outdir = os.path.join(self.eco.room['analysis'], 'do_bismark')
        self.eco.make_dir(self.outdir)

        # check_stats & run
        if not self.eco.check_stats(self.name):
            nucl_cmds = self.make_cmds(app_bsm, app_sts, self.input_dic)
            RunQsub(nucl_cmds, app_bsm.que, '8', self.eco.room['logs'], self.eco.room['scripts'], self.name)
            self.update_ecosystem(self.input_dic)
        else:
            self.update_ecosystem(self.input_dic)
        # finishing
        self.eco.sync_ecosystem()
        self.eco.add_stats(self.name)

    def update_ecosystem(self, input_dic):
        for sam_num, f_type_dic in input_dic.iteritems():
            lvs = [self.name, sam_num, 'nucleotide.report', os.path.join(self.outdir, sam_num,
                       '{0}_1.clean_bismark_bt2_pe.deduplicated.nucleotide_stats.txt'.format(sam_num))]
            self.eco.add_to_ecosystem(lvs)

    def make_cmds(self, app_bsm, app_sts, input_dic):
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
            opts.append(self.eco.eco_dic['do_bismark_dedup'][sam_num]['dedup.bam'])
            cmds.append(' '.join(opts))
        return cmds

    def select_input(self):
        input_dic = dict()
        input_dic = self.eco.eco_dic['do_filter_fastq']
        return input_dic


class Do_bismark_report:
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
            RunQsub(report_cmds, app.que, '8',
                    self.eco.room['logs'],
                    self.eco.room['scripts'],
                    '{0}_each'.format(self.name))
            self.update_ecosystem_report(self.input_dic)
            #
            summary_cmds = self.make_summary_cmds(app, self.input_dic)
            RunQsub(summary_cmds, app.que, '8',
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
        for name in ['do_bismark_map','do_bismark_dedup','do_bismark_call','do_bismark_nucl']:
            if self.eco.eco_dic.has_key(name):
                for sam_num, f_type_dic in self.eco.eco_dic[name].iteritems():
                    input_dic.setdefault(sam_num, {})
                    input_dic[sam_num].update(f_type_dic)
        return input_dic

def main():
    print ("Not supported")

if __name__=='__main__':
    main()

