import os
import sys

class MAKE_CMD:
    def __init__(self):
        self.bismark_home = '/BiO/BioPeople/siyoo/MINE/tools/Bismark_v0.19.1'
        self.bismark_map_exe = os.path.join(self.bismark_home, 'bismark')
        self.bismark_dedup_exe = os.path.join(self.bismark_home, 'deduplicate_bismark')
        self.bismark_call_exe = os.path.join(self.bismark_home, 'bismark_methylation_extractor')
        self.bismark_bam2nuc_exe = os.path.join(self.bismark_home, 'bam2nuc')
        self.bismark_cov2cyt_exe = os.path.join(self.bismark_home, 'coverage2cytosine')
        self.bismark_report_exe = os.path.join(self.bismark_home, 'bismark2report')
        self.bismark_summaryreport_exe = os.path.join(self.bismark_home, 'bismark2summary')

        self.bowtie2_home = '/BiO/BioTools/bowtie/bowtie2-2.2.3'
        self.samtools_home = '/BiO/BioTools/samtools/samtools-1.2'

    def load_conf(self, conf_fn):
        self.conf_dic = dict()
        for line in open(conf_fn):
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split('\t')
            #
            tbi_id = items[0]
            sam_id = items[1]
            read_strand = items[2]
            read_path = items[3]
            #
            self.conf_dic.setdefault(tbi_id, {}).setdefault('sam_id', sam_id)
            self.conf_dic.setdefault(tbi_id, {}).setdefault('read', {}).setdefault(read_strand, read_path)
        for tbi_id, info_dic in self.conf_dic.items():
            print(tbi_id, info_dic['sam_id'], info_dic['read']['1'], info_dic['read']['2'])

    def bismark_map(self, outdir, sample, read_1, read_2, bismark_ref_home):
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        _cmds = list()
        _cmds.append(self.bismark_map_exe)
        _cmds.append('--parallel 8')
        _cmds.append('--fastq')
        _cmds.append('--phred33-quals')
        _cmds.append('--output_dir')
        _cmds.append(os.path.join(outdir, sample))
        _cmds.append('--bowtie2')
        _cmds.append('-N 0')
        _cmds.append('-L 20')
        _cmds.append('-p 2')
        _cmds.append('--path_to_bowtie')
        _cmds.append(self.bowtie2_home)
        _cmds.append('--temp_dir')
        _cmds.append(os.path.join(outdir, sample, 'tmp'))
        _cmds.append('--unmapped')
        _cmds.append('--samtools_path')
        _cmds.append(self.samtools_home)
        _cmds.append(bismark_ref_home)
        _cmds.append('-1')
        _cmds.append(read_1)
        _cmds.append('-2')
        _cmds.append(read_2)
        cmd = ' '.join(_cmds)
        return cmd

    def bismark_dedup(self, outdir, sample):
        in_fn = os.path.join(outdir, sample, '{0}_1_bismark_bt2_pe.bam'.format(sample))
        _cmds = list()
        _cmds.append(self.bismark_dedup_exe)
        _cmds.append('--paired')
        _cmds.append('--output_dir')
        _cmds.append(os.path.join(outdir, sample))
        _cmds.append('--bam')
        _cmds.append('--samtools_path')
        _cmds.append(self.samtools_home)
        _cmds.append(in_fn)
        cmd = ' '.join(_cmds)
        return cmd

    def bismark_call(self, outdir, sample, bismark_ref_home):
        in_fn = os.path.join(outdir, sample, '{0}_1_bismark_bt2_pe.deduplicated.bam'.format(sample))
        _cmds = list()
        _cmds.append(self.bismark_call_exe)
        _cmds.append('--multicore 8')
        _cmds.append('--paired-end')
        _cmds.append('--report')
        _cmds.append('--gzip')
        _cmds.append('--output')
        _cmds.append(os.path.join(outdir, sample))
        _cmds.append('--samtools_path')
        _cmds.append(self.samtools_home)
        _cmds.append('--CX_context')
        _cmds.append('--zero_based')
        _cmds.append('--cytosine_report')
        _cmds.append('--genome_folder')
        _cmds.append(bismark_ref_home)
        _cmds.append(in_fn)
        cmd = ' '.join(_cmds)
        return cmd

    def bismark_cov2cyt(self, outdir, sample, bismark_ref_home):
        in_fn = os.path.join(outdir, sample, '{0}_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov'.format(sample))
        _cmds = list()
        _cmds.append(self.bismark_cov2cyt_exe)
        _cmds.append('--dir')
        _cmds.append(os.path.join(outdir, sample))
        _cmds.append('--genome_folder')
        _cmds.append(bismark_ref_home)
        _cmds.append('--CX_context')
        _cmds.append('-o')
        _cmds.append('{0}.CX'.format(sample))
        _cmds.append(in_fn)
        cmd = ' '.join(_cmds)
        return cmd

    def bismark_bam2nuc(self, outdir, sample, bismark_ref_home):
        in_fn = os.path.join(outdir, sample, '{0}_1_bismark_bt2_pe.deduplicated.bam'.format(sample))
        _cmds = list()
        _cmds.append(self.bismark_bam2nuc_exe)
        _cmds.append('--dir')
        _cmds.append(os.path.join(outdir, sample))
        _cmds.append('--samtools_path')
        _cmds.append(self.samtools_home)
        _cmds.append('--genome_folder')
        _cmds.append(bismark_ref_home)
        _cmds.append(in_fn)
        cmd = ' '.join(_cmds)
        return cmd

    def bismark_report(self, outdir, sample):
        _cmds = list()
        _cmds.append(self.bismark_report_exe)
        _cmds.append('--dir')
        _cmds.append(os.path.join(outdir, sample))
        _cmds.append('--output')
        _cmds.append('{0}.bismark_report.html'.format(sample))
        _cmds.append('--alignment_report')
        _cmds.append(os.path.join(outdir, sample, '{0}_1_bismark_bt2_PE_report.txt'.format(sample)))
        _cmds.append('--dedup_report')
        _cmds.append(os.path.join(outdir, sample, '{0}_1_bismark_bt2_pe.deduplication_report.txt'.format(sample)))
        _cmds.append('--splitting_report')
        _cmds.append(os.path.join(outdir, sample, '{0}_1_bismark_bt2_pe.deduplicated_splitting_report.txt'.format(sample)))
        _cmds.append('--mbias_report')
        _cmds.append(os.path.join(outdir, sample, '{0}_1_bismark_bt2_pe.deduplicated.M-bias.txt'.format(sample)))
        _cmds.append('--nucleotide_report')
        _cmds.append(os.path.join(outdir, sample, '{0}_1_bismark_bt2_pe.deduplicated.nucleotide_stats.txt'.format(sample)))
        cmd = ' '.join(_cmds)
        return cmd

    def bismark_summaryreport(self, outdir, samples):
        ln_srcs = list()
        ln_dsts = list()
        for sample in samples:
            ln_srcs.append(os.path.join(outdir, sample, '{0}_1_bismark_bt2_pe.bam'.format(sample)))
            ln_srcs.append(os.path.join(outdir, sample, '{0}_1_bismark_bt2_PE_report.txt'.format(sample)))
            ln_srcs.append(os.path.join(outdir, sample, '{0}_1_bismark_bt2_pe.deduplication_report.txt'.format(sample)))
            ln_srcs.append(os.path.join(outdir, sample, '{0}_1_bismark_bt2_pe.deduplicated_splitting_report.txt'.format(sample)))
            #
            ln_dsts.append(os.path.join(outdir, '{0}_1_bismark_bt2_pe.bam'.format(sample)))
            ln_dsts.append(os.path.join(outdir, '{0}_1_bismark_bt2_PE_report.txt'.format(sample)))
            ln_dsts.append(os.path.join(outdir, '{0}_1_bismark_bt2_pe.deduplication_report.txt'.format(sample)))
            ln_dsts.append(os.path.join(outdir, '{0}_1_bismark_bt2_pe.deduplicated_splitting_report.txt'.format(sample)))
        #
        cmds = list()
        for idx, ln_src in enumerate(ln_srcs):
            cmds.append('ln -s {0} {1}'.format(ln_src, ln_dsts[idx]))
        #
        cmds.append('cd {0}'.format(outdir))
        cmds.append(self.bismark_summaryreport_exe)
        for ln_dst in ln_dsts:
            cmds.append('unlink {0}'.format(ln_dst))
        return cmds




def main(args):
    mkcmd = MAKE_CMD()
    mkcmd.load_conf(args.conf_fn)
    #
    for tbi_id, info_dic in mkcmd.conf_dic.items():
        sample = info_dic['sam_id']
        read_1 = info_dic['read']['1']
        read_2 = info_dic['read']['2']
        #
        print(mkcmd.bismark_map(args.outdir, sample, read_1, read_2, args.bismark_ref_home))
    for tbi_id, info_dic in mkcmd.conf_dic.items():
        sample = info_dic['sam_id']
        read_1 = info_dic['read']['1']
        read_2 = info_dic['read']['2']
        #
        print(mkcmd.bismark_dedup(args.outdir, sample))
    for tbi_id, info_dic in mkcmd.conf_dic.items():
        sample = info_dic['sam_id']
        read_1 = info_dic['read']['1']
        read_2 = info_dic['read']['2']
        #
        print(mkcmd.bismark_call(args.outdir, sample, args.bismark_ref_home))
    for tbi_id, info_dic in mkcmd.conf_dic.items():
        sample = info_dic['sam_id']
        read_1 = info_dic['read']['1']
        read_2 = info_dic['read']['2']
        #
        print(mkcmd.bismark_cov2cyt(args.outdir, sample, args.bismark_ref_home))
    for tbi_id, info_dic in mkcmd.conf_dic.items():
        sample = info_dic['sam_id']
        read_1 = info_dic['read']['1']
        read_2 = info_dic['read']['2']
        #
        print(mkcmd.bismark_bam2nuc(args.outdir, sample, args.bismark_ref_home))
    for tbi_id, info_dic in mkcmd.conf_dic.items():
        sample = info_dic['sam_id']
        read_1 = info_dic['read']['1']
        read_2 = info_dic['read']['2']
        #
        print(mkcmd.bismark_report(args.outdir, sample))
    for cmd in mkcmd.bismark_summaryreport(args.outdir, [_dic['sam_id'] for _id, _dic in mkcmd.conf_dic.items()]):
        print(cmd)




if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    #parser.add_argument('conf_fn')
    #parser.add_argument('outdir')
    #parser.add_argument('bismark_ref_home')
    parser.add_argument('--conf_fn', default='bismark_runner.conf')
    parser.add_argument('--outdir', default='/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark')
    parser.add_argument('--bismark_ref_home', default='/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Ref/chr')
    args = parser.parse_args()
    main(args)
