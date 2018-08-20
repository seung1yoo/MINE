#!/usr/bin/python

import os
import sys
#
from lib.Qsub import *
from Bio import SeqIO
from subprocess import Popen, PIPE
import matplotlib as mpl
import matplotlib.pylab as plt
#
bismark_build = '/BiO/BioPeople/siyoo/MINE/tools/Bismark_v0.19.1/bismark_genome_preparation'
samtools = '/BiO/BioTools/Rine/Tools/samtools/current/samtools'
bowtie2 = '/BiO/BioTools/bowtie/bowtie2-2.2.3'
newcpgreport = '/usr/bin/newcpgreport'
bedtools = '/usr/bin/bedtools'
classify_promt = '/BiO/BioPeople/siyoo/MINE/bin/mine_build_promt.py'

def fai_to_len(fai):
    cmds = list()
    cmds.append('cut')
    cmds.append('-f-2')
    cmds.append(fai)
    proc = Popen(cmds, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate()
    if err:
        print err
        sys.exit()
    return out

def do_samtools_faidx(fa_fn):
    global samtools
    fai_fn = '{0}.fai'.format(fa_fn)
    if not os.path.exists(fai_fn):
        cmds = list()
        cmds.append(samtools)
        cmds.append('faidx')
        cmds.append(fa_fn)
        proc = Popen(cmds, stdout=PIPE, stderr=PIPE)
        out, err = proc.communicate()
        if err:
            print err
            sys.exit()
    return fai_fn

def select_chrs(ref_fa_fn):
    ref_fai_fn = do_samtools_faidx(ref_fa_fn)
    #
    chrs = list()
    out = fai_to_len(ref_fai_fn)
    for line in out.split('\n'):
        if not line:
            continue
        chr_id, chr_len = line.split('\t')
        first_chr_id = chr_id[0]
        if first_chr_id in ['1','2','3','4','5','6','7','8','9']:
            chrs.append(chr_id)
        elif chr_id in ['X','Y','MT']:
            chrs.append(chr_id)
        else:
            continue
    #
    return chrs

def fa_filter_by_chr(ref_fa_fn, chrs, home, name):
    refchr_fa_fn = '{0}/{1}.chr.fa'.format(home, name)
    if not os.path.exists(refchr_fa_fn):
        refchr_fa_fh = open(refchr_fa_fn, 'w')
        for record in SeqIO.parse(open(ref_fa_fn), 'fasta'):
            if record.id in chrs:
                print 'select chr : {0}'.format(record.id)
                SeqIO.write(record, refchr_fa_fh, 'fasta')
            else:
                print 'skip chr : {0}'.format(record.id)
        refchr_fa_fh.close()
    #
    refchr_fai_fn = do_samtools_faidx(refchr_fa_fn)
    out = fai_to_len(refchr_fai_fn)
    refchr_len_fn = '{0}/{1}.chr.len'.format(home, name)
    refchr_len_fh = open(refchr_len_fn, 'w')
    for line in out.split('\n'):
        if line:
            refchr_len_fh.write('{0}\n'.format(line))
    refchr_len_fh.close()
    #
    return refchr_fa_fn, refchr_len_fn

def fa_split(fa, chrs, home):
    home_chr = os.path.join(home, 'chr')
    if not os.path.exists(home_chr):
        os.makedirs(home_chr)
    #
    ref_fa_split_dic = dict()
    total_chr = len(chrs)
    exist_chr = 0
    for _chr in chrs:
        split_fa_fn = os.path.join(home_chr, '{0}.fa'.format(_chr))
        ref_fa_split_dic.setdefault(_chr, split_fa_fn)
        if os.path.exists(split_fa_fn):
            exist_chr += 1
    #
    if exist_chr in [total_chr]:
        pass
    else:
        for record in SeqIO.parse(open(fa), 'fasta'):
            print 'split chr : {0}'.format(record.id)
            #
            out = open('{0}/{1}.fa'.format(home_chr,record.id),'w')
            SeqIO.write(record, out, 'fasta')
            out.close()
        #
    return ref_fa_split_dic

def do_bismark_genome_preparation(refchr_fa_fn, home):
    global bismark_build
    global bowtie2
    #
    print 'bismark genome preparation...'
    home_chr = os.path.join(home, 'chr')
    if not os.path.exists('{0}/Bisulfite_Genome'.format(home_chr)):
        cmds = list()
        cmds.append(bismark_build)
        cmds.append('--path_to_bowtie')
        cmds.append(bowtie2)
        cmds.append('--bowtie2')
        cmds.append('--verbos')
        cmds.append(home_chr)
        #
        proc = Popen(cmds, stdout=PIPE, stderr=PIPE)
        out, err = proc.communicate()
    return home_chr

class ANNO_CGI:
    def __init__(self, fa_dic, len_fn, home, dir_name):
        #
        self.fa_dic = fa_dic
        self.len_dic = self.len_dic_maker(len_fn)
        #
        self.out_dir = os.path.join(home, dir_name)
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.scr_dir = os.path.join(self.out_dir, 'script')
        if not os.path.exists(self.scr_dir):
            os.makedirs(self.scr_dir)
        self.log_dir = os.path.join(self.out_dir, 'log')
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)
        #
    def len_dic_maker(self,len_fn):
        len_dic = dict()
        for line in open(len_fn):
            items = line.rstrip('\n').split('\t')
            _chr, _len = items
            len_dic.setdefault(_chr, _len)
        return len_dic
        #
    def do_newcpgreport(self):
        global newcpgreport
        self.newcpgreport = newcpgreport
        #
        self.report_dic = dict()
        cmds = list()
        for _chr, fa_fn in self.fa_dic.iteritems():
            newcpgreport_fn = os.path.join(self.out_dir, '{0}.newcpgreport'.format(_chr))
            self.report_dic.setdefault(_chr, newcpgreport_fn)
            if os.path.exists(newcpgreport_fn):
                continue
            opts = list()
            opts.append(self.newcpgreport)
            opts.append('-sequence')
            opts.append(fa_fn)
            opts.append('-window')
            opts.append('100')
            opts.append('-shift')
            opts.append('1')
            opts.append('-minlen')
            opts.append('200')
            opts.append('-minoe')
            opts.append('0.5') # default is 0.6
            opts.append('-minpc')
            opts.append('40.0') # default is 50.0
            opts.append('-outfile')
            opts.append(newcpgreport_fn)
            cmds.append(' '.join(opts))
        if cmds:
            RunQsub(cmds, 'all.q', '1', self.log_dir, self.scr_dir, 'newcpgreport')
        #
    def split_report(self, report_fh):
        lines = report_fh.readlines()
        for idx, line in enumerate(lines):
            if line.startswith('FT'):
                s_idx = idx
                break
        blocks = [lines[s_idx].rstrip('\n').split()[3]]
        for line in lines[s_idx+1:]:
            items = line.rstrip('\n').split()
            if items[0] in ['FT'] and items[1] in ['CpG'] and items[2] in ['island']:
                yield blocks
                blocks = [items[3]]
            else:
                blocks.append(' '.join(items[1:]).lstrip('/'))
        yield blocks
        #
    def filter_cgi(self):
        self.cgi_dic = dict()
        for _chr, report_fn in self.report_dic.iteritems():
            cgi_dic = dict()
            cgi_num = 0
            for blocks in self.split_report(open(report_fn)):
                if blocks[0] in ['detected']:
                    continue # FT no islands detected
                start, end = blocks[0].split('..')
                cgi_len = int(blocks[1].split('=')[1])
                gc_count = int(blocks[2].split('=')[1])
                gc_ratio = float(blocks[3].split('=')[1])
                cgi_R = float(blocks[4].split('=')[1])
                #
                # length filter
                # ref. Nucleic Acids Res. 33 (20): e176.
                # In mammalian genomes, CpG islands are typically 300-3,000 base pairs in length, 
                # and have been found in or near approximately 40% of promoters of mammalian genes.
                #
                if cgi_len > 3000 or cgi_len < 200:
                    continue
                #
                # a GC percentage greater than 50%, and an observed-to-expected CpG ratio greater than 60%. 
                #
                if gc_ratio < 50.0:
                    continue
                if cgi_R < 0.6:
                    continue
                #
                cgi_num += 1
                cgi_dic.setdefault(cgi_num, {}).setdefault('start', start)
                cgi_dic.setdefault(cgi_num, {}).setdefault('end', end)
                cgi_dic.setdefault(cgi_num, {}).setdefault('len', cgi_len)
                cgi_dic.setdefault(cgi_num, {}).setdefault('gc_count', gc_count)
                cgi_dic.setdefault(cgi_num, {}).setdefault('gc_ratio', gc_ratio)
                cgi_dic.setdefault(cgi_num, {}).setdefault('R', cgi_R)
                #
            cgi_fn  = os.path.join(self.out_dir, '{0}.CGI'.format(_chr))
            #
            self.cgi_dic.setdefault('CGI', {}).setdefault(_chr, cgi_fn)
            #
            cgi_fh = open(cgi_fn, 'w')
            for num, info_dic in sorted(cgi_dic.iteritems()):
                cgi_fh.write('{0}\n'.format('\t'.join(
                    [_chr,str(info_dic['start']),str(info_dic['end']),str(info_dic['len']),
                        str(info_dic['gc_count']),str(info_dic['gc_ratio']),str(info_dic['R'])])))
            cgi_fh.close()
            #
    def extract_around_cgi(self, name, in_bed_dic, offset):
        # name = ['Nshelf','Nshore','Sshore','Sshelf']
        # N is North, S is South
        for _chr, fn in in_bed_dic.iteritems():
            self.cgi_dic.setdefault(name, {}).setdefault(_chr, os.path.join(self.out_dir, '{0}.{1}'.format(_chr,name)))
        #
        bed_fh_dic = dict()
        for _chr, fn in self.cgi_dic[name].iteritems():
            bed_fh_dic.setdefault(_chr, open(fn,'w'))
        #
        for _chr, in_bed_fn in in_bed_dic.iteritems():
            for line in open(in_bed_fn):
                items = line.rstrip('\n').split('\t')
                _chr = items[0]
                start = items[1]
                end = items[2]
                #
                if name.startswith('N'):
                    new_start = int(start)+int(offset)
                    if new_start < 1:
                        new_start = 1
                    new_end = int(start)-1
                    if new_end < 1:
                        new_end = 1
                elif name.startswith('S'):
                    new_start = int(end)+1
                    if new_start > int(self.len_dic[_chr]):
                        new_start = int(self.len_dic[_chr])
                    new_end = int(end)+int(offset)
                    if new_end > int(self.len_dic[_chr]):
                        new_end = int(self.len_dic[_chr])
                #
                bed_fh_dic[_chr].write('{0}\n'.format('\t'.join(
                    [str(x) for x in [_chr, new_start, new_end]])))
        for _chr, fh in bed_fh_dic.iteritems():
            fh.close()


    def stats(self,_type):
        #
        count_dic = dict()
        for _chr, fn in self.cgi_dic[_type].iteritems():
            for line in open(fn):
                items = line.rstrip('\n').split()
                count_dic.setdefault('len',[]).append(int(items[3]))
                count_dic.setdefault('GC',[]).append(float(items[5]))
                count_dic.setdefault('obs/exp',[]).append(float(items[6]))
        #
        plt.clf()
        n, bins, patches = plt.hist(count_dic['len'], log=True, bins=100)
        out = open('{0}/{1}_Length.distribution.xls'.format(self.out_dir,_type), 'w')
        for idx, _n in enumerate(n):
            out.write('{0}\n'.format('\t'.join([str(x) for x in [bins[idx+1],_n]])))
        out.close()
        plt.xlim(0,)
        plt.xlabel('Length of CpG islands({0})'.format(_type))
        plt.ylabel('log(Count)')
        plt.savefig('{0}/{1}_Length.distribution.png'.format(self.out_dir,_type))
        #
        plt.clf()
        n, bins, patches = plt.hist(count_dic['GC'], bins=50)
        out = open('{0}/{1}_GCpct.distribution.xls'.format(self.out_dir,_type), 'w')
        for idx, _n in enumerate(n):
            out.write('{0}\n'.format('\t'.join([str(x) for x in [bins[idx+1],_n]])))
        out.close()
        plt.xlim(0,100)
        plt.xlabel('GC percent of CpG islands({0})'.format(_type))
        plt.ylabel('Count')
        plt.savefig('{0}/{1}_GCpct.distribution.png'.format(self.out_dir,_type))
        #
        plt.clf()
        n, bins, patches = plt.hist(count_dic['obs/exp'], bins=50)
        out = open('{0}/{1}_ObsExp.distribution.xls'.format(self.out_dir,_type), 'w')
        for idx, _n in enumerate(n):
            out.write('{0}\n'.format('\t'.join([str(x) for x in [bins[idx+1],_n]])))
        out.close()
        plt.xlim(0,)
        plt.xlabel('CpG ratio (obs/exp) of CpG islands({0})'.format(_type))
        plt.ylabel('Count')
        plt.savefig('{0}/{1}_ObsExp.distribution.png'.format(self.out_dir,_type))
        #


class ANNO_PROMOTER:
    def __init__(self, gtf_fn, ref_fa_dic, chrs, len_fn, home, dir_name):
        #
        self.lim_n = 2000
        self.lim_s = 500
        self.feature = 'transcript'
        self.attr = 'transcript_id'
        #
        self.gtf_fn = gtf_fn
        self.ref_fa_dic = ref_fa_dic
        self.chrs = chrs
        self.len_dic = self.len_dic_maker(len_fn)
        #
        self.out_dir = os.path.join(home, dir_name)
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.scr_dir = os.path.join(self.out_dir, 'script')
        if not os.path.exists(self.scr_dir):
            os.makedirs(self.scr_dir)
        self.log_dir = os.path.join(self.out_dir, 'log')
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)
        #
    def len_dic_maker(self,len_fn):
        len_dic = dict()
        for line in open(len_fn):
            items = line.rstrip('\n').split('\t')
            _chr, _len = items
            len_dic.setdefault(_chr, _len)
        return len_dic
        #
    def find_attr_id(self, string, target):
        for attr in string.split(';'):
            if not len(attr.split()) in [2]:
                print string
                print attr
                sys.exit()
            key, _value = attr.split()
            value = _value.strip('"')
            if key in [target]:
                return value
        return None
        #
    def extract_promt(self):
        bed_dic = dict()
        for _chr in self.chrs:
            bed_dic.setdefault(_chr, os.path.join(self.out_dir, '{0}.promoter'.format(_chr)))
        #
        bed_fh_dic = dict()
        for _chr, fn in bed_dic.iteritems():
            bed_fh_dic.setdefault(_chr, open(fn,'w'))
        #
        for line in open(self.gtf_fn):
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split('\t')
            if not items[0] in self.chrs:
                continue
            if not items[2] in [self.feature]:
                continue
            #
            _chr = items[0]
            start = items[3]
            end = items[4]
            _id = self.find_attr_id(items[8], self.attr)
            if _id in [None]:
                print items[8]
                sys.exit()
            #
            p_start = int(start)-int(self.lim_n)
            if p_start < 1:
                p_start = 1
            p_end = int(end)+int(self.lim_s)
            if p_end > int(self.len_dic[_chr]):
                p_end = int(self.len_dic[_chr])
            #
            bed_fh_dic[_chr].write('{0}\n'.format('\t'.join([_chr,str(p_start),str(p_end),
                '{0}_{1}_{2}_{3}'.format(_chr,str(p_start),str(p_end),_id)])))
        for _chr, fh in bed_fh_dic.iteritems():
            fh.close()
        return bed_dic
        #
    def extract_seq(self, bed_dic, name):
        global bedtools
        #
        fa_dic = dict()
        for _chr, bed_fn in bed_dic.iteritems():
            fa_fn = '{0}.fa'.format(bed_fn)
            fa_dic.setdefault(_chr, fa_fn)
        #
        cmds = list()
        for _chr, bed_fn in bed_dic.iteritems():
            if os.path.exists(fa_dic[_chr]):
                continue
            opts = list()
            opts.append(bedtools)
            opts.append('getfasta')
            opts.append('-name')
            opts.append('-fi')
            opts.append(self.ref_fa_dic[_chr])
            opts.append('-bed')
            opts.append(bed_fn)
            opts.append('-fo')
            opts.append(fa_dic[_chr])
            cmds.append(' '.join(opts))
        if cmds:
            RunQsub(cmds, 'all.q', '1', self.log_dir, self.scr_dir, 'getfasta')
        return fa_dic

    def classify_promt(self, fa_dic):
        global bedtools
        global classify_promt
        #
        cls_fn_dic = dict()
        for _chr, fa_fn in fa_dic.iteritems():
            cls_fn_dic.setdefault(_chr,{}).setdefault('HCP', os.path.join(self.out_dir, '{0}.HCP'.format(_chr)))
            cls_fn_dic.setdefault(_chr,{}).setdefault('ICP', os.path.join(self.out_dir, '{0}.ICP'.format(_chr)))
            cls_fn_dic.setdefault(_chr,{}).setdefault('LCP', os.path.join(self.out_dir, '{0}.LCP'.format(_chr)))
        #
        cmds = list()
        for _chr, fa_fn in fa_dic.iteritems():
            if os.path.exists(cls_fn_dic[_chr]['HCP']) and \
               os.path.exists(cls_fn_dic[_chr]['ICP']) and \
               os.path.exists(cls_fn_dic[_chr]['LCP']):
                continue
            opts = list()
            opts.append('/usr/bin/python')
            opts.append(classify_promt)
            opts.append('--fa-fn')
            opts.append(fa_fn)
            opts.append('--chr-id')
            opts.append(_chr)
            opts.append('--out-dir')
            opts.append(self.out_dir)
            cmds.append(' '.join(opts))
        #
        if cmds:
            RunQsub(cmds, 'all.q', '1', self.log_dir, self.scr_dir, 'cls_promt')
        #
        cmds = list()
        for _chr, fa_fn in fa_dic.iteritems():
            if os.path.exists(cls_fn_dic[_chr]['HCP']):
                continue
            opts = list()
            opts.append(bedtools)
            opts.append('sort')
            opts.append('-i')
            opts.append('{0}.raw'.format(cls_fn_dic[_chr]['HCP']))
            opts.append('>')
            opts.append(cls_fn_dic[_chr]['HCP'])
            cmds.append(' '.join(opts))
        #
        if cmds:
            RunQsub(cmds, 'all.q', '1', self.log_dir, self.scr_dir, 'bedsort_HCP_raw')
        #
        cmds = list()
        for _chr, fa_fn in fa_dic.iteritems():
            if os.path.exists(cls_fn_dic[_chr]['ICP']):
                continue
            opts = list()
            opts.append(bedtools)
            opts.append('sort')
            opts.append('-i')
            opts.append('{0}.raw'.format(cls_fn_dic[_chr]['ICP']))
            opts.append('>')
            opts.append(cls_fn_dic[_chr]['ICP'])
            cmds.append(' '.join(opts))
        #
        if cmds:
            RunQsub(cmds, 'all.q', '1', self.log_dir, self.scr_dir, 'bedsort_ICP_raw')
        #
        cmds = list()
        for _chr, fa_fn in fa_dic.iteritems():
            if os.path.exists(cls_fn_dic[_chr]['LCP']):
                continue
            opts = list()
            opts.append(bedtools)
            opts.append('sort')
            opts.append('-i')
            opts.append('{0}.raw'.format(cls_fn_dic[_chr]['LCP']))
            opts.append('>')
            opts.append(cls_fn_dic[_chr]['LCP'])
            cmds.append(' '.join(opts))
        #
        if cmds:
            RunQsub(cmds, 'all.q', '1', self.log_dir, self.scr_dir, 'bedsort_LCP_raw')
        #
        return cls_fn_dic

    def stats(self, cls_fn_dic):
        stats_fn = os.path.join(self.out_dir, 'cgp.stats.xls')
        stats_dic = dict()
        for _chr, _type_dic in cls_fn_dic.iteritems():
            for _type, cls_fn in _type_dic.iteritems():
                count = 0
                for line in open(cls_fn):
                    count += 1
                stats_dic.setdefault(_chr, {}).setdefault(_type, str(count))

        out = open(stats_fn, 'w')
        out.write('{0}\n'.format('\t'.join(['CHR','No.HCP','No.ICP','No.LCP'])))
        for _chr in self.chrs:
            new_items = list()
            new_items.append(_chr)
            for _type in ['HCP', 'ICP', 'LCP']:
                new_items.append(stats_dic[_chr][_type])
            out.write('{0}\n'.format('\t'.join(new_items)))
        out.close()

        return stats_fn



def cg_site_finder(fa_split_dic, home, name):
    out_dir = os.path.join(home, name)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #
    cg_site_dic = dict()
    for _chr, fa_fn in fa_split_dic.iteritems():
        cg_site_fn = os.path.join(out_dir, '{0}.cgsite'.format(_chr)) # cganno
        cg_site_dic.setdefault(_chr, cg_site_fn)
    #
    for _chr, fa_fn in fa_split_dic.iteritems():
        #
        if os.path.exists(cg_site_dic[_chr]):
            continue
        #
        print 'find CG site : {0}'.format(_chr)
        out = open(cg_site_dic[_chr], 'w')
        for record in SeqIO.parse(open(fa_fn), 'fasta'):
            pre_base = 'N'
            for idx, base in enumerate(str(record.seq).upper()):
                pos = idx+1
                if '{0}{1}'.format(pre_base, base) in ['CG']:
                    out.write('{0}\t{1}\t{2}\t{3}{4}\n'.format(
                        _chr, str(pos), str(pos+1), record.seq[idx-1], record.seq[idx]))
                pre_base = base
            #
        out.close()
        #
    return cg_site_dic

def cg_site_stats(cg_site_dic, len_fn, bin_size):
    len_dic = dict()
    for line in open(len_fn):
        items = line.rstrip('\n').split('\t')
        len_dic.setdefault(items[0], int(items[1]))
    #
    cg_stats_dic = dict()
    for _chr, cg_site_fn in cg_site_dic.iteritems():
        cg_stats_dic.setdefault(_chr, '{0}.dist.png'.format(cg_site_fn))
    #
    for _chr, cg_site_fn in cg_site_dic.iteritems():
        print 'call CG distribution : {0}'.format(_chr)
        if os.path.exists(cg_stats_dic[_chr]):
            continue
        #
        count_dic = dict()
        for pos in range(1,len_dic[_chr]+1):
            count_dic.setdefault(pos, 0)
        #
        for line in open(cg_site_fn):
            items = line.rstrip('\n').split('\t')
            start = int(items[1])
            count_dic[start] = 1
        #
        bin_dic = dict()
        for poss in split_bin(1,len_dic[_chr]+1,bin_size):
            counts = list()
            for pos in poss:
                counts.append(count_dic[pos])
            bin_dic.setdefault(pos, sum(counts)/float(len(counts)))
        #
        x_s = list()
        y_s = list()
        for pos, avg in sorted(bin_dic.iteritems()):
            x_s.append(pos)
            y_s.append(avg)
        #
        plt.clf()
        f1 = plt.figure(figsize=(20,5))
        plt.plot(x_s, y_s, c='k', lw=1, ls='-')
        plt.ylabel("Avg.CpG (bin size:{0})".format(bin_size))
        plt.xlabel("Genomic position of {0}".format(_chr))
        plt.title("CpG dinucleotide distribution")
        plt.savefig(cg_stats_dic[_chr])

    return cg_stats_dic

def split_bin(start, end, bin_size):
    poss = list()
    for pos in range(start, end):
        if len(poss) in [bin_size]:
            yield poss
            poss = [pos]
        else:
            poss.append(pos)
    yield poss


class ANNO_GENOMIC_RESION(ANNO_PROMOTER):
    def __init__(self, gtf_fn, chrs, len_fn, home, dir_name):
        #
        self.feature = 'transcript'
        self.attr = 'transcript_id'
        #
        self.gtf_fn = gtf_fn
        self.chrs = chrs
        self.len_dic = self.len_dic_maker(len_fn)
        #
        self.out_dir = os.path.join(home, dir_name)
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.scr_dir = os.path.join(self.out_dir, 'script')
        if not os.path.exists(self.scr_dir):
            os.makedirs(self.scr_dir)
        self.log_dir = os.path.join(self.out_dir, 'log')
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)
        #
        self.gr_bed_dic = dict()
        #
        self.extract_region('genebody') #inter_genomic region
        self.extract_region('UP1K') #inter_genomic region
        self.extract_region('DW1K') #inter_genomic region
        #
        self.extract_region('5UTR')
        self.extract_region('3UTR')
        #self.extract_region('1EXON') # becareful
        self.extract_region('EXON')
        self.extract_region('CDS')
        #

    def extract_region(self, name):
        for _chr in self.chrs:
            self.gr_bed_dic.setdefault(name, {}).setdefault(
                    _chr, os.path.join(self.out_dir, '{0}.{1}'.format(_chr,name)))
        #
        bed_fh_dic = dict()
        for _chr, fn in self.gr_bed_dic[name].iteritems():
            bed_fh_dic.setdefault(_chr, open(fn,'w'))
        #
        for line in open(self.gtf_fn):
            if line.startswith('#'):
                continue
            #
            items = line.rstrip('\n').split('\t')
            if not items[0] in self.chrs:
                continue
            #
            _chr = items[0]
            feature = items[2]
            start = items[3]
            end = items[4]
            attribute = items[8]
            #
            if name in ['genebody'] and feature in ['transcript']:
                gr_start, gr_end, gr_id = self.extract_genebody(_chr, start, end, attribute)
            elif name in ['UP1K'] and feature in ['transcript']:
                gr_start, gr_end, gr_id = self.extract_up1k(_chr, start, end, attribute)
            elif name in ['DW1K'] and feature in ['transcript']:
                gr_start, gr_end, gr_id = self.extract_dw1k(_chr, start, end, attribute)
            elif name in ['5UTR'] and feature in ['five_prime_utr']:
                gr_start, gr_end, gr_id = self.extract_5utr(_chr, start, end, attribute)
            elif name in ['3UTR'] and feature in ['three_prime_utr']:
                gr_start, gr_end, gr_id = self.extract_3utr(_chr, start, end, attribute)
            elif name in ['1EXON'] and feature in ['exon'] and self.find_attr_id(attribute, 'exon_number') in ['1']:
                gr_start, gr_end, gr_id = self.extract_1exon(_chr, start, end, attribute)
            elif name in ['EXON'] and feature in ['exon']:
                gr_start, gr_end, gr_id = self.extract_exon(_chr, start, end, attribute)
            elif name in ['CDS'] and feature in ['CDS']:
                gr_start, gr_end, gr_id = self.extract_cds(_chr, start, end, attribute)
            else:
                continue
            #
            if gr_id in [None]:
                print 'ERROR',
                print name,
                print attribute
                sys.exit()
            #
            bed_fh_dic[_chr].write('{0}\n'.format('\t'.join([_chr,str(gr_start),str(gr_end),gr_id])))
            #
        for _chr, fh in bed_fh_dic.iteritems():
            fh.close()

    def extract_genebody(self, _chr, start, end, attribute):
        return start, end, self.find_attr_id(attribute, 'transcript_id')
    def extract_up1k(self, _chr, start, end, attribute):
        gr_start = int(start)-1000
        if gr_start < 1:
            gr_start = 1
        gr_end = int(start)-1
        if gr_end < 1:
            gr_end = 1
        return gr_start, gr_end, self.find_attr_id(attribute, 'transcript_id')
    def extract_dw1k(self, _chr, start, end, attribute):
        gr_start = int(end)+1
        if gr_start > int(self.len_dic[_chr]):
            gr_start = int(self.len_dic[_chr])
        gr_end = int(end)+1000
        if gr_end > int(self.len_dic[_chr]):
            gr_end = int(self.len_dic[_chr])
        return gr_start, gr_end, self.find_attr_id(attribute, 'transcript_id')
    def extract_5utr(self, _chr, start, end, attribute):
        return start, end, self.find_attr_id(attribute, 'transcript_id')
    def extract_3utr(self, _chr, start, end, attribute):
        return start, end, self.find_attr_id(attribute, 'transcript_id')
    def extract_1exon(self, _chr, start, end, attribute):
        return start, end, self.find_attr_id(attribute, 'transcript_id')
    def extract_exon(self, _chr, start, end, attribute):
        return start, end, self.find_attr_id(attribute, 'transcript_id')
    def extract_cds(self, _chr, start, end, attribute):
        return start, end, self.find_attr_id(attribute, 'transcript_id')

def main(args):
    if not os.path.exists(args.home):
        os.makedirs(args.home)
    ##### Genome prepraration #####
    chrs = select_chrs(args.ref)
    print 'selected chromosomes : {0}'.format(chrs)
    refchr_fa_fn, refchr_len_fn = fa_filter_by_chr(args.ref, chrs, args.home, args.name)
    print 'selected chromosome fasta  : {0}'.format(refchr_fa_fn)
    print 'selected chromosome length : {0}'.format(refchr_len_fn)
    #
    bismark_dir = do_bismark_genome_preparation(refchr_fa_fn, args.home)
    print 'bismark <genome_folder> : {0}'.format(bismark_dir)
    #
    ref_fa_split_dic = fa_split(refchr_fa_fn, chrs, args.home)
    for chr_id, refchr_fa_fn in ref_fa_split_dic.iteritems():
        print 'splited fasta by selected chromosome : {0} ==> {1}'.format(chr_id, refchr_fa_fn)
    #
    cg_site_dic = cg_site_finder(ref_fa_split_dic, args.home, 'anno_CpG')
    for chr_id, cg_site_fn in cg_site_dic.iteritems():
        print 'find cgsite : {0} ==> {1}'.format(chr_id, cg_site_fn)
    cg_stats_dic = cg_site_stats(cg_site_dic, refchr_len_fn, 100000)
    for chr_id, cg_stats_fn in cg_stats_dic.iteritems():
        print 'cgsite distribution : {0} ==> {1}'.format(chr_id, cg_stats_fn)
    #
    ##### Gene coordinate prepraration #####
    # CpG islands
    anno_cgi = ANNO_CGI(ref_fa_split_dic, refchr_len_fn, args.home, 'anno_CGI')
    anno_cgi.do_newcpgreport()
    anno_cgi.filter_cgi()
    anno_cgi.extract_around_cgi('Nshore', anno_cgi.cgi_dic['CGI'], '-2000')
    anno_cgi.extract_around_cgi('Nshelf', anno_cgi.cgi_dic['Nshore'], '-2000')
    anno_cgi.extract_around_cgi('Sshore', anno_cgi.cgi_dic['CGI'], '2000')
    anno_cgi.extract_around_cgi('Sshelf', anno_cgi.cgi_dic['Sshore'], '2000')
    for cgi_type, _dic in anno_cgi.cgi_dic.iteritems():
        for chr_id, cgi_fn in _dic.iteritems():
            print 'find CpG island : {0} ==> {1},{2}'.format(cgi_type, chr_id, cgi_fn)
    #
    anno_cgi.stats('CGI')
    # promoters, HCP, ICP, LCP
    anno_promt = ANNO_PROMOTER(args.gtf, ref_fa_split_dic, chrs, refchr_len_fn, args.home, 'anno_Promoter')
    promt_bed_dic = anno_promt.extract_promt()
    for chr_id, promt_bed_fn in promt_bed_dic.iteritems():
        print 'Promoter region : {0} ==> {1}'.format(chr_id, promt_bed_fn)
    promt_fa_dic = anno_promt.extract_seq(promt_bed_dic, 'promoter')
    promt_cls_dic = anno_promt.classify_promt(promt_fa_dic)
    for chr_id, _type_dic in promt_cls_dic.iteritems():
        for _type, cls_fn in _type_dic.iteritems():
            print 'classification Promoter CpG : {0} ==> {1}, {2}'.format(chr_id, _type, cls_fn)
    promt_stats_fn = anno_promt.stats(promt_cls_dic)
    print 'Stats of classification Promoter CpG : {0}'.format(promt_stats_fn)
    # Interesting Genomic Resion
    anno_gr = ANNO_GENOMIC_RESION(args.gtf, chrs, refchr_len_fn, args.home, 'anno_GenomicRegion')
    for gr_type, _dic in anno_gr.gr_bed_dic.iteritems():
        for chr_id, gr_fn in _dic.iteritems():
            print 'Genomic region : {0} ==> {1}, {2}'.format(gr_type, chr_id, gr_fn)
    #
    ##### Merge Annotations #####










if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--home',
            default='/BiO/BioPeople/siyoo/MINE/ref/Chlorocebus_sabaeus/ENS93')
    parser.add_argument('-r', '--ref',
            default='Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa')
    parser.add_argument('-g', '--gtf',
            default='Chlorocebus_sabaeus.ChlSab1.1.93.gtf')
    parser.add_argument('-n', '--name',
            default='C_sabaeus_ENS93')
    args = parser.parse_args()
    main(args)
