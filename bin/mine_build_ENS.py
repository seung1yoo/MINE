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
import json
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
                #print 'select chr : {0}'.format(record.id)
                SeqIO.write(record, refchr_fa_fh, 'fasta')
            else:
                #print 'skip chr : {0}'.format(record.id)
                pass
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
            #print 'split chr : {0}'.format(record.id)
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
    if not os.path.exists(home_chr):
        os.makedirs(home_chr)
    if not os.path.exists('{0}/Bisulfite_Genome'.format(home_chr)):
        cmds = list()
        cmds.append(bismark_build)
        cmds.append('--path_to_bowtie')
        cmds.append(bowtie2)
        cmds.append('--bowtie2')
        cmds.append('--verbos')
        cmds.append(home_chr)
        #
        print cmds
        proc = Popen(cmds, stdout=PIPE, stderr=PIPE)
        out, err = proc.communicate()
        print out
        print err
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
                    [_chr,
                     str(info_dic['start']),
                     str(info_dic['end']),
                     '{0}_{1}_{2}_CGI_CHR{0}CGI{3:0>10}|{4}|{5}|{6}|{7}'.format(
                         _chr,
                         str(info_dic['start']),
                         str(info_dic['end']),
                         str(num),
                         str(info_dic['len']),
                         str(info_dic['gc_count']),
                         str(info_dic['gc_ratio']),
                         str(info_dic['R']))
                     ]
                     )
                     )
                     )
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
                attr = items[3]
                cgi_id = attr.split('_')[-1].split('|')[0]
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
                    [str(x) for x in [
                        _chr,
                        new_start,
                        new_end,
                        '{0}_{1}_{2}_{3}_{4}'.format(_chr,str(new_start),str(new_end),name,cgi_id)
                        ]
                    ]
                    )))
        for _chr, fh in bed_fh_dic.iteritems():
            fh.close()
        #
    def stats(self,_type):
        #
        count_dic = dict()
        for _chr, fn in self.cgi_dic[_type].iteritems():
            for line in open(fn):
                items = line.rstrip('\n').split('\t')
                units = items[3].split('_')[-1].split('|')
                count_dic.setdefault('len',[]).append(int(units[1]))
                count_dic.setdefault('GC',[]).append(float(units[3]))
                count_dic.setdefault('obs/exp',[]).append(float(units[4]))
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
        cgi_stats_dic = dict()
        cgi_stats_dic.setdefault('len','{0}/{1}_Length.distribution.png'.format(self.out_dir,_type))
        cgi_stats_dic.setdefault('GC','{0}/{1}_GCpct.distribution.png'.format(self.out_dir,_type))
        cgi_stats_dic.setdefault('obs/exp','{0}/{1}_ObsExp.distribution.png'.format(self.out_dir,_type))
        #
        return cgi_stats_dic


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
                '{0}_{1}_{2}_Promoter_{3}'.format(_chr,str(p_start),str(p_end),_id)])))
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
        self.bed_sort()

    def bed_sort(self):
        global bedtools
        #
        cmds = list()
        for name, chr_dic in self.gr_bed_dic.iteritems():
            for _chr, fn in chr_dic.iteritems():
                if os.path.exists(fn):
                    continue
                opts = [bedtools]
                opts.append('sort')
                opts.append('-i')
                opts.append('{0}.raw'.format(fn))
                opts.append('>')
                opts.append(fn)
                cmds.append(' '.join(opts))
            #
        if cmds:
            RunQsub(cmds, 'all.q', '1', self.log_dir, self.scr_dir, 'bedsort.gr')

    def extract_region(self, name):
        for _chr in self.chrs:
            self.gr_bed_dic.setdefault(name, {}).setdefault(
                    _chr, os.path.join(self.out_dir, '{0}.{1}'.format(_chr,name)))
        #
        bed_fh_dic = dict()
        for _chr, fn in self.gr_bed_dic[name].iteritems():
            bed_fh_dic.setdefault(_chr, open('{0}.raw'.format(fn),'w'))
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
            bed_fh_dic[_chr].write('{0}\n'.format('\t'.join([_chr,str(gr_start),str(gr_end),
                '{0}'.format('_'.join([_chr,str(gr_start),str(gr_end),name,gr_id]))])))
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

class ANNO_MERGE:
    def __init__(self, cg_dic, cgi_dic, promt_dic, promt_cls_dic, gr_dic, home, dir_name):
        #
        global bedtools
        self.bedtools = bedtools
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
        self.merge_dic = dict()
        #
        self.add_dic(cg_dic, 'CpG')
        self.add_dic(cgi_dic['CGI'], 'CGI')
        self.add_dic(cgi_dic['Nshore'], 'Nshore')
        self.add_dic(cgi_dic['Nshelf'], 'Nshelf')
        self.add_dic(cgi_dic['Sshore'], 'Sshore')
        self.add_dic(cgi_dic['Sshelf'], 'Sshelf')
        self.add_dic(promt_dic, 'Promoter')
        self.add_dic2(promt_cls_dic, 'Promoter')
        self.add_dic(gr_dic['genebody'], 'genebody')
        self.add_dic(gr_dic['UP1K'], 'UP1K')
        self.add_dic(gr_dic['DW1K'], 'DW1K')
        self.add_dic(gr_dic['5UTR'], '5UTR')
        self.add_dic(gr_dic['3UTR'], '3UTR')
        self.add_dic(gr_dic['EXON'], 'EXON')
        self.add_dic(gr_dic['CDS'], 'CDS')
        #
        for _chr, name_dic in self.merge_dic.iteritems():
            for name, a_fn in name_dic.iteritems():
                #print _chr, name, a_fn
                pass
        #
        self.intersect_dic = dict()
        #
        for _chr, name_dic in self.merge_dic.iteritems():
            self.exe_intersectbed(_chr, name_dic)
        #
        for _chr, pair_dic in self.intersect_dic.iteritems():
            for pair, a_fn in pair_dic.iteritems():
                #print _chr, pair, a_fn
                pass
        #
        self.annoTable_dic = dict()
        #
        self.anno_cpg()
        #
        self.merge_annoTable()

    def merge_annoTable(self):
        annoTable_fn = os.path.join(self.out_dir, 'all.annoTable')
        self.annoTable_dic.setdefault('all', annoTable_fn)
        if os.path.exists(annoTable_fn):
            pass

    def anno_cpg(self):
        #### pair list-up
        ### Interesting Genomic Regions
        ## Inter-genomic region
        #CpG_UP1K
        #CpG_genebody
        #CpG_DW1K
        ## Intra-genomic region
        #CpG_5UTR
        #CpG_CDS
        #CpG_EXON
        #CpG_3UTR
        ## Promoter region
        #CpG_Promoter
        #CpG_HCP
        #CpG_ICP
        #CpG_LCP
        ## CpG island region
        #CpG_Nshelf
        #CpG_Nshore
        #CpG_CGI
        #CpG_Sshore
        #CpG_Sshelf
        #
        regions = ['UP1K','genebody','DW1K',
                   '5UTR','CDS','EXON','3UTR',
                   'Promoter','HCP','ICP','LCP',
                   'Nshelf','Nshore','CGI','Sshore','Sshelf']
        #
        for _chr, pair_dic in self.intersect_dic.iteritems():
            print 'ANNO_MERGE-ing... {0}'.format(_chr)
            #
            annoTable_fn = os.path.join(self.out_dir, '{0}.annoTable'.format(_chr))
            self.annoTable_dic.setdefault(_chr, annoTable_fn)
            if os.path.exists(annoTable_fn):
                continue
            #
            anno_dic = self.init_anno_dic(_chr, regions)
            #
            for region in regions:
                intersect_fn = pair_dic['CpG_{0}'.format(region)]
                for line in open(intersect_fn):
                    items = line.rstrip('\n').split('\t')
                    s_pos = items[1]
                    #
                    if region in ['3UTR','5UTR','CDS','EXON','UP1K','genebody','DW1K']:
                        _id = items[7]
                    elif region in ['Promoter','HCP','ICP','LCP']:
                        _id = items[7]
                    elif region in ['CGI','Nshelf','Nshore','Sshelf','Sshore']:
                        _id = items[7]
                    #
                    if not _id in anno_dic[s_pos][region]:
                        anno_dic[s_pos][region].append(_id)
                #
            #
            out = open(annoTable_fn, 'w')
            new_titles = self.make_annoTable_titles(regions)
            out.write('#{0}\n'.format('\t'.join(new_titles)))
            for line in open(self.merge_dic[_chr]['CpG']):
                items = line.rstrip('\n').split('\t')
                s_pos = items[1]
                e_pos = items[2]
                _type = items[3]
                #
                new_items = [_chr,s_pos,e_pos,_type]
                _count = 0
                for region in regions:
                    _count += len(anno_dic[s_pos][region])
                    new_items.append(str(len(anno_dic[s_pos][region])))
                for region in regions:
                    if anno_dic[s_pos][region]:
                        new_items.append(','.join((anno_dic[s_pos][region])))
                    else:
                        new_items.append('-')
                #
                if not _count in [0]:
                    out.write('{0}\n'.format('\t'.join(new_items)))
            out.close()
            #
    def make_annoTable_titles(self, regions):
        new_titles = ['CHR','START','END','TYPE']
        for region in regions:
            new_titles.append('NUM.{0}'.format(region))
        for region in regions:
            new_titles.append('MEMBER.{0}'.format(region))
        return new_titles

    def init_anno_dic(self, _chr, regions):
        anno_dic = dict()
        for line in open(self.merge_dic[_chr]['CpG']):
            items = line.rstrip('\n').split('\t')
            s_pos = items[1]
            for region in regions:
                anno_dic.setdefault(s_pos, {}).setdefault(region, [])
        return anno_dic

    def exe_intersectbed(self, _chr, name_dic):
        #
        for name, a_fn in name_dic.iteritems():
            pair = 'CpG_{0}'.format(name)
            intersect_fn = os.path.join(self.out_dir,'{0}.intersect.{1}'.format(_chr, pair))
            self.intersect_dic.setdefault(_chr, {}).setdefault(pair, intersect_fn)
        #
        cmds = list()
        for name, a_fn in name_dic.iteritems():
            pair = 'CpG_{0}'.format(name)
            if name in ['CpG']:
                continue
            if os.path.exists(self.intersect_dic[_chr][pair]):
                continue
            opts = [self.bedtools]
            opts.append('intersect')
            opts.append('-a')
            opts.append(name_dic['CpG'])
            opts.append('-b')
            opts.append(a_fn)
            opts.append('-wa')
            opts.append('-wb')
            opts.append('-sorted')
            opts.append('>')
            opts.append(self.intersect_dic[_chr][pair])
            cmds.append(' '.join(opts))
        if cmds:
            RunQsub(cmds, 'all.q', '1', self.log_dir, self.scr_dir, 'intersect.{0}'.format(_chr))

    def add_dic(self, a_dic, name):
        for _chr, a_fn in a_dic.iteritems():
            if os.path.exists(a_fn):
                self.merge_dic.setdefault(_chr, {}).setdefault(name, a_fn)
            else:
                print 'not exists {0}'.format(a_fn)
                sys.exit()

    def add_dic2(self, a_dic, name):
        for _chr, type_dic in a_dic.iteritems():
            for _type, a_fn in type_dic.iteritems():
                if os.path.exists(a_fn):
                    self.merge_dic.setdefault(_chr, {}).setdefault(_type, a_fn)
                else:
                    print 'not exists {0}'.format(a_fn)
                    sys.exit()

class CONF_MAKER:
    def __init__(self, args, refchr_len_fn, bismark_dir, cg_stats_dic, cgi_stats_dic, promt_stats_fn, annoTable_dic):
        self.home = args.home
        self.name = args.name
        #
        self.species_name = args.species_name
        self.species_alias = args.species_alias
        self.assembly_accession = args.assembly_accession
        self.geneset_version = args.geneset_version
        #
        self.refchr_len_fn = refchr_len_fn
        self.bismark_dir = bismark_dir
        self.cg_stats_dic = cg_stats_dic # {chr : png}
        self.cgi_len_png = cgi_stats_dic['len']
        self.cgi_gc_png = cgi_stats_dic['GC']
        self.cgi_obsexp_png = cgi_stats_dic['obs/exp']
        self.promt_stats_fn = promt_stats_fn
        self.annoTable_dic = annoTable_dic
        #
        #####
        #
        self.conf_dic = dict()
        self.conf_dic.setdefault(self.name, {})
        #
        self.add_str('HOME',self.home)
        self.add_str('SPECIES_NAME',self.species_name)
        self.add_str('SPECIES_ALIAS',self.species_alias)
        self.add_str('ASSEMBLY_ACC',self.assembly_accession)
        self.add_str('GENESET_VER',self.geneset_version)
        #
        self.add_fn('CHR_LEN_FN',self.refchr_len_fn)
        self.add_dir('BISMARK_DIR',self.bismark_dir)
        self.add_dic('CGSITE_PNG',self.cg_stats_dic)
        self.add_fn('CGI_LEN_PNG',self.cgi_len_png)
        self.add_fn('CGI_GC_PNG',self.cgi_gc_png)
        self.add_fn('CGI_OBSEXP_PNG',self.cgi_obsexp_png)
        self.add_fn('PROMTER_STATS_FN',self.promt_stats_fn)
        self.add_dic('CGSITE_ANNO',self.annoTable_dic)
        #
        self.conf_fn = os.path.join(self.home,'ref_{0}.json'.format(self.name))
        conf_fh = open(self.conf_fn,'w')
        json.dump(self.conf_dic,conf_fh,indent=4)
        conf_fh.write('\n')
        conf_fh.close()
        #
    def add_str(self, key, value):
        self.conf_dic[self.name].setdefault(key, value)
    def add_fn(self, key, value):
        self.conf_dic[self.name].setdefault(key, value.replace(self.home, '[HOME]'))
    def add_dir(self, key, value):
        self.conf_dic[self.name].setdefault(key, value.replace(self.home, '[HOME]'))
    def add_dic(self, key, a_dic):
        for _chr, fn in a_dic.iteritems():
            self.conf_dic[self.name].setdefault(key, {})
            self.conf_dic[self.name][key].setdefault(_chr,fn.replace(self.home, '[HOME]'))
        #
    #



def main(args):
    ##### Make Reference HOME #####
    print 'home PATH : {0}'.format(args.home)
    if not os.path.exists(args.home):
        os.makedirs(args.home)
    ##### Genome prepraration #####
    chrs = select_chrs(args.ref)
    print 'selected chromosomes : {0}'.format(chrs)
    refchr_fa_fn, refchr_len_fn = fa_filter_by_chr(args.ref, chrs, args.home, args.name)
    print 'selected chromosome fasta  : {0}'.format(refchr_fa_fn)
    print 'selected chromosome length : {0}'.format(refchr_len_fn)
    #
    print 'split fasta...'
    ref_fa_split_dic = fa_split(refchr_fa_fn, chrs, args.home)
    for chr_id, refchr_fa_fn in ref_fa_split_dic.iteritems():
        print 'splited fasta by selected chromosome : {0} ==> {1}'.format(chr_id, refchr_fa_fn)
    #
    print 'bismark building...'
    bismark_dir = do_bismark_genome_preparation(refchr_fa_fn, args.home)
    print 'bismark <genome_folder> : {0}'.format(bismark_dir)
    #
    print 'find CpG site...'
    cg_site_dic = cg_site_finder(ref_fa_split_dic, args.home, 'anno_CpG')
    for chr_id, cg_site_fn in cg_site_dic.iteritems():
        print 'CpG site : {0} ==> {1}'.format(chr_id, cg_site_fn)
    cg_stats_dic = cg_site_stats(cg_site_dic, refchr_len_fn, 100000)
    for chr_id, cg_stats_fn in cg_stats_dic.iteritems():
        print 'CpG distribution : {0} ==> {1}'.format(chr_id, cg_stats_fn)
    #
    ##### Gene coordinate prepraration #####
    # CpG islands
    print 'annotation CpG islands... (CGI)'
    anno_cgi = ANNO_CGI(ref_fa_split_dic, refchr_len_fn, args.home, 'anno_CGI')
    anno_cgi.do_newcpgreport()
    anno_cgi.filter_cgi() # always re-run
    for chr_id, fn in anno_cgi.cgi_dic['CGI'].iteritems():
        print 'CGI : {0} ==> {1}'.format(chr_id, fn)
    #
    print 'annotation around CpG islands... (Nshore, Nshelf, Sshore, Sshelf)'
    anno_cgi.extract_around_cgi('Nshore', anno_cgi.cgi_dic['CGI'], '-2000') # always re-run
    anno_cgi.extract_around_cgi('Nshelf', anno_cgi.cgi_dic['Nshore'], '-2000') # always re-run
    anno_cgi.extract_around_cgi('Sshore', anno_cgi.cgi_dic['CGI'], '2000') # always re-run
    anno_cgi.extract_around_cgi('Sshelf', anno_cgi.cgi_dic['Sshore'], '2000') # always re-run
    #
    for chr_id, fn in anno_cgi.cgi_dic['Nshore'].iteritems():
        print 'Nshore : {0} ==> {1}'.format(chr_id, fn)
    for chr_id, fn in anno_cgi.cgi_dic['Nshelf'].iteritems():
        print 'Nshelf : {0} ==> {1}'.format(chr_id, fn)
    for chr_id, fn in anno_cgi.cgi_dic['Sshore'].iteritems():
        print 'Sshore : {0} ==> {1}'.format(chr_id, fn)
    for chr_id, fn in anno_cgi.cgi_dic['Sshelf'].iteritems():
        print 'Sshelf : {0} ==> {1}'.format(chr_id, fn)
    #
    cgi_stats_dic = anno_cgi.stats('CGI') # always re-run
    print 'Stats of CpG island : {0}'.format(cgi_stats_dic['len'])
    print 'Stats of CpG island : {0}'.format(cgi_stats_dic['GC'])
    print 'Stats of CpG island : {0}'.format(cgi_stats_dic['obs/exp'])
    # promoters, HCP, ICP, LCP
    print 'annotation promoter region...'
    anno_promt = ANNO_PROMOTER(args.gtf, ref_fa_split_dic, chrs, refchr_len_fn, args.home, 'anno_Promoter')
    print 'extract promoter region...'
    promt_bed_dic = anno_promt.extract_promt() # always re-run
    for chr_id, promt_bed_fn in promt_bed_dic.iteritems():
        print 'Promoter region : {0} ==> {1}'.format(chr_id, promt_bed_fn)
    #
    print 'extract promoter sequence...'
    promt_fa_dic = anno_promt.extract_seq(promt_bed_dic, 'promoter')
    #
    print 'classify promoter region...'
    promt_cls_dic = anno_promt.classify_promt(promt_fa_dic)
    for chr_id, _type_dic in promt_cls_dic.iteritems():
        for _type, cls_fn in _type_dic.iteritems():
            print 'classification Promoter CpG : {0} ==> {1}, {2}'.format(chr_id, _type, cls_fn)
    #
    print 'stats promoter region...'
    promt_stats_fn = anno_promt.stats(promt_cls_dic) # always re-run
    print 'Stats of classification Promoter CpG : {0}'.format(promt_stats_fn)
    #
    # Interesting Genomic Region
    print 'annotation genomic region...'
    anno_gr = ANNO_GENOMIC_RESION(args.gtf, chrs, refchr_len_fn, args.home, 'anno_GenomicRegion')
    for gr_type, _dic in anno_gr.gr_bed_dic.iteritems():
        for chr_id, gr_fn in _dic.iteritems():
            print 'Genomic region : {0} ==> {1}, {2}'.format(gr_type, chr_id, gr_fn)
    #
    ##### Merge Annotations #####
    print 'merge annotations...'
    anno_merge = ANNO_MERGE(cg_site_dic, anno_cgi.cgi_dic, promt_bed_dic, promt_cls_dic, anno_gr.gr_bed_dic, args.home, 'anno_Merge')
    for chr_id, annoTable_fn in anno_merge.annoTable_dic.iteritems():
        print chr_id, annoTable_fn
    #
    ##### Configure MINE refs.json ####
    conf_obj = CONF_MAKER(args,
                          refchr_len_fn,
                          bismark_dir,
                          cg_stats_dic,
                          cgi_stats_dic,
                          promt_stats_fn,
                          anno_merge.annoTable_dic)
    print conf_obj.conf_fn
    print '#DONE'

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
    #
    parser.add_argument('--species-name', default='Chlorocebus sabaeus', help='Chlorocebus sabaeus (green monkey)')
    parser.add_argument('--species-alias', default='ChlSab1.1', help='ChlSab1.1')
    parser.add_argument('--assembly-accession', default='GCA_000409795.2', help='NCBI(GeneBank):GCA_000409795.2')
    parser.add_argument('--geneset-version', default='Ensembl 93', help='Ensembl 93')
    args = parser.parse_args()
    main(args)

