#!/usr/bin/python

import os
import sys
#import cProfile, pstats, io
#import bz2
import gzip 


N_SAMPLE_LINE_S=50000

SANGER_QS_ENC_START=33
SCORING_TYPE_S=(
        ('Sanger',33,73,SANGER_QS_ENC_START),
        ('Illumina 1.3+/1.5+',59,104,64),
        ('Solexa',64,104,64),
#        ('Illumina 1.8+',33,74,33), # 20150916 by brandon ,, hiseq 4000 scale max=75
#        ('Illumina 1.8+',33,75,33),
        ('Illumina 1.8+',33,108,33),  #20151025 by kjh ,, data scale
        )

# Quality Score Type
TYPE_SANGER_SCORE = 1
START_SANGER_SCORE= 33
#
TYPE_ILLUMINA_1_0 = 2
START_ILLUMINA_1_0= 59+5
#
TYPE_ILLUMINA_1_3 = 3
START_ILLUMINA_1_3= 64

class Info :
    def __init__(self):
        self.score_type=None
        self.fq_1=None
        self.fq_2=None
        self.n_ratio=None
        self.q_cutoff=None
        self.q_ratio=None
        self._log=None
        self.fix_score=False
    def __del__(self):
        if self._log!=None  :
            self._log.close()
    def c_name(self) :
        buf=[]
        for i_ch in range(len(self.fq_1.out_fn)) :
            ch_1=self.fq_1.out_fn[i_ch]
            ch_2=self.fq_2.out_fn[i_ch]
            if ch_1==ch_2 :
                buf.append(ch_1)
            else :
                break
        fn="".join(buf)
        if fn.endswith("_") :
            fn=fn[:-1]
        return fn
    def report(self):
        return "%s.log"%self.c_name()
    def log_err(self,msg):
        if self._log==None :
            self._log=file("%s.err"%self.c_name(),'wt')
            self._log.write("\n")
            self._log.write("#\n")
            self._log.write("# FILTER_N_RATION\t%.2f\n"%self.n_ratio)
            self._log.write("# FILTER_Q_CUTOFF\t%d\n"%self.q_cutoff)
            self._log.write("# FILTER_Q_RATIO\t%.1f\n"%self.q_ratio)
            self._log.write("#\n")
            self._log.write("# NAME")
            self._log.write("\tN_RATIO_1\tQ_MEAN_1\tQ_RATIO_1")
            self._log.write("\tN_RATIO_2\tQ_MEAN_2\tQ_RATIO_2")
            self._log.write("\n")
        self._log.write("%s\n"%msg)
    def _collect_sample_qscore(self,fn) :
        qs_s=set()
        #
        _in=None
        if fn.endswith(".gz") :
            _in=os.popen("zcat %s"%fn)
            #_in=gzip.open('%s'%fn, 'rb') # 20140924, jinsilhanna,hmkim
#        elif fn.endswith(".bz2") or fn.endswith(".bz") :
#            _in=bz2.BZ2File(fn)
        else :
            _in=file(fn)
        #
        n_line=0
        while True :
            line=_in.readline()
            if line=='' :
                break
            _in.readline()
            _in.readline()
            q_line=_in.readline()
            for ch in q_line.strip() :
                qs_s.add(ord(ch))
            n_line+=1
            if n_line>=N_SAMPLE_LINE_S :
                break
        #
        _in.close()
        return qs_s
    def calibrate(self) :
        self.qscore_s=list(self._collect_sample_qscore(self.fq_1.in_fn))
        self.qscore_s.sort()
        #
        is_selected=False
        for score_type in SCORING_TYPE_S :
            #print self.qscore_s
            #print score_type
            name,min_qs,max_qs,zero_qs=score_type
            if len(self.qscore_s)==0 :
                continue
            if min(self.qscore_s)<min_qs :
                continue
            if max(self.qscore_s)>max_qs :
                continue
            self.score_type=score_type
            is_selected=True
            break
        #
        if is_selected :
            return True
        return False

def read_pe_read(info,fq_1,fq_2) :
    in_1=info.fq_1.get_in()
    in_2=info.fq_2.get_in()
    while True :
        read_1=SeqRead(info.score_type[3],in_1,info.fix_score)
        read_2=SeqRead(info.score_type[3],in_2,info.fix_score)
        if read_1.title=='' or  read_2.title=='' :
            break
        if not read_1.title.startswith("@") :
            continue
        yield read_1,read_2

class InFq :
    def __init__(self,in_fn,out_fn):
        self.in_fn=in_fn
        self.out_fn=out_fn
        
        if self.out_fn.endswith(".gz") :
            self._out=file(self.out_fn[:-3],'wt')
            #self._out=gzip.open(self.out_fn, 'wb') # 20140924 jinsilhanna,hmkim
        else :
            self._out=file(self.out_fn,'wt')
    def __del__(self):
        self._out.close()
        if self.out_fn.endswith(".gz") :
            os.system("gzip -f %s"%(self.out_fn[:-3]))
    def is_finished(self):
        if os.path.exists(self.out_fn) and os.path.getsize(self.out_fn)>50000 :
            return True
        return False
    def get_in(self):
        _in=None
        if self.in_fn.endswith(".gz") :
            _in=os.popen("zcat %s"%(self.in_fn))
#        elif self.in_fn.endswith(".bz2") or  self.in_fn.endswith(".bz") :
#            _in=bz2.BZ2File(self.in_fn)
        else :
            _in=file(self.in_fn)
        return _in
    def out_write(self,msg) :
        self._out.write(msg)

class SeqRead :
    def __init__(self,sc_start,_in,fix_score) :
        self.sc_start=sc_start
        self.title=_in.readline().strip()
        self.seq=_in.readline().strip()
        self.comment=_in.readline().strip()
        self.score=_in.readline().strip()
        self._q_score=[]
        self._q_score_norm=[]
        for sc in self.score :
            self._q_score.append(ord(sc)-self.sc_start)
            self._q_score_norm.append(min(40,ord(sc)-self.sc_start))
        #
        sub=[]
        for _qs in self._q_score_norm :
            sub.append(chr(_qs+self.sc_start))
        self.score_norm="".join(sub)
        self.fix_score=fix_score
        #
        self._trim_start=None
        self._trim_end=None
    def calc_n_ratio(self) :
        n_n=0
        for n in self.seq :
            if n=='N' :
                n_n+=1
        ratio=1.0*n_n/len(self.seq)
        return ratio
    def calc_q_score(self,cutoff) :
        _sum=0.0
        n_incorrect=0
        for q in self._q_score :
            if q<=cutoff :
                n_incorrect+=1
            _sum+=q
        q_mean=1.0*_sum/len(self.score)
        q_ratio=1.0*n_incorrect/len(self.score)
        return q_mean, q_ratio
    def __str__(self):
        buf=[]
        buf.append("%s\n"%self.title)
        if self._trim_start!=None :
            buf.append("%s\n"%self.seq[self._trim_start:self._trim_end+1])
        else :
            buf.append("%s\n"%self.seq)
        buf.append("%s\n"%self.comment)
        #
        if self.fix_score :
            if self._trim_end!=None :
                buf.append("%s\n"%self.score_norm[self._trim_start:self._trim_end+1])
            else :
                buf.append("%s\n"%self.score_norm)
        else :
            if self._trim_end!=None :
                buf.append("%s\n"%self.score[self._trim_start:self._trim_end+1])
            else :
                buf.append("%s\n"%self.score)
        #
        return "".join(buf)
    def __len__(self):
        if self._trim_start!=None :
            return self._trim_end-self._trim_start+1
        return len(self.seq)

def get_sc_start(score_type) :
    sc_start=None
    if score_type==TYPE_SANGER_SCORE :
        sc_start=START_SANGER_SCORE
    elif score_type==TYPE_ILLUMINA_1_0 :
        sc_start=START_ILLUMINA_1_0
    elif score_type==TYPE_ILLUMINA_1_3 :
        sc_start=START_ILLUMINA_1_3
    else :
        sys.stderr.write("ERROR: Unknown Score type: %s\n"%score_type)
        return
    return sc_start

def trimming(info,read) :
    start=0
    end=len(read.score)-1
    for idx in range(len(read.score)) :
        if (read._q_score[idx]>=info.trimming) :
            break
        start=idx
    for idx in range(len(read.score)-1,0,-1) :
        if (read._q_score[idx]>=info.trimming) :
            break
        end=idx-1
    if start==0 and (end==len(read.score)-1) :
        return
    read._trim_start=start
    read._trim_end=end
    return "TRIMMING [%d:%d]"%(start+1,end+1)

def run_pe(info) :
    removed_pair_order_s=[]
    if not (info.fq_1.is_finished() and info.fq_2.is_finished()) :
        n_all=0
        n_corr=0
        n_incorr=0
        n_trimming=0
        n_len_s=[]
        #
        for read_1,read_2 in read_pe_read(info,info.fq_1,info.fq_2) :

            if len(read_1)==0 or len(read_2)==0 :
                msg="Zero length sequence (%s)"%read_1.title.strip()
                sys.stdout.write("%s\n"%msg)
                continue

            n_ratio_1=read_1.calc_n_ratio()
            n_ratio_2=read_2.calc_n_ratio()

            q_mean_1,q_ratio_1=read_1.calc_q_score(info.q_cutoff)
            q_mean_2,q_ratio_2=read_2.calc_q_score(info.q_cutoff)

            is_correct=True
            error_type=None
            if (q_ratio_1 >= info.q_ratio) :
                is_correct=False
                error_type="Q_RATIO_1"
            elif (q_ratio_2 >= info.q_ratio) :
                is_correct=False
                error_type="Q_RATIO_2"
            elif (q_mean_1 < info.q_cutoff) :
                is_correct=False
                error_type="Q_MEAN_1"
            elif (q_mean_2 < info.q_cutoff) :
                is_correct=False
                error_type="Q_MEAN_2"
            elif (n_ratio_1>=info.n_ratio) :
                is_correct=False
                error_type="N_RATIO_1"
            elif (n_ratio_2>=info.n_ratio) :
                is_correct=False
                error_type="N_RATIO_2"
            #
            n_all+=2
            if not is_correct :
                n_incorr+=2
                buf=[]
                buf.append(read_1.title.strip())
                buf.append("%s"%error_type)
                buf.append("%8.3f\t%8.3f\t%8.3f"%(n_ratio_1,q_mean_1,q_ratio_1))
                buf.append("%8.3f\t%8.3f\t%8.3f"%(n_ratio_2,q_mean_2,q_ratio_2))
                info.log_err("\t".join(buf))
                continue
            #
            if info.trimming>0 :
                log_1=trimming(info,read_1)
                if log_1!=None :
                    n_trimming+=1
                log_2=trimming(info,read_2)
                if log_2!=None :
                    n_trimming+=1
            #
            n_corr+=2
            info.fq_1.out_write(str(read_1))
            info.fq_2.out_write(str(read_2))
            n_len_s.append(len(read_1)+len(read_2))

        out=file(info.report(),'wt')
        out.write("#\n")
        out.write("VERSION\t2\n")
        out.write("SOURCE\t1\t%s\n"%info.fq_1.in_fn)
        out.write("OUTPUT\t1\t%s\n"%info.fq_1.out_fn)
        out.write("SOURCE\t2\t%s\n"%info.fq_2.in_fn)
        out.write("OUTPUT\t2\t%s\n"%info.fq_2.out_fn)
        out.write("FILTER_N_RATION\t%g\n"%info.n_ratio)
        out.write("FILTER_Q_CUTOFF\t%d\n"%info.q_cutoff)
        out.write("FILTER_Q_RATIO\t%g\n"%info.q_ratio)
        out.write("FILTER_TRIMMING\t%g\n"%info.trimming)
        out.write("N_RAW_READ\t%s\n"%n_all)
        out.write("N_SELECT\t%s\t%.3f\n"%(n_corr,1.0*n_corr/n_all))
        out.write("N_UNSELECT\t%s\t%.3f\n"%(n_incorr,1.0*n_incorr/n_all))
        out.write("N_TRIMMING\t%s\t%.3f\n"%(n_trimming,1.0*n_trimming/n_all))
        out.write("N_BASEPAIR\t%.1f\t%s\n"%\
                (1.0*sum(n_len_s)/len(n_len_s),sum(n_len_s)))
        out.close()
    #

def usage(out) :
    out.write("python %s "%__file__)
    out.write(" [N_RATIO] [Q_CUTOFF] [Q_RATIO] [TRIMMING_Q]")
    out.write(" [IN FQ_1] [OUT FQ_1] [IN FQ_2] [OUT FQ_2]")
    out.write(" [FIX_SCORE]")
    out.write("\n")
    out.write("    N_RATIO   : (default 0.1)\n")
    out.write("    Q_CUTOFF  : (default 20)\n")
    out.write("    Q_RATIO   : (default 0.40)\n")
    out.write("    TRIMMING_Q: (default 0: No trimming)\n")
    out.write("    FIX_SCORE : T or F: Convert Illumina Phred+33-like score to  Phred+33\n")

def main() :
    if len(sys.argv)<9 :
        usage(sys.stdout)
        return

    info=Info()
    info.n_ratio=float(sys.argv[1])
    info.q_cutoff=int(sys.argv[2])
    info.q_ratio=float(sys.argv[3])
    info.trimming=int(sys.argv[4])
    
    fq_1=InFq(sys.argv[5],sys.argv[6])
    fq_2=InFq(sys.argv[7],sys.argv[8])

    if len(sys.argv)>9 :
        if sys.argv[9].upper()=='Y' :
            info.fix_score=True
        else :
            info.fix_score=False
    #
    info.fq_1=fq_1
    info.fq_2=fq_2
    if not info.calibrate() :
        sys.stderr.write("ERROR: Failed to calibrate quality scores.\n")
        return False

    run_pe(info)

if __name__=='__main__' :
    #pr = cProfile.Profile()
    #pr.enable()
    main()
    #pr.disable()
    #s = io.StringIO()
    #out=file("profile.single",'wt')
    #ps = pstats.Stats(pr, stream=out)
    #ps.print_stats()
    #out.close()

