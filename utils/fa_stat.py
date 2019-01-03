#!/usr/bin/python

import sys
import os
import getopt

STEP_START=0
STEP_END=20000
STEP_STEP=100


class Run :
    def __init__(self):
        self.in_fn_s=[]
        self._out_fn=None
        self.assembly=False
        self.step_start=STEP_START
        self.step_end=STEP_END
        self.step_step=STEP_STEP
        self.min_length=0
    def out_fn(self):
        if self._out_fn!=None :
            return self._out_fn
        if self.min_length==0 :
            return "%s.stat"%(self.in_fn_s[0])
        else :
            return "%s.%d.stat"%(self.in_fn_s[0],self.min_length)
    def raw_fn(self):
        return "%s.raw"%(self.out_fn())


class Read :
    def __init__(self):
        self.title=None
        self.seq=None
        self.cmt=None
        self.qul=None
        #
        self.nt_gc=0
        self.nt_at=0
        self.nt_nn=0
    def qual_stat(self):
        _sum=0
        for nt in self.seq.upper() :
            if nt=='A' or nt=='T' :
                self.nt_at+=1
            elif nt=='G' or nt=='C' :
                self.nt_gc+=1
            else :
                self.nt_nn+=1


class ComposeStat :
    def __init__(self):
        self.n_base=0
        self.nt_gc=0
        self.nt_at=0
        self.nt_nn=0
    def add_read(self,read):
        self.n_base+=len(read.seq)
        self.nt_at+=read.nt_at
        self.nt_gc+=read.nt_gc
        self.nt_nn+=read.nt_nn
    def to_report(self,out) :
        out.write("# GC-contents, GC-contents (excluse N)\n")
        out.write("GC\t%.2f\t%.2f\n"%\
                (100.0*self.nt_gc/(self.n_base-self.nt_nn),\
                100.0*self.nt_gc/(self.n_base)))
        out.write("# 'N'%\n")
        out.write("NN\t%d\t%.2f\n"%(self.nt_nn,100.0*self.nt_nn/self.n_base))
     

class ReadStat :
    def __init__(self):
        self.n_read=0
    def add_read(self,read):
        self.n_read+=1
    def to_report(self,out) :
        out.write("# N_READ\n")
        out.write("RD\t%d\n"%(self.n_read))


class BaseStat :
    def __init__(self):
        self.n_base=0
    def add_read(self,read):
        self.n_base+=len(read.seq)
    def to_report(self,out) :
        out.write("# N_BASE\tN_BASE_Q30\tN_BASE_Q20\n")
        out.write("BS\t%d\n"%(self.n_base))


class LengthStat :
    def __init__(self,run):
        self.n_read=0
        self.n_base=0
        self._len_s=[]
        #
        self._start=run.step_start
        self._end=run.step_end
        self._step=run.step_step
        #
        self.run=run
    def add_read(self,read):
        self.n_read+=1
        self.n_base+=len(read.seq)
        if self.run.assembly :
            self._len_s.append(len(read.seq))
    def to_report(self,out) :
        if not self.run.assembly :
            out.write("# LEN_TOTAL\tLEN_AVG\n")
            out.write("LN\t%d\t%.1f\n"%\
                    (self.n_base,self.n_base/self.n_read))
            return

        self._len_s.sort() 
        self._len_s.reverse()
        nXX=calc_Nxx(self._len_s)
        #
        out.write("# LEN_TOTAL\tLEN_MAX\tLEN_MIN\tLEN_AVG\n")
        out.write("LN\t%d\t%d\t%d\t%.1f\n"%\
                (self.n_base,max(self._len_s),min(self._len_s),\
                self.n_base/self.n_read))
        #
        out.write("# Cate\tRank\tLength\n")
        for cut,rank,_len in nXX :
            out.write("NXX\tN%d\t%d\t%d\n"%(cut,rank,_len))
        #
        out.write("# Length Distribution\n")
        for name,mrgn,freq in len_dist(self._len_s,self._start,self._end,self._step):
            out.write("DST\t%s\t%s\t%d\n"%(name,mrgn,freq))
     

def read_fasta(fn) :
    _in=None
    if fn.endswith(".gz") :
        _in=os.popen("zcat %s"%fn)
    else :
        _in=file(fn)

    title=None
    buf=[]
    for line in _in :
        if line.startswith(">") :
            if title!=None :
                f=Read()
                f.title=title
                f.seq="".join(buf)
                yield f
            title=line.strip()
            buf=[]
        else :
            buf.append(line.strip())
    if title!=None :
        f=Read()
        f.title=title
        f.seq="".join(buf)
        yield f


def select_qual_base(_fn) :
    i_seq=0
    q_value=set()
    for seq in read_fastq(_fn) :
        for q in seq.qul :
            q_value.add(ord(q))
        i_seq+=1
        if i_seq>=1000 :
            break
    
    q_value=list(q_value)
    q_value.sort()

    if q_value[0]>=64 :
        return 64
    elif q_value[0]>=59 :
        return 64
    elif q_value[0]>=33 :
        return 33
    else :
        print "ERROR: Unknown Q_score (%s)"%str(q_value)


def calc_Nxx(_len_s):
    n_total=sum(_len_s)
    nXX_s=[]
    for cut in range(10,100,10) :
        _sum=0
        for i_rank,_len in enumerate(_len_s) :
            _sum+=_len
            if 100.0*_sum/n_total >= cut :
                nXX_s.append((cut,i_rank+1,_len))
                break
    return nXX_s


def len_dist(_len_s,min_cut=0,max_cut=3000,step=200) :
    dist=[]
    for i in range(min_cut,max_cut,step) :
        dist.append( ["%d-%d"%(i+1,i+step),i+step,0] )
    dist.append( [">%d"%(max_cut),max_cut+1,0] )
    #
    for _len in _len_s :
        i_track=(_len-min_cut)/step
        if _len%step==0 :
            i_track-=1
        if _len>max_cut :
            i_track=max_cut/step
        dist[i_track][2]+=1
    return dist


def interpret(run,argv) :
    optlist, args = getopt.getopt(argv,'o:s:e:t:an:')
    run.in_fn_s=args
    for param,val in optlist :
        if param=='-o' :
            run._out_fn=val
        elif param=='-s' :
            run.step_start=int(val)
        elif param=='-e' :
            run.step_end=int(val)
        elif param=='-t' :
            run.step_step=int(val)
        elif param=='-a' :
            run.assembly=True
        elif param=='-n' :
            run.min_length=int(val)
        else :
            run.log_err("UNKONWN OPTIONS: %s / %s"%(param,val))
            return False

    return True


def usage() :
    print ""
    print "python %s [OPTIONS] [FQ_1] [[FQ_2] ...... [FQ_x]]"%__file__
    print ""
    print "\t-o\t[OUTPUT]\t<default:[FQ_1].stat"
    print ""
    print "\t-a\tAssembly Report"
    print "\t-n\tMinLength"
    print ""
    print "\t-s\t[STEP_START]\t<default:%s>"%(STEP_START)
    print "\t-e\t[STEP_END]\t<default:%s>"%(STEP_END)
    print "\t-t\t[STEP_STEP]\t<default:%s>"%(STEP_STEP)
    print ""


def main() :
    run=Run()
    if not interpret(run,sys.argv[1:]) :
        usage()
        return

    if run.in_fn_s==None or len(run.in_fn_s)==0 :
        usage()
        return

    comp_stat=ComposeStat()
    read_stat=ReadStat()
    base_stat=BaseStat()
    lngt_stat=LengthStat(run)
    
    stat_s=[]
    _len_s=[]
    for in_fn in run.in_fn_s :
        for read in read_fasta(in_fn) :
            if len(read.seq)<run.min_length :
                continue
            read.qual_stat()
            comp_stat.add_read(read)
            read_stat.add_read(read)
            base_stat.add_read(read)
            lngt_stat.add_read(read)
            if read_stat.n_read%1000==0 :
                print "Processing %s reads of %s"%(read_stat.n_read,in_fn)
            buf=[]
            buf.append(len(read.seq))
            buf.append("%d"%(read.nt_gc))
            buf.append("%d"%(read.nt_at))
            buf.append("%d"%(read.nt_nn))
            buf.append(read.title)
            stat_s.append(buf)
    #
    out=file(run.out_fn(),'wt')
    for in_fn in run.in_fn_s :
        out.write("FL\t%s\n"%(in_fn))
    base_stat.to_report(out)
    read_stat.to_report(out)
    comp_stat.to_report(out)
    lngt_stat.to_report(out)
    out.close()
    #
    stat_s.sort()
    stat_s.reverse()
    out_raw=file(run.raw_fn(),'wt')
    out_raw.write("#Name\tN_LENGTH")
    out_raw.write("\tMEAN_QUAL\tBASE(>Q30)\tBASE(>Q20)")
    out_raw.write("\tNT_GC\tNT_AT\tNT_NN")
    out_raw.write("\n")
    _sum=0
    for i_rank,record in enumerate(stat_s) :
        _sum+=record[0]
        out_raw.write("%d\t"%(i_rank+1))
        out_raw.write("%d\t"%(record[0]))
        out_raw.write("%d\t"%_sum)
        out_raw.write("\t".join(record[1:]))
        out_raw.write("\n")
    out_raw.close()
    os.system("gzip -f %s"%(run.raw_fn()))


if __name__=='__main__' :
    main()

