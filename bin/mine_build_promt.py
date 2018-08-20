import os
from Bio import SeqIO


def sliding(seq, size, offset):
    while seq:
        yield seq[:size]
        seq = seq[offset:]

def cal_gc(seq):
    n_a = seq.count('A')
    n_t = seq.count('T')
    n_g = seq.count('G')
    n_c = seq.count('C')
    n_n = seq.count('N')
    size = len(seq)
    #
    #gc_content = (n_g+n_c/size-n_n)*100
    gc_content = ((n_g+n_c)/float(size))*100
    #
    return round(gc_content,3)

def cal_obsexp(seq):
    obs = seq.count('CG')
    #
    n_g = seq.count('G')
    n_c = seq.count('C')
    size = len(seq)
    exp = (n_g*n_c)/float(size)
    #
    if exp in [0.0]:
        return 0.0
    return round(obs/exp,3)

def main(fa_fn, chr_id, out_dir):
    window_size = 500
    offset = 5
    #
    fh_dic = dict()
    fh_dic.setdefault('HCP', open(os.path.join(out_dir,'{0}.HCP.raw'.format(chr_id)),'w'))
    fh_dic.setdefault('ICP', open(os.path.join(out_dir,'{0}.ICP.raw'.format(chr_id)),'w'))
    fh_dic.setdefault('LCP', open(os.path.join(out_dir,'{0}.LCP.raw'.format(chr_id)),'w'))
    #
    for record in SeqIO.parse(open(fa_fn), 'fasta'):
        _chr, start, end, _id = record.id.split('_')
        seq = str(record.seq)
        #
        pre_type = ''
        starts = list()
        for n, window in enumerate(sliding(seq, window_size, offset)):
            if len(window) < window_size:
                continue
            w_start = int(start)+offset*n
            w_end = w_start+len(window)
            #
            gc_content = cal_gc(window.upper())
            obs_exp = cal_obsexp(window.upper())
            #
            if obs_exp > 0.75 and gc_content > 55.0:
                w_type = 'HCP'
            elif obs_exp < 0.48:
                w_type = 'LCP'
            else:
                w_type = 'ICP'
            #
            #print w_start, w_end, len(window), gc_content, obs_exp, w_type
            starts.append(w_start)
            #
            if pre_type and pre_type not in [w_type]:
                cp_start = sorted(starts)[0]
                cp_end = sorted(starts, reverse=True)[0]-1
                #print _chr, start, end, _id, pre_type, cp_start, cp_end
                fh_dic[pre_type].write('{0}\n'.format('\t'.join([str(x) for x in [_chr, cp_start, cp_end,
                    '{0}'.format('_'.join([_chr,str(cp_start),str(cp_end),pre_type,_id]))]])))
                starts = [w_start]
            #
            pre_type = w_type
            #
        cp_start = sorted(starts)[0]
        cp_end = w_end
        #print _chr, start, end, _id, w_type, cp_start, cp_end
        fh_dic[pre_type].write('{0}\n'.format('\t'.join([str(x) for x in [_chr, cp_start, cp_end,
            '{0}'.format('_'.join([_chr,str(cp_start),str(cp_end),w_type,_id]))]])))
        pre_type = ''
        starts = list()
        #
    for _type, fh in fh_dic.iteritems():
        fh.close()


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--fa-fn', default='/BiO/BioPeople/siyoo/MINE/ref/Chlorocebus_sabaeus/ENS93/anno_Promt/1.promoter.fa')
    parser.add_argument('--chr-id', default='1')
    parser.add_argument('--out-dir', default='/BiO/BioPeople/siyoo/MINE/ref/Chlorocebus_sabaeus/ENS93/anno_Promt')
    args = parser.parse_args()
    main(args.fa_fn, args.chr_id, args.out_dir)
