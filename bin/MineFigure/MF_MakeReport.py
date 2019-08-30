
# siyoo@dragon:/BiO/BioProjects/TBD190324-ALLBT-Pig-Bisulfite-20190624/MineFigure
# python /BiO/BioPeople/siyoo/MINE/bin/MineFigure/MF_MakeReport.py ./ ./Report


import os
import glob


class MF_MakeReport:
    def __init__(self):
        pass

    def load_indir(self, indir):
        self.indir = indir

    def load_outdir(self, outdir):
        self.outdir_home = outdir
        if not os.path.isdir(self.outdir_home):
            os.mkdir(self.outdir_home)
        self.outdir_rate = os.path.join(self.outdir_home, 'AllGeneMethylRate')
        if not os.path.isdir(self.outdir_rate):
            os.mkdir(self.outdir_rate)
        self.outdir_bybin = os.path.join(self.outdir_home, 'AllGeneMethylRateByBin')
        if not os.path.isdir(self.outdir_bybin):
            os.mkdir(self.outdir_bybin)
        self.outdir_dmc = os.path.join(self.outdir_home, 'AllGeneMethylRateDMC')
        if not os.path.isdir(self.outdir_dmc):
            os.mkdir(self.outdir_dmc)

    def parse_allgenemethylratedmc(self):
        pergene_fn_s = glob.glob('AllGeneMethylRateDMC/*.ByBinPerGene')
        for pergene_fn in pergene_fn_s:
            fh_dic = self.make_pergene_fh_dic(pergene_fn)
            self.write_report_pergene(pergene_fn, fh_dic)
        sum_fn_s = glob.glob('AllGeneMethylRateDMC/*.ByBinSum')
        for sum_fn in sum_fn_s:
            out_sum_fn = '{0}/{1}.xls'.format(self.outdir_dmc, sum_fn.split('/')[-1])
            _cmd_s = ['cp']
            _cmd_s.append(sum_fn)
            _cmd_s.append(out_sum_fn)
            os.system(' '.join(_cmd_s))

    def make_pergene_fh_dic(self, pergene_fn):
        fh_dic = dict()
        for line in open(pergene_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['cont_case']:
                headers = items
                continue
            cont_case = items[0].replace(':','_vs_')
            if not cont_case in fh_dic:
                fh_dic.setdefault(cont_case, open('{0}/{1}.{2}.xls'.format(self.outdir_dmc, pergene_fn.split('/')[-1], cont_case),'w'))
                fh_dic[cont_case].write('{0}\n'.format('\t'.join(headers)))
        return fh_dic

    def write_report_pergene(self, pergene_fn, fh_dic):
        for line in open(pergene_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['cont_case']:
                continue
            cont_case = items[0].replace(':','_vs_')
            fh_dic[cont_case].write('{0}\n'.format('\t'.join(items)))
        for cont_case, fh in fh_dic.items():
            fh.close()

    def parse_allgenemethylrate(self):
        rate_fn_s = glob.glob('AllGeneMethylRate/*._methyl.origin')
        for rate_fn in rate_fn_s:
            out_rate_fn = '{0}/{1}.xls'.format(self.outdir_rate, rate_fn.split('/')[-1])
            _cmd_s = ['cp']
            _cmd_s.append(rate_fn)
            _cmd_s.append(out_rate_fn)
            os.system(' '.join(_cmd_s))
        rate_fn_s = glob.glob('AllGeneMethylRate/*._methyl.extend')
        for rate_fn in rate_fn_s:
            out_rate_fn = '{0}/{1}.xls'.format(self.outdir_rate, rate_fn.split('/')[-1])
            _cmd_s = ['cp']
            _cmd_s.append(rate_fn)
            _cmd_s.append(out_rate_fn)
            os.system(' '.join(_cmd_s))

    def parse_allgenemethylratebybin(self):
        methylbybin_fn_s = glob.glob('AllGeneMethylRateByBin/*._methylByBin')
        for methylbybin_fn in methylbybin_fn_s:
            fh_dic = self.make_methylbybin_fh_dic(methylbybin_fn)
            self.write_report_methylbybin(methylbybin_fn, fh_dic)
        sum_fn_s = glob.glob('AllGeneMethylRateByBin/*._methylByBinSum')
        for sum_fn in sum_fn_s:
            out_sum_fn = '{0}/{1}.xls'.format(self.outdir_bybin, sum_fn.split('/')[-1])
            _cmd_s = ['cp']
            _cmd_s.append(sum_fn)
            _cmd_s.append(out_sum_fn)
            os.system(' '.join(_cmd_s))

    def make_methylbybin_fh_dic(self, methylbybin_fn):
        fh_dic = dict()
        for line in open(methylbybin_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['target_id']:
                headers = items
                continue
            chrom = items[1]
            if not chrom in fh_dic:
                fh_dic.setdefault(chrom, open('{0}/{1}.{2}.xls'.format(self.outdir_bybin, methylbybin_fn.split('/')[-1], chrom),'w'))
                fh_dic[chrom].write('{0}\n'.format('\t'.join(headers)))
        return fh_dic

    def write_report_methylbybin(self, methylbybin_fn, fh_dic):
        for line in open(methylbybin_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['target_id']:
                continue
            chrom = items[1]
            fh_dic[chrom].write('{0}\n'.format('\t'.join(items)))
        for chrom, fh in fh_dic.items():
            fh.close()



def main(args):
    mf_mr = MF_MakeReport()
    mf_mr.load_indir(args.indir)
    mf_mr.load_outdir(args.outdir)

    mf_mr.parse_allgenemethylrate()
    mf_mr.parse_allgenemethylratebybin()
    mf_mr.parse_allgenemethylratedmc()


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('indir')
    parser.add_argument('outdir')
    args = parser.parse_args()
    main(args)


