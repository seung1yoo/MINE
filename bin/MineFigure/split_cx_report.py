

import sys

class CxReportSpliter:
    def __init__(self, cxreport_fn):
        self.cxreport_fn = cxreport_fn

    def fh_dic_opener(self):
        cpg_oa_fn = '{0}_CpG_OA'.format(self.cxreport_fn)
        cpg_ot_fn = '{0}_CpG_OT'.format(self.cxreport_fn)
        cpg_ob_fn = '{0}_CpG_OB'.format(self.cxreport_fn)
        chg_oa_fn = '{0}_CHG_OA'.format(self.cxreport_fn)
        chg_ot_fn = '{0}_CHG_OT'.format(self.cxreport_fn)
        chg_ob_fn = '{0}_CHG_OB'.format(self.cxreport_fn)
        chh_oa_fn = '{0}_CHH_OA'.format(self.cxreport_fn)
        chh_ot_fn = '{0}_CHH_OT'.format(self.cxreport_fn)
        chh_ob_fn = '{0}_CHH_OB'.format(self.cxreport_fn)

        self.fh_dic = dict()
        self.fh_dic.setdefault('cpg', open(cpg_oa_fn, 'w'))
        self.fh_dic.setdefault('cpg_ot', open(cpg_ot_fn, 'w'))
        self.fh_dic.setdefault('cpg_ob', open(cpg_ob_fn, 'w'))
        self.fh_dic.setdefault('chg', open(chg_oa_fn, 'w'))
        self.fh_dic.setdefault('chg_ot', open(chg_ot_fn, 'w'))
        self.fh_dic.setdefault('chg_ob', open(chg_ob_fn, 'w'))
        self.fh_dic.setdefault('chh', open(chh_oa_fn, 'w'))
        self.fh_dic.setdefault('chh_ot', open(chh_ot_fn, 'w'))
        self.fh_dic.setdefault('chh_ob', open(chh_ob_fn, 'w'))

    def fh_dic_closer(self):
        for cx_type, fh in self.fh_dic.items():
            fh.close()

    def spliter(self):
        for line in open(self.cxreport_fn):
            items = line.rstrip('\n').split('\t')
            #[#CHROM] [POS] [STRAND] [#METHYL] [#UNMETHYL] [CONTEXT] [BASES]
            _cx_type = []
            if items[5] in ['CG']:
                _cx_type.append('cpg')
            elif items[5] in ['CHG']:
                _cx_type.append('chg')
            elif items[5] in ['CHH']:
                _cx_type.append('chh')
            else:
                print('check the context')
                print(items)
                sys.exit()
            #
            if items[2] in ['+']:
                _cx_type.append('ot')
            elif items[2] in ['-']:
                _cx_type.append('ob')
            else:
                print('check the strand')
                print(items)
                sys.exit()
            #
            cx_type = '_'.join(_cx_type)
            self.fh_dic[cx_type].write(line)
            self.fh_dic[_cx_type[0]].write(line)


def main(args):
    crs = CxReportSpliter(args.cxreport)
    crs.fh_dic_opener()
    crs.spliter()
    crs.fh_dic_closer()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('cxreport')
    args = parser.parse_args()
    main(args)
