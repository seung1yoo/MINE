
import os


class DMCbyBin:
    def __init__(self, in_fn, in_context, out_prefix):
        self.in_fn = in_fn
        self.in_context = in_context
        self.out_prefix = out_prefix
        self.interest_member_dic = dict()

    def setup_bin_order(self, bin_size_dic):
        self.bin_order_dic = dict()
        n = 0
        for region in ['up', 'origin', 'down']:
            bin_size = bin_size_dic[region]
            for bin_num in range(bin_size):
                n += 1
                self.bin_order_dic.setdefault('{0}_{1}'.format(region, str(bin_num)), n)

    def load_interest(self, interest_fn):
        self.interest_member_dic = dict()
        for line in open(interest_fn):
            items = line.rstrip('\n').split()
            gene_id = items[0]
            function = items[1]
            self.interest_member_dic.setdefault(gene_id, [])
            if not function in self.interest_member_dic[gene_id]:
                self.interest_member_dic[gene_id].append(function)
        return self.interest_member_dic

    def table2dic(self):
        self.all_delta_dic = dict()
        self.interest_delta_dic = dict()

        self.all_density_dic = dict()
        self.interest_density_dic = dict()

        for line in open(self.in_fn):

            items = line.rstrip('\n').split('\t')

            if items[0] in ['target_id']:
                idx_dic = dict()
                diffdelta_idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                    if item.endswith('avg_diffdelta'):
                        diffdelta_idx_dic.setdefault(item, idx)
                continue

            target_id = items[idx_dic['target_id']]
            #
            region_bin_type = items[idx_dic['region_type']]
            region_type = region_bin_type.split('_')[0]
            #
            diffdelta_dic = self.make_value_dic(items, diffdelta_idx_dic)
            #
            for compare_id, avg_diffdelta in diffdelta_dic.items():
                self.all_delta_dic.setdefault(compare_id, {}).setdefault(region_bin_type, [])
                self.all_delta_dic[compare_id][region_bin_type].append(avg_diffdelta)
                #
                self.all_density_dic.setdefault(compare_id, {}).setdefault(region_bin_type, 0)
                self.all_density_dic.setdefault(compare_id, {}).setdefault(region_type, 0)
                n_bin_dmc = int(items[idx_dic['{0}:n_bin_dmc'.format(compare_id)]])
                n_region_dmc = int(items[idx_dic['{0}:n_region_dmc'.format(compare_id)]])
                self.all_density_dic[compare_id][region_bin_type] += n_bin_dmc
                self.all_density_dic[compare_id][region_type] += n_region_dmc
                #
                if target_id in self.interest_member_dic:
                    self.interest_delta_dic.setdefault(compare_id, {}).setdefault(region_bin_type, [])
                    self.interest_delta_dic[compare_id][region_bin_type].append(avg_diffdelta)
                    #
                    self.interest_density_dic.setdefault(compare_id, {}).setdefault(region_bin_type, 0)
                    self.interest_density_dic.setdefault(compare_id, {}).setdefault(region_type, 0)
                    n_bin_dmc = int(items[idx_dic['{0}:n_bin_dmc'.format(compare_id)]])
                    n_region_dmc = int(items[idx_dic['{0}:n_region_dmc'.format(compare_id)]])
                    self.interest_density_dic[compare_id][region_bin_type] += n_bin_dmc
                    self.interest_density_dic[compare_id][region_type] += n_region_dmc

    def dic2dataframe(self):
        # DiffDelta
        out_fh = open('{0}.DiffDelta'.format(self.out_prefix), 'w')
        out_fh.write('{0}\n'.format('\t'.join(['category','compare','bin_order','avg_diffdelta'])))

        for compare_id, region_bin_type_dic in self.all_delta_dic.items():
            for region_bin_type, avg_diffdelta_s in region_bin_type_dic.items():
                region_bin_order = str(self.bin_order_dic[region_bin_type])
                new_avg_diffdelta = sum(avg_diffdelta_s)/len(avg_diffdelta_s)
                out_fh.write('{0}\n'.format('\t'.join(['AllGene',compare_id,region_bin_order,str(new_avg_diffdelta)])))
        for compare_id, region_bin_type_dic in self.interest_delta_dic.items():
            for region_bin_type, avg_diffdelta_s in region_bin_type_dic.items():
                region_bin_order = str(self.bin_order_dic[region_bin_type])
                new_avg_diffdelta = sum(avg_diffdelta_s)/len(avg_diffdelta_s)
                out_fh.write('{0}\n'.format('\t'.join(['Interested',compare_id,region_bin_order,str(new_avg_diffdelta)])))
        out_fh.close()

        # Density # Not yet
        out_fh = open('{0}.Density'.format(self.out_prefix), 'w')
        out_fh.close()

    def make_value_dic(self, items, _idx_dic):
        _dic = dict()
        for colname, idx in _idx_dic.items():
            compare_id = colname.split(':')[0]
            value = float(items[idx])
            _dic.setdefault(compare_id, value)
        return _dic



def load_interest(interest_fn):
    interest_member_dic = dict()
    for line in open(interest_fn):
        items = line.rstrip('\n').split()
        gene_id = items[0]
        function = items[1]
        interest_member_dic.setdefault(gene_id, [])
        if not function in interest_member_dic[gene_id]:
            interest_member_dic[gene_id].append(function)
    return interest_member_dic

def setup_bin_order(bin_size_dic):
    bin_order_dic = dict()
    n = 0
    for region in ['up', 'origin', 'down']:
        bin_size = bin_size_dic[region]
        for bin_num in range(bin_size):
            n += 1
            bin_order_dic.setdefault('{0}_{1}'.format(region, str(bin_num)), n)
    return bin_order_dic

class MethylRateByBin:
    def load_in_fn(self, in_fn):
        self.in_fn = in_fn

    def load_in_context(self, in_context):
        self.in_context = in_context

    def load_interest_member_dic(self, interest_member_dic):
        self.int_mem_dic = interest_member_dic

    def load_bin_order_dic(self, bin_order_dic):
        self.bin_order_dic = bin_order_dic

    def load_out_prefix(self, out_prefix):
        self.out_prefix = out_prefix

    def table2dic(self):
        for line in open(self.in_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['target_id']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            #
            target_id = items[idx_dic['target_id']]
            region_bin_id = items[idx_dic['region_type']]
            region_bin_order = self.bin_order_dic[region_bin_id]
            #
            # ING.......

def main(args):
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    out_prefix = os.path.join(args.out_dir, '{0}'.format(args.out_prefix))

    if args.interest_fn:
        interest_member_dic = load_interest(args.interest_fn)
    else:
        interest_member_dic = dict()

    bin_order_dic = setup_bin_order({'origin':50,'up':10,'down':10})

    if args.kind in ['DMCbyBin']:
        mf = DMCbyBin(args.in_fn, args.in_context, out_prefix)
        mf.setup_bin_order({'origin':50,'up':10,'down':10})
        if args.interest_fn:
            mf.load_interest(args.interest_fn)
        mf.table2dic()
        mf.dic2dataframe()
    elif args.kind in ['MethylRateByBin']:
        mf = MethylRateByBin()
        mf.load_in_fn(args.in_fn)
        mf.load_in_context(args.in_context)
        mf.load_interest_member_dic(interest_member_dic)
        mf.load_bin_order_dic(bin_order_dic)
        mf.load_out_prefix(out_prefix)
        mf.table2dic()



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--kind', choices=(['DMCbyBin', 'MethylRateByBin']),
            default='MethylRateByBin')
    parser.add_argument('--in-fn',
            default='./AllGeneMethylRateByBin/AllGeneMethylRateByBin.CpG.OA._methylByBin')
    parser.add_argument('--in-context',
            default='CpG')
    parser.add_argument('--out-dir',
            default='MF_Table2DataFrame4Draw')
    parser.add_argument('--out-prefix',
            default='MethylRateByBin.CpG')
    parser.add_argument('--interest-fn',
            default='GeneList')
    args = parser.parse_args()
    main(args)

