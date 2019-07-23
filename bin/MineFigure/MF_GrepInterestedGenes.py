

def grep_interested(_fn, interested_s, interested_name):
    out_fh = open('{0}.Interested.{1}'.format(_fn, interested_name),'w')
    for line in open(_fn):
        items = line.rstrip('\n').split('\t')
        if items[0] in ['target_id']:
            out_fh.write(line)
            continue
        if items[0] in interested_s:
            out_fh.write(line)
    out_fh.close()
    return '{0}.Interested.{1}'.format(_fn, interested_name)

def make_region_bin_s(bin_size_dic={'origin':50,'up':10,'down':10}):
    region_bin_id_s = list()
    for region in ['up', 'origin', 'down']:
        bin_size = bin_size_dic[region]
        for bin_num in range(bin_size):
            region_bin_id_s.append('{0}_{1}'.format(region, str(bin_num)))
    return region_bin_id_s

def summer(_interested_fn, region_bin_s, sample_orders):
    sum_dic = dict()
    for line in open(_interested_fn):
        items = line.rstrip('\n').split('\t')
        if items[0] in ['target_id']:
            idx_dic = dict()
            for idx, item in enumerate(items):
                idx_dic.setdefault(item, idx)
            target_colnames = items[6:]
            continue
        bin_id = items[idx_dic['region_type']]
        #
        for colname in target_colnames:
            sample, _type = colname.split(':')
            if _type in ['n_methyl']:
                sum_dic.setdefault(bin_id, {}).setdefault(sample, {}).setdefault('n_met', 0)
                sum_dic[bin_id][sample]['n_met'] += int(items[idx_dic[colname]])
            if _type in ['n_unmethyl']:
                sum_dic.setdefault(bin_id, {}).setdefault(sample, {}).setdefault('n_unmet', 0)
                sum_dic[bin_id][sample]['n_unmet'] += int(items[idx_dic[colname]])
    for region_bin, sample_dic in sum_dic.items():
        for sample, info_dic in sample_dic.items():
            n_met = info_dic['n_met']
            n_unmet = info_dic['n_unmet']
            depth = n_met + n_unmet
            if depth:
                r_met = n_met/float(depth)*100
            else:
                r_met = 0.0
            sum_dic.setdefault(region_bin, {}).setdefault(sample, {}).setdefault('r_met', r_met)

    binsum_fn = '{0}.Sum'.format(_interested_fn)
    out_fh = open(binsum_fn, 'w')
    headers_1 = ['region_type']
    headers_2 = ['{0}:{1}'.format(sample, tag) for sample in sample_orders \
                 for tag in ['n_methyl','n_unmethyl','r_methyl']]
    out_fh.write('{0}\t{1}\n'.format('\t'.join(headers_1), '\t'.join(headers_2)))
    for region_bin in region_bin_s:
        items = [region_bin]
        for sample in sample_orders:
            info_dic = sum_dic[region_bin][sample]
            items.append(str(info_dic['n_met']))
            items.append(str(info_dic['n_unmet']))
            items.append(str(info_dic['r_met']))
        out_fh.write('{0}\n'.format('\t'.join(items)))
    out_fh.close()
    return binsum_fn

def make_sample_orders():
    return ['S1', 'S2', 'PIG_3', 'PIG_4']

def main(args):
    interested_s = list()
    for line in open(args.genelist_fn):
        interested_s.append(line.rstrip('\n').split()[0])
    interested_name = args.interested_name
    sample_orders = args.samples
    region_bin_s = make_region_bin_s(bin_size_dic={'origin':50,'up':10,'down':10})
    #
    #cpg_fn = 'AllGeneMethylRateByBin/AllGeneMethylRateByBin.CpG.OA._methylByBin'
    #chg_fn = 'AllGeneMethylRateByBin/AllGeneMethylRateByBin.CHG.OA._methylByBin'
    #chh_fn = 'AllGeneMethylRateByBin/AllGeneMethylRateByBin.CHH.OA._methylByBin'
    #
    #cpg_interested_fn = grep_interested(cpg_fn, interested_s, interested_name)
    #chg_interested_fn = grep_interested(chg_fn, interested_s, interested_name)
    #chh_interested_fn = grep_interested(chh_fn, interested_s, interested_name)
    #
    #summer(cpg_interested_fn, region_bin_s, sample_orders)
    #summer(chg_interested_fn, region_bin_s, sample_orders)
    #summer(chh_interested_fn, region_bin_s, sample_orders)
    #
    interested_fn = grep_interested(args.methylByBin_fn, interested_s, interested_name)
    summer(interested_fn, region_bin_s, sample_orders)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('genelist_fn')
    parser.add_argument('interested_name')
    parser.add_argument('methylByBin_fn')
    parser.add_argument('--samples', nargs='+')
    args = parser.parse_args()
    main(args)
