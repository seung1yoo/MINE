


def make_region_bin_s(bin_size_dic={'origin':50,'up':10,'down':10}):
    region_bin_id_s = list()
    for region in ['up', 'origin', 'down']:
        bin_size = bin_size_dic[region]
        for bin_num in range(bin_size):
            region_bin_id_s.append('{0}_{1}'.format(region, str(bin_num)))
    return region_bin_id_s

def family_dic_maker():
    family_dic = dict()
    for line in open('GeneList.v2.ALL'):
        items = line.rstrip('\n').split()
        family_dic.setdefault(items[0], items[1])
    return family_dic


def main(args):
    region_bin_s = make_region_bin_s(bin_size_dic={'origin':50,'up':10,'down':10})
    family_dic = family_dic_maker()


    print(args)
    info_dic = dict()
    value_dic = dict()
    for line in open(args.infn):
        items = line.rstrip('\n').split('\t')
        if items[0] in ['target_id']:
            idx_dic = dict()
            for idx, item in enumerate(items):
                idx_dic.setdefault(item, idx)
            continue
        target_id = items[idx_dic['target_id']]
        info_s = items[1:5]
        info_s.append(family_dic[target_id])
        info_dic.setdefault(target_id, info_s) # chrom start end strand family(function)
        region_type = items[idx_dic['region_type']]
        for item, idx in idx_dic.items():
            if item.endswith(':r_methyl'):
                sample_id = item.split(':')[0]
                value_dic.setdefault(target_id, {}).setdefault(sample_id, {}).setdefault(region_type, items[idx])
    #
    out = open('{0}.PerGene'.format(args.infn), 'w')
    header_1 = ['target_id','chrom','start','end','strand','function']
    header_2 = ['{0}_{1}'.format(sample, region_bin_id) for sample in ['CT12h','MJ12h','MJ24h','MJ48h'] for region_bin_id in region_bin_s]
    out.write('{0}\t{1}\n'.format('\t'.join(header_1), '\t'.join(header_2)))
    for target_id, sample_dic in sorted(value_dic.items()):
        items = [target_id]
        items.extend(info_dic[target_id])
        for sample in ['CT12h','MJ12h','MJ24h','MJ48h']:
            for region_bin_id in region_bin_s:
                items.append(sample_dic[sample][region_bin_id])
        out.write('{0}\n'.format('\t'.join(items)))
    out.close()



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--infn',
            default='AllGeneMethylRateByBin/AllGeneMethylRateByBin.CHH.OA._methylByBin.Interested.v2ALL.GroupAvg')
    args = parser.parse_args()
    main(args)
