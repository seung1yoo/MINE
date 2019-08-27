








def parse_sample_group_pairs(pairs):
    sam_grp_dic = dict()
    for pair in pairs:
        sample, group = pair.split(':')
        sam_grp_dic.setdefault(sample, group)

    _groups = [pair.split(':')[1] for pair in pairs]
    groups = list()
    for group in _groups:
        if not group in groups:
            groups.append(group)

    return sam_grp_dic, groups



def main(args):
    sam_grp_dic, groups = parse_sample_group_pairs(args.sample_group_pairs)
    print(groups)

    if args.intype in ['AllGeneMethylRate','AllGeneMethylRateByBinSum']:
        groupAvg_dic = dict()
        for line in open(args.infn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['region_type', 'target_id']:
                header_key = items[0]
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            if header_key in ['region_type']:
                region_type = items[idx_dic['region_type']]
                _id = region_type.split('_')[0]
                num = int(region_type.split('_')[1])
            elif header_key in ['target_id']:
                target_id  = items[idx_dic['target_id']]
                _id = target_id
                num = 0
            for sample, group in sam_grp_dic.items():
                colname_methyl = '{0}:{1}'.format(sample, 'n_methyl')
                colname_unmethyl = '{0}:{1}'.format(sample, 'n_unmethyl')
                n_methyl = int(items[idx_dic[colname_methyl]])
                n_unmethyl = int(items[idx_dic[colname_unmethyl]])
                #
                groupAvg_dic.setdefault(_id, {}).setdefault(num, {})
                groupAvg_dic[_id][num].setdefault(group, {})
                groupAvg_dic[_id][num][group].setdefault('n_methyl', 0)
                groupAvg_dic[_id][num][group]['n_methyl'] += n_methyl
                groupAvg_dic[_id][num][group].setdefault('n_unmethyl', 0)
                groupAvg_dic[_id][num][group]['n_unmethyl'] += n_unmethyl
        for _id, num_dic in groupAvg_dic.items():
            for num, grp_dic in num_dic.items():
                for group in groups:
                    info_dic = grp_dic[group]
                    n_methyl = int(info_dic['n_methyl'])
                    n_unmethyl = int(info_dic['n_unmethyl'])
                    depth = n_methyl + n_unmethyl
                    if depth in [0]:
                        groupAvg_dic[_id][num][group].setdefault('r_methyl', 0.0)
                    else:
                        r_methyl = n_methyl/float(depth)*100
                        groupAvg_dic[_id][num][group].setdefault('r_methyl', r_methyl)
        out_fh = open('{0}.GroupAvg'.format(args.infn),'w')
        headers_1 = ['tracking_id','num']
        headers_2 = ['{0}:{1}'.format(group, tag) for group in groups \
                     for tag in ['n_methyl','n_unmethyl','r_methyl']]
        out_fh.write('{0}\t{1}\n'.format('\t'.join(headers_1), '\t'.join(headers_2)))
        for _id, num_dic in sorted(groupAvg_dic.items(), reverse=True):
            for num, grp_dic in sorted(num_dic.items()):
                items = [_id, str(num)]
                for group in groups:
                    info_dic = grp_dic[group]
                    items.append(str(info_dic['n_methyl']))
                    items.append(str(info_dic['n_unmethyl']))
                    items.append(str(info_dic['r_methyl']))
                out_fh.write('{0}\n'.format('\t'.join(items)))
        out_fh.close()
    elif args.intype in ['AllGeneMethylRateByBin']:
        out_fh = open('{0}.GroupAvg'.format(args.infn),'w')

        for line in open(args.infn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['target_id']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                #
                headers_1 = items[:6]
                headers_2 = ['{0}:{1}'.format(group, tag) for group in groups \
                             for tag in ['n_methyl','n_unmethyl','r_methyl']]
                out_fh.write('{0}\t{1}\n'.format('\t'.join(headers_1), '\t'.join(headers_2)))
                continue
            groupAvg_dic = dict()
            for sample, group in sam_grp_dic.items():
                colname_methyl = '{0}:{1}'.format(sample, 'n_methyl')
                colname_unmethyl = '{0}:{1}'.format(sample, 'n_unmethyl')
                n_methyl = int(items[idx_dic[colname_methyl]])
                n_unmethyl = int(items[idx_dic[colname_unmethyl]])
                #
                groupAvg_dic.setdefault(group, {})
                groupAvg_dic[group].setdefault('n_methyl', 0)
                groupAvg_dic[group]['n_methyl'] += n_methyl
                groupAvg_dic[group].setdefault('n_unmethyl', 0)
                groupAvg_dic[group]['n_unmethyl'] += n_unmethyl
                #
            for group, info_dic in groupAvg_dic.items():
                n_methyl = int(info_dic['n_methyl'])
                n_unmethyl = int(info_dic['n_unmethyl'])
                depth = n_methyl + n_unmethyl
                if depth in [0]:
                    groupAvg_dic[group].setdefault('r_methyl', 0.0)
                else:
                    r_methyl = n_methyl/float(depth)*100
                    groupAvg_dic[group].setdefault('r_methyl', r_methyl)
            #
            new_items = items[:6]
            for group in groups:
                info_dic = groupAvg_dic[group]
                new_items.append(str(info_dic['n_methyl']))
                new_items.append(str(info_dic['n_unmethyl']))
                new_items.append(str(info_dic['r_methyl']))
            out_fh.write('{0}\n'.format('\t'.join(new_items)))
        out_fh.close()






if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--infn', help='MF_methylByBin or SUM',
            default='AllGeneMethylRateByBin.CpG.OA.Scaffold1._methylByBinSum')
    parser.add_argument('--intype', choices=('AllGeneMethylRate','AllGeneMethylRateByBin','AllGeneMethylRateByBinSum'),
            default='AllGeneMethylRateByBin')
    parser.add_argument('--sample-group-pairs', nargs="+",
            default= ["CT12h_1:CT12h", "CT12h_2:CT12h", "CT12h_3:CT12h",
                      "MJ12h_1:MJ12h", "MJ12h_2:MJ12h", "MJ12h_3:MJ12h",
                      "MJ24h_1:MJ24h", "MJ24h_2:MJ24h", "MJ24h_3:MJ24h",
                      "MJ48h_1:MJ48h", "MJ48h_2:MJ48h", "MJ48h_3:MJ48h"])
    args = parser.parse_args()
    main(args)
