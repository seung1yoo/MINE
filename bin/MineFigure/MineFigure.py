
import sys
import os

class FromFai:
    def __inti__(self):
        pass

    def load_fai(self, fai_fn):
        self.fai_fn = fai_fn
        if not os.path.isfile(self.fai_fn):
            print('can not find fai file : {0}'.format(self.fai_fn))
            sys.exit()
        self.make_chrom_len_dic()
        self.make_chrom_bin_dic()

    def make_chrom_len_dic(self):
        self.chrom_len_dic = dict()
        for line in open(self.fai_fn):
            items = line.rstrip('\n').split('\t')
            chrom = items[0]
            length = int(items[1])
            self.chrom_len_dic.setdefault(chrom, length)

    def make_chrom_bin_dic(self, bin_size=100):
        print('chromosome indexing... bin_size={0}'.format(bin_size))
        self.chrom_bin_dic = dict()
        for chrom, length in self.chrom_len_dic.items():
            bp_per_step = length/float(bin_size)
            for bin_num in range(bin_size):
                bin_start = int(round(bin_num*bp_per_step+1,0))
                bin_end   = int(round((bin_num+1)*bp_per_step,0))
                self.chrom_bin_dic.setdefault(chrom, {})
                self.chrom_bin_dic[chrom].setdefault('_'.join([str(bin_start),str(bin_end)]), bin_num)
        out_fh = open('{0}.chromidx'.format(self.prefix), 'w')
        for chrom, boundary_dic in self.chrom_bin_dic.items():
            bp_per_step = self.chrom_len_dic[chrom]/float(bin_size)
            for boundary, bin_num in boundary_dic.items():
                items = [chrom, str(bp_per_step)]
                items.extend(boundary.split('_'))
                items.append(str(bin_num))
                out_fh.write('{0}\n'.format('\t'.join(items)))
        out_fh.close()

    def find_chrom_bin(self, chrom, pos, bin_size=100):
        for boundary, bin_num in self.chrom_bin_dic[chrom].items():
            bin_start, bin_end = boundary.split('_')
            if int(bin_start) <= int(pos) <= int(bin_end):
                return bin_num
        return bin_size



class FromGTF:
    def __init__(self):
        pass

    def load_gtf(self, gtf_fn):
        self.gtf_fn = gtf_fn
        if not os.path.isfile(self.gtf_fn):
            print('can not find gtf file : {0}'.format(self.gtf_fn))
            sys.exit()

    def extract_boundary_from_gtf(self, target_key, extend_bp=2000):
        self.note = '# boundary ==> upstream {0}bp ~ start of first exon ~ end of last exon ~ downstream {0}bp'.format(extend_bp)
        self._boundary_dic = dict()
        for line in open(self.gtf_fn):
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split('\t')
            chrom = items[0]
            feature = items[2]
            start, end = sorted([int(pos) for pos in items[3:5]])
            strand = items[6]
            attr_dic = self.attr_parser(items[-1])
            _id = attr_dic[target_key]
            #
            if not attr_dic['gene_biotype'] in ['protein_coding']:
                continue
            if not feature in ['exon']:
                continue
            #
            if _id in self._boundary_dic:
                if not chrom in [self._boundary_dic[_id]['origin']['chrom']]:
                    print('wrong gtf : {0}'.format('\t'.join(items)))
                    sys.exit()
                if start < self._boundary_dic[_id]['origin']['start']:
                    self._boundary_dic[_id]['origin']['start'] = start
                if self._boundary_dic[_id]['origin']['end'] < end:
                    self._boundary_dic[_id]['origin']['end'] = end
            else:
                self._boundary_dic.setdefault(_id, {}).setdefault('origin', {})
                self._boundary_dic[_id]['origin'].setdefault('chrom', chrom)
                self._boundary_dic[_id]['origin'].setdefault('start', start)
                self._boundary_dic[_id]['origin'].setdefault('end', end)
                self._boundary_dic[_id]['origin'].setdefault('strand', strand)

        for _id, region_dic in self._boundary_dic.items():
            info_dic = region_dic['origin']
            chrom = info_dic['chrom']
            start = int(info_dic['start'])
            end = int(info_dic['end'])
            strand = info_dic['strand']
            #
            extend_start = start-extend_bp
            if extend_start < 1:
                extend_start = 1
            extend_end = end+extend_bp
            if extend_end > self.chrom_len_dic[chrom]:
                extend_end = self.chrom_len_dic[chrom]
            #
            if strand in ['+'] and not extend_bp in [0]:
                self._boundary_dic.setdefault(_id, {}).setdefault('up', {})
                self._boundary_dic[_id]['up'].setdefault('chrom', chrom)
                self._boundary_dic[_id]['up'].setdefault('start', extend_start)
                self._boundary_dic[_id]['up'].setdefault('end', start-1)
                self._boundary_dic[_id]['up'].setdefault('strand', strand)
                self._boundary_dic.setdefault(_id, {}).setdefault('down', {})
                self._boundary_dic[_id]['down'].setdefault('chrom', chrom)
                self._boundary_dic[_id]['down'].setdefault('start', end+1)
                self._boundary_dic[_id]['down'].setdefault('end', extend_end)
                self._boundary_dic[_id]['down'].setdefault('strand', strand)
            elif strand in ['-'] and not extend_bp in [0]:
                self._boundary_dic.setdefault(_id, {}).setdefault('up', {})
                self._boundary_dic[_id]['up'].setdefault('chrom', chrom)
                self._boundary_dic[_id]['up'].setdefault('start', end+1)
                self._boundary_dic[_id]['up'].setdefault('end', extend_end)
                self._boundary_dic[_id]['up'].setdefault('strand', strand)
                self._boundary_dic.setdefault(_id, {}).setdefault('down', {})
                self._boundary_dic[_id]['down'].setdefault('chrom', chrom)
                self._boundary_dic[_id]['down'].setdefault('start', extend_start)
                self._boundary_dic[_id]['down'].setdefault('end', start-1)
                self._boundary_dic[_id]['down'].setdefault('strand', strand)
            elif strand in ['+'] and extend_bp in [0]:
                self._boundary_dic.setdefault(_id, {}).setdefault('up', {})
                self._boundary_dic[_id]['up'].setdefault('chrom', chrom)
                self._boundary_dic[_id]['up'].setdefault('start', start)
                self._boundary_dic[_id]['up'].setdefault('end', start)
                self._boundary_dic[_id]['up'].setdefault('strand', strand)
                self._boundary_dic.setdefault(_id, {}).setdefault('down', {})
                self._boundary_dic[_id]['down'].setdefault('chrom', chrom)
                self._boundary_dic[_id]['down'].setdefault('start', end)
                self._boundary_dic[_id]['down'].setdefault('end', end)
                self._boundary_dic[_id]['down'].setdefault('strand', strand)
            elif strand in ['-'] and extend_bp in [0]:
                self._boundary_dic.setdefault(_id, {}).setdefault('up', {})
                self._boundary_dic[_id]['up'].setdefault('chrom', chrom)
                self._boundary_dic[_id]['up'].setdefault('start', end)
                self._boundary_dic[_id]['up'].setdefault('end', end)
                self._boundary_dic[_id]['up'].setdefault('strand', strand)
                self._boundary_dic.setdefault(_id, {}).setdefault('down', {})
                self._boundary_dic[_id]['down'].setdefault('chrom', chrom)
                self._boundary_dic[_id]['down'].setdefault('start', start)
                self._boundary_dic[_id]['down'].setdefault('end', start)
                self._boundary_dic[_id]['down'].setdefault('strand', strand)
            else:
                print('check the strand')
                sys.exit()
            self._boundary_dic.setdefault(_id, {}).setdefault('extend', {})
            self._boundary_dic[_id]['extend'].setdefault('chrom', chrom)
            self._boundary_dic[_id]['extend'].setdefault('start', extend_start)
            self._boundary_dic[_id]['extend'].setdefault('end', extend_end)
            self._boundary_dic[_id]['extend'].setdefault('strand', strand)
            #
        return self._boundary_dic

    def attr_parser(self, attr):
        attr_dic = dict()
        for k_v in attr.split(';'):
            if not k_v.strip():
                continue
            _k = k_v.strip().split(' ')[0]
            _v = ' '.join(k_v.strip().split()[1:])
            k = _k.strip()
            v = _v.strip().strip('"')
            attr_dic.setdefault(k, v)
        return attr_dic


class FromBED:
    def __init__(self):
        pass

    def load_region_from_bed(self, bed_fn):
        self.region_dic = dict()
        for line in open(bed_fn):
            items = line.rstrip('\n').split('\t')
            self.region_dic.setdefault(items[0], [])
            chrom = items[0]
            start, stop = sorted([int(pos) for pos in items[1:3]])
            for pos in range(start, stop+1):
                if not pos in self.region_dic[chrom]:
                    self.region_dic[chrom].append(pos)
        print('Set up complete for region dic...')


class AllGeneMethylRate(FromGTF, FromFai):
    def __init__(self):
        print('Start AllGeneMethylRate')

    def load_prefix(self, prefix):
        self.prefix = prefix

    def add_methyl_from_cx(self, cx_fn_dic, depth_cutoff=5, regions=['origin','up','down','extend']):
        b_dic = self.reform_boundary_dic(regions)
        #
        self._methyl_dic = dict()
        for sample, cx_fn in cx_fn_dic.items():
            print('add methyl from... {0}'.format(sample))
            for _id, region_dic in self._boundary_dic.items():
                for region in regions:
                    self._methyl_dic.setdefault(_id, {}).setdefault(region, {})
                    self._methyl_dic[_id][region].setdefault(sample, {})
                    self._methyl_dic[_id][region][sample].setdefault('n_met', 0)
                    self._methyl_dic[_id][region][sample].setdefault('n_unmet', 0)
            self.update_methyl_dic_from_cx(sample, cx_fn, b_dic, depth_cutoff)
        self.calculate_methylrate()
        return self._methyl_dic

    def reform_boundary_dic(self, regions):
        b_dic = dict()
        for _id, region_dic in sorted(self._boundary_dic.items()):
            print('reform boundaries...', _id)
            for region in regions:
                info_dic = region_dic[region]
                chrom = info_dic['chrom']
                start = int(info_dic['start'])
                end = int(info_dic['end'])
                strand = info_dic['strand']
                #
                boundary = '{0}_{1}_{2}'.format(str(start), str(end), strand)
                chrom_bin_s = self.find_chrom_bin(chrom, start)
                chrom_bin_e = self.find_chrom_bin(chrom, end)
                chrom_bins = sorted([chrom_bin_s, chrom_bin_e])
                for chrom_bin in range(chrom_bins[0], chrom_bins[1]+1):
                    b_dic.setdefault(chrom, {}).setdefault(chrom_bin, {})
                    b_dic[chrom][chrom_bin].setdefault(region, {}).setdefault(boundary, _id)
        return b_dic

    def update_methyl_dic_from_cx(self, sample, cx_fn, b_dic, depth_cutoff):
        for line in open(cx_fn):
            items = line.rstrip('\n').split('\t')
            chrom = items[0]
            pos = int(items[1])
            strand = items[2]
            n_met = int(items[3])
            n_unmet = int(items[4])
            #
            if n_met in [0] and n_unmet in [0]:
                continue
            depth = n_met + n_unmet
            if depth < depth_cutoff:
                continue
            #
            if not chrom in b_dic:
                continue
            #
            chrom_bin = self.find_chrom_bin(chrom, pos)
            if chrom_bin in b_dic[chrom]:
                for region, boundary_dic  in b_dic[chrom][chrom_bin].items():
                    for boundary, _id in boundary_dic.items():
                        start = int(boundary.split('_')[0])
                        end = int(boundary.split('_')[1])
                        b_strand = boundary.split('_')[2]
                        #if start <= pos <= end and strand in [b_strand]:
                        if start <= pos <= end:
                            self._methyl_dic[_id][region][sample]['n_met'] += n_met
                            self._methyl_dic[_id][region][sample]['n_unmet'] += n_unmet
            else:
                continue

    def calculate_methylrate(self):
        for _id, region_dic in self._methyl_dic.items():
            for region, sample_dic in region_dic.items():
                for sample, info_dic in sample_dic.items():
                    n_methyl = info_dic['n_met']
                    n_unmethyl = info_dic['n_unmet']
                    n_depth = n_methyl + n_unmethyl
                    if n_depth in [0]:
                        methyl_rate = 0.0
                    else:
                        methyl_rate = (n_methyl/float(n_depth))*100
                    self._methyl_dic[_id][region][sample].setdefault('r_met', methyl_rate)

    def write_boundary_table(self, region='extend'):
        boundary_fn = '{0}._boundary.{1}'.format(self.prefix, region)
        out_fh = open(boundary_fn, 'w')
        out_fh.write('{0}\n'.format(self.note))
        for _id, region_dic in sorted(self._boundary_dic.items()):
            info_dic = region_dic[region]
            chrom = info_dic['chrom']
            start = str(info_dic['start'])
            end   = str(info_dic['end'])
            strand = info_dic['strand']
            out_fh.write('{0}\n'.format('\t'.join([_id, chrom, start, end, strand, region])))
        out_fh.close()
        return boundary_fn

    def write_methyl_table(self, sample_group_s, region='extend'):
        sample_orders = [x.split(':')[0] for x in sample_group_s]
        methyl_fn = '{0}._methyl.{1}'.format(self.prefix, region)
        out_fh = open(methyl_fn, 'w')
        headers_1 = ['target_id','chrom','start','end','strand','region_type']
        headers_2 = ['{0}:{1}'.format(sample, tag) for sample in sample_orders \
                     for tag in ['n_methyl','n_unmethyl','r_methyl']]
        out_fh.write('{0}\t{1}\n'.format('\t'.join(headers_1), '\t'.join(headers_2)))
        for _id, region_dic in sorted(self._boundary_dic.items()):
            info_dic = region_dic[region]
            chrom = info_dic['chrom']
            start = str(info_dic['start'])
            end = str(info_dic['end'])
            strand = info_dic['strand']
            items = [_id, chrom, start, end, strand, region]
            for sample in sample_orders:
                items.append(str(self._methyl_dic[_id][region][sample]['n_met']))
                items.append(str(self._methyl_dic[_id][region][sample]['n_unmet']))
                items.append(str(self._methyl_dic[_id][region][sample]['r_met']))
            out_fh.write('{0}\n'.format('\t'.join(items)))
        out_fh.close()
        return methyl_fn


class AllGeneMethylRateByBin(AllGeneMethylRate):
    def __init__(self):
        print('Start AllGeneMethylRateByBin')

    def define_region_of_bin(self, bin_size_dic):
        print("define bin of region")
        for _id, region_dic in self._boundary_dic.items():
            for region, bin_size in bin_size_dic.items():
                info_dic = region_dic[region]
                chrom = info_dic['chrom']
                start = int(info_dic['start'])
                end = int(info_dic['end'])
                strand = info_dic['strand']
                #
                bp_len = end-start+1
                bp_per_step = bp_len/float(bin_size)
                for bin_num in range(bin_size):
                    if strand in ['+']:
                        region_bin_id = '{0}_{1}'.format(region, bin_num)
                    elif strand in ['-']:
                        region_bin_id = '{0}_{1}'.format(region, bin_size-bin_num-1)
                    bin_start = int(round(bin_num*bp_per_step+start,0))
                    bin_end   = int(round((bin_num+1)*bp_per_step+start-1,0))
                    #
                    bin_boundary = '_'.join([str(bin_start), str(bin_end), strand])
                    self._boundary_dic[_id].setdefault(region_bin_id, {})
                    self._boundary_dic[_id][region_bin_id].setdefault('chrom', chrom)
                    self._boundary_dic[_id][region_bin_id].setdefault('start', bin_start)
                    self._boundary_dic[_id][region_bin_id].setdefault('end', bin_end)
                    self._boundary_dic[_id][region_bin_id].setdefault('strand', strand)
                self._boundary_dic[_id][region_bin_id].setdefault('end', end) # important

        region_bin_id_s = list()
        for region in ['up', 'origin', 'down']:
            bin_size = bin_size_dic[region]
            for bin_num in range(bin_size):
                region_bin_id_s.append('{0}_{1}'.format(region, str(bin_num)))

        out = open('{0}._region_bin_boundary'.format(self.prefix), 'w')
        for _id, rb_dic in sorted(self._boundary_dic.items()):
            for region_bin_id in sorted(region_bin_id_s):
                items = [_id, region_bin_id]
                info_dic = rb_dic[region_bin_id]
                items.append(info_dic['chrom'])
                items.append(str(info_dic['start']))
                items.append(str(info_dic['end']))
                items.append(info_dic['strand'])
                out.write('{0}\n'.format('\t'.join(items)))
        out.close()

        return region_bin_id_s

    def write_methylByBin_table(self, sample_group_s, region_bin_s):
        sample_orders = [x.split(':')[0] for x in sample_group_s]
        bin_fn = '{0}._methylByBin'.format(self.prefix)
        out_fh = open(bin_fn, 'w')
        headers_1 = ['target_id','chrom','start','end','strand','region_type']
        headers_2 = ['{0}:{1}'.format(sample, tag) for sample in sample_orders \
                     for tag in ['n_methyl','n_unmethyl','r_methyl']]
        out_fh.write('{0}\t{1}\n'.format('\t'.join(headers_1), '\t'.join(headers_2)))
        for _id, region_dic in sorted(self._boundary_dic.items()):
            for region_bin in region_bin_s:
                info_dic = region_dic[region_bin]
                chrom = info_dic['chrom']
                start = str(info_dic['start'])
                end = str(info_dic['end'])
                strand = info_dic['strand']
                items = [_id, chrom, start, end, strand, region_bin]
                for sample in sample_orders:
                    items.append(str(self._methyl_dic[_id][region_bin][sample]['n_met']))
                    items.append(str(self._methyl_dic[_id][region_bin][sample]['n_unmet']))
                    items.append(str(self._methyl_dic[_id][region_bin][sample]['r_met']))
                out_fh.write('{0}\n'.format('\t'.join(items)))
        out_fh.close()
        return bin_fn

    def write_methylByBinSum_table(self, sample_group_s, region_bin_s):
        sample_orders = [x.split(':')[0] for x in sample_group_s]
        sum_dic = dict()
        for _id, region_dic in sorted(self._boundary_dic.items()):
            for region_bin in region_bin_s:
                for sample in sample_orders:
                    sum_dic.setdefault(region_bin, {}).setdefault(sample, {}).setdefault('n_met', 0)
                    sum_dic[region_bin][sample]['n_met'] += int(self._methyl_dic[_id][region_bin][sample]['n_met'])
                    sum_dic.setdefault(region_bin, {}).setdefault(sample, {}).setdefault('n_unmet', 0)
                    sum_dic[region_bin][sample]['n_unmet'] += int(self._methyl_dic[_id][region_bin][sample]['n_unmet'])
        for region_bin, sample_dic in sum_dic.items():
            for sample, info_dic in sample_dic.items():
                n_met = info_dic['n_met']
                n_unmet = info_dic['n_unmet']
                depth = n_met + n_unmet
                r_met = n_met/float(depth)*100
                sum_dic.setdefault(region_bin, {}).setdefault(sample, {}).setdefault('r_met', r_met)

        binsum_fn = '{0}._methylByBinSum'.format(self.prefix)
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



class AllGeneMethylRateDMC(AllGeneMethylRateByBin):
    def __init__(self):
        print('Start AllGeneMethylRateDMC')

    def add_methyl_from_cx_for_dmc(self, cx_fn_dic, group_dic, sample_group_s, depth_cutoff=5, selected_region='extend'):
        self._pos_dic = dict()
        for _id, region_dic in self._boundary_dic.items():
            info_dic = region_dic[selected_region]
            chrom = info_dic['chrom']
            start = int(info_dic['start'])
            end = int(info_dic['end'])
            strand = info_dic['strand']
            #
            self._pos_dic.setdefault(chrom, {}).setdefault(_id, [start, end, strand])

        self.sample_methyl_dic = dict()
        self.group_methyl_dic = dict()
        for sample, cx_fn in cx_fn_dic.items():
            out_fh = open('{0}.{1}'.format(self.prefix, sample), 'w')
            for line in open(cx_fn):
                items = line.rstrip('\n').split('\t')
                chrom = items[0]
                pos = int(items[1])
                strand = items[2]
                n_met = int(items[3])
                n_unmet = int(items[4])
                depth = n_met + n_unmet
                if depth < depth_cutoff:
                    continue
                if not chrom in self._pos_dic:
                    continue
                for _id, pos_s in self._pos_dic[chrom].items():
                    if pos_s[0] <= pos <= pos_s[1]:
                        items.append(_id)
                        items.append(str(pos_s[0]))
                        items.append(str(pos_s[1]))
                        items.append(str(pos_s[2]))
                        out_fh.write('{0}\n'.format('\t'.join(items)))
                        #
                        self.sample_methyl_dic.setdefault(chrom, {}).setdefault(pos, {})
                        self.sample_methyl_dic[chrom][pos].setdefault(sample, {}).setdefault('n_met', 0)
                        self.sample_methyl_dic[chrom][pos][sample]['n_met'] += n_met
                        self.sample_methyl_dic[chrom][pos].setdefault(sample, {}).setdefault('n_unmet', 0)
                        self.sample_methyl_dic[chrom][pos][sample]['n_unmet'] += n_unmet
                        #
                        self.group_methyl_dic.setdefault(chrom, {}).setdefault(pos, {})
                        self.group_methyl_dic[chrom][pos].setdefault(group_dic[sample], {}).setdefault('n_met', 0)
                        self.group_methyl_dic[chrom][pos][group_dic[sample]]['n_met'] += n_met
                        self.group_methyl_dic[chrom][pos].setdefault(group_dic[sample], {}).setdefault('n_unmet', 0)
                        self.group_methyl_dic[chrom][pos][group_dic[sample]]['n_unmet'] += n_unmet
                        #
            out_fh.close()

        sample_orders = [x.split(':')[0] for x in sample_group_s]
        out_fh = open('{0}.Samples'.format(self.prefix), 'w')
        headers_1 = ['chrom', 'pos']
        headers_2 = ['{0}:{1}'.format(sample, tag) for sample in sample_orders \
                     for tag in ['n_methyl','n_unmethyl','r_methyl']]
        out_fh.write('{0}\t{1}\tannoGenes\n'.format('\t'.join(headers_1), '\t'.join(headers_2)))
        for chrom, pos_dic in sorted(self.sample_methyl_dic.items()):
            for pos, sample_dic in sorted(pos_dic.items()):
                items = [chrom, str(pos)]
                for sample in sample_orders:
                    if sample in sample_dic:
                        info_dic = sample_dic[sample]
                        n_met = info_dic['n_met']
                        n_unmet = info_dic['n_unmet']
                        depth = n_met + n_unmet
                        r_met = n_met/float(depth)*100
                    else:
                        n_met = 0
                        n_unmet = 0
                        r_met = 0.0
                    #
                    items.append(str(n_met))
                    items.append(str(n_unmet))
                    items.append(str(r_met))
                items.append(self.annoGenes_at_pos(chrom, int(pos)))
                out_fh.write('{0}\n'.format('\t'.join(items)))
        out_fh.close()

        group_orders = list()
        _group_orders = [x.split(':')[1] for x in sample_group_s]
        for group in _group_orders:
            if not group in group_orders:
                group_orders.append(group)
        out_fh = open('{0}.Groups'.format(self.prefix), 'w')
        headers_1 = ['chrom', 'pos']
        headers_2 = ['{0}:{1}'.format(group, tag) for group in group_orders \
                     for tag in ['n_methyl','n_unmethyl','r_methyl']]
        out_fh.write('{0}\t{1}\tannoGenes\n'.format('\t'.join(headers_1), '\t'.join(headers_2)))
        for chrom, pos_dic in sorted(self.group_methyl_dic.items()):
            for pos, group_dic in sorted(pos_dic.items()):
                items = [chrom, str(pos)]
                for group in group_orders:
                    if group in group_dic:
                        info_dic = group_dic[group]
                        n_met = info_dic['n_met']
                        n_unmet = info_dic['n_unmet']
                        depth = n_met + n_unmet
                        r_met = n_met/float(depth)*100
                    else:
                        n_met = 0
                        n_unmet = 0
                        r_met = 0.0
                    #
                    items.append(str(n_met))
                    items.append(str(n_unmet))
                    items.append(str(r_met))
                items.append(self.annoGenes_at_pos(chrom, int(pos)))
                out_fh.write('{0}\n'.format('\t'.join(items)))
        out_fh.close()

    def annoGenes_at_pos(self, chrom, at_pos):
        annoGenes_dic = dict()
        _id_dic = self._pos_dic[chrom]
        for _id, pos_s in _id_dic.items():
            start = int(pos_s[0])
            end = int(pos_s[1])
            strand = pos_s[2]
            if start <= int(at_pos) <= end:
                annoGenes_dic.setdefault(_id, {}).setdefault('start', str(start))
                annoGenes_dic.setdefault(_id, {}).setdefault('end', str(end))
                annoGenes_dic.setdefault(_id, {}).setdefault('strand', strand)
        annoGenes = list()
        for _id, info_dic in annoGenes_dic.items():
            annoGenes.append('{0}:{1}:{2}..{3}'.format(_id, info_dic['strand'], info_dic['start'], info_dic['end']))
        return ','.join(annoGenes)




    def find_dmc_for_group(self, delta_cutoff=10.0, cont_case_s=['CT12h:MJ12h', 'CT12h:MJ24h', 'CT12h:MJ48h']):
        self._dmr_dic = dict()
        for cont_case in cont_case_s:
            cont = cont_case.split(':')[0]
            case = cont_case.split(':')[1]

            out_fh = open('{0}.DMC.{1}_vs_{2}'.format(self.prefix, cont, case), 'w')
            headers_1 = ['chrom', 'pos']
            headers_2 = ['{0}:{1}'.format(group, tag) for group in [cont, case] \
                         for tag in ['n_methyl','n_unmethyl','r_methyl']]
            headers_3 = ['diff_delta', 'diff_type', 'abs(delta)>={0}'.format(str(delta_cutoff)), 'annoGenes']
            out_fh.write('{0}\t{1}\t{2}\n'.format('\t'.join(headers_1), '\t'.join(headers_2), '\t'.join(headers_3)))

            for chrom, pos_dic in sorted(self.group_methyl_dic.items()):
                for pos, group_dic in sorted(pos_dic.items()):
                    items = [chrom, str(pos)]
                    if cont in group_dic:
                        cont_info_dic = group_dic[cont]
                        cont_n_met = cont_info_dic['n_met']
                        cont_n_unmet = cont_info_dic['n_unmet']
                        cont_depth = cont_n_met + cont_n_unmet
                        cont_r_met = cont_n_met/float(cont_depth)*100
                    else:
                        cont_n_met = 0
                        cont_n_unmet = 0
                        cont_r_met = 0.0
                        cont_depth = 0
                        #continue
                    items.append(str(cont_n_met))
                    items.append(str(cont_n_unmet))
                    items.append(str(cont_r_met))
                    if case in group_dic:
                        case_info_dic = group_dic[case]
                        case_n_met = case_info_dic['n_met']
                        case_n_unmet = case_info_dic['n_unmet']
                        case_depth = case_n_met + case_n_unmet
                        case_r_met = case_n_met/float(case_depth)*100
                    else:
                        case_n_met = 0
                        case_n_unmet = 0
                        case_r_met = 0.0
                        case_depth = 0
                        #continue
                    items.append(str(case_n_met))
                    items.append(str(case_n_unmet))
                    items.append(str(case_r_met))
                    #
                    diff_delta = case_r_met - cont_r_met
                    if diff_delta == 0.0:
                        diff_type = '-'
                    elif diff_delta > 0.0:
                        diff_type = 'hyper'
                    elif diff_delta < 0.0:
                        diff_type = 'hypo'
                    items.append(str(diff_delta))
                    items.append(diff_type)
                    if abs(diff_delta) >= delta_cutoff and not cont_depth in [0] and not case_depth in [0]:
                        items.append('Y')
                        self._dmr_dic.setdefault(chrom, {}).setdefault(pos, {})
                        self._dmr_dic[chrom][pos].setdefault(cont_case, diff_delta)
                    else:
                        items.append('N')
                    #
                    items.append(self.annoGenes_at_pos(chrom, int(pos)))
                    out_fh.write('{0}\n'.format('\t'.join(items)))
                    #
            out_fh.close()

    def cal_dmc_by_bin_for_group(self, cont_case_s, region_bin_s):
        out_fh = open('{0}.DMC.ByBin'.format(self.prefix), 'w')
        headers_1 = ['target_id','chrom','start','end','strand','region_type']
        headers_2 = ['{0}:{1}'.format(cont_case.replace(':', '_vs_'), tag) for cont_case in cont_case_s \
                     for tag in ['n_bin_dmc','n_region_dmc','d_bin_dmc','avg_diffdelta']]
        out_fh.write('{0}\t{1}\n'.format('\t'.join(headers_1), '\t'.join(headers_2)))

        byBinSum_dic = dict()
        gene_mtx_dic = dict()
        for _id, region_dic in sorted(self._boundary_dic.items()):
            _density_dic = dict()
            _diffdelta_dic = dict()
            for region, info_dic in region_dic.items():
                for cont_case in cont_case_s:
                    _density_dic.setdefault(region, {})
                    _density_dic[region].setdefault(cont_case, 0)
                    _diffdelta_dic.setdefault(region, {})
                    _diffdelta_dic[region].setdefault(cont_case, [])
                chrom = info_dic['chrom']
                start = info_dic['start']
                end = info_dic['end']
                strand = info_dic['strand']
                for pos in range(start, end+1):
                    for cont_case in cont_case_s:
                        if not chrom in self._dmr_dic:
                            continue
                        if not pos in self._dmr_dic[chrom]:
                            continue
                        if not cont_case in self._dmr_dic[chrom][pos]:
                            continue
                        _density_dic[region][cont_case] += 1
                        diffdelta = self._dmr_dic[chrom][pos][cont_case]
                        _diffdelta_dic[region][cont_case].append(diffdelta)

            for region_bin in region_bin_s:
                info_dic = region_dic[region_bin]
                chrom = info_dic['chrom']
                start = info_dic['start']
                end = info_dic['end']
                strand = info_dic['strand']
                items = [_id, chrom, str(start), str(end), strand, region_bin]

                region = region_bin.split('_')[0]
                for cont_case in cont_case_s:
                    region_n_dmc = _density_dic[region][cont_case]
                    bin_n_dmc = _density_dic[region_bin][cont_case]
                    if region_n_dmc:
                        bin_d_dmc = bin_n_dmc/float(region_n_dmc)*100
                    else:
                        bin_d_dmc = 0.0
                    items.append(str(bin_n_dmc))
                    items.append(str(region_n_dmc))
                    items.append(str(bin_d_dmc))

                    diffdelta_s = _diffdelta_dic[region_bin][cont_case]
                    if len(diffdelta_s):
                        diffdelta_avg = sum(diffdelta_s)/len(diffdelta_s)
                    else:
                        diffdelta_avg = 0.0
                    items.append(str(diffdelta_avg))

                    byBinSum_dic.setdefault(region_bin, {}).setdefault(cont_case, {}).setdefault('bin_n_dmc', 0)
                    byBinSum_dic.setdefault(region_bin, {}).setdefault(cont_case, {}).setdefault('region_n_dmc', 0)
                    byBinSum_dic.setdefault(region_bin, {}).setdefault(cont_case, {}).setdefault('diffdelta_s', [])
                    byBinSum_dic[region_bin][cont_case]['bin_n_dmc'] += bin_n_dmc
                    byBinSum_dic[region_bin][cont_case]['region_n_dmc'] += region_n_dmc
                    byBinSum_dic[region_bin][cont_case]['diffdelta_s'].extend(diffdelta_s)

                    gene_mtx_dic.setdefault(cont_case, {}).setdefault(_id, {}).setdefault(region_bin, {}).setdefault('diffdelta_s', [])
                    gene_mtx_dic[cont_case][_id][region_bin]['diffdelta_s'].extend(diffdelta_s)

                out_fh.write('{0}\n'.format('\t'.join(items)))

        out_fh.close()

        #ByBinSum
        out_fh = open('{0}.DMC.ByBinSum'.format(self.prefix), 'w')
        headers_1 = ['region_type']
        headers_2 = ['{0}:{1}'.format(cont_case.replace(':', '_vs_'), tag) for cont_case in cont_case_s \
                     for tag in ['n_bin_dmc','n_region_dmc','d_bin_dmc','avg_diffdelta']]
        out_fh.write('{0}\t{1}\n'.format('\t'.join(headers_1), '\t'.join(headers_2)))

        for region_bin in region_bin_s:
            items = [region_bin]
            for cont_case in cont_case_s:
                bin_n_dmc = byBinSum_dic[region_bin][cont_case]['bin_n_dmc']
                region_n_dmc = byBinSum_dic[region_bin][cont_case]['region_n_dmc']
                bin_d_dmc = bin_n_dmc/float(region_n_dmc)*100
                items.append(str(bin_n_dmc))
                items.append(str(region_n_dmc))
                items.append(str(bin_d_dmc))

                diffdelta_s = byBinSum_dic[region_bin][cont_case]['diffdelta_s']
                if len(diffdelta_s):
                    diffdelta_avg = sum(diffdelta_s)/len(diffdelta_s)
                else:
                    diffdelta_avg = 0.0
                items.append(str(diffdelta_avg))

            out_fh.write('{0}\n'.format('\t'.join(items)))

        out_fh.close()

        #Delta Diff per gene/region_bin
        out_fh = open('{0}.DMC.ByBinPerGene'.format(self.prefix), 'w')
        headers_1 = ['cont_case','tracking_id','value_type']
        headers_2 = region_bin_s
        out_fh.write('{0}\t{1}\n'.format('\t'.join(headers_1), '\t'.join(headers_2)))

        for cont_case in cont_case_s:
            for _id, region_bin_dic in sorted(gene_mtx_dic[cont_case].items()):
                items = [cont_case]
                items.append(_id)
                items.append('diffdelta_avg')
                for region_bin in region_bin_s:
                    diffdelta_s = region_bin_dic[region_bin]['diffdelta_s']
                    if len(diffdelta_s):
                        diffdelta_avg = sum(diffdelta_s)/len(diffdelta_s)
                    else:
                        diffdelta_avg = 0.0
                    items.append(str(diffdelta_avg))
                out_fh.write('{0}\n'.format('\t'.join(items)))
        out_fh.close()


def load_fn(sample_group_s, fns):
    fn_dic = dict()
    group_dic = dict()
    for idx, sample_group in enumerate(sample_group_s):
        print('LOAD_FN', 'INPUT', sample_group, fns[idx])
        sample = sample_group.split(':')[0]
        group  = sample_group.split(':')[1]
        if os.path.isfile(fns[idx]):
            fn_dic.setdefault(sample, fns[idx])
            group_dic.setdefault(sample, group)
        else:
            print('Can not find {0}'.format(fns[dix]))
            sys.exit()
    return fn_dic, group_dic


def main(args):
    cx_fn_dic, group_dic = load_fn(args.sample_group_s, args.cxreports)

    if args.kind in ['AllGeneMethylRate']:
        mf = AllGeneMethylRate()
        mf.load_prefix(args.prefix)
        mf.load_gtf(args.gtf)
        mf.load_fai(args.fai)
        #
        boundary_dic = mf.extract_boundary_from_gtf('gene_id', extend_bp=2000)
        _boundary_fn = mf.write_boundary_table(region='origin')
        _boundary_fn = mf.write_boundary_table(region='up')
        _boundary_fn = mf.write_boundary_table(region='down')
        _boundary_fn = mf.write_boundary_table(region='extend')
        #
        regions = ['origin','up','down','extend']
        methyl_dic = mf.add_methyl_from_cx(cx_fn_dic, depth_cutoff=5, regions=regions)
        _methyl_fn = mf.write_methyl_table(args.sample_group_s, region='origin')
        _methyl_fn = mf.write_methyl_table(args.sample_group_s, region='up')
        _methyl_fn = mf.write_methyl_table(args.sample_group_s, region='down')
        _methyl_fn = mf.write_methyl_table(args.sample_group_s, region='extend')
        #
    elif args.kind in ['AllGeneMethylRateByBin']:
        mf = AllGeneMethylRateByBin()
        mf.load_prefix(args.prefix)
        mf.load_gtf(args.gtf)
        mf.load_fai(args.fai)
        #
        boundary_dic = mf.extract_boundary_from_gtf('gene_id', extend_bp=2000)
        _boundary_fn = mf.write_boundary_table(region='origin')
        _boundary_fn = mf.write_boundary_table(region='up')
        _boundary_fn = mf.write_boundary_table(region='down')
        _boundary_fn = mf.write_boundary_table(region='extend')
        #
        region_bin_s = mf.define_region_of_bin(bin_size_dic={'origin':50,'up':10,'down':10})
        methylByBin_dic = mf.add_methyl_from_cx(cx_fn_dic, depth_cutoff=5, regions=region_bin_s)
        methylByBin_fn = mf.write_methylByBin_table(args.sample_group_s, region_bin_s)
        methylByBinSum_fn = mf.write_methylByBinSum_table(args.sample_group_s, region_bin_s)
        #
    elif args.kind in ['AllGeneMethylRateDMC']:
        mf = AllGeneMethylRateDMC()
        mf.load_prefix(args.prefix)
        mf.load_gtf(args.gtf)
        mf.load_fai(args.fai)
        #
        boundary_dic = mf.extract_boundary_from_gtf('gene_id', extend_bp=2000)
        mf.add_methyl_from_cx_for_dmc(cx_fn_dic, group_dic, args.sample_group_s)
        mf.find_dmc_for_group(delta_cutoff=10.0, cont_case_s=args.cont_case_s)
        #
        region_bin_s = mf.define_region_of_bin(bin_size_dic={'origin':50,'up':10,'down':10})
        mf.cal_dmc_by_bin_for_group(args.cont_case_s, region_bin_s)






if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--kind', choices=(['AllGeneMethylRate', 'AllGeneMethylRateByBin', 'AllGeneMethylRateDMC']),
            default='AllGeneMethylRateByBin')
    parser.add_argument('--prefix', help='requied input (AllGeneMethylRate)',
            default='AllGeneMethylRateByBin.CpG.OA.Scaffold1')
    parser.add_argument('--sample-group-s', nargs='+', help='requied input (AllGeneMethylRate)',
            default=['CT12h_1:CT12h', 'CT12h_2:CT12h', 'CT12h_3:CT12h',
                     'MJ12h_1:MJ12h', 'MJ12h_2:MJ12h', 'MJ12h_3:MJ12h',
                     'MJ24h_1:MJ24h', 'MJ24h_2:MJ24h', 'MJ24h_3:MJ24h',
                     'MJ48h_1:MJ48h', 'MJ48h_2:MJ48h', 'MJ48h_3:MJ48h'])
    parser.add_argument('--cont-case-s', nargs='+', help='requied input (AllGeneMethylRateDMC)',
            default=['CT12h:MJ12h', 'CT12h:MJ24h', 'CT12h:MJ48h'])
    parser.add_argument('--cxreports', nargs='+', help='requied input (AllGeneMethylRate)',
            default=['/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/12h_control_1/12h_control_1.CX.CX_report.txt_CpG_OA.Scaffold1',
                     '/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/12h_control_2/12h_control_2.CX.CX_report.txt_CpG_OA.Scaffold1',
                     '/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/12h_control_3/12h_control_3.CX.CX_report.txt_CpG_OA.Scaffold1',
                     '/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/12h_MJ_1/12h_MJ_1.CX.CX_report.txt_CpG_OA.Scaffold1',
                     '/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/12h_MJ_2/12h_MJ_2.CX.CX_report.txt_CpG_OA.Scaffold1',
                     '/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/12h_MJ_3/12h_MJ_3.CX.CX_report.txt_CpG_OA.Scaffold1',
                     '/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/24h_MJ_1/24h_MJ_1.CX.CX_report.txt_CpG_OA.Scaffold1',
                     '/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/24h_MJ_2/24h_MJ_2.CX.CX_report.txt_CpG_OA.Scaffold1',
                     '/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/24h_MJ_3/24h_MJ_3.CX.CX_report.txt_CpG_OA.Scaffold1',
                     '/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/48h_MJ_1/48h_MJ_1.CX.CX_report.txt_CpG_OA.Scaffold1',
                     '/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/48h_MJ_2/48h_MJ_2.CX.CX_report.txt_CpG_OA.Scaffold1',
                     '/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Bismark/48h_MJ_3/48h_MJ_3.CX.CX_report.txt_CpG_OA.Scaffold1'])
    parser.add_argument('--gtf', help='requied input (AllGeneMethylRate)',
            default='/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Ref/gff2gtf/PGJ_1.0.geneset.gtf')
    parser.add_argument('--fai', help='requied input (AllGeneMethylRate)',
            default='/BiO/BioPeople/siyoo/TBD160182-TBI-Bellflower-Bisulfite-20190318/Ref/PGJ_1.0.genome.fa.fai')
    args = parser.parse_args()
    main(args)
