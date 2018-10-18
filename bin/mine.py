#!/usr/bin/python

import sys
import os
import logging
import json
import glob
#
from lib.Qsub import *
from lib.Tools import *
from subprocess import Popen, PIPE
#
from do_link_raw import *
from do_fastqc import *
from do_stat_fastq import *
from do_filter_fastq import *
from do_bismark import *
from do_target_check import *
from do_metilene import *
#
#from do_metilene_input import *
#from do_metilene_parse import *



class Mine:
    def __init__(self, config):
        self.setup_logging()
        #
        self.config_fn = config
        self.config_parsing()
        self.config_checking()
        self.idmatch_dic = self.idmatch()

    def setup_logging(self):
        self.logger = logging.getLogger(self.__class__.__name__)
        ## self.logger.setLevel(logging.[*])
        ## [DEBUG<INFO<WARNING<ERROR<CRITICAL]
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter(\
            '%(asctime)s : '
            '%(name)s : '
            '%(levelname)s : '
            '%(message)s')
        stream_hander = logging.StreamHandler()
        stream_hander.setFormatter(formatter)
        self.logger.addHandler(stream_hander)
        #
        self.log_fn = '{0}.log'.format(self.__class__.__name__)
        file_handler = logging.FileHandler(self.log_fn)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def config_parsing(self):
        self.species = None
        self.home = None
        self.output = None
        self.sample_dic = dict()
        self.dmr_dic = dict()
        self.targetbed = None
        for line in open(self.config_fn):
            items = line.rstrip('\n').split()
            if items[0].startswith('#'):
                continue
            elif items[0] in ['SPECIES']:
                self.species = items[1].strip()
            elif items[0] in ['HOME']:
                self.home = items[1].strip()
            elif items[0] in ['OUTPUT']:
                self.output = items[1].strip().replace('[HOME]',self.home)
            elif items[0] in ['INPUT']:
                num = 'S{0:0>4}'.format(items[1].strip())
                pair = items[2].strip()
                library = self.convert_library(items[3].strip())
                path = items[4].strip().replace('[HOME]',self.home)
                label = items[5].strip()
                self.sample_dic.setdefault(num, {}).setdefault(
                        pair, {}).setdefault('library', library)
                self.sample_dic.setdefault(num, {}).setdefault(
                        pair, {}).setdefault('path', path)
                self.sample_dic.setdefault(num, {}).setdefault(
                        pair, {}).setdefault('label', label)
            elif items[0] in ['DMR']:
                num = 'DMR{0:0>4}'.format(items[1].strip())
                pair = items[2].strip()
                groups = pair.split('_')
                if not len(groups) in [2]:
                    self.logger.error('Wrong DMR pair')
                    sys.exit()
                group1s = ['S{0:0>4}'.format(x) for x in groups[0].split('-')]
                group2s = ['S{0:0>4}'.format(x) for x in groups[1].split('-')]
                #
                self.dmr_dic.setdefault(num, {}).setdefault('G1', group1s)
                self.dmr_dic.setdefault(num, {}).setdefault('G2', group2s)
            elif items[0] in ['TARGET']:
                if items[1] in ['ON']:
                    self.targetbed = items[2]
                    if not os.path.exists(self.targetbed):
                        self.logger.error("NOT EXIST : {0}".format(items[2]))
                        sys.exit()
                elif items[1] in ['OFF']:
                    self.targetbed = None
            else:
                self.logger.error("UNKNOWN LINE IN CONFIG='{0}'".format(line))
                sys.exit()

    def convert_library(self, libraryType):
        lib_kit_dic = {'1':'unstrand specific',
                       '2':'first strand (TruSeq , SureSelect)',
                       '3':'second strand'}
        return lib_kit_dic[libraryType]

    def config_checking(self):
        if self.species==None:
            self.logger.error("NONE PARAMETER IN CONFIG = SPECIES")
            sys.exit()
        if self.home==None:
            self.logger.error("NONE PARAMETER IN CONFIG = HOME")
            sys.exit()
        if self.output==None:
            self.logger.error("NONE PARAMETER IN CONFIG = OUTPUT")
            sys.exit()
        if len(self.sample_dic) in [0]:
            self.logger.error("NONE PARAMETER IN CONFIG = INPUT")
            sys.exit()
        if len(self.dmr_dic) in [0]:
            self.logger.error("NONE PARAMETER IN CONFIG = DMR")
            sys.exit()

    def idmatch(self):
        idmatch_dic = dict()
        for num, pair_dic in self.sample_dic.iteritems():
            for pair, info_dic in pair_dic.iteritems():
                idmatch_dic.setdefault(num, info_dic['label'])
        return idmatch_dic

    def add_tools(self, tools_fn):
        self.tools = json.load(open(tools_fn))

    def add_refs(self, refs_fn):
        self.refs = json.load(open(refs_fn))[self.species]
        self.refs['BISMARK_DIR'] = self.refs['BISMARK_DIR'].replace('[HOME]', self.refs['HOME'])

class Eco(Mine):
    def __init__(self):
        self.setup_logging()

    def build_home(self, home_path):
        self.home_path = home_path
        self.make_dir(self.home_path)
        #
        self.room = dict()
        self.add_room('analysis')
        self.add_room('scripts')
        self.add_room('logs')
        self.add_room('stats')
        self.add_room('report')

    def add_room(self, name):
        dir_path = os.path.join(self.home_path,name)
        self.make_dir(dir_path)
        self.room.setdefault(name, dir_path)

    def init_ecosystem(self):
        self.eco_dic = dict()
        self.eco_json = ''
        _eco_fn = 'ecosystem.json'
        self.eco_fn = os.path.join(self.home_path, _eco_fn)
        self.sync_ecosystem()

    def sync_ecosystem(self):
        self.eco_json = json.dumps(self.eco_dic)
        eco_fh = open(self.eco_fn, 'w')
        json.dump(self.eco_dic, eco_fh, indent=4)
        eco_fh.write('\n')
        eco_fh.close()

    def add_to_ecosystem(self, lvs):
        try:
            lv1, lv2, lv3, lv4 = lvs
        except:
            self.logger.error('add_to_ecosystem input error : {0}'.format(lvs))
            sys.exit()
        if os.path.exists(lv4):
            self.eco_dic.setdefault(lv1, {})
            self.eco_dic[lv1].setdefault(lv2, {})
            self.eco_dic[lv1][lv2].setdefault(lv3, lv4)
        else:
            self.logger.error('{0} is NOT exist {0}'.format(lv4))
            sys.exit()

    def add_stats(self, lv1):
        if lv1 in self.eco_dic:
            cmds = ['touch']
            cmds.append(os.path.join(self.room['stats'],'{0}.done'.format(lv1)))
            proc = Popen(cmds)
            proc.communicate()
        else:
            self.logger.error("{0} is not in ecosystem. stop add_stats.".format(lv1))
            sys.exit()

    def check_stats(self, lv1):
        _stats = [x.split('/')[-1] for x in glob.glob('{0}/*'.format(self.room['stats']))]
        stats = [x.split('.done')[0] for x in _stats]
        if lv1 in stats:
            self.logger.info("Found pre-run {0}".format(lv1))
            return 1
        return 0

    def make_dir(self, dir_path):
        if os.path.exists(dir_path):
            pass
        else:
            os.makedirs(dir_path)

    def find_file_extension(self, fn):
        if fn.endswith('fq.gz'):
            return 'fq.gz'
        elif fn.endswith('fq'):
            return 'fq'
        elif fn.endswith('fastq.gz'):
            return 'fastq.gz'
        elif fn.endswith('fastq'):
            return 'fastq'
        else:
            self.logger.error("Cannot find file extension")
            sys.exit()

    def join_fn_by_dot(self, fns):
        return '.'.join(fns)

def main(mode, config, tools_fn, refs_fn):
    #
    mine = Mine(config)
    mine.add_tools(tools_fn)
    mine.add_refs(refs_fn)
    mine.logger.info("WELCOME MINE :)")
    #
    mine.logger.debug(mine.config_fn)
    mine.logger.debug(mine.species)
    mine.logger.debug(mine.home)
    mine.logger.debug(mine.output)
    mine.logger.debug(mine.tools)
    mine.logger.debug(json.dumps(mine.sample_dic,indent=4))
    mine.logger.debug(json.dumps(mine.idmatch_dic,indent=4))
    mine.logger.debug(json.dumps(mine.dmr_dic,indent=4))
    #
    eco = Eco()
    eco.build_home(mine.output)
    eco.init_ecosystem()
    #
    eco.logger.debug(eco.eco_fn)
    eco.logger.debug(json.dumps(eco.room,indent=4))
    #
    pipes = []
    if mode in ['basic']:
        pipes.append('do_link_raw')
        pipes.append('do_fastqc_raw')
        pipes.append('do_stat_fastq_raw')
        #
        pipes.append('do_filter_fastq')
        pipes.append('do_fastqc_clean')
        pipes.append('do_stat_fastq_clean')
        #
        pipes.append('do_bismark_map')
        pipes.append('do_bismark_dedup')
        pipes.append('do_bismark_call')
        pipes.append('do_bismark_nucl')
        pipes.append('do_bismark_report')
        #
        pipes.append('do_targetcheck_sort')
        pipes.append('do_targetcheck_grep')
        #
        pipes.append('do_metilene_prepare')
        pipes.append('do_metilene_exe')
        pipes.append('do_metilene_parse')
        #
        ## Report
        #
        #pipes.append('report_stat_fastq_raw')
        #pipes.append('report_fastqc_raw')
        #pipes.append('report_stat_fastq_clean')
        #pipes.append('report_fastqc_clean')
        #
    #
    for idx, pipe in enumerate(pipes):
        mine.logger.info('PIPE START ==> {0:0>4}:{1}'.format(str(idx),pipe))
        #
        if pipe in ['do_link_raw']:
            Do_link_raw(pipe, mine, eco)
        elif pipe in ['do_fastqc_raw']:
            Do_fastqc(pipe, mine, eco)
        elif pipe in ['do_stat_fastq_raw']:
            Do_stat_fastq(pipe, mine, eco)
        elif pipe in ['do_filter_fastq']:
            Do_filter_fastq(pipe, mine, eco)
        elif pipe in ['do_fastqc_clean']:
            Do_fastqc(pipe, mine, eco)
        elif pipe in ['do_stat_fastq_clean']:
            Do_stat_fastq(pipe, mine, eco)
        elif pipe in ['do_bismark_map']:
            Do_bismark_map(pipe, mine, eco)
        elif pipe in ['do_bismark_dedup']:
            Do_bismark_dedup(pipe, mine, eco)
        elif pipe in ['do_bismark_call']:
            Do_bismark_call(pipe, mine, eco)
        elif pipe in ['do_bismark_nucl']:
            Do_bismark_nucl(pipe, mine, eco)
        elif pipe in ['do_bismark_report']:
            Do_bismark_report(pipe, mine, eco)
        elif pipe in ['do_targetcheck_sort']:
            Do_targetcheck_sort(pipe, mine, eco)
        elif pipe in ['do_targetcheck_grep']:
            Do_targetcheck_grep(pipe, mine, eco)

        elif pipe in ['do_metilene_prepare']:
            Do_metilene_prepare(pipe, mine, eco)
        elif pipe in ['do_metilene_exe']:
            Do_metilene_exe(pipe, mine, eco)
        elif pipe in ['do_metilene_parse']:
            Do_metilene_parse(pipe, mine, eco)
        else:
            mine.logger.error('Unknown PIPE ==> {0:0>4}:{1}'.format(str(idx),pipe))
            sys.exit()
        #
        mine.logger.info('PIPE END   ==> {0:0>4}:{1}'.format(str(idx),pipe))



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description="MINE [Methylation(for bisulfite-Seq) Analysis Pipeline]")
    parser.add_argument('mode', choices=['build','basic'], default='basic')
    parser.add_argument('config', default='config.in')
    parser.add_argument('--tools-fn', default='/BiO/BioPeople/siyoo/MINE/bin/lib/tools.json')
    parser.add_argument('--refs-fn', default='/BiO/BioPeople/siyoo/MINE/ref/refs.json')
    parser.add_argument('--lib-path', default='/BiO/BioPeople/siyoo/MINE/bin/lib')
    args = parser.parse_args()
    #
    #sys.path.insert(0, args.lib_path)
    #
    main(args.mode, args.config, args.tools_fn, args.refs_fn)
