#!/usr/bin/python

import json

def main(args):

    eco_dic = json.load(open(args.ecosystem))
    #
    if args.name:
        print args.name
        print json.dumps(eco_dic[args.name], indent=4)
    else:
        for key in sorted(eco_dic.keys()):
            print key
    #

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--ecosystem')
    parser.add_argument('--name')
    args = parser.parse_args()
    main(args)
