#!/usr/bin/python

import sys
import os
import commands
import time
import glob

def manageQsub(_qlist, _qname, _logdir, _thread):
    jobs = []
    for qsub in _qlist :
        if _thread == '1' :
            Jname = "{0}/{1}".format(_logdir, qsub.split("/")[-1])
            qsubOutput = commands.getoutput("qsub -q %s -cwd -o %s.stout -e %s.error -S /bin/bash %s" %(_qname, Jname, Jname, qsub))
            job = qsubOutput.split(" ")[2]
            jobs.append(job)
        else :
            Jname = "{0}/{1}".format(_logdir, qsub.split("/")[-1])
            qsubOutput = commands.getoutput("qsub -q %s -cwd -pe smp %s -o %s.stout -e %s.error -S /bin/bash %s" %(_qname, _thread, Jname, Jname, qsub))
            job = qsubOutput.split(" ")[2]
            jobs.append(job)

    while True:
        qstats = commands.getoutput("qstat")
        qstatsList = qstats.split("\n")
        qstatLen = len(qstatsList)
        if qstatLen == 0:
            break ## "no qstat"

        qjobs = []
        for idx in range(2, qstatLen):
            qjob = qstatsList[idx].strip().split(" ")[0]
            qjobs.append(qjob)

        for i in jobs:
            if not i in qjobs:
                jobs.remove(i) ## "%s remove" % i

        if len(jobs) == 0:
            break ## "job end"
        time.sleep(5)

def RunQsub(_oriqlist, _qname, _thread, _logdir, _scdir, _runName) :
    qlist = []
    for num in range(0, len(_oriqlist)) :
        _w = open("{0}/{1}.{2}.sh".format(_scdir, _runName, str(num)), 'w')
        _w.write("#!/bin/sh\n\n# Start\ndate\n{0}\n# End\ndate".format(_oriqlist[num]))
        _w.close()
        qlist.append("{0}/{1}.{2}.sh".format(_scdir, _runName, str(num)))
    manageQsub(qlist, _qname, _logdir, _thread)

