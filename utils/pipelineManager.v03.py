#!/usr/bin/python
import os
import sys
import subprocess
import commands
import getopt
import time
import threading
import smtplib
from email.MIMEText import MIMEText


CMD_DIC = {"cd":"os.chdir('%s')",
           "mkdir":"os.mkdir('%s')",
}


IDX = 0
SCTRIP_LIST = []
MAX_THREADS = None
E_MAIL = None
SUBJECT = "[report] Pipeline Execution Report"
ENGINE = None
MSG = ""

def init():
	global SCRIPT_LIST, MAX_THREADS, E_MAIL, SUBJECT, ENGINE
	options, args = getopt.getopt(sys.argv[1:], "t:e:s:")
	for op, p in options:
		if op == "-t": MAX_THREADS = p
		elif op == "-e": E_MAIL = p
		elif op == "-s": SUBJECT = p
		else: print "Unknown option : ", op
	ENGINE = args[0]
	SCRIPT_LIST = sys.stdin.read().split("\n")

	if MAX_THREADS == None:
		if ENGINE == "SGE":
			MAX_THREADS = getSGACpuCnt()
		elif ENGINE == "SA":
			MAX_THREADS = getSACpuCnt()

def getSGACpuCnt():
	qstatFC = commands.getoutput("qstat -f")
	qstatFList = qstatFC.split("\n")
	cpuCnt = 0
	for qstatF in qstatFList:
		if qstatF.find("all.q") != -1:
			words = qstatF.split()
			cpuCnt += int(words[2].split("/")[-1])
	return cpuCnt

def getSACpuCnt():
	cpuCnt = int(commands.getoutput("grep -c processor /proc/cpuinfo").strip())
	return cpuCnt

def fork():
	global IDX, ENGINE
	#print "#fork s"
	scriptList = []
	while True:
		script = getLine()
		if not script: break
		if "#join" in script: break
		elif "#fork_" in script:
			engine = script.split()[0].split("_")[-1]
			selectFork(engine)
		elif "#begin_"in script:
			engine = script.split()[0].split("_")[-1]
			selectBegin(engine)
		elif "#fork" in script: fork()
		elif "#begin" in script: begin()
		else: 
			#print script + " &&&&&&& " + str(IDX)
			scriptList.append(script)
	#print "fork " ,scriptList
	if ENGINE == "SGE": runSunGridEngine(scriptList)
	elif ENGINE == "SA": runStandAlone(scriptList)
	#print "#join e"

def begin():
	global IDX, ENGINE
	#print "#begin s"
	while True:
		scriptList = []
		script = getLine()
		if not script: break
		if "#end" in script: break
		elif "#fork_" in script:
			engine = script.split()[0].split("_")[-1]
			selectFork(engine)
		elif "#begin_"in script:
			engine = script.split()[0].split("_")[-1]
			selectBegin(engine)
		elif "#fork" in script: fork()
		elif "#begin"in script: begin()
		else:
			#print script + " ------- " + str(IDX)
			scriptList.append(script)
			if ENGINE == "SGE": runSunGridEngine(scriptList)
			elif ENGINE == "SA": runStandAlone(scriptList)
	#print "#end e"

def selectFork(_engine):
	global IDX
	#print "#fork s"
	scriptList = []
	while True:
		script = getLine()
		if not script: break
		if "#join" in script: break
		elif "#fork_" in script:
			engine = script.split()[0].split("_")[-1]
			selectFork(engine)
		elif "#begin_"in script:
			engine = script.split()[0].split("_")[-1]
			selectBegin(engine)
		elif "#fork" in script: fork()
		elif "#begin" in script: begin()
		else:
			#print script + " &&&&&&& " + str(IDX)
			scriptList.append(script)
	#print "fork " ,scriptList
	if _engine == "SGE": runSunGridEngine(scriptList)
	elif _engine == "SA": runStandAlone(scriptList)
	#print "#join e"

def selectBegin(_engine):
	global IDX
	#print "#begin s"
	while True:
		scriptList = []
		script = getLine()
		if not script: break
		if "#end" in script: break
		elif "#fork_" in script:
			engine = script.split()[0].split("_")[-1]
			selectFork(engine)
		elif "#begin_"in script:
			engine = script.split()[0].split("_")[-1]
			selectBegin(engine)
		elif "#fork" in script: fork()
		elif "#begin"in script: begin()
		else:
			#print script + " ------- " + str(IDX)
			scriptList.append(script)
			if _engine == "SGE": runSunGridEngine(scriptList)
			elif _engine == "SA": runStandAlone(scriptList)

def getLine():
	global IDX, SCRIPT_LIST, MSG
	line = None
	if IDX >= len(SCRIPT_LIST): line = None
	else:
		print SCRIPT_LIST[IDX]
		MSG += SCRIPT_LIST[IDX] + "\n"
		line = SCRIPT_LIST[IDX].strip()
	IDX += 1
	return line

def sendMail(_msg):
	global E_MAIL, SUBJECT
	smtpServerAddr = 'smtp.gmail.com'
	smtpServerPort = 587
	gmailUserId    = 'pipelinesmanager'
	gmailUserPw    = 'pipelinesmanager'
	fromMailAddr   = 'pipelinesmanager@gmail.com'

	msg = _msg + "\n\n\npipelines manager v 0.3\nmade by JongSoo-Kim"

	fromMailAddr = fromMailAddr
	utf8Msg = unicode(msg).encode('utf-8')

	mimeText = MIMEText(utf8Msg, _charset = 'utf-8')
	mimeText['Subject'] = SUBJECT
	mimeText['From'] = fromMailAddr
	mimeText['To'] = E_MAIL

	server = smtplib.SMTP(smtpServerAddr, smtpServerPort)
	server.ehlo()
	server.starttls()
	server.ehlo()
	server.login(gmailUserId, gmailUserPw)
	server.sendmail(fromMailAddr, E_MAIL.split(","), mimeText.as_string())
	server.close()

def runThread(_script):
	commands.getoutput("%s" % (_script))
	time.sleep(0)

def runStandAlone(_scriptList):
	global MAX_THREADS, CMD_DIC
	maxCpuCnt = int(MAX_THREADS)
	threads = []
	for script in _scriptList:
		cmd = script.split()[0]
		cmdArg = script.split()[-1]
		if cmd in CMD_DIC:
			eval(CMD_DIC[cmd] % cmdArg)
			continue
		while True:
			curCpuCnt = len(threads)
			#print curCpuCnt, maxCpuCnt
			if curCpuCnt < maxCpuCnt:
				th = threading.Thread(target=runThread, args=(script,))
				th.start()
				threads.append(th)
				break
			else:
				for th in threads:
					if th.isAlive() == False: threads.remove(th)
				time.sleep(5)

	for th in threads:
		th.join()

def runSunGridEngine(_scriptList):
	global MAX_THREADS, CMD_DIC
	compilerDic = {"java" : "/usr/java/latest/bin/",
					"python" : "/usr/bin/",
					"perl" : "/usr/bin/",
	}
	maxCpuCnt = int(MAX_THREADS)
	jobs = []
	for script in _scriptList:
		cmd = script.split()[0]
		cmdArg = script.split()[-1]
		if cmd in CMD_DIC:
			eval(CMD_DIC[cmd] % cmdArg)
			continue
		while True:
			curCpuCnt = len(jobs)
			#print curCpuCnt, maxCpuCnt
			if curCpuCnt < maxCpuCnt:
				compiler = script.split()[0]
				if compiler in compilerDic:
					compilerPath = compilerDic[compiler]
					qsubOutput = commands.getoutput("qsub -cwd -S %s%s" % (compilerPath, script))
				else:
					qsubOutput = commands.getoutput("qsub -cwd -b y '%s'" % (script))
				job = qsubOutput.split(" ")[2]
				jobs.append(job)
				break
			else:
				qjobs = []
				qstat = commands.getoutput("qstat")
				qstatList = qstat.split("\n")
				switch = False
				for line in qstatList:
					if line[:5] == "-----":
						switch = True
						continue
					if switch == True:
						words = line.strip().split()
						qjob = words[0]
						qjobs.append(qjob)

				for job in jobs:
					if not job in qjobs:
						jobs.remove(job)
				time.sleep(5)

	while True:
		qstats = commands.getoutput("qstat")
		qstatsList = qstats.split("\n")
		qstatLen = len(qstatsList)
		if qstatLen == 0:
			#print "no qstat"
			break

		qjobs = []
		for idx in range(2, qstatLen):
			qjob = qstatsList[idx].strip().split(" ")[0]
			qjobs.append(qjob)
		for i in jobs:
			if not i in qjobs:
				jobs.remove(i)

		if len(jobs) == 0:
			#print "job end"
			break
		time.sleep(5)

if __name__ == "__main__":
	if len(sys.argv) == 1:
		print "Usage : python %s [-t nThrds] [-e emailAddr ex)jinsilhanna@gmail.com,jongsoo.kim@therabio.kr] [-s emailSubject] <in.SGE|SA> <-.standard input>" % (sys.argv[0])
		sys.exit()
	else:
		init()

	try:
		MSG += "----------------------\nStart : %s\n----------------------\n" % time.strftime("%Y/%m/%d/%H:%M:%S")
		print "----------------------\nStart : %s\n----------------------\n" % time.strftime("%Y/%m/%d/%H:%M:%S"),
		if E_MAIL != None: sendMail(MSG)
		while True:
			scriptList = []
			script = getLine()
			if not script: break
			if "#fork_" in script:
				engine = script.split()[0].split("_")[-1]
				selectFork(engine)
			elif "#begin_"in script:
				print "------------b", script
				engine = script.split()[0].split("_")[-1]
				selectBegin(engine)
			elif "#fork" in script: fork()
			elif "#begin" in script: begin()
			else:
				#print script + " ------- " + str(IDX)
				scriptList.append(script)
				if ENGINE == "SGE": runSunGridEngine(scriptList)
				elif ENGINE == "SA": runStandAlone(scriptList)
		MSG += "----------------------\nEnd : %s\n----------------------\n" % time.strftime("%Y/%m/%d/%H:%M:%S")
		print "----------------------\nEnd : %s\n----------------------\n" % time.strftime("%Y/%m/%d/%H:%M:%S"),
		if E_MAIL != None: sendMail(MSG)
	except Exception, msg:
		print msg
		MSG += "Error\n-----------------------------------\n" + msg + "\n-----------------------------------\n"
		if E_MAIL != None: sendMail(MSG)
		sys.exit()
