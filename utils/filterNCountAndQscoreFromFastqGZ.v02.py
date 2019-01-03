#!/usr/bin/python
import sys, gzip, subprocess, getopt, os, commands


#variables
BUFFER_SIZE = 10000
PROPER_N_MIN_RATE = 10.0
PROPER_AVERAGE_SCORE = 20.0
PROPER_QUALITY_2_DOWN_RATE = 5.0
PLATFORM_ZERO_DIC = {"Sanger" : -33, "Solexa" : -64, "Illumina 1.3+" : -64, "Illumina 1.5+" : -64, "Illumina 1.8+" : -33}
SEED_LENGTH = -10

FIRST_IN_FASTQ = None
SECOND_IN_FASTQ = None
FIRST_OUT_FASTQ = None
SECOND_OUT_FASTQ = None
RESULT_OUTPUT = None

def init():
	global FIRST_IN_FASTQ, SECOND_IN_FASTQ, FIRST_OUT_FASTQ, SECOND_OUT_FASTQ, RESULT_OUTPUT
	options, args = getopt.getopt(sys.argv[1:], "")
#	for op, p in options:
#		print op, p
#		if op == "-t": THREADS = p
#		elif op == "-l": SEED_LEN = p
#		elif op == "-k": MAX_SEED_DIFF = p
#		elif op == "-B": BWA = p
#		else: print "Unknown option : ", op
	FIRST_IN_FASTQ = args[0]
	SECOND_IN_FASTQ = args[1]
	FIRST_OUT_FASTQ = args[2]
	SECOND_OUT_FASTQ = args[3]
	RESULT_OUTPUT = args[4]
	
def checkFlatform():
	# quality #
	qualityIdcDic = {}
	firReadFile = gzip.open(FIRST_IN_FASTQ, "rb")
	lineIdx = -1
	for line in firReadFile :
		lineIdx += 1
		if lineIdx == 100000 : break
		if lineIdx % 4 != 3 : continue
		for idc in line.strip() : qualityIdcDic[idc] = None
	firReadFile.close()

	secReadFile = gzip.open(SECOND_IN_FASTQ, "rb")
	lineIdx = -1
	for line in secReadFile :
		lineIdx += 1
		if lineIdx == 100000 : break
		if lineIdx % 4 != 3 : continue
		for idc in line.strip() : qualityIdcDic[idc] = None
	secReadFile.close()
	
	platform = None
	qualityIdcDicKeys = qualityIdcDic.keys()
	minQualityIdc = min(qualityIdcDicKeys)
	maxQualityIdc = max(qualityIdcDicKeys)
	if   minQualityIdc == "!" and maxQualityIdc == "I" :
		platform = "Sanger"
	elif minQualityIdc == ";" and maxQualityIdc == "h" :
		platform = "Solexa"
	elif minQualityIdc == "@" and maxQualityIdc == "h" :
		platform = "Illumina 1.3+"
	elif minQualityIdc == 'B' and maxQualityIdc == "h" :
		platform = "Illumina 1.5+"
	elif minQualityIdc == "!" and maxQualityIdc == "J" :
		platform = "Illumina 1.8+"
	else :
		if   minQualityIdc < ";" :
			platform = "Illumina 1.8+"
		elif minQualityIdc < "@" :
			platform = "Solexa"
		elif minQualityIdc < "B" :
			platform = "Illumina 1.3+"
		else :
			platform = "Illumina 1.5+"
	return platform
	
# N count
def filterCountN(seq):
	count = 0
	countLen = len(seq)
	for s in seq:
		if s == "N":
			count += 1
		else:
			pass
	filter = None
	if count * 100.0 / countLen >= PROPER_N_MIN_RATE:
		filter = True
	else:
		filter = False
	return filter

# qulaity filter
def filterQuality(_qualities, _flatform):
	convertValue = PLATFORM_ZERO_DIC[_flatform]
	qualitiesLen = len(_qualities)
	qualityTotal = 0
	quality2Down = 0
	for quality in _qualities:
		qValue = ord(quality) + convertValue
		qualityTotal += qValue
		if qValue <= 2:
			quality2Down += 1
	
	filter = None
	if (qualityTotal * 1.0 / qualitiesLen <= PROPER_AVERAGE_SCORE) or (quality2Down * 100.0 / qualitiesLen >= PROPER_QUALITY_2_DOWN_RATE):
		filter = True
	else:
		filter = False
	return filter

def writeBuffer(fo, bufferList, forcedWrite):
	if forcedWrite == False:
		if len(bufferList) >= BUFFER_SIZE:
			for element in bufferList:
				fo.write("%s\n" % element)
			del bufferList[:]
	else:
		for element in bufferList:
			fo.write("%s\n" % element)
		del bufferList[:]

if __name__ == "__main__":
	if len(sys.argv) != 6:
		print "Usage : python %s <in.first.fq.gz> <in.second.fq.gz> <out.first.fq.gz> <out.second.fq.gz> <out.result.txt>" % sys.argv[0]
		sys.exit()
	else:
		init()
		
	flatform = checkFlatform()
	'''
	fiFirst = gzip.open(FIRST_IN_FASTQ, "rb")
	fiSecond = gzip.open(SECOND_IN_FASTQ, "rb")
	
	foFirst = gzip.open(FIRST_OUT_FASTQ, "wb")
	foSecond = gzip.open(SECOND_OUT_FASTQ, "wb")

	fiFirst = subprocess.Popen("zcat %s" % FIRST_IN_FASTQ, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
	fiSecond = subprocess.Popen("zcat %s" % SECOND_IN_FASTQ, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout

	foFirst = subprocess.Popen("gzip > %s" % FIRST_OUT_FASTQ, shell=True, stdin=subprocess.PIPE).stdin
	foSecond = subprocess.Popen("gzip > %s" % SECOND_OUT_FASTQ, shell=True, stdin=subprocess.PIPE).stdin
	'''

	fiFirst = os.popen("zcat %s" % FIRST_IN_FASTQ)
	fiSecond = os.popen("zcat %s" % SECOND_IN_FASTQ)

	foFirst = os.popen("gzip > %s" % FIRST_OUT_FASTQ, "w")
	foSecond = os.popen("gzip > %s" % SECOND_OUT_FASTQ, "w")

	bufFirstList = []
	bufSecondList = []
	
	totalReadCnt = 0
	properReadCnt = 0
	filterReadCnt = 0

	while True:
		line1_id = fiFirst.readline().rstrip()
		line1_seq = fiFirst.readline().rstrip()[:SEED_LENGTH]
		line1_strand = fiFirst.readline().rstrip()
		line1_quality = fiFirst.readline().rstrip()[:SEED_LENGTH]
		line2_id = fiSecond.readline().rstrip()
		line2_seq = fiSecond.readline().rstrip()[:SEED_LENGTH]
		line2_strand = fiSecond.readline().rstrip()
		line2_quality = fiSecond.readline().rstrip()[:SEED_LENGTH]
		
		if not line1_id or not line2_id:
			break

		totalReadCnt += 1
		if line1_seq == line2_seq:
			filterReadCnt += 1
			continue
		
		if (filterCountN(line1_seq) == True) or (filterCountN(line2_seq) == True):
			filterReadCnt += 1
			continue

		if (filterQuality(line1_quality, flatform) == True) or (filterQuality(line2_quality, flatform) == True):
			filterReadCnt += 1
			continue
		
		properReadCnt += 1
		bufFirstList.append(line1_id)
		bufFirstList.append(line1_seq)
		bufFirstList.append(line1_strand)
		bufFirstList.append(line1_quality)
		bufSecondList.append(line2_id)
		bufSecondList.append(line2_seq)
		bufSecondList.append(line2_strand)
		bufSecondList.append(line2_quality)

		writeBuffer(foFirst, bufFirstList, False)
		writeBuffer(foSecond, bufSecondList, False)
		
		if totalReadCnt % BUFFER_SIZE == 0:
			print str(totalReadCnt)
		
	writeBuffer(foFirst, bufFirstList, True)
	writeBuffer(foSecond, bufSecondList, True)

	fiSecond.close()
	fiFirst.close()
	foFirst.close()
	foSecond.close()
	
	fo = open(RESULT_OUTPUT, "w")
	headBufList = []
	headBufList.append("totalReadCnt")
	headBufList.append("properReadCnt")
	headBufList.append("filterReadCnt")
	headBufList.append("properReadRate")
	fo.write("#%s\n" % "\t".join(headBufList))
	
	valueBufList = []
	valueBufList.append(str(totalReadCnt * 2))
	valueBufList.append(str(properReadCnt * 2))
	valueBufList.append(str(filterReadCnt * 2))
	valueBufList.append("%.2f" % (properReadCnt * 100.0 / totalReadCnt))
	fo.write("%s\n" % "\t".join(valueBufList))
	fo.close()
