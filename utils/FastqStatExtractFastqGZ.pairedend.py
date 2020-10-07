#!/usr/bin/python
import sys, gzip, subprocess, getopt, os, commands

PLATFORM_ZERO_DIC = {"Sanger" : -33, "Solexa" : -64, "Illumina 1.3+" : -64, "Illumina 1.5+" : -64, "Illumina 1.8+" : -33}

FIRST_IN_FASTQ = None
SECOND_IN_FASTQ = None
RESULT_OUTPUT = None

def init():
	global FIRST_IN_FASTQ, SECOND_IN_FASTQ, RESULT_OUTPUT
	options, args = getopt.getopt(sys.argv[1:],"")
	
	FIRST_IN_FASTQ = args[0]
	SECOND_IN_FASTQ = args[1]
	RESULT_OUTPUT = args[2]

def checkFlatform():
        # quality #
        qualityIdcDic = {}
        firReadFile = gzip.open(FIRST_IN_FASTQ, 'rb') 
        lineIdx = -1
        for line in firReadFile :
                lineIdx += 1
                if lineIdx == 100000 : break
                if lineIdx % 4 != 3 : continue
                for idc in line.strip() : qualityIdcDic[idc] = None

	firReadFile.close()

        secReadFile = gzip.open(SECOND_IN_FASTQ, 'rb')
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

# N and CG count
def NandGCCount(_seq) :
	N_count = 0
	GC_count = 0
	countLen = len(_seq)
	for s in _seq :
		if s == "N" :
			N_count += 1
		elif s == "G" or s == "C" :
			GC_count += 1
		else :
			pass
	return N_count, GC_count

# qulity count
def QualityCount(_qualities, _flatform) :
	convertValue = PLATFORM_ZERO_DIC[_flatform]
	qualitiesLen = len(_qualities)
	qualityQ20MoreTotal = 0
	qualityQ30MoreTotal = 0
	qualityQ20MoreRead = 0
	qualityQ30MoreRead = 0
	TotalQValue = 0
	for quality in _qualities :
		qValue = ord(quality) + convertValue
		TotalQValue += qValue
		if qValue >= 20 :
			qualityQ20MoreTotal += 1
		if qValue >= 30 :
			qualityQ30MoreTotal += 1
	if 30 <= (TotalQValue * 1.0 / qualitiesLen) :
		qualityQ30MoreRead += 1
	if 20 <= (TotalQValue * 1.0 / qualitiesLen) :
		qualityQ20MoreRead += 1
	return qualityQ20MoreTotal, qualityQ30MoreTotal, qualityQ20MoreRead, qualityQ30MoreRead

if __name__ == "__main__" :
	if len(sys.argv) != 4 :
		print "Usage : python %s <in.first.fq.gz> <in.second.fq.gz> <out.result.txt>" % sys.argv[0]
		sys.exit()
	else :
		init()

	flatform = checkFlatform()

#	fiFirst = os.popen("zcat %s" % FIRST_IN_FASTQ)
	fiFirst = gzip.open(FIRST_IN_FASTQ, 'rb')
#	fiSecond = os.popen("zcat %s" % SECOND_IN_FASTQ)
	fiSecond = gzip.open(SECOND_IN_FASTQ, 'rb')

	totalReadCnt = 0
	totalLen = 0
	totalGCCnt = 0
	totalQ30AboveReadCntR1 = 0
	totalQ30AboveReadCntR2 = 0
	totalQ20AboveReadCntR1 = 0
	totalQ20AboveReadCntR2 = 0
	totalQ30AboveBaseCntR1 = 0
	totalQ30AboveBaseCntR2 = 0
	totalQ20AboveBaseCntR1 = 0
	totalQ20AboveBaseCntR2 = 0
	totalNCnt = 0
	totalNZeroCnt = 0
	totalN5BelowCnt = 0
	read1Lengh = 0
	read2Lengh = 0

	while True :
		line1_id = fiFirst.readline().rstrip()
		line1_seq = fiFirst.readline().rstrip()
		line1_strand = fiFirst.readline().rstrip()
		line1_quality = fiFirst.readline().rstrip()
		line2_id = fiSecond.readline().rstrip()
		line2_seq = fiSecond.readline().rstrip()
		line2_strand = fiSecond.readline().rstrip()
		line2_quality = fiSecond.readline().rstrip()

		if not line1_id or not line2_id :
			break

		totalReadCnt += 2
		totalLen += len(line1_seq)
		read1Lengh += len(line1_seq)
		totalLen += len(line2_seq)
		read2Lengh += len(line2_seq)
		nCount1, gcCount1 = NandGCCount(line1_seq)
		nCount2, gcCount2 = NandGCCount(line2_seq)
		if nCount1 + nCount2 == 0 :
			totalNZeroCnt += 2
		if ((nCount1+nCount2) * 1.0 / (len(line1_seq)+len(line2_seq))) < 0.05 :
			totalN5BelowCnt += 2
		totalNCnt += nCount1 + nCount2
		totalGCCnt += gcCount1 + gcCount2
		
		check_QB201, check_QB301, check_QR201, check_QR301  = QualityCount(line1_quality, flatform)
		check_QB202, check_QB302, check_QR202, check_QR302 = QualityCount(line2_quality, flatform)
		totalQ30AboveReadCntR1 += check_QR301
		totalQ30AboveReadCntR2 += check_QR302
		totalQ20AboveReadCntR1 += check_QR201
		totalQ20AboveReadCntR2 += check_QR202
		totalQ30AboveBaseCntR1 += check_QB301
		totalQ30AboveBaseCntR2 += check_QB302
		totalQ20AboveBaseCntR1 += check_QB201
		totalQ20AboveBaseCntR2 += check_QB202

	fo = open(RESULT_OUTPUT, "w")
	headBufList = []
	headBufList.append("totalReadCnt")
	headBufList.append("totalLength")
	headBufList.append("totalGCCnt")
	headBufList.append("NZeroIncludeReadCnt")
	headBufList.append("N5IncludeReadCnt")
	headBufList.append("totalNCnt")
	headBufList.append("totalQ30ReadCntR1")
	headBufList.append("totalQ30ReadCntR2")
	headBufList.append("totalQ30BaseCntR1")
	headBufList.append("totalQ30BaseCntR2")
	headBufList.append("totalQ20ReadCntR1")
	headBufList.append("totalQ20ReadCntR2")
	headBufList.append("totalQ20BaseCntR1")
	headBufList.append("totalQ20BaseCntR2")
	headBufList.append("Read1Length")
	headBufList.append("Read2Length")
	fo.write("#%s\n" % "\t".join(headBufList))

	valueBufList = []
	valueBufList.append(str(totalReadCnt))
	valueBufList.append(str(totalLen))
	valueBufList.append(str(totalGCCnt))
	valueBufList.append(str(totalNZeroCnt))
	valueBufList.append(str(totalN5BelowCnt))
	valueBufList.append(str(totalNCnt))
	valueBufList.append(str(totalQ30AboveReadCntR1))
	valueBufList.append(str(totalQ30AboveReadCntR2))
	valueBufList.append(str(totalQ30AboveBaseCntR1))
	valueBufList.append(str(totalQ30AboveBaseCntR2))
	valueBufList.append(str(totalQ20AboveReadCntR1))
	valueBufList.append(str(totalQ20AboveReadCntR2))
	valueBufList.append(str(totalQ20AboveBaseCntR1))
	valueBufList.append(str(totalQ20AboveBaseCntR2))
	valueBufList.append(str(read1Lengh))
	valueBufList.append(str(read2Lengh))
	fo.write("%s\n" % "\t".join(valueBufList))
	fo.close()

