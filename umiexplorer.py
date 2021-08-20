'''
Created on Mar 29, 2021

@author: S. Makhortov, A. Maslov
'''

import os
#import sys
import datetime
import pysam

from params import Params
from utils import Utils
import family

class UMIExplorer(object):

    def __init__(self):
        self.refData = None
        self.alnGData = None
        self.vcfData = None

        self.btResData = None
        self.vcfResData = None
        self.alnResData = None
        
        self.families = []   # Dynamic family list

        self.nTotalBases = 0 # positions covered
        self.nQualBases = 0  # positions qualified
        self.sStatTemplate = "##Stat:TB=%d;QB=%d"

    @staticmethod
    def umiFromName(readName, lenCalc = False):
        sUMI = readName.split(Params.NAME_DELIM)[-1]
        if lenCalc:
            Params.umiLen = (len(sUMI) - 1) / 2
        return sUMI
    
    # Base search at reference position (0-based) in the fullformed families       
    def testFamPosBase(self, chBase, nRefStart, famBesides = None):
        myFam = None
        for nFam in range(0, len(self.families)):
            fam = self.families[nFam] # Current family
            if fam == famBesides: continue # Besides this fam 
            if not fam.isClosing: break # End of closing
            chFamBase = fam.getBaseAt(nRefStart)
            if chFamBase == chBase:
                myFam = fam
                break
        return myFam

    # Back track fullformed families at the reference position (0-based)       
    def btFamPosition(self, nRefStart, famBesides = None):
        for nFam in range(0, len(self.families)):
            fam = self.families[nFam] # Current family
            if fam == famBesides: continue # Besides this fam 
            if not fam.isClosing: break # End of closing
            fam.backTrackPos(self.btResData, nRefStart)

    # Close down old families        
    def umiFinishFamilies(self, nCurrSBamPos = None, bLog = False):
        posCoverSet = set() # Positions covered
        if bLog: 
            Utils.logOut("FamCount = {}, nCurrSBamPos = {}\n".format(len(self.families), nCurrSBamPos), True)
        nCloseStop = len(self.families)
        for nFam in range(0, nCloseStop):
            fam = self.families[nFam] # Current family
            if  (nCurrSBamPos) and (nCurrSBamPos <= fam.nEnd):
                nCloseStop = nFam
                break # Not old Family
                '''
                # Search in the next
                continue
                '''
            if bLog: 
                Utils.logOut("CurrFam = {} ReadCount = {} fam.nStart = {} fam.nEnd = {}\
                              \n".format(fam.sUMI, fam.cntPositive + fam.cntNegative, fam.nStart, fam.nEnd), True)
            self.families[nFam].isClosing = True # Mark for closing
            #self.families.pop(nFamDel)
            
            if Params.command == Params.cmdUMI: 
                fam.umiElection()
                fam.printOut()
            elif Params.command == Params.cmdRMDUP:
                for myRead in fam.myReads:
                    self.alnResData.write(myRead.origRead)

            if Params.command != Params.cmdRMDUP:
                posCoverSet.update(range(fam.nStart, fam.nEnd + 1)) # Positions covered by family

            if Params.logFamilies: fam.logOut()

        if Params.command == Params.cmdVC:
            for nFam in range(0, nCloseStop):
                fam = self.families[nFam]
                if fam.isQualified(): # Family is qualified
                    self.nQualBases += fam.varCalling(self) #.refData, self.alnGData, self.vcfData, self.vcfResData, btResData)

        for nFam in reversed(range(0, nCloseStop)):
            fam = self.families.pop(nFam)
            del fam

        self.nTotalBases += len(posCoverSet) # Positions covered count
        posCoverSet.clear()
        
    def umiAnalyze(self):
#        origStdOut = sys.stdout
#        fstd = open('out.txt', 'w')
#        sys.stdout = fstd
        #print datetime.datetime.now()
        alnSData = pysam.AlignmentFile(Params.bamSFile, 'rb')  # @UndefinedVariable
        
        if Params.logFamilies:
            Utils.logFile = open(Params.bamSFile + '.log', 'w')

        if Params.command == Params.cmdVC:
            fPath = os.path.dirname(Params.bamSFile)
            fName = os.path.basename(Params.bamSFile)
            fName = fName.rsplit('.', 1)[0] + '.' + Params.chrName # remove ".bam", add Chr

            btFile = os.path.join(fPath, fName + '.backtrack')
            vcfResFile = os.path.join(fPath, fName + '.vcf.tmp')
            
            self.refData = pysam.FastaFile(Params.refFile) # @UndefinedVariable
            self.alnGData = pysam.AlignmentFile(Params.bamGFile, 'rb')  # @UndefinedVariable
            self.vcfAllData = pysam.VariantFile(Params.vcfAllFile, 'r')  # @UndefinedVariable , index_filename = Params.vcfAllFile + '.tbi'
            self.vcfPersData = pysam.VariantFile(Params.vcfPersFile, 'r')  # @UndefinedVariable , index_filename = Params.vcfPersFile + '.tbi'

            self.btResData = open(btFile, "w+")
            self.vcfResData = open(vcfResFile, "w+")
            strOut = "##Params:BQ=%d;MapQ=%d;UC=%d;SC=%d;FC=%d;MC=%d;BD=%d\n" % \
                              (Params.minBaseQuality, Params.minMapQuality, Params.minUmiComplexity, \
                               Params.minStrandCount, Params.minFamilyCount, \
                               Params.maxCoincidence, Params.minBreakDistance)
            self.vcfResData.write(strOut)
            nHeaderPos = self.vcfResData.tell()
            strOut = self.sStatTemplate % (self.nTotalBases, self.nQualBases) + ' ' * 20 + '\n'
            self.vcfResData.write(strOut)
            self.vcfResData.flush()
            self.lstGerml = []
        elif Params.command == Params.cmdRMDUP:
          
            fName = os.path.splitext(Params.bamSFile)[0]
            if Params.chrName:
                fName = fName + '.' + Params.chrName
            bamResFile = fName + '.rmdup.bam' 

            self.alnResData = pysam.AlignmentFile(bamResFile, 'wb', template = alnSData)  # @UndefinedVariable

        currChrom, currFam = Params.chrName, None
        for read in alnSData.fetch(Params.chrName, until_eof = True):
            #if read.reference_start > 1000000: break
            # Testing for good read
            if read.is_unmapped: continue
            if (Params.command == Params.cmdVC) and \
                not (read.is_proper_pair and 'M' in read.cigarstring and not read.is_secondary): 
                continue
            if (Params.command != Params.cmdRMDUP) and (read.mapping_quality < Params.minMapQuality): #   
                continue

            if (currChrom == None) or (read.reference_name != currChrom):
                self.umiFinishFamilies() # Finish chromosome families rest
                currFam = None
                currChrom = read.reference_name
                if Params.logFamilies:
                    myStr = "*** currChrom = {}\n".format(currChrom) 
                    Utils.logOut(str(datetime.datetime.now()) + ": " + myStr, True)
                    print myStr

            sUMI = UMIExplorer.umiFromName(read.query_name, True)
            nStart, nEnd = read.reference_start, read.reference_end
            bLog = False #140319680 in range(nStart, nEnd)
            if bLog: 
                Utils.logOut("*** currRead = {} nStart = {} nEnd = {}\
                             \n".format(read.query_name, nStart, nEnd), True)
            if  (currFam == None) or (not currFam.isMyRead(sUMI, nStart, nEnd)):
                self.umiFinishFamilies(nStart, bLog) # Finish old families
                currFam = None
                for fam in reversed(self.families): # Family search
                    if fam.isMyRead(sUMI, nStart, nEnd):
                        currFam = fam
                        break

                if currFam == None: # Not found - create new 
                    bCreateOk = True
                    if Params.command == Params.cmdRMDUP:
                        bCreateOk = self.testFamCreating(nStart, nEnd)
                                
                    if bCreateOk:
                        currFam = family.Family(currChrom, sUMI, nStart, nEnd)
                        self.families.append(currFam)
 
            if Params.command == Params.cmdRMDUP:
                if currFam: 
                    currFam.replRead(read, sUMI) # nPosChange = 
                    #self.testSorted(currFam, nPosChange) # Start pos may be changed

            else: currFam.addRead(read, sUMI)

        self.umiFinishFamilies() # Finish families rest

        if Params.command == Params.cmdUMI:
            print "TotalBases = ", self.nTotalBases
        elif Params.command == Params.cmdRMDUP:
            self.alnResData.close()
        elif Params.command == Params.cmdVC:   
            # Stat writing
            strOut = self.sStatTemplate % (self.nTotalBases, self.nQualBases)
            self.vcfResData.seek(nHeaderPos)
            self.vcfResData.write(strOut)
            self.vcfResData.close()
            del self.lstGerml
            # Rename res file
            os.rename(vcfResFile, vcfResFile[0:len(vcfResFile)-4])
            
            self.btResData.close()

            self.vcfPersData.close()
            self.vcfAllData.close()
            self.alnGData.close()
            self.refData.close()

        if Params.logFamilies: Utils.logFile.close() 
        alnSData.close()
        #print datetime.datetime.now()
#        if origStdOut:
#            sys.stdout = origStdOut
#            fstd.close()

    def testFamCreating(self, nStart, nEnd):
        bCreateOk = True
        if len(self.families) < Params.maxPosFamilies:
            return bCreateOk

        nStartCover, nEndCover = 0, 0
        for fam in self.families: # Family iterate
            if nStart in range(fam.nStart, fam.nEnd): nStartCover += 1
            if nEnd - 1 in range(fam.nStart, fam.nEnd): nEndCover += 1
        
        bCreateOk = (nStartCover < Params.maxPosFamilies) or (nStartCover < Params.maxPosFamilies)
        return bCreateOk

    def testSorted(self, fam, nChangeSign):
        nIdx = self.families.index(fam)
        if nChangeSign > 0: # Is higher
            for i in range(nIdx, len(self.families) - 1):
                if self.families[i].nStart <= self.families[i+1].nStart: break
                tmp = self.families[i]
                self.families[i] = self.families[i+1]
                self.families[i+1] = tmp
        elif nChangeSign < 0: # Is lower
            for i in reversed(range(1, nIdx + 1)):
                if self.families[i].nStart >= self.families[i-1].nStart: break
                tmp = self.families[i]
                self.families[i] = self.families[i-1]
                self.families[i-1] = tmp
    
    def bamTrunc(self, nPosMin, nPosMax):
#        origStdOut = None #sys.stdout
#        fstd = open('out.txt', 'w')
#        sys.stdout = fstd
        print datetime.datetime.now()
        print "Chr = {} nPosStart = {} nPosStop = {}".format(Params.chrName, nPosMin, nPosMax)
        alnSData = pysam.AlignmentFile(Params.bamSFile, 'rb')  # @UndefinedVariable
        
        fName = os.path.splitext(Params.bamSFile)[0]
        bamResFile = fName + '.trunc.bam' 
#            print bamResFile
        alnResData = pysam.AlignmentFile(bamResFile, 'wb', template = alnSData)  # @UndefinedVariable

        for read in alnSData.fetch(Params.chrName, nPosMin, nPosMax, until_eof = True):
            alnResData.write(read)
#           if nPosLeave:
#               if read.reference_start > nPosLeave: break
            
        alnSData.close()
        alnResData.close()
        pysam.index(bamResFile)   # @UndefinedVariable
        print datetime.datetime.now()
