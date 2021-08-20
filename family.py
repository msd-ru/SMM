'''
Created on Mar 30, 2021

@author: S. Makhortov, A. Maslov
'''

import sys
import datetime

from params import Params
from utils import Utils

class Read(object):

    def __init__(self, aFam, read, sUMI, isPositive, isOverlapped = False):
        '''
        Constructor
        '''
        self.fam = aFam # My Family
        self.origRead = read
        self.sUMI = sUMI
        self.isPositive = isPositive
        self.isOverlapped = isOverlapped # True if is second in the pair
        
        self.myBrother = None # Overlapped read index

        if isPositive: self.sUMIMod = sUMI
        else: self.sUMIMod = Utils.umiSwap(sUMI)
        
        self.baseList = None
            
    def sameOrientation(self, read):
        res = (self.origRead.is_reverse == read.is_reverse) and (self.origRead.is_read1 == read.is_read1)
        res = res or (self.origRead.is_reverse != read.is_reverse) and (self.origRead.is_read1 != read.is_read1)
        return res
    
    def calcBaseList(self):
        resList = [] # List of string
        read = self.origRead
        
        nReadPos = self.fam.nStart - read.reference_start
        if nReadPos < 0:
            resList = [None] * (-nReadPos)
            nReadPos = 0

        flagM = 0  # Indicates that M section was not processed yet
        for i in range(len(read.cigartuples)):
            nSect, nCount = read.cigartuples[i][0], read.cigartuples[i][1]

            if nSect == Params.SectM: # Put bases (M section)
                for j in range(nReadPos, nReadPos + nCount):
                    if read.query_qualities[j] < Params.minBaseQuality:
                        resList += Params.chNoQual
                    else:    
                        resList += read.query_sequence[j]
                nReadPos += nCount
                flagM = 1
            elif nSect == Params.SectI: # Put inserted bases (I section)
                if flagM: # if I is after M then
                    resList[-1] += read.query_sequence[nReadPos : nReadPos + nCount] # Add insertion to last base
                # Otherwise treat it as soft-clipped
                nReadPos += nCount

            elif nSect == Params.SectD: # Put 'D' if deleted (D section)
                resList += [Params.chDel] * nCount

            elif nSect == Params.SectS: # Skip soft-clipped part (S section)
                nReadPos += nCount
                
            elif nSect == Params.SectH: # Ignore hard-clipped part (H section)
                pass
            else:
                sys.exit('Unrecognized section code in Cigar tuple:\t%s' %nSect)
                
        nRefLen = self.fam.nEnd - self.fam.nStart # 0-based, nEnd - behind the  
        nMyLen = len(resList)
        
        # Set length
        if nMyLen < nRefLen: resList += [None] * (nRefLen - nMyLen)
        else: del resList[nRefLen:]
       
        self.baseList = resList

    def printOut(self):
        read = self.origRead
        myStr = "{}\t{}\t{}".format(read.query_name, read.reference_start, read.reference_end)
        print myStr

    def logOut(self):
        read = self.origRead
        myStr = "{}\t{}\t{}".format(read.query_name, read.reference_start, read.reference_end)
        Utils.logOut(myStr, True)

    def getBaseAt(self, nRefStart):
        sBase = None
        if nRefStart in range(self.fam.nStart, self.fam.nEnd):
            if not self.baseList: self.calcBaseList() 
            sBase = self.baseList[nRefStart - self.fam.nStart]
        return sBase

    def backTrackPos(self, btResData, nRefStart):
        chReadBase = self.getBaseAt(nRefStart)
        rd = self.origRead
        # Read info string
        outString = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chReadBase, int(self.isPositive), rd.query_name, \
                                            rd.reference_start, rd.reference_end, rd.next_reference_start, rd.cigarstring, rd.mapping_quality) 
        btResData.write(outString) 

class Family(object):

    def __init__(self, aChrom, sUMI, nStart, nEnd):
        '''
        Constructor
        '''
        self.currChrom = aChrom
        self.sUMI = sUMI
        self.nStart = nStart
        self.nEnd = nEnd

        self.myReads = []     # Records of all family reads
        self.cntPositive = 0
        self.cntNegative = 0

        self.isClosing = False
        self.consList = None # Consensus list
        self.consListP = None # Consensus list +
        self.consListN = None # Consensus list -
        
#        self.posQualList = [] # Positions qualified
        self.allNames = []  # All read names - for rmdup command
        
    def isQualified(self):
        return ((Params.minStrandCount == 0) or (self.cntPositive >= Params.minStrandCount) and (self.cntNegative >= Params.minStrandCount)) and \
                ((Params.minFamilyCount == 0) or (self.cntPositive + self.cntNegative >= Params.minFamilyCount))
            
    def isMyRead(self, sUMI, nStart, nEnd):
#        print 'Fam: ', self.sUMI, self.nStart, self.nEnd
#        print 'New: ', sUMI, nStart, nEnd

        isMy = (Utils.umiDistance(sUMI, self.sUMI) <= Params.maxUMIDist)
        if isMy:
            isMy = isMy and (Utils.segOverlap(nStart, nEnd, self.nStart, self.nEnd, Params.maxPosDist))
                    
        return isMy 

    def replRead(self, read, sUMI): # For rmdup command
        newRead = read
        newName = Utils.pureName(newRead.query_name) # Name without UMI
        
        # Test for brother read overlap (1,2)
        isOverlap = False
        nCount = len(self.myReads)
        if nCount == 1:
            currRead = self.myReads[0].origRead
            if newName == Utils.pureName(currRead.query_name):
                isOverlap = True #Utils.segOverlap(newRead.reference_start, newRead.reference_end, currRead.reference_start, currRead.reference_end)
                if isOverlap: # 2 reads are overlapping - store both brothers
                    self.myReads[0].myBrother = nCount # Store new index to the brother
        
        if (not isOverlap) and nCount: # Competitor testing
            bNewBetter = False
            
            if not newName in self.allNames: # Is it not second brother for the lost read 
                for myRead in self.myReads:
                    if newRead.mapping_quality > myRead.origRead.mapping_quality:
                        bNewBetter = True
                        break
                    
            if bNewBetter: # New read is better - remove all older
                while len(self.myReads): 
                    dr = self.myReads.pop(0)
                    del dr
            else: newRead = None # New read is not better
        
        if newRead:
            myRead = Read(self, newRead, sUMI, isOverlap)
            self.myReads.append(myRead)
        
        self.allNames.append(newName) # For the next brother testing

        # Family pos correction
        nReadStart, nReadEnd = read.reference_start, read.reference_end
        if (nReadStart < self.nStart): self.nStart = nReadStart
        if (nReadEnd > self.nEnd): self.nEnd = nReadEnd

    def addRead(self, read, sUMI):
        if not len(self.myReads): isPositive = True 
        else: isPositive = self.myReads[0].sameOrientation(read)
        
        nReadStart, nReadEnd = read.reference_start, read.reference_end
        
        # Find read overlap (1,2)
        myName = Utils.pureName(read.query_name) # Name without UMI
        nCount = len(self.myReads)
        iBrother1 = None
        if Params.testOverlapPair: # Need to test for overlapping pairs ### 
            for i in range(0, nCount):
                currRead = self.myReads[i].origRead
                if Utils.pureName(currRead.query_name) == myName:
                    # isOverlap = Utils.segOverlap(nReadStart, nReadEnd, currRead.reference_start, currRead.reference_end)
                    #if isOverlap: # Reads overlapping
                    iBrother1 = i
                    self.myReads[iBrother1].myBrother = nCount # Store my index to the brother
                    break

        myRead = Read(self, read, sUMI, isPositive, iBrother1 != None)
        if iBrother1 == None: # Not overlapped pair in the one family
            if isPositive: self.cntPositive += 1 
            else: self.cntNegative += 1
        else: myRead.myBrother = iBrother1
        self.myReads.append(myRead)
        
        # Family pos correction
        if (nReadStart < self.nStart): self.nStart = nReadStart
        if (nReadEnd > self.nEnd): self.nEnd = nReadEnd

    def umiElection(self):
        # UMI election
        umiList = []
        nCount = len(self.myReads)
        for i in range(0, nCount):
            umiList.append(self.myReads[i].sUMIMod)

        iMax = Utils.mostCommon(umiList)
        self.sUMI = umiList[iMax] # more often met UMI
        
        del umiList[:]

    def logOut(self, bReads = False):
        myStr = "Close family: UMI={} nStart={} nEnd={}\n".\
                 format(self.sUMI, self.nStart, self.nEnd)
        Utils.logOut(myStr, True)

        if bReads:
            for myRead in self.myReads:
                myRead.logOut()
    
    # Positions qualified calculation
    def posQualified(self, alnGData):
        self.posQualList = []

        for plpCol in alnGData.pileup(Params.chrName, self.nStart, self.nEnd):
            if plpCol.nsegments >= Params.minUmiComplexity:
                self.posQualList.append(plpCol.reference_pos)

        return len(self.posQualList)
        
    def calcConsensusList(self):

        nQualCount = 0 # Min count for pos qualified 
        if Params.minStrandCount > 0:
            nQualCount = Params.minStrandCount + Params.minStrandCount
        if Params.minFamilyCount > 0:
            nQualCount = max(nQualCount, Params.minFamilyCount)
        
        ll = [] # List of base list
        nReadCount = len(self.myReads)
        for i in range(0, nReadCount):
            myRead = self.myReads[i]
            if not myRead.baseList: myRead.calcBaseList()
            ll.append(myRead.baseList) 

        self.consList = []
        self.consListP, self.consListN = None, None
        
        nReadLen = self.nEnd - self.nStart
        strListCol = [""] * nReadCount
        #strListColP, strListColN = [], []
        for k in range(0, nReadLen):
            for i in range(0, nReadCount):
                currBase = ll[i][k]
                myRead = self.myReads[i]
                if myRead.isOverlapped: # Is the secondary read
                    brother = self.myReads[myRead.myBrother]
                    if brother.baseList[k] != None:
                        currBase = None # The second is not taken into account 
                strListCol[i] = currBase

            # Common consensus
            iMax = Utils.mostCommon(strListCol, nMinCount = nQualCount)
            if iMax < 0: sCons = Params.chNoQual 
            else: sCons = strListCol[iMax]  # more often met base
            self.consList.append(sCons)

        del strListCol[:]
 
    def calcConsensusListPN(self):
        ll = [] # List of base list
        nReadCount = len(self.myReads)
        for i in range(0, nReadCount):
            myRead = self.myReads[i]
            if not myRead.baseList: myRead.calcBaseList()
            ll.append(myRead.baseList) 

        self.consListP, self.consListN = [], []
        
        nReadLen = self.nEnd - self.nStart
        strListColP, strListColN = [], []
        for k in range(0, nReadLen):
            for i in range(0, nReadCount):
                currBase = ll[i][k]
                myRead = self.myReads[i]
                if myRead.isOverlapped: # Is the secondary read
                    brother = self.myReads[myRead.myBrother]
                    if brother.baseList[k] != None:
                        currBase = None # The second is not taken into account 
                if myRead.isPositive: strListColP.append(currBase)
                else: strListColN.append(currBase)
            # Consensus +
            iMax = Utils.mostCommon(strListColP, nMinCount = Params.minStrandCount)
            if iMax < 0: sCons = Params.chNoQual 
            else: sCons = strListColP[iMax]  # more often met base
            self.consListP.append(sCons)
            # Consensus -
            iMax = Utils.mostCommon(strListColN, nMinCount = Params.minStrandCount)
            if iMax < 0: sCons = Params.chNoQual 
            else: sCons = strListColN[iMax]  # more often met base
            self.consListN.append(sCons)
        del strListColP[:]
        del strListColN[:]

    # Base in consensus list, refpos 0-based
    def getBaseAt(self, nRefStart):
        sBase = None
        if nRefStart in range(self.nStart, self.nEnd):
            if not self.consList: self.calcConsensusList()
            sBase = self.consList[nRefStart - self.nStart]
        return sBase

    # True if it is break distance
    def testBreakDist(self, nRefStart):
        isBreak = (nRefStart - self.nStart) <= Params.minBreakDistance
        isBreak = isBreak or (self.nEnd - nRefStart <= Params.minBreakDistance) 
        if not isBreak:
            if not self.consList: self.calcConsensusList()
            k = nRefStart - self.nStart
            for i in range(1, Params.minBreakDistance + 1):
                if (self.consList[k-i] == Params.chDel) or \
                    len(self.consList[k-i]) > 1: # Left test: DEL/INS
                    isBreak = True
                    break
                if (self.consList[k+i] == Params.chDel) or \
                    len(self.consList[k+i]) > 1: # Right test: DEL/INS
                    isBreak = True
                    break
        
        return isBreak

    # Calc the consensus char break distance in the Family 
    # (True if it is really break distance)
    def getBreakDist(self, nRefStart, chBase):
        nToLeft, nToRight = None, None 
        k = nRefStart - self.nStart # Relatively pos in all my reads base lists
        nReadCount = len(self.myReads)

        for i in range(0, nReadCount):
            myRead = self.myReads[i]
            if myRead.baseList[k] != chBase: continue # Is not my base char
            if myRead.isOverlapped: # Is the secondary read
                brother = self.myReads[myRead.myBrother]
                if brother.baseList[k] != None: continue # The second is not taken into account
            # Left dist calc
            nDist = nRefStart - myRead.origRead.reference_start
            if (nDist > nToLeft): nToLeft = nDist 
            # Right dist calc
            nDist = myRead.origRead.reference_end - nRefStart
            if (nDist > nToRight): nToRight = nDist 

        nDist = min(nToLeft, nToRight)
        isBreak = nDist <= Params.minBreakDistance 
        if not isBreak:
            if not self.consList: self.calcConsensusList()
            for i in range(1, Params.minBreakDistance + 1):
                if (self.consList[k-i] == Params.chDel) or \
                    len(self.consList[k-i]) > 1: # Left test: DEL/INS
                    isBreak = True
                    break
                if (self.consList[k+i] == Params.chDel) or \
                    len(self.consList[k+i]) > 1: # Right test: DEL/INS
                    isBreak = True
                    break
        
        return isBreak, nDist

    # True if variant position is strange
    def isStrangeVar(self, nVarPos):
        bStrange = False

        for myRead in self.myReads:
            read = myRead.origRead
            if nVarPos == read.next_reference_start - 1 or nVarPos == read.reference_end:
                bStrange = True
                break
        
        return bStrange


    def backTrackPos(self, btResData, nRefStart, isMarked = False):
        if not nRefStart in range(self.nStart, self.nEnd):
            return
        # Consensus base calc
        sCons = self.getBaseAt(nRefStart)
        # Family title
        outString = ""
        #if isMarked: outString += "*** "
        
        outString += "{}\t{} ({}-{})".format(self.sUMI, self.cntPositive + self.cntNegative, self.cntPositive, self.cntNegative)
        outString += "\t{}\t{}\t\n".format(self.cntPositive > 0 and self.cntNegative > 0, sCons) 
        btResData.write(outString) 
        # Family reads
        for myRead in self.myReads:
            myRead.backTrackPos(btResData, nRefStart)

    # Base test at reference position (0-based) in my reads (all identical?)       
    def testReadsPosBase(self, nRefStart, isPositive = None):
        chTestBase = None
        for myRead in self.myReads:
            chReadBase = myRead.getBaseAt(nRefStart)
            if chReadBase == Params.chNoQual: continue # No quality
            if (isPositive != None) and (myRead.isPositive != isPositive):
                continue # Is other strand
            if chTestBase == None: # The first calc
                chTestBase = chReadBase
                continue 
            if chReadBase != chTestBase: # Not equal founded
                chTestBase = None
                break
        return chTestBase # 100% Consensus or None

    # As string list
    def getPileupColumn(self, plpColumn):
        strListCol = []
        for plpRead in plpColumn.pileups:
            if plpRead.is_refskip: continue # query position is None if is_refskip is set
            if plpRead.is_del: # Coming soon
                continue            
            chBase = plpRead.alignment.query_sequence[plpRead.query_position]
            strListCol.append(chBase)
        return strListCol
    
    def varCalling(self, umiExp): #refData, alnGData, vcfAllData, vcfPersData, vcfResData):
        nQualBases = 0 # return value
        strRef = umiExp.refData.fetch(Params.chrName, self.nStart, self.nEnd).upper()

        vcfAllRecs = []
        for dbRec in umiExp.vcfAllData.fetch(Params.chrName, self.nStart, self.nEnd):
            vcfAllRecs.append(dbRec)
        vcfPersRecs = []
        for dbRec in umiExp.vcfPersData.fetch(Params.chrName, self.nStart, self.nEnd):
            vcfPersRecs.append(dbRec)
        if not self.consList: self.calcConsensusList()
        varRecs, nSomas = [], 0
        for plpCol in umiExp.alnGData.pileup(Params.chrName, self.nStart, self.nEnd, ignore_overlaps=True, truncate=True):
            #if plpCol.reference_pos < self.nStart: continue
            #if plpCol.reference_pos >= self.nEnd: break
            if plpCol.nsegments < Params.minUmiComplexity: continue # is not qualified

            nQualBases += 1
            k = plpCol.reference_pos - self.nStart # 0-based!
            
            baseS = self.consList[k] # [0]
            if baseS == Params.chNoQual: continue # No quality/qualified
            if len(baseS) != 1: continue # Coming soon (InDels)
            if baseS == Params.chDel: continue # Coming soon
            
            if baseS == strRef[k]: # same as in Reference
                # Test for presence in G-VCF (Person) file
                varID, altBases = self.getDBVariant(vcfPersRecs, plpCol.reference_pos)
                if varID != None: # Is found
                    altBases = ", ".join(altBases) # Tuple to comma string
                    keyGerml = [Params.chrName, plpCol.reference_pos, varID, strRef[k]] 
                    if keyGerml in umiExp.lstGerml: continue # Is duplicate
                    umiExp.lstGerml.append(keyGerml)
                    varRecs.append([plpCol.reference_pos, varID, strRef[k], altBases, Params.vkSNV, None, Params.vtGerm]) # varType must be the last

            else: # differs from the Reference
                btAddInfo = ";GC={};SCP={};SCN={};".format(plpCol.nsegments, self.cntPositive, self.cntNegative)  # Info to the backtracking
                # a)
                varID, altBases = self.getDBVariant(vcfAllRecs, plpCol.reference_pos, baseS) # , Params.vkSNV
                fName = Params.vcfAllFile
                isGermline = (varID != None)
                if not isGermline:  # b)
                    varID, altBases = self.getDBVariant(vcfPersRecs, plpCol.reference_pos, baseS) # , Params.vkSNV
                    fName = Params.vcfPersFile
                    isGermline = (varID != None)
                if isGermline:
                    keyGerml = [Params.chrName, plpCol.reference_pos, varID, strRef[k]] 
                    if keyGerml in umiExp.lstGerml: continue # Is duplicate
                    umiExp.lstGerml.append(keyGerml)

                    btAddInfo += "VCF=" + fName
                    varType = Params.vtGerm
                else:
                    # Somatic variant testing
                    nSomas += 1
                    if nSomas > 1: continue # 1. prohibited  

                    isBreak, nEdgeDist = self.getBreakDist(plpCol.reference_pos, baseS)
                    if isBreak: continue # 2. Break founded
                    #if self.testBreakDist(plpCol.reference_pos): continue # 2. Break founded
                    if self.isStrangeVar(plpCol.reference_pos): continue  # Is Strange Variant (by Maslov)

                    btAddInfo += "BD={};".format(nEdgeDist) 
                    varType = Params.vtUndef  # light soma
                    isIdent = False
                    if Params.maxCoincidence > 0: # 3.
                        strListCol = self.getPileupColumn(plpCol) # G-column as string list
                        if strListCol.count(baseS) >= Params.maxCoincidence: 
                            isIdent = True
                            btAddInfo += "RR=#3" 
                    if not isIdent:
                        famOther = umiExp.testFamPosBase(baseS, plpCol.reference_pos, famBesides=self)
                        if famOther != None:
                            isIdent = True
                            btAddInfo += "RR=#4" 
                    baseP, baseN = None, None
                    if not isIdent: # 5. +
                        baseP = self.testReadsPosBase(plpCol.reference_pos, isPositive = True)
                        if baseP == None: # 5. Not 100% identical in the (+) strand
                            isIdent = True
                            btAddInfo += "RR=#5+" 
                    if not isIdent: # 5. -
                        baseN = self.testReadsPosBase(plpCol.reference_pos, isPositive = False)
                        if baseN == None: # 5. Not 100% identical in the (-) strand
                            isIdent = True
                            btAddInfo += "RR=#5-" 
                    if not isIdent:
                        if baseP != baseN: # 6. Two strand bases are not identical
                            isIdent = True
                            btAddInfo += "RR=#6" 
                    if not isIdent:
                        varType = Params.vtSoma  # hard soma
                
                varRecs.append([plpCol.reference_pos, varID, strRef[k], baseS, Params.vkSNV, btAddInfo, varType]) # varType must be the last  

        # Calc for isAgreedPN
        isAgreedPN = True
        if nSomas:
            if not self.consListP: self.calcConsensusListPN()
            for k in range(0, self.nEnd - self.nStart):
                if (self.consListP[k] != Params.chNoQual) and (self.consListN[k] != Params.chNoQual) and \
                    (self.consListP[k] != self.consListN[k]):
                    isAgreedPN = False
                    break  
        
        # Put info to files
        for rec in varRecs:
            if nSomas > 1 and rec[-1] in (Params.vtUndef, Params.vtSoma):
                continue  # > 1 - prohibited
            if rec[-1] == Params.vtSoma and not isAgreedPN:
                rec[-1] = Params.vtUndef
                rec[-2] += "RR=#7" # 7. Not agreed strands in family

            varAddInfo = "VC={};SAO={};FM={}".format(rec[4], rec[6], self.sUMI) # vcKind, vcType, UMI
            self.putVariant(umiExp.vcfResData, rec[0], rec[1], rec[2], rec[3], varAddInfo) # Variant to file
            if rec[5] != None:
                varAddInfo += rec[5]
                self.putBackTrack(umiExp, rec[0], rec[1], rec[2], rec[3], varAddInfo) # Backtrack to file
        del varRecs[:]
           
        return nQualBases
                
    def putVariant(self, vcfData, nRefStart, varID, sRef, sAlter, varAddInfo):
        nRefPos = nRefStart + 1 # 1-based
        if varID == None: varID = Params.chNoInfo 
        outString = "{}\t{}\t{}\t".format(Params.chrName, nRefPos, varID)
        outString += "{}\t{}\t{}\t{}\t".format(sRef, sAlter, Params.chNoInfo, Params.chNoInfo)
        outString += varAddInfo
        outString += '\n' 
        #outString = str(datetime.datetime.now()) + ": " + outString
        vcfData.write(outString)
        vcfData.flush()

    def putBackTrack(self, umiExp, nRefStart, varID, sRef, sAlter, btAddInfo):
        umiExp.btResData.write("\n===>\t")
        
        self.putVariant(umiExp.btResData, nRefStart, varID, sRef, sAlter, btAddInfo)

        self.backTrackPos(umiExp.btResData, nRefStart, True) # Track this family
        umiExp.btFamPosition(nRefStart, self) # Track all fullformed families, besides self

        umiExp.btResData.flush()

    # Variant ID from VCF DB records
    def getDBVariant(self, vcfRecs, nRefStart, sAlter = None, varKind = None):
        varID, sResAlts = None, None
        
        while len(vcfRecs) and (vcfRecs[0].start < nRefStart):
            dbRec = vcfRecs.pop(0) # Delete smaller
            del dbRec
             
        for dbRec in vcfRecs:
            #if dbRec.start < nRefStart: continue # Skip smaller
            if dbRec.start > nRefStart: break  # Higher - end of work

            if (sAlter == None or sAlter in dbRec.alts) and (varKind == None or dbRec.info['VC'] == varKind):
                varID = dbRec.id
                if sAlter == None: sResAlts = dbRec.alts
                else: sResAlts = sAlter 
                break

        return varID, sResAlts
    
    def printOut(self, bReads = False):
        myStr = "{}\t{}:{}-{}\t{}\t{}".format(self.sUMI, self.currChrom, self.nStart, \
                          self.nEnd, self.cntPositive, self.cntNegative)
        if bReads:
            for myRead in self.myReads:
                myRead.printOut()
                
        print myStr

    def clear(self):
        del self.myReads[:]
