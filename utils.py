'''
Created on Apr 16, 2021

@author: S. Makhortov, A. Maslov
'''

from params import Params
import datetime

class Utils(object):
    '''
    classdocs
    '''
    logFile = None
    
    def __init__(self, params):
        '''
        Constructor
        '''

    @staticmethod
    def logOut(myStr, putTime = False):
        logStr = myStr
        if putTime:
            logStr = str(datetime.datetime.now()) + ": " + logStr
        Utils.logFile.write(logStr)
        Utils.logFile.flush()

    @staticmethod
    def posDistance(nStart1, nEnd1, nStart2, nEnd2):
        return max(abs(nStart1 - nStart2), abs(nEnd1 - nEnd2))        
    
    @staticmethod
    def segOverlap(nStart1, nEnd1, nStart2, nEnd2, minOver = 0):
        isOverlap = nStart1 in range(nStart2, nEnd2 - minOver)
        isOverlap = isOverlap or (nStart2 in range(nStart1, nEnd1 - minOver))
        return isOverlap        

    @staticmethod
    def umiSwap(umi):
        uL, uR = umi.split(Params.UMI_DELIM)
        res = uR + Params.UMI_DELIM + uL
        return res        

    @staticmethod
    def umiDistance(u1, u2):
        # print '0: ', u1, u2
        u1a, u1b = u1.split(Params.UMI_DELIM)
        u2a, u2b = u2.split(Params.UMI_DELIM)
        d1, d2 = 0, 0
        for i in range(0, Params.umiLen):
            if u1a[i] != u2a[i]: d1 += 1
            if u1b[i] != u2b[i]: d1 += 1
            if u1a[i] != u2b[i]: d2 += 1
            if u1b[i] != u2a[i]: d2 += 1    
        return min(d1, d2)        

    @staticmethod
    def pureName(readName): # Name without UMI
        res = readName
        nLastDelim = res.rindex(Params.NAME_DELIM)
        if nLastDelim > 0:
            res = res[:nLastDelim] 
        return res

    @staticmethod
    def mostCommon(strList, nMinCount = 0, fMinPart = 0, sBesides = None):
        nStrs = len(strList)
        if not nStrs: return -1
         
        counts = []
        for i in range(0, nStrs):
            counts.append(0)
            if strList[i] == sBesides: continue # Not counted by condition
            if strList[i] == None: continue   # Empty - not counted
            counts[i] += 1
            for j in range(i+1, nStrs):
                if strList[j] == strList[i]:
                    strList[j] = None
                    counts[i] += 1
        iMax = 0; # Search for max count
        for i in range(1, nStrs):
            if counts[i] > counts[iMax]: iMax = i
        
        # Requirements test
        if counts[iMax] < nMinCount:
            iMax = -1
        elif fMinPart and ((float(counts[iMax]) / nStrs) < fMinPart):
            iMax = -1
        
        del counts[:]
        return iMax # The index of
    
    '''
    @staticmethod
    def umiDistanceLight(u1, u2):
        d = 0
        for i in range(0, Params.umiLen + Params.umiLen + 1):
            if u1[i] != u2[i]: d += 1
        return d        
    '''
    '''      
    def isUMIinList(self, u0, umiList): #Compare one UMI with every UMI in the list
        res = True    
        for umi in umiList:
            if umiDistance(umi, u0) > Params.maxUMIDist:    #if distance with any member is more than MinUmiDistance, then return False
                res = False
                break

        return res
    '''

        