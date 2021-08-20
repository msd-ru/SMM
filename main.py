'''
Created on Mar 29, 2021

@author: S. Makhortov, A. Maslov
'''

import sys, os

from umiexplorer import UMIExplorer
from params import Params

def main():
    print sys.argv
    
    myIniFile = os.path.join(os.curdir, "params.ini")
    Params.init(myIniFile)
    
    umiExplorer = UMIExplorer()
    nParams = 0
    
    nParams += 1
    if len(sys.argv) > nParams:
        Params.command = sys.argv[nParams] # Command

    nParams += 1
    if len(sys.argv) > nParams: 
        Params.bamSFile = sys.argv[nParams] # BAM-file
        
    nParams += 1
    if len(sys.argv) > nParams: 
        Params.chrName = sys.argv[nParams]  # Chromosome
        if Params.chrName == "-": # No Chr 
            Params.chrName = None

    '''
    print Params.command 
    print Params.bamSFile 
    print ' Chr =', Params.chrName
    print os.getcwd()
    '''
    if (Params.command != Params.cmdRMDUP):
        
        if (Params.command == Params.cmdTR): # BAM truncation
            nPosMin = None
            nParams += 1
            if len(sys.argv) > nParams: 
                nPosMin = int(sys.argv[nParams])
                
            nPosMax = None
            nParams += 1
            if len(sys.argv) > nParams: 
                nPosMax = int(sys.argv[nParams])

            umiExplorer.bamTrunc(nPosMin, nPosMax)
            return
 
        if Params.command == Params.cmdVC:
            nParams += 1
            if len(sys.argv) > nParams: 
                Params.bamGFile = sys.argv[nParams] # G-file
            nParams += 1
            if len(sys.argv) > nParams: 
                Params.refFile = sys.argv[nParams]  # Reference-file
            nParams += 1
            if len(sys.argv) > nParams: 
                Params.vcfAllFile = sys.argv[nParams]  # VCF All input file
            nParams += 1
            if len(sys.argv) > nParams: 
                Params.vcfPersFile = sys.argv[nParams]  # VCF Person input file
            
    umiExplorer.umiAnalyze()

if __name__ == '__main__':
    main()