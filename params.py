'''
Created on Apr 16, 2021

@author: S. Makhortov, A. Maslov
'''
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

class Params(object):
    '''
    classdocs
    '''
    # commands
    command = "???"
    cmdUMI = "umi"
    cmdRMDUP = "rmdup"
    cmdVC = "vc"
    cmdTR = "trunc"

    # File names
    bamSFile = ""
    bamGFile = ""
    refFile = ""
    vcfAllFile = ""
    vcfPersFile = ""
    
    chrName = None
    
    NAME_DELIM = ':'
    UMI_DELIM = '+'

    # Parameters
    umiLen = 9
    maxUMIDist = 2
    maxPosDist = 5
    minUmiComplexity = 1

    minBaseQuality = 13
    minMapQuality = 60

    minStrandCount = 7
    minFamilyCount = 0
    
    maxPosFamilies = 100 # For rmdup command
    
    logFamilies = 0
    testOverlapPair = 1 # Test for overlapped pair in the one family

    maxCoincidence = 1
    minBreakDistance = 5
    
    # Special chars
    chNoQual = 'Q'
    chDel = 'D'
    chNoBase = 'N'
    chNoInfo = '.'
    
    # Cigar sections
    SectM = 0 
    SectI = 1
    SectD = 2
    SectS = 4
    SectH = 5
    
    # Variant kinds (VC)
    vkSNV = "SNV"
    vkDel = "DEL" 
    vkIns = "INS"  

    # Variant types (SOA)
    vtUndef = 0
    vtGerm = 1
    vtSoma = 2
    
    #minDomine = 0.3
    #minCount = 3

    def __init__(self, cfgFile):
        '''
        Constructor
        '''
    
    @classmethod
    def init(cls, cfgFile):
        cfg = configparser.ConfigParser()
        if cfg.read(cfgFile):
            Params.bamSFile = cls.getStrParam(cfg, "File", "bamSFile", Params.bamSFile)
            Params.bamGFile = cls.getStrParam(cfg, "File", "bamGFile", Params.bamGFile)
            Params.refFile = cls.getStrParam(cfg, "File", "refFile", Params.refFile)
            Params.vcfAllFile = cls.getStrParam(cfg, "File", "vcfAllFile", Params.vcfAllFile)
            Params.vcfPersFile = cls.getStrParam(cfg, "File", "vcfPersFile", Params.vcfPersFile)

            Params.chrName = cls.getStrParam(cfg, "File", "chrName", Params.chrName)
            if (Params.chrName): Params.chrName = Params.chrName.strip()
            
            Params.NAME_DELIM = cls.getStrParam(cfg, "File", "NAME_DELIM", Params.NAME_DELIM).strip()

            Params.minBaseQuality = cls.getIntParam(cfg, "File", "minBaseQuality", Params.minBaseQuality)
            Params.minMapQuality = cls.getIntParam(cfg, "File", "minMapQuality", Params.minMapQuality)

            Params.UMI_DELIM = cls.getStrParam(cfg, "UMI", "UMI_DELIM", Params.UMI_DELIM).strip() #':'
            Params.umiLen = cls.getIntParam(cfg, "UMI", "umiLen", Params.umiLen)
            Params.maxUMIDist = cls.getIntParam(cfg, "UMI", "maxUMIDist", Params.maxUMIDist)
            Params.minUmiComplexity = cls.getIntParam(cfg, "UMI", "minUmiComplexity", Params.minUmiComplexity)

            Params.maxPosDist = cls.getIntParam(cfg, "Family", "maxPosDist", Params.maxPosDist)
            
            Params.minStrandCount = cls.getIntParam(cfg, "Family", "minStrandCount", Params.minStrandCount)
            Params.minFamilyCount = cls.getIntParam(cfg, "Family", "minFamilyCount", Params.minFamilyCount)
            Params.testOverlapPair = cls.getIntParam(cfg, "Family", "testOverlapPair", Params.testOverlapPair)
            
            if (Params.minFamilyCount > 0) and (Params.minFamilyCount < Params.minStrandCount + Params.minStrandCount):
                Params.minFamilyCount = Params.minStrandCount + Params.minStrandCount
                print "minFamilyCount is automatically adjusted to 2*minStrandCount"

            Params.maxPosFamilies = cls.getIntParam(cfg, "Family", "maxPosFamilies", Params.maxPosFamilies)
            
            Params.logFamilies = cls.getIntParam(cfg, "Family", "logFamilies", Params.logFamilies)

            Params.maxCoincidence = cls.getIntParam(cfg, "Variant", "maxCoincidence", Params.maxCoincidence)
            Params.minBreakDistance = cls.getIntParam(cfg, "Variant", "minBreakDistance", Params.minBreakDistance)
        
        else: print "File not open: " + cfgFile
             
    @staticmethod
    def getStrParam(cfg, sect, opt, defValue = None):
        try: retValue = cfg.get(sect, opt)
        except: retValue = defValue
        return retValue
            
    @staticmethod
    def getIntParam(cfg, sect, opt, defValue = None):
        try: retValue = cfg.getint(sect, opt)
        except: retValue = defValue
        return retValue
        