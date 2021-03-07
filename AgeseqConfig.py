# -*- coding: utf-8 -*-


class ASConfig(object):
    """
    This is a AGESeqConfig class for saving AGEseq Configuration
    """

    def __init__(self, mismatch_cutoff, min_cutoff,
                 wt_like_report, indel_report, remove_files):
        self.mismatch_cutoff = mismatch_cutoff
        self.min_cutoff = min_cutoff
        self.wt_like_report = wt_like_report
        self.indel_report = indel_report
        self.remove_files = remove_files

    def returnConfig(self):
        print("Mismatch cutoff:", self.mismatch_cutoff)
        print("Minimal cutoff:", self.min_cutoff)
        print("Top WT reported:", self.wt_like_report)
        print("Top Indel reported:", self.indel_report)
        print("Remove intermediate/BLAT files:", self.remove_files)


class BLATConfig(object):
    """
    This is a BLATConfig class for saving BLAT Configuration
    """

    def __init__(self, tileSize, oneOff,
                 maxGap, minIdentity, minScore):
        self.tileSize = tileSize
        self.oneOff = oneOff
        self.maxGap = maxGap
        self.minIdentity = minIdentity
        self.minScore = minScore

    def returnBLATConfig(self):
        """Report BLATConfiguration"""
        print("BLAT is set as following:")
        print("tileSize:", self.tileSize)
        print("oneOff:", self.oneOff)
        print("maxGap:", self.maxGap)
        print("minIdentity:", self.minIdentity)
        print("minScore:", self.minScore)


# Preset:
ASPRESET = ASConfig(mismatch_cutoff=0.1,
                    min_cutoff=0,
                    wt_like_report=20,
                    indel_report=50,
                    remove_files=True)

BLATPRESET = BLATConfig(tileSize=7,
                        oneOff=1,
                        maxGap=20,
                        minIdentity=70,
                        minScore=20)

PRESET_PAM = {
    "SpCas9": "NGG",
    "SaCas9": "NGRRN",
    "NmeCas9": "NNNNGATT",
    "CjCas9": "NNNNRYAC",
    "StCas9": "NNAGAAW",
    "LbCpf1": "TTTV",
    "AsCpf1": "TTTV"
}


def readConfigFile(configFile):
    """
    Reading configuration from the configFile
    Return an instance of ASConfig and a BLATConfig instance
    via a tuple.
    """
    cfile = open(configFile, "r")
    paradict = dict()
    print("Reading Configuration as following!")
    for cfileline in cfile:
        if cfileline.startswith("#") or not cfileline.strip():
            pass
        else:
            spline = cfileline.split()
            paradict[spline[0]] = spline[2]
            print(spline[0]+"\t"+spline[2])

    USER_ASCONF = ASConfig(mismatch_cutoff=paradict['mismatch_cutoff'],
                           min_cutoff=paradict['min_cutoff'],
                           wt_like_report=paradict['wt_like_report'],
                           indel_report=paradict['indel_report'],
                           remove_files=paradict['remove_files'])
    USER_BLCONF = BLATConfig(tileSize=paradict['blat_tileSize'],
                             oneOff=paradict['blat_oneOff'],
                             maxGap=paradict['blat_maxGap'],
                             minIdentity=paradict['blat_minIdentity'],
                             minScore=paradict['blat_minScore'])
    return (USER_ASCONF, USER_BLCONF)


if __name__ == '__main__':
    print("Default presets for AGEseq are:")
    ASPRESET.returnConfig()
    BLATPRESET.returnBLATConfig()
