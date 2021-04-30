from collections import Counter
from glob import glob
from Bio import SeqIO
from Bio import SearchIO
import os
import sys
import subprocess
import gzip
from Bio import AlignIO
import pathlib
from pathlib import Path
from modules.reads_handler import ReadsHandler

pwd = pathlib.Path.cwd()
user_os = sys.platform

READS_PATH = pwd / "reads"

TEMP_TARGET_FILE = "TEMP_TARGET.fa"
ALIGNED_TARGET_FILE = "TEMP_TARGET_aligned.fa"
prew_path = os.getcwd()
READS_HANDLER = ReadsHandler()

# Ensures only one ReadHandler object is created to handle all files
READS_HANDLER = ReadsHandler()

class ASTarget(object):
    """Ageseq target class for storaging target sequence infomation
    """

    def __init__(self, seq_dict):
        self.seq_dict = seq_dict
        self.seqid_list = list()
        self.tar_aln = ""
        self.target_features = dict()
        for keys in self.seq_dict:
            self.seqid_list.append(keys)

    def showAsFasta(self):
        """Print out all target sequences in a fasta style.
        """
        for x in self.seq_dict:
            print(">"+x+"\n"+self.seq_dict[x]+"\n")

    def toFastaFile(self, out_file_name):
        """Write all target sequences into the given out_file_name file
        for later running blat_win_x84.
        """
        ofile = open(out_file_name, 'w')
        for x in self.seq_dict:
            ofile.write(">"+x+"\n"+self.seq_dict[x]+"\n")
        ofile.close()

    def getTargetSeq(self, tid):
        "Return the target sequence based on the ID provided."
        return self.seq_dict[tid]

    def gettarlist(self):
        "Return IDs of stored target sequences"
        return self.seqid_list

    def gettardict(self):
        "Return the dictionary which stores all targets"
        return self.seq_dict

    def add_target_alignment(self, alignment):
        "Adding an alignment to the ASTarget class"
        self.tar_aln = alignment

    def add_target_features(self, feature):
        "Adding variants/features to the current ASTarget class"
        self.target_features = feature

    def align_target(self, aligner):
        # issue: cross platform, assume user add to path
        #issue: temp_target_file
        """Aligning sequences in the target file with the aligner provided.
           Write the aligned targets into a file and add alignment to the
           ASTarget object. Then return the alignment.
        """

        muscle_executable = "muscle"

        muscle_cline = f'{muscle_executable}' \
            f' -in {TEMP_TARGET_FILE}' \
            f' -out {ALIGNED_TARGET_FILE}'
        print("Running muscle as:\n"+str(muscle_cline))
        subprocess.call(muscle_cline, shell=True)

        target_alignment = AlignIO.read(ALIGNED_TARGET_FILE, "fasta")
        self.add_target_alignment(target_alignment)
        return target_alignment

    def target_polymorphism(self):
        """This will call the align_target() function to directly, and
           analyze the alignment return a dictionary with identified variants
        """
        aligned_target = self.align_target("muscle")
        aligned_length = aligned_target.get_alignment_length()
        target_feature = dict()
        for i in range(0, aligned_length):
            aln_col = aligned_target[:, i]
            aln_col_base = "".join(dict.fromkeys(aln_col))
            if len(aln_col_base) > 1:
                target_feature[i] = aln_col
        self.add_target_features(target_feature)
        return target_feature


class ASSample(object):
    """Ageseq sample class for storaging reads infomation
    """

    def __init__(self, bioiter):
        self.bioiter = bioiter
        self.SampleSize = len(self.bioiter)
        self.SampleList = list(self.bioiter.keys())

    def toFastaFile(self, out_file_name):
        ofile = open(out_file_name, 'w')
        count = 0
        for id in self.bioiter:
            rec = self.bioiter[id]
            ofile.write(">"+str(id)+"\n"+str(rec.seq)+"\n")
            count += 1
        print("Writing %i reads into file\n" % count)

    def getallreadsids(self):
        return self.SampleList

    def getReadSeq(self, qid):
        return self.bioiter[qid].seq


class ASCore(object):
    """ASCore class functions as a core agent during the analysis of the
       relationship bewteen one read and target as well as the main place to
       save analysis related data. It is essential a combination of HSP parsed
       from PSL alignment and sequences with addtional analysis functions.

       :param qid: QueryID or the Read ID, is required during the initialization.
       The rest attributes are added or updated with the class function accordingly.
    """

    def __init__(self, qid):
        """QueryID or the Read ID, is required during the initialization.
        The rest attributes are added or updated with the class function accordingly.
        """
        self.read_id = qid
        "QueryID is required for an instance of ASCore object"
        self.tdict = dict()
        self.besthit_id = ""
        self.besthit_score = int()
        self.hsp = ""
        self.hspstr = ""
        self.editPattern = ""
        self.qseq = ""
        self.bhseq = ""
        self.mismatch = ""
        self.mismatch_count = 0
        self.indel = list()
        self.indel_count = 0
        self.var_mask = ""

    def assignEditPat(self):
        """Does not return, but updates the editPattern of the instance based on
        the mismatch and indels information in the instance.
        """
        old_mismatch = self.getPSLMismatch()
        old_indel = self.indel
        new_indel = old_indel.copy()
        if self.mismatch_count == 0 and self.indel_count == 0:
            self.editPattern = "WT"
        elif self.indel_count == 0:
            new_mismatch = old_mismatch.copy()
            for emis in old_mismatch:
                if emis in self.var_mask:
                    new_mismatch.remove(emis)
                    # logfile.write(str(old_mismatch[emis])+"\n")
            if len(new_mismatch) == 0:
                self.editPattern = "WT"
            else:
                self.editPattern = new_mismatch
        else:
            new_indel = old_indel.copy()
            for eind in old_indel:
                if eind in self.var_mask:
                    new_indel.remove(eind)
            if self.mismatch_count == 0:
                old_mismatch = []
            new_mismatch = old_mismatch.copy()
            for emis in old_mismatch:
                if emis in self.var_mask:
                    new_mismatch.remove(emis)
            if len(new_mismatch) == 0 and len(new_indel) == 0:
                self.editPattern = "WT"
            elif len(new_mismatch) == 0:
                self.editPattern = new_indel
            else:
                self.editPattern = new_indel + new_mismatch

    def getEditPat(self):
        """Returns the editPattern"""
        return self.editPattern

    def updateHSP(self, hsp):
        """Adding/Update HSP from BLAT to the instance. To use hsp, you can
        follow the document on BlatIO in Biopython
        :param hsp: The hsp object from BlatIO
        """
        self.hsp = hsp

    def updateHSPstr(self, hspstr):
        """Adding/Update parsed HSP strings. This is a more convient way to save
        parsed HSP results. In the psl_parse function, this is used for saving
        the parsed HSP in a concatenated string format.

        :param str hspstr: a hsp string separated by space.
        This string starts with number of gapopenings by
        number of "*", no "*" if no gap openings. Then, followed by:
        all query start positions
        all query span lengths
        all hit start positions
        all hit span lengths
        number of gap opening in the query
        number of gap opening in the hits
        strand of hit
        number of mismatches
        """
        self.hspstr = hspstr

    def updateBH(self, newBH):
        """Adds/Updates BestHit

        :param newBH: the ID of best hit
        """
        self.besthit_id = newBH

    def addINDEL(self, newINDEL):
        """Adds snewINDEL into the indel list through append()
        and update the indel_count attribute.

        :param newINDEL: String of coded INDEL
        e.g. 16D2:GT this encodes 2 bp deletion happened on the 16 target
        position (0-pos), the corresponding nucleotides are GT in the target.
        This string is automatically generated in the main function during
        psl_parse().
        """

        self.indel.append(newINDEL)
        self.indel_count = len(self.indel)

    def updateBHScore(self, newBHscore):
        """Adding/Update BestHit score

        :param newBHscre: integral score
        """
        self.besthit_score = newBHscore

    def updateEP(self, newEP):
        """Adding/Update editPattern

        :param newEP: string of editPattern
        """
        self.editPattern = newEP

    def updateQSEQ(self, qseq):
        """Adding/Update query sequence

        :param qseq: string of query sequence
        """
        self.qseq = qseq

    def updateBHSEQ(self, bhseq):
        """Adding/Update best hit sequence

        :param bhseq: string of best hit sequence
        """
        self.bhseq = bhseq

    def updateMismatch(self, new_mismatch):
        """Adding/Update mismatches between a read and its best hit and
        updating number of mismatches automatically.

        :param new_mismatch: dict e.g. {36: "G->T"} target 36 (0-pos) bp G has
        a mismatch T in query
        """
        self.mismatch = new_mismatch
        self.mismatch_count = len(new_mismatch)

    def updateMismatchCount(self, mcount):
        """Updates mimatch count

        :param mcount: int, mismatch count
        """
        self.mismatch_count = mcount

    def updatePSLTargets(self, psl_target_dict):
        """Updates/Adds PSL object from BlatIO

        :param psl_target_dict: dict, dictionary generatd by SearchIO to_dict()
        """
        self.tdict = psl_target_dict

    def getPSLMismatch(self):
        """Returns mismatches from the alignment. A list will be returned if
        there are multiple mismatches
        """
        if self.mismatch_count > 0:
            return list(self.mismatch.values())
        else:
            return self.mismatch

    def getPSLTargets(self):
        """Returns the dictionary with the targetID as keys, and corresponding
        alignment scores as values.
        """
        return str(self.tdict)

    def getreadID(self):
        """Returns the reads ID as strings
        """
        return str(self.read_id)

    def getBH(self):
        """Returns the best Hit ID as strings
        """
        return str(self.besthit_id)

    def getHitScore(self, targetID):
        """Returns alignment score of given target, integral.
        :param targetID: str, the ID of target
        """
        try:
            return self.tdict[targetID]
        except KeyError:
            print(str(targetID)+" does not exists in targets!\n")

    def getBHScore(self):
        """Returns alignment score of best hit, integral.
        """
        return self.getHitScore(self.besthit_id)

    def printmismatch(self):
        """Return number of mismatches plus the mimatch dictionary as one line
        string (used for printing log information).
        """
        ostr = str(self.mismatch_count)+" Mismatch(es)\t" + \
            str(self.getPSLMismatch())+"\n"
        return ostr

    def getINDEL(self):
        """Returns the list of indels
        """
        if len(self.indel) == 0:
            print("No INDEL")
            return []
        else:
            return self.indel

    def getMismatchCount(self):
        """Returns the number of mismatches
        """
        return self.mismatch_count

    def updateMask(self, mask_dict):
        """Updates/Adds the dictionary of masked variants.

        :param mask_dict: dict, the dictionary used for tracking variants should
        be masked (Variants like wobble base, low frequency errors).
        """
        self.var_mask = mask_dict

    def getHSPrange(self):
        """Returns the start and end (full span) in a list of the target aligned
        by read.
        """
        start = 0
        end = 0
        for i in self.hsp:
            start = i.hit_start_all[0]
            end = i.hit_start_all[-1]+i.hit_span_all[-1]
            break
        return [start, end]


def readv1TargetFile(v1_target):
    """Read the specified target file and returns a target object

    :param v1_target: an input file formatted as AGEseq v1.
    """
    cfile = open(v1_target, "r")
    tseq_dict = dict()
    print("Reading Target File!")
    seq_count = 0
    for cfileline in cfile:
        if cfileline.startswith("#") or not cfileline.strip():
            pass
        else:
            tfline = cfileline.split()
            if tfline[0] in tseq_dict.keys():
                print("Duplicated seqID are found in the target, please modify it!")
                exit
            else:
                tseq_dict[tfline[0]] = tfline[1]
                seq_count += 1
    print(str(seq_count)+" sequences are loaded from target file!")
    USER_TARGETFILE = ASTarget(tseq_dict)
    return USER_TARGETFILE


def readsToSample(reads_file_id):
    """Read the specified read file and returns a sample object

    :param reads_file_id: str, a read file name
    """
    return ASSample(READS_HANDLER.handle(reads_file_id))

    # rf = reads_file_id
    # print(reads_file_id)
    # if reads_file_id.endswith("gz"):
    #     file_handler = gzip.gzopen(rf, "rt")
    # else:
    #     file_handler = open(rf, "r")
    #
    # # issue use try to catch both fasta and fastq
    # # to do if else try
    # try:
    #     seq_parser = SeqIO.to_dict(SeqIO.parse(file_handler, "fastq"))
    # except:
    #     seq_parser = SeqIO.to_dict(SeqIO.parse(file_handler, "fasta"))
    # print("fasta/fastq file is opened!")
    # return ASSample(seq_parser)


def parseBLATlegacymod(psl_file):
    """Wrapper of psl to psl_dictionary.Returns a dictionary generated by
    SearchIO.to_dict()

    :param psl_file: str, psl_file name
    """
    psl_parse = SearchIO.parse(psl_file, 'blat-psl')
    psl_dict = SearchIO.to_dict(psl_parse)
    return psl_dict
    # Only Matched bases > 0.2*query(reads)length will be processed


# def getquerybesthit(qid, qblatres):
#     """Returns the best hitID and the score of besthit from
#     """
#     max_score = 0
#     bh_tid = ""
#     for hit in qblatres:
#         if hit.score > max_score:
#             max_score = hit.hsp.score
#             bh_tid = hit.id
#         else:
#             next
#     return bh_tid, max_score


if __name__ == '__main__':
    pass
