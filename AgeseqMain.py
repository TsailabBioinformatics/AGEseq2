# -*- coding: utf-8 -*-
"""This is the main entry of AGEseq.
By running this script, AGEseq will automatically run with default setting.
"""

from __future__ import division
from datetime import datetime
import os
import sys
import AgeseqConfig as ageconf
import AgeseqIO
import subprocess
from collections import defaultdict
from Bio import Seq
import pathlib
from pathlib import Path
pwd = pathlib.Path.cwd()
# pwd = os.path.dirname(__file__)
user_os = sys.platform
prew_path = os.getcwd()
os.chdir(pwd)
READS_PATH = pwd / "reads"
TEMP_TARGET_FILE = "TEMP_TARGET.fa"
READS_FILE_LIST = list()
DEF_BLAT_PATH = pwd / "blat"
#record for script running
sysdatetime = datetime.now()
dt_string = sysdatetime.strftime("%Y%m%d_%H_%M_%S")
logfile_name = "AGESeq_run_"+dt_string+".log"
logfile = open(logfile_name, "w")
WELCOM_MESSAGE = "Thank you for using AGEseq!\n"
mod = "legacy"
logfile.write(WELCOM_MESSAGE)

READ_SNP_MINR_FREQ = 0.05  #Target specific, not overall frequency
READ_INDEL_MINR_FREQ = 0.01 #Target specific, not overall frequency
READ_SNP_MINIMAL_SUPPORT = 3
READ_INDEL_MINIMAL_SUPPORT = 5
WOBBLE_BASE = True
WOBBLE_FREQ_RANGE = (0.35, 0.65)

def main():
    # checking path/dir/
    pconf = check_conf()
    USER_ASCONF = pconf[0]
    USER_BLCONF = pconf[1]
    USER_BLCONF.returnBLATConfig()
    DEF_BLAT_PATH = check_blat()
    USER_TARGET = load_target()
    target_feature = USER_TARGET.target_polymorphism()
    logfile.write("Polymorphic sites in targets:"+"\n")
    logfile.write(str(target_feature)+"\n")
    USER_READSFILE_LIST = get_reads_file()
    # for loop for what #issue
    for urfile in USER_READSFILE_LIST:
        read_file = AgeseqIO.readsToSample(urfile)
        read2fas_file = str(urfile)+".fa"
        read_file.toFastaFile(read2fas_file)
        # Run blat
        blat_outfmt="psl"
        os_blat = "blat"
        file_blat_in = read2fas_file
        blat_out = f'{file_blat_in}_blat_crispr.psl'
        blat_cmd = f'{os_blat} {TEMP_TARGET_FILE} {file_blat_in} ' \
                   f' -tileSize={USER_BLCONF.tileSize}' \
                   f' -oneOff={USER_BLCONF.oneOff}' \
                   f' -maxGap={USER_BLCONF.maxGap}' \
                   f' -minIdentity={USER_BLCONF.minIdentity}' \
                   f' -out={blat_outfmt}' \
                   f' -minScore={USER_BLCONF.minScore} {blat_out}'
        logfile.write("Running blat as:\n"+str(blat_cmd)+"\n")
        subprocess.call(blat_cmd, shell=True)
        naive_assign_rst = psl_parse(read_file, USER_TARGET, blat_out)
        assigned_core = naive_assign_rst[0]
        mismatch_collection = naive_assign_rst[1]
        indel_collection = naive_assign_rst[2]
        target_assigned_count = naive_assign_rst[3]
        target_var_freq = cal_target_var_freq(mismatch_collection,
                                              indel_collection,
                                              target_assigned_count)
        logfile.write("Assignment begins!"+"\n")
        file_assign_sum = assign_mask(assigned_core, mismatch_collection,
                    indel_collection, target_assigned_count, WOBBLE_BASE)
        for t in file_assign_sum:
            for ep in file_assign_sum[t]:
                logfile.write(t+"\t"+ep+"\t"+str(file_assign_sum[t][ep])+"\n")
        for eachtarget in sorted(target_assigned_count):
            target_assigned_prop = round(target_assigned_count[eachtarget]/len(assigned_core)*100, 1)
            logfile.write("\n"+eachtarget+"\t"+str(target_assigned_count[eachtarget])+"\t"+
                          str(target_assigned_prop)+"% of assigned reads"+"\n")

def cal_target_var_freq(mismatch_collection, indel_collection,
                        target_assigned_count):
    target_var_freq = defaultdict(dict)
    for target in target_assigned_count:
        for mismatch in mismatch_collection[target]:
            totalmismatch = mismatch_collection[target][mismatch]
            target_var_freq[target][mismatch] = totalmismatch
        for mismatch in indel_collection[target]:
            totalmismatch = indel_collection[target][mismatch]
            target_var_freq[target][mismatch] = totalmismatch
    return target_var_freq


def assign_mask(assigned_core, mismatch_collection,
                indel_collection, target_assigned_count, WOBBLE_BASE):
    file_summary = defaultdict(dict)
    for qid in assigned_core:
        as_core = assigned_core[qid]
        coreID = as_core.getreadID()
        coreBH = as_core.getBH()
        coreBHscore = as_core.getBHScore()
        coremismatch = as_core.getPSLMismatch()
        coremismatch_count = as_core.getMismatchCount()
        coremismatchstr = as_core.printmismatch()
        coreINDEL = as_core.getINDEL()
        logfile.write(str(coreID)+"\t"+coreBH+"\t"+str(coreBHscore)+"\n")
        logfile.write(coremismatchstr)
        var_mask = dict()
        if coremismatch_count > 0:
            for mismatch in coremismatch:
                str_mismatch = mismatch
                totalmismatch = mismatch_collection[coreBH][str_mismatch]
                misfreq = round(totalmismatch/target_assigned_count[coreBH], 3)
                if misfreq < READ_SNP_MINR_FREQ or totalmismatch < READ_SNP_MINIMAL_SUPPORT:
                    var_mask[mismatch] = "masked"
                    logfile.write(str_mismatch+"\t"+str(misfreq)+"\t"+
                                  " removed!\n")
                else:
                    if WOBBLE_BASE:
                        target_range = as_core.getHSPrange()
                        mispos = int(mismatch.split(":")[0])
                        if (mispos -target_range[0] < 12 or abs(mispos -target_range[1]) <12) and (misfreq < WOBBLE_FREQ_RANGE[1] and misfreq > WOBBLE_FREQ_RANGE[0]):
                            var_mask[mismatch] = "masked"
                            logfile.write(str_mismatch+"\t"+str(misfreq)+" masked as a wobblebase!\n")
                        else:
                            logfile.write(str_mismatch+"\t"+str(misfreq)+" kept!\n")
                    else:
                        logfile.write(str_mismatch+"\t"+str(misfreq)+" kept!\n")
        if len(coreINDEL) > 0:
            logfile.write(str(len(coreINDEL))+" INDEL(s) is/are found!\n")
            for indel in coreINDEL:
                totalindel = indel_collection[coreBH][indel]
                indelfreq = round(totalindel/target_assigned_count[coreBH], 3)
                if indelfreq < READ_INDEL_MINR_FREQ or totalindel < READ_INDEL_MINIMAL_SUPPORT:
                    var_mask[indel] = "masked"
                    logfile.write(indel+"\t"+str(indelfreq)+"\t"+
                                  " removed!\n")
                else:
                    logfile.write(indel+"\t"+str(indelfreq)+" kept!\n\n")
        else:
            logfile.write("No INDEL(s)!"+"\n\n")
        as_core.updateMask(var_mask)
        as_core.assignEditPat()
        coreEP = as_core.getEditPat()
        if type(coreEP) is list:
            coreEP = coreEP[0]
        if coreBH in file_summary:
            if coreEP in file_summary[coreBH]:
                file_summary[coreBH][coreEP] += 1
            else:
                file_summary[coreBH][coreEP] = 1
        else:
            file_summary[coreBH][coreEP] = 1
        logfile.write(str(coreEP)+" assigned!\n\n")
    return file_summary

def masked_edit_pattern(assigned_core):
    assigned_core

def psl_parse(read_file, target_file, blat_out):
    sample_readslist = read_file.getallreadsids()
    number_reads_in_sample = len(sample_readslist)
    ASCore_collection = dict()
    psl_dict = AgeseqIO.parseBLATlegacymod(blat_out)
    q_besthit_score = dict()
    q_besthit = dict()
    tier_count = 0
    number_q_in_psl = len(psl_dict)
    perc_hit = number_q_in_psl/number_reads_in_sample*100
    #to2file = open("test_pa.out", "w")
    logfile.write("%i reads have hits to the target sequences provided (%d%% of the sample)!" % (number_q_in_psl, perc_hit)+"\n")
    for qid in psl_dict:
        ASCore_collection[qid] = AgeseqIO.ASCore(qid)
        psl_t_dict = dict()
        for hit in psl_dict[qid]:
            for hsp in hit:
                psl_t_dict[hit.id] = hsp.score
                if qid in q_besthit.keys():
                    if hsp.score > q_besthit_score[qid]:
                        q_besthit_score[qid] = hsp.score
                        q_besthit[qid] = hit.id
                    elif hsp.score == q_besthit_score[qid]:
                        logfile.write(qid+" Tie Hit found!"+str(hit.id)+"\n")
                        logfile.write(qid+" Tie Hit found!"+q_besthit[qid]+"\n")
                        tier_count += 1
                    else:
                        continue
                else:
                    q_besthit_score[qid] = hsp.score
                    q_besthit[qid] = hit.id
        ASCore_collection[qid].updatePSLTargets(psl_t_dict)
    #for eachread in sample_readslist:
    print("Found %i reads have tie hit" % tier_count)
    PATTERN_COLLECTION = defaultdict(dict)
    MISMATCH_COLLECTION = defaultdict(dict)
    target_assigned_count = dict()

    for query in q_besthit:
        ASCore_collection[query].updateBH(q_besthit[query])
        q_hsp = psl_dict[query][q_besthit[query]]
        ASCore_collection[query].updateHSP(q_hsp)
        q_seq = read_file.getReadSeq(query)
        t_seq = target_file.getTargetSeq(q_besthit[query])
        if q_besthit[query] in target_assigned_count:
            target_assigned_count[q_besthit[query]] += 1
        else:
            target_assigned_count[q_besthit[query]] = 1
        #logfile.write(str(query)+"\t"+q_besthit[query]+
        #              "\t"+str(q_besthit_score[query])+"\n")
        #to2file.write(str(ASCore_collection[query].getPSLTargets())+"\n")
        for ehsp in q_hsp:
            q_start_all = ehsp.query_start_all
            q_span_all = ehsp.query_span_all
            h_start_all = ehsp.hit_start_all
            h_span_all = ehsp.hit_span_all
            q_gap = ehsp.hit_gapopen_num
            h_gap = ehsp.query_gapopen_num
            mismatch = ehsp.mismatch_num
            hit_frag_seq = []
            query_frag_seq = []
            a = len(ehsp)
            strand_info = []
            for i in range(0, a):
                if ehsp[i].query_strand == -1:
                    strand_info.append("-")
                else:
                    strand_info.append("+")
            for i in range(0, len(h_start_all)):
                s = t_seq[h_start_all[i]:h_start_all[i]+h_span_all[i]]
                hit_frag_seq.append(s)
                fq_seq = q_seq[q_start_all[i]:q_start_all[i] + q_span_all[i]]
                query_frag_seq.append(fq_seq)
            if strand_info[0] == "-":
                ucount = 0
                for seqs in query_frag_seq:
                    query_frag_seq[ucount]=Seq.reverse_complement(query_frag_seq[ucount])
                    ucount += 1
        psl_str = str("*"*ehsp.gapopen_num+str(q_start_all)+" "+str(q_span_all)
                      + " "+str(h_start_all)+" "+str(h_span_all)+" "+str(q_gap)
                      + " "+str(h_gap)+" "+str(strand_info)+" "
                      + str(mismatch)+"\n")
        ASCore_collection[qid].updateHSPstr(psl_str)
        MISMATCH = dict()
        mismatch_count = 0
        for i in range(0, len(hit_frag_seq)):
            blast_link = []
            for j in range(0, len(hit_frag_seq[i])):
                if hit_frag_seq[i][j] == query_frag_seq[i][j]:
                    blast_link.append("*")
                else:
                    blast_link.append(".")
                    mismatch_pos = h_start_all[i] + j + 1
                    MISMATCH[mismatch_count] = str(mismatch_pos) + ":" + str(hit_frag_seq[i][j]) + "->"+str(query_frag_seq[i][j])
                    if MISMATCH[mismatch_count] in MISMATCH_COLLECTION[q_besthit[query]]:
                        MISMATCH_COLLECTION[q_besthit[query]][MISMATCH[mismatch_count]] += 1
                    else:
                        MISMATCH_COLLECTION[q_besthit[query]][MISMATCH[mismatch_count]] = 1
                    mismatch_count += 1
            #logfile.write(str(hit_frag_seq[i])+"\n")
            #logfile.write(str(query_frag_seq[i]) +"\n")
            #logfile.write(str("".join(blast_link)) + "\n\n")
        #to2file.write(str(MISMATCH) + "\n\n")

        if mismatch == 0:
            pass
        else:
            ASCore_collection[query].updateMismatch(MISMATCH)
            ASCore_collection[query].updateMismatchCount(mismatch)
        #logfile.write(str(ASCore_collection[query].mismatch_count)+" Mismatch(es)\t"+str(ASCore_collection[query].getPSLMismatch())+"\n")
        for i in range(0, q_gap):
            if h_start_all[i] + h_span_all[i] < h_start_all[i+1]:
                del_start = h_start_all[i] + h_span_all[i]
                del_end = h_start_all[i+1] - 1
                del_size = del_end - del_start + 1
                del_seq = t_seq[del_start: del_start+del_size]
                del_pat = str(del_start) + "D" + str(del_size) + ":" + str(del_seq)
                ASCore_collection[query].addINDEL(del_pat)
                if del_pat in PATTERN_COLLECTION[q_besthit[query]]:
                    PATTERN_COLLECTION[q_besthit[query]][del_pat] += 1
                else:
                    PATTERN_COLLECTION[q_besthit[query]][del_pat] = 1
                #logfile.write(del_pat + "\n\n")
        for i in range(0, h_gap):
            if h_start_all[i] + h_span_all[i] == h_start_all[i+1]:
                ins_start_ontarget = h_start_all[i+1]
                ins_start_onquery = q_start_all[i] + q_span_all[i]
                ins_end_onquery = q_start_all[i+1] - 1
                if strand_info[0] == "-":
                    ins_start_onquery = q_start_all[i+1] + q_span_all[i+1]
                    ins_end_onquery = q_start_all[i] - 1
                ins_size = ins_end_onquery - ins_start_onquery +1
                ins_seq = q_seq[ins_start_onquery: ins_start_onquery+ins_size]
            ins_pat = str(ins_start_ontarget) + "I" + str(ins_size) + ":" + str(ins_seq)
            ASCore_collection[query].addINDEL(ins_pat)
            if ins_pat in PATTERN_COLLECTION[q_besthit[query]]:
                PATTERN_COLLECTION[q_besthit[query]][ins_pat] += 1
            else:
                PATTERN_COLLECTION[q_besthit[query]][ins_pat] = 1
            #logfile.write(ins_pat + "\n\n")
    return ASCore_collection, MISMATCH_COLLECTION, PATTERN_COLLECTION, target_assigned_count
'''
    for eachtarget in sorted(target_assigned_count):
        target_assigned_prop = round(target_assigned_count[eachtarget]/number_q_in_psl*100, 1)
        to2file.write("\n"+eachtarget+"\t"+str(target_assigned_count[eachtarget])+"\t"+
                      str(target_assigned_prop)+"% of assigned reads"+"\n")
    to2file.write("\n"+str(len(ASCore_collection))+"\n")
    to2file.write("Below are counts of INDELs\n")
    to2file.write("TARGET\tEDITPOS:BP\tCount\n")
    for hit in sorted(PATTERN_COLLECTION):
        for pat in sorted(PATTERN_COLLECTION[hit], key=PATTERN_COLLECTION[hit].get ,reverse=True):
            to2file.write(hit+"\t"+pat+"\t"+str(PATTERN_COLLECTION[hit][pat])+"\n")
    for hit in sorted(MISMATCH_COLLECTION):
        for pat in sorted(MISMATCH_COLLECTION[hit], key=MISMATCH_COLLECTION[hit].get ,reverse=True):
            mis_freq = str(round(MISMATCH_COLLECTION[hit][pat]/target_assigned_count[hit],3))
            to2file.write(hit+"\t"+pat+"\t"+str(MISMATCH_COLLECTION[hit][pat])+
                          "\t"+mis_freq+"\n")
'''
    # Write output
    # Clean up


def check_conf():
    # checking the configuration file
    fp_AGEseqConf = pwd / "AGEseq.conf"
    print(fp_AGEseqConf)
    if os.path.exists(fp_AGEseqConf):
        print("Configuration File is found:"+str(fp_AGEseqConf))
        param = ageconf.readConfigFile(fp_AGEseqConf)
        USER_ASCONF = param[0]
        USER_BLCONF = param[1]
        USER_ASCONF.returnConfig()
        USER_BLCONF.returnBLATConfig()
        return param
    else:
        print("Couldn't find AGEseq.conf in your directory: "
              + str(pwd))
        print("Make sure you have AGEseq.conf with working folder: "
              + str(pwd))
        print("For details, please refer to the manual!")
        print("If you want to run AGEseq with the default please use -d mode")
        print("Ageseq has stopped!")
        exit


def check_blat():
    # checking environment blat

    if user_os == "win":
        DEF_BLAT_PATH = pwd / "blat.exe"
    elif "linux" in user_os or "darwin" in user_os:
        DEF_BLAT_PATH = pwd / "blat"
    else:
        print("OS couldn't be determined!")
        exit
    if not os.path.exists(DEF_BLAT_PATH):
        print("BLAT couldn't be located at " + str(DEF_BLAT_PATH))
        exit
    else:
        print("BLAT is located at " + str(DEF_BLAT_PATH))
    print("AGEseq will run on a "+user_os+" machine!")
    return DEF_BLAT_PATH

def load_target():
    #Load target
    v1_targetfile= pwd / "targets.txt"
    if os.path.exists(v1_targetfile):
        print("Target File is found:" + str(v1_targetfile))
        UTarget = AgeseqIO.readv1TargetFile(v1_targetfile)
        #UTarget.showAsFasta()
        UTarget.toFastaFile(TEMP_TARGET_FILE)
        return UTarget
    else:
        print("Couldn't locate target file:"+str(v1_targetfile))
        exit

def get_reads_file():
    if os.path.isdir(READS_PATH):
        print("Found Reads Directory: " + str(READS_PATH) + " successfully!")
        READS_FILE_LIST = os.listdir(READS_PATH)
        return READS_FILE_LIST
    else:
        print(str(READS_PATH) + "doesn't exists")
        exit

HELP = ""
import argparse
if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='AGEseq2')
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                        const=sum, default=max,
                        help='sum the integers (default: find the max)')

    args = parser.parse_args()
    # print(args.accumulate(args.integers))


    main()
    print(HELP + "\n")
    # tutorial()
