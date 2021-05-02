# -*- coding: utf-8 -*-
"""This is the main entry of AGEseq.
By running this script, AGEseq will automatically run with default setting.
"""

from __future__ import division
from handlers.indels_snp_store import IndelSNPStore
import argparse
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
import platform
import pandas as pd

pwd = Path.cwd()
user_os = platform.system()  # "Windows" "Linux" "Darwin"
prew_path = os.getcwd()
os.chdir(pwd)
READS_FILE_LIST = list()
DEF_BLAT_PATH = pwd / "blat"
""" record for script running """
sysdatetime = datetime.now()
dt_string = sysdatetime.strftime("%Y%m%d_%H_%M_%S")
logfile_name = "AGESeq_run_"+dt_string+".log"
logfile = open(logfile_name, "w")
WELCOM_MESSAGE = "Thank you for using AGEseq!\n"
mod = "legacy"
logfile.write(WELCOM_MESSAGE)
""" analysis result writing """
# raw_output = open("AGESeq_raw_"+dt_string+".csv", "w")
# summary_output = open("AGESeq_output_"+dt_string+".csv", "w")


""" intermediate files list """
intermediate_file_list = []


def main():
    """ checking path/dir/ """
    [USER_ASCONF, USER_BLCONF] = check_conf()

    """ load configuration """
    USER_BLCONF.returnBLATConfig()
    DEF_BLAT_PATH = check_blat()
    USER_TARGET = load_target()

    """ using MUSCLE to check input polymorphism """
    # target_feature = USER_TARGET.target_polymorphism()
    # logfile.write("Polymorphic sites in targets:"+"\n")
    # logfile.write(str(target_feature)+"\n")

    USER_READSFILE_LIST = get_reads_file()

    # init lists for DataFrame
    input_file_name_list = []
    target_list = []
    aligned_target_list = []
    aligned_read_list = []
    file_hits_num_list = []
    target_hits_num_list = []
    editing_pattern_hits_num_list = []
    editing_pattern_list = []

    """ for each reads file, get their read assigned to targets """
    for urfile in USER_READSFILE_LIST:
        read_file = AgeseqIO.readsToSample(urfile)
        read2fas_file = str(urfile)+".fa"
        intermediate_file_list.append(read2fas_file)
        read_file.toFastaFile(read2fas_file)

        """ blat cmd generation and execution """
        blat_out = f'{read2fas_file}_blat_crispr.psl'
        intermediate_file_list.append(blat_out)
        blat_cmd = f'{DEF_BLAT_PATH} {target_file}.fa {read2fas_file} ' \
                   f' -tileSize={USER_BLCONF.tileSize}' \
                   f' -oneOff={USER_BLCONF.oneOff}' \
                   f' -maxGap={USER_BLCONF.maxGap}' \
                   f' -minIdentity={USER_BLCONF.minIdentity}' \
                   f' -out="psl"' \
                   f' -minScore={USER_BLCONF.minScore} {blat_out}'
        logfile.write("Running blat as:\n"+str(blat_cmd)+"\n")
        subprocess.call(blat_cmd, shell=True)

        """ parsing blat result """
        [assigned_core, mismatch_collection, indel_collection,
            target_assigned_count] = psl_parse(read_file, USER_TARGET, blat_out)
        # target_var_freq = cal_target_var_freq(mismatch_collection,
        #                                       indel_collection,
        #                                       target_assigned_count)
        logfile.write("Assignment begins!"+"\n")

        # filter out unqualified mismatch and only retain editing pattern
        file_assign_sum = assign_mask(assigned_core, mismatch_collection,
                                                   indel_collection, target_assigned_count, WOBBLE_BASE)
        """ register each editing pattern to list for fil writing """
        for target in file_assign_sum:
            # target: 'AMD1a'
            for ep in file_assign_sum[target]:
                input_file_name_list.append(urfile)
                target_list.append(target)
                aligned_target_list.append(USER_TARGET.getTargetSeq(target))
                # I cannot figure out where to get the original reads sequence
                aligned_read_list.append("aligned reads #issue")
                file_hits_num_list.append(
                                    "file wise hits number #issue")
                target_hits_num_list.append(target_assigned_count[target])
                editing_pattern_hits_num_list.append(
                    file_assign_sum[target][ep])
                editing_pattern_list.append(ep)
                logfile.write(target+"\t"+ep+"\t" +
                              str(file_assign_sum[target][ep])+"\n")
        # percentage of target
        for eachtarget in sorted(target_assigned_count):
            target_assigned_prop=round(
                target_assigned_count[eachtarget]/len(assigned_core)*100, 1)
            logfile.write("\n"+eachtarget+"\t"+str(target_assigned_count[eachtarget])+"\t" +
                          str(target_assigned_prop)+"% of assigned reads"+"\n")

    """ write dataframe and to csv """
    summary_table=pd.DataFrame.from_dict({
        "input_file": input_file_name_list,
        "target": target_list,
        "aligned_target": aligned_target_list,
        "aligned_read": aligned_read_list,
        "Total hits": file_hits_num_list,
        "sub hits": target_hits_num_list,
        "editing pattern hits": editing_pattern_hits_num_list,
        "editing pattern": editing_pattern_list
    })
    summary_table.to_csv("AGESeq_summary_"+dt_string+".csv")
    """ clean up intermediate files """
    for intermediate_file in intermediate_file_list:
        os.remove(intermediate_file)


def cal_target_var_freq(mismatch_collection, indel_collection,
                        target_assigned_count):
    target_var_freq=defaultdict(dict)
    for target in target_assigned_count:
        for mismatch in mismatch_collection[target]:
            totalmismatch=mismatch_collection[target][mismatch]
            target_var_freq[target][mismatch]=totalmismatch
        for mismatch in indel_collection[target]:
            totalmismatch=indel_collection[target][mismatch]
            target_var_freq[target][mismatch]=totalmismatch
    return target_var_freq


def assign_mask(assigned_core, mismatch_collection,
                indel_collection, target_assigned_count, WOBBLE_BASE):
    """ filter out those unqualified mismatch and retain only editing pattern """
    file_summary=defaultdict(dict)

    for reads, as_core in assigned_core.items():
        coreBH=as_core.getBH()  # target
        coreBHscore=as_core.getBHScore()  # int
        coremismatch=as_core.getPSLMismatch()  # ['146:A->G']
        coremismatch_count=as_core.getMismatchCount()  # int
        # "1 Mismatch(es)\t['146:A->G']\n"
        coremismatchstr=as_core.printmismatch()
        coreINDEL=as_core.getINDEL()
        logfile.write(str(reads)+"\t"+coreBH+"\t"+str(coreBHscore)+"\n")
        logfile.write(coremismatchstr)
        var_mask=dict()
        """ go through each mismatch """
        if coremismatch_count > 0:
            for str_mismatch in coremismatch:
                totalmismatch=mismatch_collection[coreBH][str_mismatch]
                misfreq=round(totalmismatch/target_assigned_count[coreBH], 3)
                """ too few mapped reads and wobble base case """
                if misfreq < READ_SNP_MINR_FREQ or totalmismatch < READ_SNP_MINIMAL_SUPPORT:
                    var_mask[str_mismatch]="masked"
                    logfile.write(str_mismatch+"\t"+str(misfreq)+"\t" +
                                  " removed!\n")
                else:
                    if WOBBLE_BASE:
                        target_range=as_core.getHSPrange()
                        mispos=int(str_mismatch.split(":")[0])
                        if (mispos - target_range[0] < 12 or abs(mispos - target_range[1]) < 12) and (misfreq < WOBBLE_FREQ_RANGE[1] and misfreq > WOBBLE_FREQ_RANGE[0]):
                            var_mask[str_mismatch]="masked"
                            logfile.write(
                                str_mismatch+"\t"+str(misfreq)+" masked as a wobblebase!\n")
                        else:
                            logfile.write(str_mismatch+"\t" +
                                          str(misfreq)+" kept!\n")
                    else:
                        logfile.write(str_mismatch+"\t" +
                                      str(misfreq)+" kept!\n")
        """ go through each INDEL """
        if len(coreINDEL) > 0:
            logfile.write(str(len(coreINDEL))+" INDEL(s) is/are found!\n")
            for indel in coreINDEL:
                totalindel=indel_collection[coreBH][indel]
                indelfreq=round(totalindel/target_assigned_count[coreBH], 3)
                if indelfreq < READ_INDEL_MINR_FREQ or totalindel < READ_INDEL_MINIMAL_SUPPORT:
                    var_mask[indel]="masked"
                    logfile.write(indel+"\t"+str(indelfreq)+"\t" +
                                  " removed!\n")
                else:
                    logfile.write(indel+"\t"+str(indelfreq)+" kept!\n\n")
        else:
            logfile.write("No INDEL(s)!"+"\n\n")
        as_core.updateMask(var_mask)
        as_core.assignEditPat()
        coreEP_list=as_core.getEditPat()
        
        # ['221:C->A', '227:T->A', '229:A->C', '233:C->G', '241:C->T', '244:A->T', '250:C->T', '254:C->A', '263:C->T', '267:G->T', '280:C->A', '285:G->A', '288:C->T', '298:C->G']
        # ['181D7:CAAAGAG', '193I1:A', '209:T->C', '211:C->T', '220:C->T', '240:G->A', '243:G->A', '245:C->A', '259:G->T', '268:G->T', '273:G->T', '278:C->A', '285:G->A', '286:A->T']
        # EP Editing pattern
        # nested double dictionary
        
        # data=IndelSNPStore.parseCoreEP(coreEP_list)
        
        # BH best hit: reads hitting target
        # EP editing pattern:
        for coreEP in coreEP_list:
            if coreBH in file_summary:
                if coreEP in file_summary[coreBH]:
                    file_summary[coreBH][coreEP] += 1
                else:
                    file_summary[coreBH][coreEP]=1
            else:
                file_summary[coreBH][coreEP]=1
            logfile.write(str(coreEP)+" assigned!\n\n")

    return file_summary


# def masked_edit_pattern(assigned_core):
#     assigned_core


def psl_parse(read_file, target_file, blat_out):
    """
        psl file
        to
        assigned_core, mismatch_collection, indel_collection, target_assigned_count

    """
    sample_readslist=read_file.getallreadsids()
    number_reads_in_sample=len(sample_readslist)
    ASCore_collection=dict()
    psl_dict=AgeseqIO.parseBLATlegacymod(blat_out)
    """
    {
        reads_id: QueryResult object
    }
    """
    q_besthit_score=dict()
    q_besthit=dict()
    tier_count=0
    number_q_in_psl=len(psl_dict)
    perc_hit=number_q_in_psl/number_reads_in_sample*100

    logfile.write(
        f"{number_q_in_psl} reads have hits to the target sequences provided ({perc_hit}% of the sample)!\n")

    # find the best hit?
    for qid in psl_dict:
        # qid: reads id
        ASCore_collection[qid]=AgeseqIO.ASCore(qid)
        psl_t_dict=dict()
        for hit in psl_dict[qid]:
            for hsp in hit:
                psl_t_dict[hit.id]=hsp.score
                if qid in q_besthit.keys():
                    if hsp.score > q_besthit_score[qid]:
                        q_besthit_score[qid]=hsp.score
                        q_besthit[qid]=hit.id
                    elif hsp.score == q_besthit_score[qid]:
                        logfile.write(qid+" Tie Hit found!"+str(hit.id)+"\n")
                        logfile.write(qid+" Tie Hit found!" +
                                      q_besthit[qid]+"\n")
                        tier_count += 1
                    else:
                        continue
                else:
                    q_besthit_score[qid]=hsp.score
                    q_besthit[qid]=hit.id
        ASCore_collection[qid].updatePSLTargets(psl_t_dict)

    # q_besthit_score[read] = hsp.score
    # q_besthit[read] = target
    # psl_t_dict
    # {   target: score   }

    print("Found %i reads have tie hit" % tier_count)

    PATTERN_COLLECTION=defaultdict(dict)
    MISMATCH_COLLECTION=defaultdict(dict)
    target_assigned_count=dict()

    for query in q_besthit:
        ASCore_collection[query].updateBH(q_besthit[query])
        q_hsp=psl_dict[query][q_besthit[query]]
        ASCore_collection[query].updateHSP(q_hsp)
        q_seq=read_file.getReadSeq(query)
        t_seq=target_file.getTargetSeq(q_besthit[query])
        if q_besthit[query] in target_assigned_count:
            target_assigned_count[q_besthit[query]] += 1
        else:
            target_assigned_count[q_besthit[query]]=1
        # logfile.write(str(query)+"\t"+q_besthit[query]+
        #              "\t"+str(q_besthit_score[query])+"\n")
        # to2file.write(str(ASCore_collection[query].getPSLTargets())+"\n")
        for ehsp in q_hsp:
            q_start_all=ehsp.query_start_all
            q_span_all=ehsp.query_span_all
            h_start_all=ehsp.hit_start_all
            h_span_all=ehsp.hit_span_all
            q_gap=ehsp.hit_gapopen_num
            h_gap=ehsp.query_gapopen_num
            mismatch=ehsp.mismatch_num
            hit_frag_seq=[]
            query_frag_seq=[]
            a=len(ehsp)
            strand_info=[]
            for i in range(0, a):
                if ehsp[i].query_strand == -1:
                    strand_info.append("-")
                else:
                    strand_info.append("+")
            for i in range(0, len(h_start_all)):
                s=t_seq[h_start_all[i]:h_start_all[i]+h_span_all[i]]
                hit_frag_seq.append(s)
                fq_seq=q_seq[q_start_all[i]:q_start_all[i] + q_span_all[i]]
                query_frag_seq.append(fq_seq)
            if strand_info[0] == "-":
                ucount=0
                for seqs in query_frag_seq:
                    query_frag_seq[ucount]=Seq.reverse_complement(
                        query_frag_seq[ucount])
                    ucount += 1
        psl_str=str("*"*ehsp.gapopen_num+str(q_start_all)+" "+str(q_span_all)
                      + " "+str(h_start_all)+" "+str(h_span_all)+" "+str(q_gap)
                      + " "+str(h_gap)+" "+str(strand_info)+" "
                      + str(mismatch)+"\n")
        ASCore_collection[qid].updateHSPstr(psl_str)
        MISMATCH=dict()
        mismatch_count=0
        for i in range(0, len(hit_frag_seq)):
            blast_link=[]
            for j in range(0, len(hit_frag_seq[i])):
                if hit_frag_seq[i][j] == query_frag_seq[i][j]:
                    blast_link.append("*")
                else:
                    blast_link.append(".")
                    mismatch_pos=h_start_all[i] + j + 1
                    MISMATCH[mismatch_count]=str(
                        mismatch_pos) + ":" + str(hit_frag_seq[i][j]) + "->"+str(query_frag_seq[i][j])
                    if MISMATCH[mismatch_count] in MISMATCH_COLLECTION[q_besthit[query]]:
                        MISMATCH_COLLECTION[q_besthit[query]
                                            ][MISMATCH[mismatch_count]] += 1
                    else:
                        MISMATCH_COLLECTION[q_besthit[query]
                                            ][MISMATCH[mismatch_count]]=1
                    mismatch_count += 1
            # logfile.write(str(hit_frag_seq[i])+"\n")
            # logfile.write(str(query_frag_seq[i]) +"\n")
            # logfile.write(str("".join(blast_link)) + "\n\n")
        # to2file.write(str(MISMATCH) + "\n\n")

        if mismatch == 0:
            pass
        else:
            ASCore_collection[query].updateMismatch(MISMATCH)
            ASCore_collection[query].updateMismatchCount(mismatch)
        # logfile.write(str(ASCore_collection[query].mismatch_count)+" Mismatch(es)\t"+str(ASCore_collection[query].getPSLMismatch())+"\n")
        for i in range(0, q_gap):
            if h_start_all[i] + h_span_all[i] < h_start_all[i+1]:
                del_start=h_start_all[i] + h_span_all[i]
                del_end=h_start_all[i+1] - 1
                del_size=del_end - del_start + 1
                del_seq=t_seq[del_start: del_start+del_size]
                del_pat=str(del_start) + "D" + \
                    str(del_size) + ":" + str(del_seq)
                ASCore_collection[query].addINDEL(del_pat)
                if del_pat in PATTERN_COLLECTION[q_besthit[query]]:
                    PATTERN_COLLECTION[q_besthit[query]][del_pat] += 1
                else:
                    PATTERN_COLLECTION[q_besthit[query]][del_pat]=1
                # logfile.write(del_pat + "\n\n")
        # insert
        for i in range(0, h_gap):
            if h_start_all[i] + h_span_all[i] == h_start_all[i+1]:
                ins_start_ontarget=h_start_all[i+1]
                ins_start_onquery=q_start_all[i] + q_span_all[i]
                ins_end_onquery=q_start_all[i+1] - 1
                if strand_info[0] == "-":
                    ins_start_onquery=q_start_all[i+1] + q_span_all[i+1]
                    ins_end_onquery=q_start_all[i] - 1
                ins_size=ins_end_onquery - ins_start_onquery + 1
                ins_seq=q_seq[ins_start_onquery: ins_start_onquery+ins_size]
                ins_pat=str(ins_start_ontarget) + "I" + \
                    str(ins_size) + ":" + str(ins_seq)
                ASCore_collection[query].addINDEL(ins_pat)
                if ins_pat in PATTERN_COLLECTION[q_besthit[query]]:
                    PATTERN_COLLECTION[q_besthit[query]][ins_pat] += 1
                else:
                    PATTERN_COLLECTION[q_besthit[query]][ins_pat]=1
                # logfile.write(ins_pat + "\n\n")

    return ASCore_collection, MISMATCH_COLLECTION, PATTERN_COLLECTION, target_assigned_count


'''
    for eachtarget in sorted(target_assigned_count):
        target_assigned_prop = round(
            target_assigned_count[eachtarget]/number_q_in_psl*100, 1)
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
            mis_freq = str(
                round(MISMATCH_COLLECTION[hit][pat]/target_assigned_count[hit],3))
            to2file.write(hit+"\t"+pat+"\t"+str(MISMATCH_COLLECTION[hit][pat])+
                          "\t"+mis_freq+"\n")
'''
# Write output
# Clean up


def check_conf():
    """ checking the configuration file """
    fp_AGEseqConf=pwd / "AGEseq.conf"
    print(fp_AGEseqConf)
    if os.path.exists(fp_AGEseqConf):
        print("Configuration File is found:"+str(fp_AGEseqConf))
        param=ageconf.readConfigFile(fp_AGEseqConf)
        USER_ASCONF=param[0]
        USER_BLCONF=param[1]
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
        exit()


def check_blat():
    """ Checks the user's system to determine which BLAT binary to use.
    """
    if "win" in user_os.lower():
        DEF_BLAT_PATH=pwd / "blat_binaries" / "blat.exe"
    elif "linux" in user_os.lower():
        DEF_BLAT_PATH=pwd / "blat_binaries" / "blat_linux"
    elif "darwin" in user_os:
        DEF_BLAT_PATH=pwd / "blat_binaries" / "blat_macos"
    else:
        print("OS couldn't be determined!")
        exit()
    if not os.path.exists(DEF_BLAT_PATH):
        print("BLAT couldn't be located at " + str(DEF_BLAT_PATH))
        exit()
    else:
        print("BLAT is located at " + str(DEF_BLAT_PATH))
    print("AGEseq will run on a "+user_os+" machine!")
    return DEF_BLAT_PATH


def load_target():
    # Load target
    v1_targetfile=pwd / "targets.txt"
    if os.path.exists(v1_targetfile):
        print("Target File is found:" + str(v1_targetfile))
        UTarget=AgeseqIO.readv1TargetFile(v1_targetfile)
        # UTarget.showAsFasta()
        UTarget.toFastaFile(target_file+'.fa')
        return UTarget
    else:
        print("Couldn't locate target file:"+str(v1_targetfile))
        exit()


def get_reads_file():
    if os.path.isdir(READS_PATH):
        print("Found Reads Directory: " + str(READS_PATH) + " successfully!")
        import glob
        READS_FILE_LIST=glob.glob(os.path.join(READS_PATH, '*'))

        return READS_FILE_LIST
    else:
        print(str(READS_PATH) + "doesn't exists")
        exit()


HELP=""
if __name__ == '__main__':

    # ageseq.py \
    # -t target.fa \
    # -r /dir/of/reads/ \
    # --READ_SNP_MINR_FREQ 0.05 \
    # --READ_INDEL_MINR_FREQ 0.01 \
    # --READ_SNP_MINIMAL_SUPPORT 3 \
    # --READ_INDEL_MINIMAL_SUPPORT 5 \
    # --WOBBLE_BASE "true" \
    # --WOBBLE_FREQ_RANGE 0.35 0.65
    parser=argparse.ArgumentParser(
        description='AGEseq2 introduction words #issue')
    parser.add_argument('-t', type=str,
                        help='fasta format target file')
    parser.add_argument('-r', type=str,
                        help='fasta/fastq format reads dir')
    parser.add_argument('--READ_SNP_MINR_FREQ', type=float,
                        help='')
    parser.add_argument('--READ_INDEL_MINR_FREQ', type=float,
                        help='')
    parser.add_argument('--READ_SNP_MINIMAL_SUPPORT', type=float,
                        help='')
    parser.add_argument('--READ_INDEL_MINIMAL_SUPPORT', type=float,
                        help='')
    parser.add_argument("-w", "--WOBBLE_BASE", help="remove wobble base",
                        action="store_true")
    parser.add_argument('--WOBBLE_FREQ_RANGE', nargs=2, type=float,
                        help='')
    args=parser.parse_args()

    target_file=args.t
    READS_PATH=args.r

    # issue assign value if specified in command
    READ_SNP_MINR_FREQ=0.05  # Target specific, not overall frequency
    READ_INDEL_MINR_FREQ=0.01  # Target specific, not overall frequency
    READ_SNP_MINIMAL_SUPPORT=3
    READ_INDEL_MINIMAL_SUPPORT=5
    WOBBLE_BASE=True
    WOBBLE_FREQ_RANGE=(0.35, 0.65)

    print(READS_PATH)
    main()
