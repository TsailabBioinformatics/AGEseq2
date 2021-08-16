# -*- coding: utf-8 -*-
"""This is the main entry of AGEseq.
By running this script, AGEseq will automatically run with default setting.
"""

from __future__ import division
from handlers.indels_snp_store import IndelSNPStore
import argparse
from datetime import datetime
from copy import deepcopy
import os
import sys
import re
import AgeseqConfig as ageconf
import AgeseqIO
import subprocess
from collections import defaultdict
from Bio import Seq
import pathlib
from pathlib import Path
import platform
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('-t', '--target', type=str,
                    help='target file[targets.txt]', default="targets.txt")
parser.add_argument('-r', '--read', type=str,
                    help='fasta/fastq format reads dir', default="./reads")
parser.add_argument('-sa', '--skip_aln', type=int, help='set 0 to show alignment in the log file.',
                    default=1)
args = parser.parse_args()

TARGET_FILE = args.target
READS_PATH = args.read
SKIP_ALN = args.skip_aln

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
mod = "legacy"
logfile.write("Thank you for using AGEseq!\n")

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
    USER_TARGET = load_target(TARGET_FILE)
    fasta_TARGET_FILE = USER_TARGET.fastaname
    """ using MUSCLE to check input polymorphism """
    # target_feature = USER_TARGET.target_polymorphism()
    # logfile.write("Polymorphic sites in targets:"+"\n")
    # logfile.write(str(target_feature)+"\n")

    USER_READSFILE_LIST = get_reads_file()

    # initiation of lists for output DataFrame
    input_file_name_list = []
    target_list = []
    aligned_target_list = []
    aligned_read_list = []
    #file_hits_num_list = []
    target_hits_num_list = []
    editing_pattern_hits_num_list = []
    editing_pattern_list = []
    g_editPat_collapse = defaultdict(dict)
    """ for each reads file, get their read assigned to targets """
    for urfile in USER_READSFILE_LIST:
        print("Readingfile"+str(urfile))
        read_file = AgeseqIO.readsToSample(urfile)

        read2fas_file = str(urfile)+".fa"
        intermediate_file_list.append(read2fas_file)
        read_file.toFastaFile(read2fas_file)
        sub_USER_TARGET = load_target(TARGET_FILE)
        """ blat cmd generation and execution """
        blat_out = f'{read2fas_file}_blat_crispr.psl'
        intermediate_file_list.append(blat_out)
        blat_cmd = f'{DEF_BLAT_PATH} {fasta_TARGET_FILE} {read2fas_file} ' \
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
        [file_assign_sum, target_consensus_stat, aligned_read_consensus] = assign_mask(assigned_core, mismatch_collection,
                                      indel_collection, target_assigned_count,
                                      USER_ASCONF.WOBBLE_BASE,
                                      USER_ASCONF.READ_SNP_MINR_FREQ,
                                      USER_ASCONF.READ_SNP_MINIMAL_SUPPORT,
                                      USER_ASCONF.READ_INDEL_MINR_FREQ,
                                      USER_ASCONF.READ_INDEL_MINIMAL_SUPPORT,
                                      USER_ASCONF.WOBBLE_FREQ_HIGH,
                                      USER_ASCONF.WOBBLE_FREQ_LOW,
                                      sub_USER_TARGET, urfile)
        #logfile.write("#####DEBUG AREA START####\n\n")
        logfile.write("\n#Following Eidted Pattern needs to be merged#\n\n")
        merged_file_assign_sum = deepcopy(file_assign_sum)
        for cseq in target_consensus_stat:
            #logfile.write(cseq+"\n")
            for bh in target_consensus_stat[cseq]:
                #logfile.write("+"+str(bh)+"\n")
                if len(target_consensus_stat[cseq][bh]) > 1:
                    top_ep_val = 0
                    top_ep_id = ""
                    for ep in target_consensus_stat[cseq][bh]:
                        if target_consensus_stat[cseq][bh][ep] > top_ep_val:
                            top_ep_id = ep
                            top_ep_val = target_consensus_stat[cseq][bh][ep]
                    for ep in target_consensus_stat[cseq][bh]:
                        if ep != top_ep_id:
                            if ep in g_editPat_collapse[bh]:
                                continue
                            else:
                                g_editPat_collapse[bh][ep] = top_ep_id
                                logfile.write("Edited Pattern: "+ep+"\t"+bh+" will be merged with "+top_ep_id+"\n")
        #logfile.write("\n#####AREA END####\n\n")
        for bh in file_assign_sum:
            for ep in file_assign_sum[bh]:
                if ep in g_editPat_collapse[bh]:
                    original_ep_count = merged_file_assign_sum[bh][ep]
                    collapse_ep = g_editPat_collapse[bh][ep]
                    merged_file_assign_sum[bh][collapse_ep] = merged_file_assign_sum[bh][collapse_ep] + original_ep_count
                    del merged_file_assign_sum[bh][ep]

        # percentage of target
        for eachtarget in sorted(target_assigned_count):
            target_assigned_prop = round(
                target_assigned_count[eachtarget]/len(assigned_core)*100, 1)
            logfile.write("\n"+eachtarget+"\t"+str(target_assigned_count[eachtarget])+"\t"
                          + str(target_assigned_prop)+"% of assigned reads"+"\n")

        """ register each editing pattern to list for file writing """
        for bh in merged_file_assign_sum:
            for ep in merged_file_assign_sum[bh]:
                urfilename=urfile.split("/")
                if "windows" in platform.system().lower():
                    urfilename=urfile.split("\\")
                input_file_name_list.append(urfilename[-1])
                target_list.append(bh)
                aligned_target_list.append(sub_USER_TARGET.getAlignedTarget(bh, ep))
                aligned_read_list.append(aligned_read_consensus[bh][ep])
                target_hits_num_list.append(target_assigned_count[bh])
                ep_hitnum = 0
                if ep in g_editPat_collapse[bh]:
                    ep_hitnum = ep_hitnum + merged_file_assign_sum[bh][g_editPat_collapse[bh][ep]]
                else:
                    ep_hitnum = merged_file_assign_sum[bh][ep]
                editing_pattern_hits_num_list.append(ep_hitnum)
                editing_pattern_list.append(ep)
                logfile.write(bh+"\t"+ep+"\t"
                              + str(merged_file_assign_sum[bh][ep])+"\n")
    """ write dataframe and to csv """
    summary_table = pd.DataFrame.from_dict({
        "input_file": input_file_name_list,
        "target": target_list,
        "aligned_target": aligned_target_list,
        "aligned_consensus": aligned_read_list,
        "sub hits": target_hits_num_list,
        "editing pattern hits": editing_pattern_hits_num_list,
        "editing pattern": editing_pattern_list
    })
    summary_table.to_csv("AGESeq_summary_"+dt_string+".csv")
    logfile.write("\nSummary File has been written to AGESeq_summary_"+dt_string+".csv"+"\n")
    logfile.write("\n"+"-"*60+"\n")
    """ clean up intermediate files """
    for intermediate_file in intermediate_file_list:
        os.remove(intermediate_file)


def cal_target_var_freq(mismatch_collection, indel_collection,
                        target_assigned_count):
    """Count the variants from mismatches and indels collected per target
    :param mismatch_collection: dictionary, mismatch_collection
    :param indel_collection: dictionary, indel_collection
    :param target_assigned_count: dictionary, target_assigned_count
    """
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
                indel_collection, target_assigned_count,
                WOBBLE_BASE, READ_SNP_MINR_FREQ, READ_SNP_MINIMAL_SUPPORT,
                READ_INDEL_MINR_FREQ, READ_INDEL_MINIMAL_SUPPORT,
                WOBBLE_FREQ_HIGH, WOBBLE_FREQ_LOW, sub_USER_TARGET, urfile):
    """Assign editted pattern based on the mismatch and indels identified from
    alignment.
    :param mismatch_collection: dictionary, mismatch_collection
    :param indel_collection: dictionary, indel_collection
    :param target_assigned_count: dictionary, target_assigned_count
    :return:
         - file_summary - original summary of edited patterns identified per target
         - target_consensus_stat - a dictionary that contains detailed counts of consensus sequences per file
         - aligned_read_consensus - a dictionary that contains consensus sequences indexed by their target and edited pattern
    """
    file_summary = defaultdict(dict)
    target_consensus_stat = defaultdict(dict)
    aligned_read_consensus = defaultdict(dict)
    for reads, as_core in assigned_core.items():
        coreBH = as_core.getBH()  # target
        coreBHscore = as_core.getBHScore()  # int
        coremismatch = as_core.getPSLMismatch()  # ['146:A->G']
        coremismatch_count = as_core.getMismatchCount()  # int
        # "1 Mismatch(es)\t['146:A->G']\n"
        coremismatchstr = as_core.printmismatch()
        coreINDEL = as_core.getINDEL()
        logfile.write(str(reads)+"\t"+coreBH+"\t"+str(coreBHscore)+"\n")
        logfile.write(coremismatchstr)
        var_mask = dict()
        """ go through each mismatch """
        if coremismatch_count > 0:
            pos_mask_list = []
            for str_mismatch in coremismatch:
                totalmismatch = mismatch_collection[coreBH][str_mismatch]
                misfreq = round(totalmismatch/target_assigned_count[coreBH], 3)
                misstr = str_mismatch.split(':')
                pos_x = misstr[0]
                """ too few mapped reads and wobble base case """
                if misfreq < READ_SNP_MINR_FREQ or totalmismatch < READ_SNP_MINIMAL_SUPPORT:
                    var_mask[str_mismatch] = "masked"
                    pos_mask_list.append(pos_x)
                    logfile.write(str_mismatch+"\t"+str(misfreq)+"\t"
                                  + " removed!\n")
                else:
                    if WOBBLE_BASE:
                        target_range = as_core.getHSPrange()
                        mispos = int(str_mismatch.split(":")[0])
                        if (mispos - target_range[0] < 12 or abs(mispos - target_range[1]) < 12) and (misfreq < WOBBLE_FREQ_HIGH and misfreq > WOBBLE_FREQ_LOW):
                            var_mask[str_mismatch] = "masked"
                            pos_mask_list.append(pos_x)
                            logfile.write(
                                str_mismatch+"\t"+str(misfreq)+" masked as a wobblebase!\n")
                        else:
                            logfile.write(str_mismatch+"\t"
                                          + str(misfreq)+" kept!\n")
                    else:
                        logfile.write(str_mismatch+"\t"
                                      + str(misfreq)+" kept!\n")
                sub_USER_TARGET.maskTargetSNP(coreBH, pos_mask_list)
        """ go through each INDEL """
        if len(coreINDEL) > 0:
            logfile.write(str(len(coreINDEL))+" INDEL(s) is/are found!\n")
            for indel in coreINDEL:
                totalindel = indel_collection[coreBH][indel]
                indelfreq = round(totalindel/target_assigned_count[coreBH], 3)
                if indelfreq < READ_INDEL_MINR_FREQ or totalindel < READ_INDEL_MINIMAL_SUPPORT:
                    var_mask[indel] = "masked"
                    logfile.write(indel+"\t"+str(indelfreq)+"\t("+str(totalindel)
                                  + ") removed!\n")
                else:
                    logfile.write(indel+"\t"+str(indelfreq)+" kept!\n\n")
        else:
            logfile.write("No INDEL(s)!"+"\n\n")
        as_core.updateMask(var_mask)
        as_core.assignEditPat()
        coreEP = str(as_core.getEditPat()[0])
        aln_str = as_core.aln2str()
        aln_feature_str = as_core.aln2feature(sub_USER_TARGET)
        #logfile.write("INDEL:"+str(as_core.getINDEL())+"\n")
        if SKIP_ALN == 0:
            logfile.write("Alignment:"+"\n"+str(aln_str[0])+"\t (T)\n")
            logfile.write(str(aln_str[1])+"\t (R)\n"+aln_feature_str[0]+"\t (F)\n")
            logfile.write(as_core.printConsensus()+"\t (C)\n\n")
            #logfile.write(striped_consenseq+"\t (SC)\n\n")

        #create a record for aligned target sequence

        if coreBH in aligned_read_consensus:
            if coreEP in aligned_read_consensus[coreBH]:
                if abs(len(aln_str[0])-len(as_core.printConsensus())) < abs(len(sub_USER_TARGET.aligned_target[coreBH][coreEP])-len(aligned_read_consensus[coreBH][coreEP])):
                    aligned_read_consensus[coreBH][coreEP] = as_core.printConsensus()
                    sub_USER_TARGET.aligned_target[coreBH][coreEP] = aln_str[0]
            else:
                aligned_read_consensus[coreBH][coreEP] = as_core.printConsensus()
                sub_USER_TARGET.aligned_target[coreBH][coreEP]= aln_str[0]
        else:
            aligned_read_consensus[coreBH] = dict()
            aligned_read_consensus[coreBH][coreEP] = as_core.printConsensus()
            sub_USER_TARGET.aligned_target[coreBH] = dict()
            sub_USER_TARGET.aligned_target[coreBH][coreEP] = aln_str[0]

        # Separated INDEL patterns caused by BLAT needs to be collapsed before being assigned to final stat
        striped_consenseq = re.sub(r'-', '', as_core.printConsensus())
        # logfile.write("***\t"+coreBH+"\t"+coreEP+"\n")
        #if striped_consenseq == "ttacgcggcaattataagtgtaggttagcactgactggattgcctgacgaagtttatgataaggaatgggatttgataatgattgatgcaagagggtacttcccggaggcaccagggaggatggcggcgatattttcagcggcggtgatgg":
        #    if coreEP == "193D3:GCC" and coreBH == "AMD2a":
        #        logfile.write("$$$$$$$$$$$$$$ECX\n\n")
        #    if coreEP == "WT" and coreBH == "AMD2a":
        #        logfile.write("$$$$$$$$$$$$$$WCX\n\n")

        if striped_consenseq in target_consensus_stat:
            if coreBH in target_consensus_stat[striped_consenseq]:
                if coreEP in target_consensus_stat[striped_consenseq][coreBH]:
                    target_consensus_stat[striped_consenseq][coreBH][coreEP] += 1
                else:
                    target_consensus_stat[striped_consenseq][coreBH][coreEP] = 1
            else:
                target_consensus_stat[striped_consenseq][coreBH][coreEP] = 1
        else:
            target_consensus_stat[striped_consenseq]=dict()
            target_consensus_stat[striped_consenseq][coreBH]=dict()
            target_consensus_stat[striped_consenseq][coreBH][coreEP] = 1

        if coreBH in file_summary:
            if coreEP in file_summary[coreBH]:
                file_summary[coreBH][coreEP] += 1
            else:
                file_summary[coreBH][coreEP] = 1
        else:
            file_summary[coreBH][coreEP] = 1
        logfile.write(str(coreEP)+" assigned!\n\n")

    return file_summary, target_consensus_stat, aligned_read_consensus


def psl_parse(read_file, target_file, blat_out):
    """Parse the output psl file from BLAT to multiple dictionaries for reads (with hit only),
       all mismatches, all indels, and number of reads assigned to each target
       :param read_file: reads file
       :param target_file: target file
       :param blat_out: the PSL output from BLAT
        assigned_core, mismatch_collection, indel_collection, target_assigned_count
       :return:
            - ASCore_collection - a dictionary that contains the all assigned read IDs as keys mapped to ASCore instances
            - MISMATCH_COLLECTION - a default dictionary that contains the all identified mismatches (indexed by target, and mismatch string).
            - PATTERN_COLLECTION - Similar to MISMATCH_COLLECTION except being used for INDELS
            - target_assigned_count - a dictionary hoding target keys and number of reads hitting to the target
    """
    sample_readslist = read_file.getallreadsids()
    number_reads_in_sample = len(sample_readslist)
    ASCore_collection = dict()
    psl_dict = AgeseqIO.parseBLATlegacymod(blat_out)
    """
    {
        reads_id: QueryResult object
    }
    """
    q_besthit_score = dict()
    q_besthit = dict()
    tier_count = 0
    number_q_in_psl = len(psl_dict)
    perc_hit = number_q_in_psl/number_reads_in_sample*100
    logfile.write(
        f"{number_q_in_psl} reads have hits to the target sequences provided ({perc_hit}% of the sample)!\n")
    # find the best hit?
    for qid in psl_dict:
        # qid: reads id
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
                        logfile.write(qid+" Tie Hit found!"+str(hit.id)+" & "
                                      + q_besthit[qid]+"\n")
                        tier_count += 1
                    else:
                        continue
                else:
                    q_besthit_score[qid] = hsp.score
                    q_besthit[qid] = hit.id
        ASCore_collection[qid].updatePSLTargets(psl_t_dict)

    # q_besthit_score[read] = hsp.score
    # q_besthit[read] = target
    # psl_t_dict
    # {   target: score   }

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
        # logfile.write(str(query)+"\t"+q_besthit[query]+
        #              "\t"+str(q_besthit_score[query])+"\n")
        # to2file.write(str(ASCore_collection[query].getPSLTargets())+"\n")
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
            strand_info = []
            for i in range(0, len(ehsp)):
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
                    query_frag_seq[ucount] = Seq.reverse_complement(
                        query_frag_seq[ucount])
                    ucount += 1
            ASCore_collection[query].hit_aligned_block = hit_frag_seq
            ASCore_collection[query].query_aligned_block = query_frag_seq
        psl_str = ["*"*ehsp.gapopen_num, q_start_all, q_span_all,
                   h_start_all, h_span_all, q_gap,
                   h_gap, strand_info,
                   mismatch]
        ASCore_collection[query].updateHSPstr(psl_str)
        MISMATCH = dict()
        mismatch_count = 0
        #feature_align_block = []
        for i in range(0, len(hit_frag_seq)):
            blast_link = []
            for j in range(0, len(hit_frag_seq[i])):
                if hit_frag_seq[i][j] == query_frag_seq[i][j]:
                    blast_link.append("*")
                else:
                    mismatch_pos = h_start_all[i] + j + 1
                    MISMATCH[mismatch_count] = str(
                        mismatch_pos) + ":" + str(hit_frag_seq[i][j]) + "->"+str(query_frag_seq[i][j])
                    if MISMATCH[mismatch_count] in MISMATCH_COLLECTION[q_besthit[query]]:
                        MISMATCH_COLLECTION[q_besthit[query]
                                            ][MISMATCH[mismatch_count]] += 1
                    else:
                        MISMATCH_COLLECTION[q_besthit[query]
                                            ][MISMATCH[mismatch_count]] = 1
                    mismatch_count += 1
            #feature_align_block.append(str("".join(blast_link)))
        #ASCore_collection[query].updatebesthitalign(align_block)

        if mismatch != 0:
            ASCore_collection[query].updateMismatch(MISMATCH)
            ASCore_collection[query].updateMismatchCount(mismatch)
        # logfile.write(str(ASCore_collection[query].mismatch_count)+" Mismatch(es)\t"+str(ASCore_collection[query].getPSLMismatch())+"\n")
        for i in range(0, q_gap):
            if h_start_all[i] + h_span_all[i] < h_start_all[i+1]:
                del_start = h_start_all[i] + h_span_all[i]
                del_end = h_start_all[i+1] - 1
                del_size = del_end - del_start + 1
                del_seq = t_seq[del_start: del_start+del_size]
                del_pat = str(del_start) + "D" + \
                    str(del_size) + ":" + str(del_seq)
                ASCore_collection[query].addINDEL(del_pat)
                if del_pat in PATTERN_COLLECTION[q_besthit[query]]:
                    PATTERN_COLLECTION[q_besthit[query]][del_pat] += 1
                else:
                    PATTERN_COLLECTION[q_besthit[query]][del_pat] = 1
                # logfile.write(del_pat + "\n\n")
        # insert
        for i in range(0, h_gap):
            if h_start_all[i] + h_span_all[i] == h_start_all[i+1]:
                ins_start_ontarget = h_start_all[i+1]
                ins_start_onquery = q_start_all[i] + q_span_all[i]
                ins_end_onquery = q_start_all[i+1] - 1
                if strand_info[0] == "-":
                    ins_start_onquery = q_start_all[i+1] + q_span_all[i+1]
                    ins_end_onquery = q_start_all[i] - 1
                ins_size = ins_end_onquery - ins_start_onquery + 1
                ins_seq = q_seq[ins_start_onquery: ins_start_onquery+ins_size]
                ins_pat = str(ins_start_ontarget) + "I" + \
                    str(ins_size) + ":" + str(ins_seq)
                ASCore_collection[query].addINDEL(ins_pat)
                if ins_pat in PATTERN_COLLECTION[q_besthit[query]]:
                    PATTERN_COLLECTION[q_besthit[query]][ins_pat] += 1
                else:
                    PATTERN_COLLECTION[q_besthit[query]][ins_pat] = 1
                # logfile.write(ins_pat + "\n\n")

    return ASCore_collection, MISMATCH_COLLECTION, PATTERN_COLLECTION, target_assigned_count


def check_conf():
    """ checking the configuration file
    """
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
        exit()


def check_blat():
    """ Checks the user's system to determine which BLAT binary to use.
    """
    if "windows" in user_os.lower():
        DEF_BLAT_PATH = ".\\blat_binaries\\blat.exe"
    elif "linux" in user_os.lower():
        DEF_BLAT_PATH = pwd / "blat_binaries" / "blat_linux"
    elif "darwin" in user_os.lower():
        DEF_BLAT_PATH = pwd / "blat_binaries" / "blat_macos"
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


def load_target(TARGET_FILE):
    # Load target
    v1_targetfile = pwd / TARGET_FILE
    if os.path.exists(v1_targetfile):
        print("Target File is found:" + str(v1_targetfile))
        UTarget = AgeseqIO.readv1TargetFile(v1_targetfile)
        fasta_TARGET_FILE = "AGESeq_run_"+dt_string+".target.fa"
        UTarget.toFastaFile(fasta_TARGET_FILE)
        return UTarget
    else:
        print("Couldn't locate target file:"+str(v1_targetfile))
        exit()


def get_reads_file():
    if os.path.isdir(READS_PATH):
        print("Found Reads Directory: " + str(READS_PATH) + " successfully!")
        import glob
        READS_FILE_LIST = glob.glob(os.path.join(READS_PATH, '*'))

        return READS_FILE_LIST
    else:
        print(str(READS_PATH) + "doesn't exists")
        exit()


HELP = "Running as: AgeseqMain.py -t [target_file] -r [reads_path] -sa [0|1]\n"
if __name__ == '__main__':

    # ageseq.py \
    # -t target.fa \
    # -r /dir/of/reads/ \
    print(HELP)
    main()
