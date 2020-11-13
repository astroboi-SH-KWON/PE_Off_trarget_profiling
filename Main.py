import time
import os
import multiprocessing as mp
import platform

import Util
import Logic
import LogicPrep
############### start to set env ################
WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
SYSTEM_NM = platform.system()

if SYSTEM_NM == 'Linux':
    # REAL
    REF_DIR = "../reference/hg38/"
else:
    # DEV
    REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"

TRGT_FILE = "trgt_list.txt"

PAM = ['NGG', 'NAG', 'NGA']
LEN_PAM = 3

LEN_RTT_f_PAM = 3
LEN_f_PAM = 20

MIS_MTCH_WIN = [0, 8]

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
############### end setting env #################


def main(st_time):
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    trgt_full_list = util.read_tsv_ignore_N_line('./input/' + TRGT_FILE)

    chr_arr = ['chrX', 'chrY']
    for i in range(1, 23):
        chr_arr.append('chr' + str(i))

    manager = mp.Manager()
    result_list = manager.list()
    jobs = []
    for chr_nm in chr_arr:
        # proc = mp.Process(target=count_mismatch, args=(trgt_full_list, chr_nm, st_time))
        proc = mp.Process(target=count_mismatch_by_guide_PAM_RTTonly, args=(result_list, trgt_full_list, chr_nm, st_time))
        jobs.append(proc)
        proc.start()

    for ret_proc in jobs:
        ret_proc.join()

    filtered_result_list = logic_prep.filter_out_N_seq(result_list)
    # clear ListProxy object
    result_list[:] = []
    srted_result_list = logic_prep.sort_list_by_ele(filtered_result_list, 0)
    filtered_result_list.clear()

    header = ['Type', 'Chr', 'Location', 'Off-target sequence', 'mismatch_tot_cnt', 'mis_mtch_guide_cnt(1~20)',
              'mis_mtch_PAM_cnt(21~23)', 'mis_mtch_RTTonly_cnt(24~)', 'strand']
    try:
        util.make_excel('./output/tot_result', header, srted_result_list)
    except Exception as err:
        util.make_tsv('./output/tot_result', header, srted_result_list)

    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - st_time))


def count_mismatch_by_guide_PAM_RTTonly(result_list, trgt_full_list, chr_nm, st_time):
    print(chr_nm, 'count_mismatch_by_guide_PAM_RTTonly ::::::::::::::::::')
    util = Util.Utils()
    logic = Logic.Logics()

    fa_seq, comp_fa_seq = util.read_file_by_biopython(REF_DIR + chr_nm + '.fa', 'fasta')
    len_fa_seq = len(fa_seq)

    for i in range(len_fa_seq):
        if fa_seq[i] == 'N':
            continue

        fa_p_pam = fa_seq[i: i + LEN_PAM]
        fa_m_pam = comp_fa_seq[i - LEN_PAM: i]

        for trgt_inf_arr in trgt_full_list:
            trgt_type = trgt_inf_arr[0]
            len_pbs = int(trgt_inf_arr[5])
            len_rtt = int(trgt_inf_arr[6])
            rule_seq = trgt_inf_arr[7]

            if logic.is_fa_pam_in_rule(fa_p_pam, PAM):
                p_fa_seq = fa_seq[i - LEN_f_PAM: i + len_rtt - LEN_RTT_f_PAM]
                mis_mtch_guide_cnt = logic.cnt_mismatch(p_fa_seq[:LEN_f_PAM], rule_seq[:LEN_f_PAM])
                mis_mtch_PAM_cnt = logic.cnt_mismatch(p_fa_seq[LEN_f_PAM: LEN_f_PAM + LEN_PAM],
                                                      rule_seq[LEN_f_PAM: LEN_f_PAM + LEN_PAM])
                mis_mtch_RTTonly_cnt = logic.cnt_mismatch(p_fa_seq[LEN_f_PAM + LEN_PAM:],
                                                          rule_seq[LEN_f_PAM + LEN_PAM:])
                mis_mtch_cnt = mis_mtch_guide_cnt + mis_mtch_PAM_cnt + mis_mtch_RTTonly_cnt
                if MIS_MTCH_WIN[0] <= mis_mtch_cnt <= MIS_MTCH_WIN[1]:
                    result_list.append(
                        [trgt_type, chr_nm, i - LEN_f_PAM + 1, p_fa_seq, mis_mtch_cnt, mis_mtch_guide_cnt,
                         mis_mtch_PAM_cnt, mis_mtch_RTTonly_cnt, '+'])

            if logic.is_fa_pam_in_rule(fa_m_pam[::-1], PAM):
                m_fa_seq = comp_fa_seq[i - len_rtt + LEN_RTT_f_PAM: i + LEN_f_PAM][::-1]
                mis_mtch_guide_cnt = logic.cnt_mismatch(m_fa_seq[:LEN_f_PAM], rule_seq[:LEN_f_PAM])
                mis_mtch_PAM_cnt = logic.cnt_mismatch(m_fa_seq[LEN_f_PAM: LEN_f_PAM + LEN_PAM],
                                                      rule_seq[LEN_f_PAM: LEN_f_PAM + LEN_PAM])
                mis_mtch_RTTonly_cnt = logic.cnt_mismatch(m_fa_seq[LEN_f_PAM + LEN_PAM:],
                                                          rule_seq[LEN_f_PAM + LEN_PAM:])
                mis_mtch_cnt = mis_mtch_guide_cnt + mis_mtch_PAM_cnt + mis_mtch_RTTonly_cnt

                if MIS_MTCH_WIN[0] <= mis_mtch_cnt <= MIS_MTCH_WIN[1]:
                    result_list.append(
                        [trgt_type, chr_nm, len_fa_seq - i - LEN_f_PAM, m_fa_seq, mis_mtch_cnt, mis_mtch_guide_cnt,
                         mis_mtch_PAM_cnt, mis_mtch_RTTonly_cnt, '-'])

    print(chr_nm, "count_mismatch_by_guide_PAM_RTTonly DONE ::: %.2f seconds :::" % (time.perf_counter() - st_time))


def count_mismatch(trgt_full_list, chr_nm, st_time):
    print(chr_nm, 'count_mismatch ::::::::::::::::::')
    util = Util.Utils()
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()

    fa_seq, comp_fa_seq = util.read_file_by_biopython(REF_DIR + chr_nm + '.fa', 'fasta')
    len_fa_seq = len(fa_seq)

    result_list = []
    for i in range(len_fa_seq):
        if fa_seq[i] == 'N':
            continue

        fa_p_pam = fa_seq[i: i + LEN_PAM]
        fa_m_pam = comp_fa_seq[i - LEN_PAM: i]

        for trgt_inf_arr in trgt_full_list:
            trgt_type = trgt_inf_arr[0]
            len_pbs = int(trgt_inf_arr[5])
            len_rtt = int(trgt_inf_arr[6])
            rule_seq = trgt_inf_arr[7]

            if logic.is_fa_pam_in_rule(fa_p_pam, PAM):
                p_fa_seq = fa_seq[i - LEN_f_PAM: i + len_rtt - LEN_RTT_f_PAM]
                mis_mtch_cnt = logic.cnt_mismatch(p_fa_seq, rule_seq)
                if MIS_MTCH_WIN[0] <= mis_mtch_cnt <= MIS_MTCH_WIN[1]:
                    result_list.append([trgt_type, chr_nm, i - LEN_f_PAM + 1, p_fa_seq, mis_mtch_cnt, '+'])

            if logic.is_fa_pam_in_rule(fa_m_pam[::-1], PAM):
                m_fa_seq = comp_fa_seq[i - len_rtt + LEN_RTT_f_PAM: i + LEN_f_PAM][::-1]
                mis_mtch_cnt = logic.cnt_mismatch(m_fa_seq, rule_seq)
                if MIS_MTCH_WIN[0] <= mis_mtch_cnt <= MIS_MTCH_WIN[1]:
                    result_list.append([trgt_type, chr_nm, len_fa_seq - i - LEN_f_PAM, m_fa_seq, mis_mtch_cnt, '-'])

    filtered_result_list = logic_prep.filter_out_N_seq(result_list)
    result_list.clear()
    srted_result_list = logic_prep.sort_list_by_ele(filtered_result_list, 0)
    filtered_result_list.clear()

    header = ['Type', 'Chr', 'Location', 'Off-target sequence', 'mismatch_cnt', 'strand']
    try:
        util.make_excel('./output/' + chr_nm, header, srted_result_list)
    except Exception as err:
        util.make_tsv('./output/' + chr_nm, header, srted_result_list)

    print(chr_nm, "count_mismatch DONE ::: %.2f seconds :::" % (time.perf_counter() - st_time))


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    main(start_time)
