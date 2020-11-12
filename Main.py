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

MIS_MTCH_WIN = [2, 10]

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
############### end setting env #################


def main():
    util = Util.Utils()
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()

    trgt_full_list = util.read_tsv_ignore_N_line('./input/' + TRGT_FILE)
    chr_arr = ['chrX', 'chrY']
    for i in range(1, 23):
        chr_arr.append('chr' + str(i))

    for chr_nm in chr_arr:
        proc = mp.Process(target=multi_processing, args=(trgt_full_list, chr_nm))
        proc.start()


def multi_processing(trgt_full_list, chr_nm):
    print(chr_nm, 'multi_processing ::::::::::::::::::')
    util = Util.Utils()
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()

    fa_seq, comp_fa_seq = util.read_file_by_biopython(REF_DIR + chr_nm + '.fa', 'fasta')
    len_fa_seq = len(fa_seq)

    result_list = []
    for trgt_inf_arr in trgt_full_list:
        trgt_type = trgt_inf_arr[0]
        len_pbs = int(trgt_inf_arr[5])
        len_rtt = int(trgt_inf_arr[6])
        rule_seq = trgt_inf_arr[7]

        for i in range(len_fa_seq):
            if fa_seq[i] == 'N':
                continue

            fa_p_pam = fa_seq[i: i + LEN_PAM]
            fa_m_pam = comp_fa_seq[i - LEN_PAM: i]

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


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    main()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))