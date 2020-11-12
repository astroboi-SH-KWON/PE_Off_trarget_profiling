
from astroboi_bio_tools.ToolLogic import ToolLogics
class Logics(ToolLogics):
    def is_fa_pam_in_rule(self, fa_p_pam, pam_arr):
        for pam_rule in pam_arr:
            if self.match(0, fa_p_pam, pam_rule):
                return True
        return False
    """
    same length
    """
    def cnt_mismatch(self, trgt_seq, rule_seq):
        cnt = 0
        for i in range(len(rule_seq)):
            if not self.checkSeqByChar(trgt_seq[i], rule_seq[i]):
                cnt += 1

        return cnt


