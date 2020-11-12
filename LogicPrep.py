
from astroboi_bio_tools.ToolLogicPrep import ToolLogicPreps
class LogicPreps(ToolLogicPreps):
    def filter_out_N_seq(self, result_list):
        return [tmp for tmp in result_list if 'N' not in tmp[3]]