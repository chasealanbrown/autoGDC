import rpy3.robjects as ro
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.conversion import localconverter

pandas2ri.activate()

ro.r("source('DE_analysis.R')")
deseq2_r = ro.globalenv["deseq2_r"]

def pydeseq(counts: pd.DataFrame,
            design: pd.DataFrame,
            formula: str):

    with localconverter(ro.default_converter + pandas2ri.converter):
        formula_r = Formula(formula)
        counts_r = ro.conversion.py2rpy(counts)
        design_r = ro.conversion.py2rpy((design)

        deg_r = deseq2_r(counts_df = counts_r,
                         design_matrix = design_r,
                         design_formula = formula_r)

        deg_df =  ro.conversion.rpy2py(deg_r).set_index("gene")
    return deg_df
