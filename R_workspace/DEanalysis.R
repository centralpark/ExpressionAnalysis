args <- commandArgs(trailingOnly = T)
gct_file = args[1]
col_file = args[2]
data_file = args[3]
result_file = args[4]

source('/home/staff/hes14/Roche/R_workspace/DifferentialExpression.R')

sanity.check(gct_file,col_file)
res <- differentialExpressionAnalysis(gct_file,col_file,save_data_file=data_file,
                                      result_file=result_file)