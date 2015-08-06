source('/Users/HSH/Roche/R_workspace/DifferentialExpression.R')

args <- commandArgs(trailingOnly = T)
inputFile <- args[1]
outputFileName <- args[2]

load(inputFile)
pathwayAnalysis(res, result_file=outputFileName)
