rm(list = ls())

library(RGCCA)
library(ggplot2)
library(plotly)


# functions
# ----------------------------------
loadData <- function(targFile, targSep = '\t', 
                     dataFile, dataSep = '\t'){
    targ = read.table(targFile, header = TRUE, sep = targSep, row.names = 1)
    data = as.matrix(read.table(dataFile, header = TRUE, sep = dataSep, row.names = 1, check.names = FALSE))
    mids = match(row.names(targ), colnames(data))
    data = t(data[,mids[!is.na(mids)]])
    data = data[,apply(data, 2, var) > 1e-10]
    return(list(target = targ, data = data))
}

factorPlot <- function(y1, y2, conds){
    df = data.frame(comp1 = y1, comp2 = y2, colFactor = conds)
    p <- ggplot(df, aes(comp1, comp2)) + 
         geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
         ggtitle('Factor plot') + 
         geom_text(aes(colour = colFactor, label = rownames(df)), vjust = 0,nudge_y = 0.03, size = 3) + 
         theme(legend.position = 'bottom', legend.box = 'horizontal', legend.title = element_blank())
    return(p)
}
# ----------------------------------


# load rna-seq data
res = loadData(
    '../../data/sample_sheet.csv', targSep = ',',
    '../../data/rna_vst_counts.csv', dataSep = ','
)
ttarg = res$target
tdata = res$data

# load metabolite data
res = loadData(
    '../../data/sample_sheet.csv', targSep = ',',
    '../../data/polar_pos_pqn_imputed_glog.csv', dataSep = ','
)
ptarg = res$target
pdata = res$data


# take samples available on both data
sample.ids = intersect(row.names(tdata), row.names(pdata))
tdata = tdata[sample.ids,]
pdata = pdata[sample.ids,]
A = list(tdata, pdata)


# apply RGCCA/SGCCA
A = lapply(A, function(x) scale2(x, bias = TRUE))
C = matrix(c(0, 1, 1, 0), 2, 2)

cca.res = rgcca(A, C, tau = c(1,1), ncomp = c(2,2), scale = FALSE, verbose = FALSE)
cca.res$AVE

cca.res = sgcca(A, C, c1 = c(0.2,0.2), ncomp = c(2,2), scale = FALSE, verbose = FALSE)
cca.res$AVE
# number of features selected
colSums(cca.res$a[[1]] != 0)
colSums(cca.res$a[[2]] != 0)


# CCA plots
p = factorPlot(cca.res$Y[[1]][,1], cca.res$Y[[2]][,1], ttarg[row.names(tdata),]$REF)
ggplotly(p)

p = factorPlot(cca.res$Y[[1]][,2], cca.res$Y[[2]][,2], ttarg[row.names(tdata),]$Site)
ggplotly(p)
