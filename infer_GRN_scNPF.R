# BEELINE数据，scNPF.pro（计算dropout）+rf方法
rm(list=ls())
source('GRACE_function.R')

#一、加载表达数据----
dataset=c('mHSC-L','mHSC-GM','mHSC-E')
spec=c('mouse')
dset=dataset[1]

rdir=paste('Result',dset,'TFRN',sep = '/')
system(paste("mkdir", rdir, sep=" "))

lst.data=load_data(dset,spec)
v.tfs<-lst.data$tf
m.expression<-lst.data$exp[v.tfs,]
v.tgs<-v.tfs  #TF-TF regulatory network
write.table(v.tfs,file = paste(rdir,'tfs.txt',sep='/'),sep = '\t',row.names = FALSE,col.names = FALSE)

#原始rf和xgboost，不加scNPF步骤
#rf
n.cpus <- 2 # number of available cpus for parallel random forest regression
m.rf <- compute_randomforest_based_GRN(m.expression, k="sqrt", nb.trees=1000, set.regulators=v.tfs, set.genes = v.tfs, importance.measure="IncNodePurity", seed=1234, n.cpus=n.cpus)
saveRDS(m.rf, paste(rdir,'m.grn_rf.rds',sep='/'))
link.rf=get_link(m.rf)
saveRDS(link.rf,file = paste(rdir,'link_rf.rds',sep='/'))
metric.rf=cal_AUC3(link.rf,spec,dset)
write.table(metric.rf,file =paste(rdir,'metric_rf.txt',sep='/'),sep='\t',row.names = TRUE,col.names = TRUE)

#rf换成xgboost
m.xgb <- compute_xgboost_based_GRN(m.expression,  nb.trees=500, set.regulators=v.tfs, set.genes = rownames(m.expression))
saveRDS(m.xgb, paste(rdir,'m.grn_xgboost.rds',sep='/'))
link.xgb=get_link(m.xgb)
saveRDS(link.xgb,file = paste(rdir,'link_xgboost.rds',sep='/'))
metric.xgb=cal_AUC3(link.xgb,spec,dset)
write.table(metric.xgb,file =paste(rdir,'metric_xgboost.txt',sep='/'),sep='\t',row.names = TRUE,col.names = TRUE)




# 二、scNPF-pro的多种方式（humannet、STRING、context-mode)----
source('scNPF/scNPF_pro.R')
# 1.humannet
load(file=paste(rdir,'humannet.Rdata',sep='/'))
humannet.pro.data <- scNPF.pro(x=as.matrix(m.expression),network=humannet,nThreads=8)

# 2.STRING-network(human)
load(file =paste(rdir,'string.Rdata',sep='/'))
string.pro.data <- scNPF.pro(x=as.matrix(m.expression),network=string,nThreads=8)

# 3.context-mode
source('scNPF/context_mode.R')
library(WGCNA)
context.pro.data <- scNPF.pro(x=as.matrix(m.expression),network='context',nThreads=8)

#三、推断GRN----
#1.grrf推断GRN
n.cpus <- 2 # number of available cpus for parallel random forest regression
m.grrf <- compute_grrf_based_GRN(m.expression,mat.prior, k="sqrt", nb.trees=1000, set.regulators=v.tfs, set.genes = v.tgs, importance.measure="IncNodePurity", seed=1234, n.cpus=n.cpus)
link.rf=get_link(m.grrf)

#2.rf推断GRN
n.cpus <- 2 # number of available cpus for parallel random forest regression
m.rf.context <- compute_randomforest_based_GRN(context.pro.data, k="sqrt", nb.trees=1000, set.regulators=v.tfs, set.genes = v.tgs, importance.measure="IncNodePurity", seed=1234, n.cpus=2)
saveRDS(m.rf.context, paste(rdir,'m.grn_rf_context.rds',sep='/'))
link.rf.context=get_link(m.rf.context)
saveRDS(link.rf.context,file = paste(rdir,'link_rf_context.rds',sep='/'))

#四、评估----
metric.rf.context=cal_AUC3(link.rf.context,spec,dset)
write.table(metric.rf.context,file =paste(rdir,'metric_rf_context.txt',sep='/'),sep='\t',row.names = TRUE,col.names = TRUE)

#humannet
m.rf.humannet <- compute_randomforest_based_GRN(humannet.pro.data, k="sqrt", nb.trees=1000, set.regulators=v.tfs, set.genes = v.tgs, importance.measure="IncNodePurity", seed=1234, n.cpus=2)
saveRDS(m.rf.humannet, paste(rdir,'m.grn_rf_humannet.rds',sep='/'))
link.rf.humannet=get_link(m.rf.humannet)
saveRDS(link.rf.humannet,file = paste(rdir,'link_rf_humannet.rds',sep='/'))
metric.rf.humannet=cal_AUC3(link.rf.humannet,spec,dset)
write.table(metric.rf.humannet,file =paste(rdir,'metric_rf_humannet.txt',sep='/'),sep='\t',row.names = TRUE,col.names = TRUE)


#STRING
m.rf.string <- compute_randomforest_based_GRN(string.pro.data, k="sqrt", nb.trees=1000, set.regulators=v.tfs, set.genes = v.tgs, importance.measure="IncNodePurity", seed=1234, n.cpus=2)
saveRDS(m.rf.string, paste(rdir,'m.grn_rf_string.rds',sep='/'))
link.rf.string=get_link(m.rf.string)
saveRDS(link.rf.string,file = paste(rdir,'link_rf_string.rds',sep='/'))
metric.rf.string=cal_AUC3(link.rf.string,spec,dset)
write.table(metric.rf.string,file =paste(rdir,'metric_rf_string.txt',sep='/'),sep='\t',row.names = TRUE,col.names = TRUE)

