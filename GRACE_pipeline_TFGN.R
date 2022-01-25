#新的TF-TF regulatory network的GRACE流程

rm(list=ls())

source("GRACE_function.R")
cat("\014")  
print("request and load libraries..")
request_libraries()

#数据加载----
dataset=c('mHSC-L','mHSC-GM','mHSC-E')
idx=1
l.samples=c(847,889,1071)
spec=c('mouse')
dset=dataset[idx]
n.samples=l.samples[idx]

# dir=paste('SCODE-dataset',dset,sep = '/')
# dir=getwd()
rdir=paste('Result',dset,'TFRN',sep = '/')
lst.data <- request_dataset() ### load datasets - make sure to configure GRACE_load_datasets.R first
saveRDS(lst.data,file =paste(rdir,'lst.data.rds',sep = '/'))
# lst.data <- readRDS(file =paste(rdir,'lst.data.rds',sep = '/'))
# regulators and targets
v.tfs.global <- lst.data$v.tfs
v.tgs.global <- lst.data$v.tgs

# initial high confidence network
mat.gene_regulatory_network <- as.matrix(lst.data$mat.conserved_motif_plus_expression[v.tfs.global,v.tgs.global])
# mat.gene_regulatory_network <- mapping_intervals(mat.gene_regulatory_network, min(mat.gene_regulatory_network), max(mat.gene_regulatory_network), 0, 1) 

# co-functional network
mat.cofunctional_network <- lst.data$mat.FunctionalAssociation[v.tgs.global,v.tgs.global]

# Training, testing and validation sets

# regulatory evidence based gold standard
mat.regulatory_evidence <- lst.data$mat.GS

# # co-functional evidence
mat.cofunctional_evidence <- lst.data$mat.GO.tgs[v.tgs.global,v.tgs.global]

#****
# lst.validation <- list(lst.data$mat.Chip, lst.data$mat.tissue_localization.tgs, lst.data$mat.HiC)
# names(lst.validation) <- c("n.gs.grn.Chip","n.coreg.localization","n.coreg.HiC")
# code.validation <- c("regulatory_evidence", "cofunctional_evidence", "cofunctional_evidence")

lst.validation <- list(lst.data$mat.chip, lst.data$mat.string)
names(lst.validation) <- c("n.gs.grn.Chip","n.coreg.String")
code.validation <- c("regulatory_evidence", "cofunctional_evidence")

#f-score----
# global constraints
n.cpus <- 2
beta <- 1 # equal emphasis on precision (atrm) and recall (bp coreg pairs) - introducing (novel) combined f-measure (physical + functional evidence, singularity handling, minimum size)
b.normalize_precision <- TRUE
b.jaccard = TRUE

n.sets <- 20
lambda.gridSearch <- c(0.01,seq(0.5,2.5,0.5))
#lambda.gridSearch <- seq(1.5,2,0.5)
#lambda.gridSearch <- seq(1.5,2.5,0.5)


# set.seed(1234)
# model learning related 
th.percentage.hyperparameters = 0.95 # grid search
max.call =  200 # simulated annealing 
n.models <- 100 # 100 models 
n.sample_size <- 0.632 # as in traditional bootrapping approaches

lst.benchmarks <- list(mat.regulatory_evidence, mat.cofunctional_evidence)
names(lst.benchmarks) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")

# lst.benchmarks <- list(mat.regulatory_evidence)
# names(lst.benchmarks) <- c("mat.regulatory_evidence")


if(TRUE){ # needs to be recomputed for first run
  
  # needed for GRACE initiation - using random forest # 
  df.rate_density <- compute_fmeasures_regulatory_network(mat.grn=mat.gene_regulatory_network, 
                                                          lst.benchmarks=lst.benchmarks,
                                                          n.samples = n.samples,
                                                          n.cpus=n.cpus,
                                                          lst.validation=lst.validation,
                                                          code.validation = code.validation)
  
  saveRDS(df.rate_density, paste(rdir,'df.rate_density.rds',sep = '/')) 
  
}else{
  df.rate_density <- readRDS(paste(rdir,'df.rate_density.rds',sep = '/'))
}

lst.fscore_results <- plot_fmeasure_scores(df.rate_density)

# cairo_pdf(paste(rdir,'f-score.pdf',sep = '/'), width = 7.5, height = 6, family= "Helvetica")#"DejaVu Sans")
lst.fscore_results$plot.PR 
# dev.off()

# constants needed globally
ratio.gs.grn_vs_bp.normalization <- lst.fscore_results$ratio.gs.grn_vs_bp.normalization
max.coreg.bp <- lst.fscore_results$max.coreg.bp
df.rate_density <- lst.fscore_results$df.rate_density

### variations ### 

# define limits 
gamma.upper <- alpha.upper <- quantile(mat.gene_regulatory_network[mat.gene_regulatory_network>0], 0.99999999)
#gamma.lower <- alpha.lower <- quantile(mat.gene_regulatory_network[mat.gene_regulatory_network>0], 0.9999)
####gamma.upper <- alpha.upper <- quantile(mat.gene_regulatory_network[mat.gene_regulatory_network>0], 0.9)
gamma.lower <- alpha.lower <- min(mat.gene_regulatory_network[mat.gene_regulatory_network>0]) #estimate_grn_cutoff(mat.grn=mat.grn.cns, v.cut=alpha.lower.cut)

gamma.gridSearch <- as.numeric(sapply(seq(1:n.sets), function(m) gamma.upper - (gamma.upper - gamma.lower)/n.sets * m))

th.preselection <- min(mat.gene_regulatory_network[mat.gene_regulatory_network>0]) #estimate_grn_cutoff(mat.grn=mat.grn.cns,v.cut=v.preselection)

#创建meta GRN----
source("GRACE_function.R")
# compute connectivity matrix per preslection (once)
print("construct meta co-regulation network")
lst.grace_support <- contruct_grace_meta_connectivity_network(mat.grn=mat.gene_regulatory_network, mat.FunctionalAssociation=mat.cofunctional_network, 
                                                              th.preselection = th.preselection, n.cores = n.cpus)

print("extract meta modules")
lst.modules <- prepare_meta_modules(adj.links=lst.grace_support$adj.links, df.grn = lst.grace_support$df.grn)

mat.grn.preselection <- acast(lst.grace_support$df.grn, TF~TG, value.var = "v.grn")
mat.grn.preselection[is.na(mat.grn.preselection)] <- 0
class(mat.grn.preselection) <- "numeric"   

v.tfs <- rownames(mat.grn.preselection)
v.tgs <- colnames(mat.grn.preselection)

mat.cofunctional_network <- mat.cofunctional_network[v.tgs,v.tgs]
mat.cofunctional_evidence <- mat.cofunctional_evidence[v.tgs,v.tgs]

#学习MRF模型参数----
####

tf <- intersect(rownames(mat.grn.preselection),rownames(mat.regulatory_evidence)) 
tgs <- intersect(colnames(mat.grn.preselection), colnames(mat.regulatory_evidence))

truereg.active <- mat.regulatory_evidence
truereg.active <- mat.grn.preselection[tf,tgs] * truereg.active[tf,tgs]
# truereg <- truereg.active
truereg.active[truereg.active < alpha.lower] <- 0
truereg.active[truereg.active > 0] <- 1
idx.tf <- names(which(rowSums(truereg.active) > 0))
idx.tg <- names(which(colSums(truereg.active) > 0))
truereg.active <- truereg.active[idx.tf,idx.tg]

truereg <- Matrix(0, nrow = length(rownames(mat.grn.preselection)), ncol = length(colnames(mat.grn.preselection)),
                  dimnames = list(rownames(mat.grn.preselection), colnames(mat.grn.preselection)))
tfs <- intersect(rownames(truereg.active), rownames(truereg))
tgs <- intersect(colnames(truereg.active), colnames(truereg))
truereg[tfs, tgs] <- truereg.active[tfs, tgs]

idx.gs_grn <- which(truereg == 1)
idx.tg_bp <- which(mat.cofunctional_evidence == 1)

source("GRACE_function.R")
if(TRUE){
  #### multiple rounds of bootstrap aggregation 
  strt<-Sys.time() 
  cl<-makeCluster(min(n.models, n.cpus))
  registerDoParallel(cl)
  
  ### running multiple bootstrap samples in parallel
  lst.grace_models <- foreach(i = 1:n.models, .packages=c("Matrix", "reshape2", "GenSA", "CRF", "igraph")) %dopar% {    
    # source("GRACE_optimization_routines.R")
    # source("GRACE.R")
    # source("GRACE_helperfunctions.R")
    source("GRACE_function.R")
    # use different seeds per run (simulated annealing)
    v.seed = 1234 + 158 * i
    
    lst.grace_models <- evaluate_mode_parameter(lst.grace_support=lst.grace_support, lst.modules=lst.modules, truereg=truereg, mat.cofunctional_evidence=mat.cofunctional_evidence, 
                                                
                                                gamma.gridSearch=gamma.gridSearch,
                                                lambda.gridSearch=lambda.gridSearch,
                                                th.percentage.hyperparameters=th.percentage.hyperparameters,
                                                
                                                ratio.gs.grn_vs_bp.normalization = ratio.gs.grn_vs_bp.normalization,
                                                max.coreg.bp = max.coreg.bp,
                                                
                                                n.sample_size=n.sample_size, max.call=max.call, v.seed = v.seed, beta = beta)
    
    return(lst.grace_models)
    
  }
  
  stopCluster(cl)
  print(Sys.time()-strt)
  
  saveRDS(lst.grace_models,file =paste(rdir,'lst.grace_models.rds',sep = '/'))
  
}else{
  lst.grace_models <- readRDS(file =paste(rdir,'lst.grace_models.rds',sep = '/')) # load precomputed grace models
}  

# lst.grace_models <- readRDS(file =paste( "Result/mHSC-L/TFRN",'lst.grace_models.rds',sep = '/')) # load precomputed grace models
# alpha_avg <- 0
# gamma_avg<- 0
# lambda_avg <- 0
# for (i in 1:n.models) {
#   alpha_avg  <- alpha_avg+lst.grace_models[[i]]$weights[1]
#   gamma_avg  <- gamma_avg+lst.grace_models[[i]]$weights[2]
#   lambda_avg <- lambda_avg+lst.grace_models[[i]]$weights[3]
# }
# alpha_avg=alpha_avg/n.models#0.035828
# gamma_avg=gamma_avg/n.models#0.035828
# lambda_avg=lambda_avg/n.models# 1.26744


lst.benchmarks <- list(mat.regulatory_evidence, mat.cofunctional_evidence)
names(lst.benchmarks) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")


#Bootstrapping集成final GRN----
#验证数据集
print("Computing averaged bootstrap error validation - using same size network")
### compute the bootstrap averaged errors 
validationSet = c("Gene regulatory links (Non-Specific-Chip)", 
                  "Gene regulatory links (ChIP-seq-network)", 
                  "Co-functional pairs (MouseNetv2_GS_symbol)",
                  "Protein-protein interaction (STRING)") 

v.sets_map <-  c('n.gs.grn','n.gs.grn.Chip','n.coreg.bp_tg_pairs','n.coreg.String')#df.rate_density的(2,7,3,8)列
names(v.sets_map) <- validationSet

methods_name<- c("RF_row", "xgboost",'SCODE')
l.foldchanges <- vector(mode = "list", length = length(methods_name))
for(m in 1:length(methods_name)){
  l.foldchanges[[m]] <- vector(mode = "list", length = 2)
  names(l.foldchanges[[m]]) <- c("GRACE", "Original")
  for(j in 1:2){
    l.foldchanges[[m]][[j]] <- vector(mode = "list", length = length(validationSet))
    for(s in 1:length(validationSet)){
      l.foldchanges[[m]][[j]][[s]] <- numeric(n.models)
    }
  }
}
#初始化的结构相同
l.foldchanges_comparative <-l.pvalues_comparative <-  l.pvalues <-l.foldchanges
names(l.foldchanges) <- methods_name
names(l.pvalues) <- methods_name
names(l.pvalues_comparative) <- methods_name
names(l.foldchanges_comparative) <- methods_name


v.background = df.rate_density[nrow(df.rate_density),]
n.pred.background <- df.rate_density$n.pred[nrow(df.rate_density)]
n.coreg.background <- df.rate_density$n.coreg.tg_pairs[nrow(df.rate_density)]
v.backgroundSet <- as.numeric(df.rate_density[nrow(df.rate_density),v.sets_map]) # needs to be adjusted
names(v.backgroundSet) <- validationSet

# generate training and test bootstrap prediction errors

#读取各个对比方法推断的GRN
# GENIE3, GGM, CLR, BC3NET, wGLASSO, iRafNet
l.methods <- vector(mode = "list", length = length(methods_name))
names(l.methods) <- methods_name
l.methods[[1]] <- readRDS(paste(rdir,'m.grn_rf.rds',sep = '/'))[v.tfs,v.tgs]
l.methods[[2]] <- readRDS(paste(rdir,'m.grn_xgboost.rds',sep = '/'))[v.tfs,v.tgs]

SCODE<- read.table(paste('Result',dset,'SCODE','A.txt',sep = '/'), header=TRUE, sep="\t")
SCODE <- as.matrix(SCODE)
SCODE <- mapping_intervals(SCODE, min(SCODE), max(SCODE), 0, 1) #标准化到0-1
l.methods[[3]] <- SCODE
# l.methods[[3]] <- l.methods[[4]] <- l.methods[[1]]
# l.methods[[2]] <- readRDS("drosophila_datasets/mat.grn.clr_drosophila.rds")[v.tfs,v.tgs]
# #l.methods[[3]] <- readRDS("drosophila_datasets/mat.grn.aracne_drosophila.rds")[v.tfs,v.tgs]
# # l.methods[[4]] <- readRDS("bc3net_drosophila.rds")[v.tfs,v.tgs]
# l.methods[[3]] <- readRDS("drosophila_datasets/ggm_drosophila.rds")[v.tfs,v.tgs]
# #l.methods[[5]] <- readRDS("glasso_drosophila.rds")
# l.methods[[4]] <- abs(readRDS("drosophila_datasets/wGlasso_drosophila.rds")[v.tfs,v.tgs])

#rownames(l.methods[[5]]$wi) <- colnames(l.methods[[5]]$wi) <- v.gns
#saveRDS(l.methods[[5]]$wi[v.tfs,v.tgs], "drosophila_datasets/wGlasso_drosophila.rds")
# l.methods[[5]]

# readRDS("iRafNet_drosophila")

#背景
n.links_background <- length(v.tfs) * length(v.tgs)
n.co_pairs <- length(v.tgs) * length(v.tgs)

n.gs_background <- sum(mat.regulatory_evidence[intersect(rownames(mat.regulatory_evidence), v.tfs), intersect(colnames(mat.regulatory_evidence), v.tgs)]) # 386 
n.chip_background <- sum(lst.validation$n.gs.grn.Chip[intersect(rownames(lst.validation$n.gs.grn.Chip), v.tfs), intersect(colnames(lst.validation$n.gs.grn.Chip), v.tgs)])

v.gs.grn <- c(n.gs_background, n.chip_background)
v.gs.copairs <- c(sum(mat.cofunctional_evidence[v.tgs,v.tgs]), sum(lst.validation$n.coreg.String[v.tgs,v.tgs]))

####
# t区分训练集和测试集，i遍历n.models个模型，m遍历各个对比方法，j区分GRACE和原始方法，s遍历各个验证集（benchmark+validation）
for(t in 1:2){
  
  if(t == 1){ #b.validation=FALSE训练集
    b.validation <- FALSE
  }else{
    b.validation <- TRUE
  }
  
  for(i in 1:n.models){
    
    cat("Processing... ", round(i/n.models * 100, digits = 2) , "%", "\r"); flush.console() 
    
    # use different seeds per run (simulated annealing)
    v.seed = 1234 + 158 * i
    set.seed(v.seed)
    
    #### Bootstrapping (based on regulators based meta module sampling (0.632))
    tfs.mrf <- unique(lst.grace_support$df.grn$TF) 
    tgs.mrf <- unique(lst.grace_support$df.grn$TG)
    
    lst.bootstrap_sample <- generate_bootstrapping_sample(lst.grace_support=lst.grace_support,lst.modules=lst.modules, truereg=as.matrix(truereg), mat.cofunctional_evidence=mat.cofunctional_evidence, n.sample_size = n.sample_size, b.validation = b.validation)
    lst.grace_support.bootstrap_sample <- lst.bootstrap_sample$lst.grace_support.bootstrap_sample
    lst.modules.bootstrap_sample <- lst.bootstrap_sample$lst.modules.bootstrap_sample
    
    lst.benchmarks <- list(lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.GS.bootstrap_sample, lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.cofunctional_evidence.bootstrap_sample)    #list(mat.GS.train, mat.BP.coreg.train)
    names(lst.benchmarks) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")
    
    lst.validation <- list(lst.data$mat.chip,  lst.data$mat.string)
    names(lst.validation) <- c("n.gs.grn.Chip","n.coreg.String")
    code.validation <- c("regulatory_evidence", "cofunctional_evidence")
    
    ####
    
    alpha  <- lst.grace_models[[i]]$weights[1]
    gamma  <- lst.grace_models[[i]]$weights[2]
    lambda <- lst.grace_models[[i]]$weights[3]
    lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields(lst.grace_support=lst.grace_support, lst.modules=lst.modules, 
                                                              v.alpha = alpha, v.gamma=gamma, v.lambda=lambda, b.prim=TRUE) 
    mat.p.crf <- lst.mat.p.crf$mat.p.crf
    
    ### 
    
    l.eval <- vector(mode = "list", length = 2)
    #l.eval[[1]]表示GRACE，l.eval[[2]]表示对比方法,
    # 因此l.foldchanges和l.pvalues的[[m]][[1]]都是相等的，且表示GRACE，[[m]][[2]]即original表示对比方法的结果
    l.eval[[1]] <- compute_grace_complete_assessment(mat.grace=mat.p.crf, 
                                                     lst.benchmarks=lst.benchmarks,
                                                     lst.validation=lst.validation,
                                                     code.validation = code.validation)
    
    
    n.links <- sum(mat.p.crf)
    
    #n.links_background <- length(v.tfs) * length(v.tgs)
    #n.gs_background <- sum(mat.regulatory_evidence)
    #n.co_pairs <- length(v.tgs) * length(v.tgs)
    # all methods as comparison models here #
    for(m in 1:length(l.methods)){
      
      mat.gene_regulatory_network <- as.matrix(l.methods[[m]])
      
      mat.grn.comparison <- mat.gene_regulatory_network[v.tfs,v.tgs]
      # mat.grn.comparison <- mat.grn.comparison[rownames(lst.benchmarks$mat.regulatory_evidence),]
      
      df.grn.comparison <- as.data.frame(as.table(mat.grn.comparison), stringsAsFactors = FALSE) 
      names(df.grn.comparison)[1:3] <- c("TF","TG","v.grn")
      
      df.grn.comparison <- subset(df.grn.comparison, df.grn.comparison$v.grn > 0)
      
      df.grn.comparison <- df.grn.comparison[order(-df.grn.comparison$v.grn),] 
      df.grn.comparison <- df.grn.comparison[1:n.links,]#取相同数量的top link
      
      mat.grn.comparison <- acast(df.grn.comparison, TF~TG, value.var="v.grn")
      mat.grn.comparison[is.na(mat.grn.comparison)] <- 0; 
      mat.grn.comparison[mat.grn.comparison > 0] <- 1
      
      l.eval[[2]] <- compute_grace_complete_assessment(mat.grace=mat.grn.comparison, 
                                                       lst.benchmarks=lst.benchmarks,
                                                       lst.validation=lst.validation,
                                                       code.validation = code.validation)
      
      for(j in 1:2){
        
        # print(j)
        for(s in 1:length(v.sets_map)){
          
          hitInSample <- l.eval[[j]][[v.sets_map[s]]]
          
          if(s <= 2){
            hitInPop <- v.gs.grn[s]
            sampleSize <- n.links
            failInPop <- (n.links_background - hitInPop)
            foldchange <- (hitInSample/sampleSize)/(hitInPop/n.links_background)
          }else{
            hitInPop <- v.gs.copairs[s - 2]
            sampleSize <- l.eval[[j]]$n.coreg.tg_pairs
            failInPop <- (n.co_pairs - hitInPop)
            foldchange <- (hitInSample/sampleSize)/(hitInPop/n.co_pairs)
          }
          
          #print(foldchange)
          
          p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
          #l.foldchanges表示FC=GRACE或对比方法/大背景
          l.foldchanges[[m]][[j]][[s]][[i]] <- foldchange
          l.pvalues[[m]][[j]][[s]][[i]] <- p.val
          
        }
        #print("")
      }
      
      ### comparative enrichment 
      
      for(s in 1:length(v.sets_map)){
        
        hitInSample <- l.eval[[1]][[v.sets_map[s]]]
        hitInPop <- l.eval[[2]][[v.sets_map[s]]]
        
        if(s <= 2){
          sampleSize <- n.links
          failInPop <- (n.links - hitInPop)
          foldchange <- (hitInSample/n.links)/(hitInPop/n.links)
        }else{
          sampleSize <- l.eval[[1]]$n.coreg.tg_pairs
          failInPop <- (l.eval[[2]]$n.coreg.tg_pairs - hitInPop)
          foldchange <- (hitInSample/sampleSize)/(hitInPop/l.eval[[2]]$n.coreg.tg_pairs)
        }
        
        tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, failInPop), nrow = 2)
        p.val <- fisher.test(tab)$p.value
        
        # l.foldchanges_comparative表示FC=GRACE/对比方法
        l.foldchanges_comparative[[m]][[j]][[s]][[i]] <- foldchange
        l.pvalues_comparative[[m]][[j]][[s]][[i]] <- p.val
        
      }
    }
  }
  
  if(t == 1){
    l.foldchanges.training <- l.foldchanges
    l.pvalues.training <- l.pvalues
    l.foldchanges_comparative.training <- l.foldchanges_comparative
    l.pvalues_comparative.training <- l.pvalues_comparative
  }else{
    l.foldchanges.testing <- l.foldchanges
    l.pvalues.testing <- l.pvalues
    l.foldchanges_comparative.testing <- l.foldchanges_comparative
    l.pvalues_comparative.testing <- l.pvalues_comparative
  }
  
  # compute the actual values instead of fold changes
  # individual comparative enrichment tests
  
}
### 

#GRACE方法在4个评估集下的结果（m可以为任何，因为在几种对比方法下的GRACE都是相同的）
print(lapply(l.foldchanges.testing[[m]][[1]], mean))
print(lapply(l.pvalues.testing[[m]][[1]], mean))
# print(lapply(l.foldchanges.testing[[m]][[2]], mean))
print(lapply(l.foldchanges.training[[m]][[1]], mean))

#对比方法在4个评估集下的结果
for(m in 1:length(methods_name)){
  
  print(names(l.methods)[m])
  
  # print(lapply(l.foldchanges_comparative.testing[[m]][[2]], mean))
  # print(lapply(l.pvalues_comparative.testing[[m]][[2]], mean))
  
  # print(lapply(l.foldchanges.training[[m]][[2]], mean))
  
  # print(lapply(l.foldchanges.testing[[m]][[2]], mean))
  print(lapply(l.pvalues.testing[[m]][[2]], mean))
  
}

#GRACE比上对比方法在4个评估集下的结果
for(m in 1:length(methods_name)){
  
  print(names(l.methods)[m])
  
  print(lapply(l.foldchanges_comparative.testing[[m]][[2]], mean))
  print('p-values')
  print(lapply(l.pvalues_comparative.testing[[m]][[2]], mean))
  
}

saveRDS(l.foldchanges,file =paste(rdir,'l.foldchanges.rds',sep = '/'))
saveRDS(l.foldchanges.training,file =paste(rdir,'l.foldchanges.training.rds',sep = '/'))
saveRDS(l.foldchanges.testing,file =paste(rdir,'l.foldchanges.testing.rds',sep = '/'))
saveRDS(l.foldchanges_comparative.training,file =paste(rdir,'l.foldchanges_comparative.training.rds',sep = '/'))
saveRDS(l.foldchanges_comparative.testing,file =paste(rdir,'l.foldchanges_comparative.testing.rds',sep = '/'))

saveRDS(l.pvalues,file =paste(rdir,'l.pvalues.rds',sep = '/'))
saveRDS(l.pvalues.training,file =paste(rdir,'l.pvalues.training.rds',sep = '/'))
saveRDS(l.pvalues.testing,file =paste(rdir,'l.pvalues.testing.rds',sep = '/'))
saveRDS(l.pvalues_comparative.training,file =paste(rdir,'l.pvalues_comparative.training.rds',sep = '/'))
saveRDS(l.pvalues_comparative.testing,file =paste(rdir,'l.pvalues_comparative.testing.rds',sep = '/'))

# dset= "mHSC-L"
# rdir=paste('Result',dset,'TFRN',sep = '/')
# l.foldchanges.testing <- readRDS(file=paste(rdir,'l.foldchanges.testing.rds',sep = '/'))
# l.pvalues.testing <- readRDS(file =paste(rdir,'l.pvalues.testing.rds',sep = '/'))
# l.foldchanges_comparative.testing <- readRDS(file =paste(rdir,'l.foldchanges_comparative.testing.rds',sep = '/'))
# l.pvalues_comparative.testing <- readRDS(file =paste(rdir,'l.pvalues_comparative.testing.rds',sep = '/'))
# 
