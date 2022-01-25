
# 加载BEELINE数据
load_data<-function(dset,spec){
  # browser()
  #read data
  exp=read.table(file =paste('scRNA-Seq',dset,'ExpressionData.csv',sep='/') ,sep = ',',header = T ,stringsAsFactors = FALSE)
  rownames(exp)=exp$X
  exp=subset(exp,select=-X)
  
  ord=read.table(file =paste('scRNA-Seq',dset,'GeneOrdering.csv',sep='/') ,sep = ',',header = T ,stringsAsFactors = FALSE)
  pst=read.table(file =paste('scRNA-Seq',dset,'PseudoTime.csv',sep='/') ,sep = ',',header = T ,stringsAsFactors = FALSE)
  spe_tf=read.csv(file = paste('BEELINE-Networks/',spec,'-tfs.csv',sep = ''),header = T,stringsAsFactors = FALSE)
  
  #filter genes
  ord$p_adjust=p.adjust(p=ord$VGAMpValue,method = 'bonferroni')
  tgs=ord$X[which(ord$p_adjust<0.01)]
  tfs=intersect(tgs,spe_tf$TF)
  
  #sort cells
  pst=pst[order(pst$PseudoTime),]
  exp=exp[tgs,pst$X]
  
  return(list(exp=exp,pst=pst,tfs=tfs,tgs=tgs))
}

#评估相关函数
get_link<-function(W){
  library(reshape2)
  #W_tg*tf
  link = melt(t(W), na.rm=TRUE)
  colnames(link) <- c("regulator", "target", "weight")
  link <- link[order(link$weight, decreasing=TRUE),]
  link$weight<-round(link$weight,4)
  return(link)
}
#三个评估网络整合为一个gold standard network
cal_AUC3<-function(link,species,dset){
  library(ROCR)
  library(modEvA)
  
  #Non-specific-chip-seq-network
  nonspe_chip_net=read.table(file = paste('BEELINE-Networks/Networks',species,'Non-Specific-Chip-seq-network.csv',sep='/'),
                             stringsAsFactors = FALSE,header = TRUE,sep=',',col.names = c('regulator','target'))
  nonspe_chip_net$nonSpecific=rep(1,dim(nonspe_chip_net)[1])
  #STRING-NETWORK
  STRING_net=read.table(file = paste('BEELINE-Networks/Networks',species,'STRING-network.csv',sep='/'),
                        stringsAsFactors = FALSE,header = TRUE,sep=',',col.names = c('regulator','target'))
  STRING_net$String=rep(1,dim(STRING_net)[1])
  # browser()
  #CHIP-SEQ-NETWORK
  chip_dset=as.character(unlist(strsplit(dset, split = "-")))[1]
  
  chip_net=read.table(file = paste('BEELINE-Networks/Networks',species,paste(chip_dset,'ChIP-seq-network.csv',sep='-'),sep='/'),
                      stringsAsFactors = FALSE,header = TRUE,sep=',')
  if(ncol(chip_net)==2){
    colnames(chip_net) = c('regulator','target')
  }else if(ncol(chip_net)==3){
    colnames(chip_net) = c('regulator','target','chip_score')
  }
  chip_net$Specific=rep(1,dim(chip_net)[1])
  
  score=merge(link,nonspe_chip_net,by=c('regulator','target'),all.x = TRUE)
  score=merge(score,chip_net,by=c('regulator','target'),all.x = TRUE)
  score=merge(score,STRING_net,by=c('regulator','target'),all.x = TRUE)
  
  score$nonSpecific[is.na(score$nonSpecific)]=0
  score$Specific[is.na(score$Specific)]=0
  score$String[is.na(score$String)]=0
  score$Total= score$nonSpecific|score$Specific|score$String
  score$Total=ifelse(score$Total == "TRUE", 1, 0)
  
  eval_metric=data.frame(AUC=rep(0.0,4),AUPR=rep(0.0,4))
  rownames(eval_metric)=c('STRING-network','Non-Specific-Chip',paste(dset,'Chip',sep='-'),'comprehensive')
  
  eval_metric[1,'AUC']<-roc_curve(obs=score$String,pred=score$weight)
  eval_metric[1,'AUPR']<-AUC(obs=score$String,pred=score$weight,curve = "PR", simplif=TRUE, main = "PR curve")
  
  eval_metric[2,'AUC']<-roc_curve(obs=score$nonSpecific,pred=score$weight)
  eval_metric[2,'AUPR']<-AUC(obs=score$nonSpecific,pred=score$weight,curve = "PR", simplif=TRUE, main = "PR curve")
  
  eval_metric[3,'AUC']<-roc_curve(obs=score$Specific,pred=score$weight)
  eval_metric[3,'AUPR']<-AUC(obs=score$Specific,pred=score$weight,curve = "PR", simplif=TRUE, main = "PR curve")
  
  eval_metric[4,'AUC']<-roc_curve(obs=score$Total,pred=score$weight)
  eval_metric[4,'AUPR']<-AUC(obs=score$Total,pred=score$weight,curve = "PR", simplif=TRUE, main = "PR curve")
  
  return(eval_metric)
}

#绘制ROC曲线,返回AUC
roc_curve<-  function(pred, obs) {
  predict<-prediction(pred,obs)
  curve<-performance(predict,"tpr","fpr")
  auc<-performance(predict,"auc")@y.values[[1]]
  # png(filename =file )
  plot(curve,main=paste("ROC Curve\nAUC=",auc),xlab="False Positive Rate",ylab="True Positive Rate",lwd=3)
  # dev.off()
  return(auc)
}

compute_randomforest_based_GRN <- function(mat.expression=mat.expression, k="sqrt", nb.trees=1000, set.regulators = NULL, set.genes = NULL, seed=1234, importance.measure = "impurity", n.cpus = 2){
  
  library(ranger)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  mat.expression.norm <- t(mat.expression)
  mat.expression.norm <- apply(mat.expression.norm, 2, function(x) { (x - mean(x)) / sd(x) } )
  n.genes <- dim(mat.expression.norm)[2]
  genes<- colnames(mat.expression.norm)
  n.samples <- dim(mat.expression.norm)[1]
  
  if(is.null(set.genes)){
    n.genes <- n.genes
    genes <- genes
    
  }else{
    n.genes <- length(set.genes)
    genes <- set.genes
    
  }
  
  #mat.expression.norm <- mat.expression.norm[,genes]
  
  if (is.null(set.regulators)) {
    n.regulators <- n.genes
    regulators <- genes
  } else {
    n.regulators <- length(set.regulators)
    # regulators provided by names
    if(is.character(set.regulators)){
      regulators <- set.regulators
      genes.undefined <- setdiff(regulators, genes)
      if (length(genes.undefined) != 0) {
        stop(paste("Error: genes: ", paste(genes.undefined, collapse = ",")," not represented in gene expression matrix \n", sep=""))
      }
      # regulators provided by indices
    }else if (is.numeric(set.regulators)) {
      regulators <- genes[set.regulators]
    }else{
      stop("Error: invalid regulator format")
    }
  }
  # set mtry
  if (class(k) == "numeric") {
    mtry <- K
  } else if (k == "sqrt") {
    mtry <- round(sqrt(n.regulators))
  } else if (k == "all") {
    mtry <- n.regulators-1
  } else {
    stop("Error: invalid parameter k, options: \"sqrt\", \"all\", or integer")
  }
  
  print(paste("Performing random-forest regression based gene regulatory network inference (# of decision trees per gene is: ", nb.trees, ", # of regulators per decision tree node: ", mtry, sep=""))
  mat.rapid <- matrix(0.0, nrow=n.regulators, ncol=n.genes, dimnames = list(regulators, genes))
  # browser()
  strt<-Sys.time()
  for(i in 1:n.genes){
    cat("Processing... ", round(i/n.genes * 100, digits = 2) , "%", "\r"); flush.console() 
    gn.i <- genes[i]
    # remove target gene from set of regulators
    regulators.i <- setdiff(regulators, gn.i)
    x <- mat.expression.norm[,regulators.i] # avoid coercion to numeric
    y <- mat.expression.norm[,gn.i]
    # df.xy <- cbind(as.data.frame(x),y)
    
    rf.model <- ranger(x=x,y=y, mtry=mtry, num.trees=nb.trees, importance = "impurity", num.threads = n.cpus)
    imp_scores <- rf.model$variable.importance  
    imp_scores.names <- names(imp_scores)
    mat.rapid[imp_scores.names, gn.i] <- as.numeric(imp_scores)
  }
  print(Sys.time()-strt)
  print("..finished.")
  
  return(mat.rapid / n.samples)
}       

compute_xgboost_based_GRN <- function(mat.expression=mat.expression, nb.trees=500, set.regulators = NULL, set.genes = NULL){
  
  library(xgboost)
  
  
  mat.expression.norm <- t(mat.expression)
  mat.expression.norm <- apply(mat.expression.norm, 2, function(x) { (x - mean(x)) / sd(x) } )
  n.genes <- dim(mat.expression.norm)[2]
  genes<- colnames(mat.expression.norm)
  n.samples <- dim(mat.expression.norm)[1]
  
  if(is.null(set.genes)){
    n.genes <- n.genes
    genes <- genes
    
  }else{
    n.genes <- length(set.genes)
    genes <- set.genes
    
  }
  
  #mat.expression.norm <- mat.expression.norm[,genes]
  
  if (is.null(set.regulators)) {
    n.regulators <- n.genes
    regulators <- genes
  } else {
    n.regulators <- length(set.regulators)
    # regulators provided by names
    if(is.character(set.regulators)){
      regulators <- set.regulators
      genes.undefined <- setdiff(regulators, genes)
      if (length(genes.undefined) != 0) {
        stop(paste("Error: genes: ", paste(genes.undefined, collapse = ",")," not represented in gene expression matrix \n", sep=""))
      }
      # regulators provided by indices
    }else if (is.numeric(set.regulators)) {
      regulators <- genes[set.regulators]
    }else{
      stop("Error: invalid regulator format")
    }
  }
  
  # print(paste("Performing random-forest regression based gene regulatory network inference (# of decision trees per gene is: ", nb.trees, ", # of regulators per decision tree node: ", mtry, sep=""))
  mat.rapid <- matrix(0.0, nrow=n.regulators, ncol=n.genes, dimnames = list(regulators, genes))
  
  strt<-Sys.time()
  for(i in 1:n.genes){
    cat("Processing... ", round(i/n.genes * 100, digits = 2) , "%", "\r"); flush.console() 
    gn.i <- genes[i]
    # remove target gene from set of regulators
    regulators.i <- setdiff(regulators, gn.i)
    x <- mat.expression.norm[,regulators.i, drop=FALSE] # avoid coercion to numeric
    y <- mat.expression.norm[,gn.i]
    # df.xy <- cbind(as.data.frame(x),y)
    # browser()
    # rf.model <- ranger(x=x,y=y, mtry=mtry, num.trees=nb.trees, importance = "impurity", num.threads = n.cpus)
    #boosting训练
    xgb.model <- xgboost(data = x, 
                         label = y, 
                         eta = 0.1,
                         max_depth =3 , 
                         nround=nb.trees, 
                         subsample = 0.8,
                         colsample_bytree = 0.6,
                         objective = "reg:squarederror",
                         nthread = 3,
                         verbose=0
                         # ,
                         # print_every_n=150
    )
    #调用特征重要性
    imp <- xgb.importance(model=xgb.model)
    xgb.plot.importance(imp)
    mat.rapid[imp$Feature, gn.i] <- as.numeric(imp$Importance)
  }
  print(Sys.time()-strt)
  print("..finished.")
  
  return(mat.rapid / n.samples)
}       




#GRACE的函数----
#GRACE_pipeline中用到的函数

#GRACE_helperfunctions.R
request_libraries <- function(){
  
  # is.installed("graphics")
  # is.installed("Matrix")
  # is.installed("ggplot2")
  # is.installed("ROCR")
  # is.installed("reshape2")
  # is.installed("CRF")
  # is.installed("caTools")
  # is.installed("plyr")
  # is.installed("foreach")
  # is.installed("doParallel")
  # is.installed("igraph")
  # is.installed("GenSA")
  # is.installed("ranger")
  # 
  # #detach("package:TFBSTools", unload=TRUE)
  # 
  # if(!require(netbenchmark)){
  #   source("https://bioconductor.org/biocLite.R")
  #   biocLite("netbenchmark")
  # }
  
  
  library(graphics)
  library(Matrix)
  library(ggplot2)
  library(ROCR)
  library(reshape2)
  library(CRF)
  library(caTools)
  library(plyr)
  library(foreach)
  library(doParallel)
  library(igraph)
  library(GenSA)
  library(ranger)
  # library(netbenchmark)
  
}

#改写
request_dataset1 <- function(){
  # browser()
  # set paths to novel datasets 
  lst.data <- vector(mode = "list", length = 8)
  
  message("loading datasets")
  
  ### 二、load DNA binding set 
  # the filtered regulatory network based on integrating binding with developmental gene expression ( as computed in GRACE_initieal_GRN_inference.R )
  m.initial_grn <- readRDS(paste(rdir,'m.grn_rf_string.rds',sep = '/'))
  v.tfs <- rownames(m.initial_grn)
  v.tgs <- colnames(m.initial_grn)
  v.gns <- unique(c(v.tfs, v.tgs))
  
  #一、 load regulatory evidence (2 column dataframe with experimental regulatory links)
  #reference_TFTF_network.txt for SCODE-datasets
  # df.GS <- read.table(paste(dir,'reference_TFTF_network.txt',sep = '/'), header = FALSE,col.names = c('TF','Target','TFid','Targetid'), sep = "\t", quote = "", stringsAsFactors = FALSE)
  df.GS <- read.table(paste('BEELINE-Networks/Networks',spec,'Non-Specific-Chip-seq-network.csv',sep='/'),
                      stringsAsFactors = FALSE,header = TRUE,sep=',',col.names = c("gn.a", "gn.b"))
  df.GS["value"] <- 1
  m.GS <- df_link_to_mat2(df.GS,v.tfs,v.tgs)
  # m.GS <- acast(df.GS, gn.a~gn.b)
  # m.GS[is.na(m.GS)] <- 0
  # m.GS[m.GS != 0] <- 1
  # class(m.GS) <- "numeric"   
  
  # 三、load function evidence  (2 column dataframe with experimental gene ontology annotations)
  df.cofunctional_benchmark <- read.table('co-function/MouseNetV2_GS_symbol.txt', header = FALSE,col.names = c("gn.a", "gn.b","value"), sep = "\t", stringsAsFactors = FALSE)
  df.cofunctional_benchmark["value"] <- 1
  mat.tgs.cofunctional_benchmark <- df_link_to_mat2(df.cofunctional_benchmark,v.tgs,v.tgs)
  # names(df.cofunctional_benchmark) <- c("gn.a", "gn.b")
  # df.cofunctional_benchmark$gn.a <- toupper(df.cofunctional_benchmark$gn.a)
  # df.cofunctional_benchmark$gn.b <- toupper(df.cofunctional_benchmark$gn.b)
  # df.cofunctional_benchmark <- subset(df.cofunctional_benchmark, df.cofunctional_benchmark$gn.a %in% v.gns)
  # df.cofunctional_benchmark <- subset(df.cofunctional_benchmark, df.cofunctional_benchmark$gn.b %in% v.gns)
  # df.cofunctional_benchmark["val"] <- 1
  # m.cofunctional_benchmark <- acast(df.cofunctional_benchmark, gn.a~gn.b, value.var = "val")
  # m.cofunctional_benchmark[is.na(m.cofunctional_benchmark)] <- 0
  # class(m.cofunctional_benchmark) <- "numeric"
  # m.cofunctional_benchmark <- as(m.cofunctional_benchmark, "CsparseMatrix")
  # 
  # mat.tgs.cofunctional_benchmark <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
  # gns.rows <- intersect(rownames(mat.tgs.cofunctional_benchmark), rownames(m.cofunctional_benchmark))
  # gns.cols <- intersect(colnames(mat.tgs.cofunctional_benchmark), colnames(m.cofunctional_benchmark))
  # mat.tgs.cofunctional_benchmark[gns.rows, gns.cols] <- m.cofunctional_benchmark[gns.rows, gns.cols]
  # 
  
  # 四、load cofunctional network (2 column dataframe)
  df.cofunctional <- read.table('co-function/MouseNetV2_symbol.txt', header  = FALSE,col.names = c("gn.a", "gn.b", "value"), sep = "\t", stringsAsFactors = FALSE)
  mat.RN <- df_link_to_mat2(df.cofunctional,v.tgs,v.tgs)
  # names(df.cofunctional) <- c("gn.a", "gn.b", "value")
  # df.cofunctional$gn.a <- toupper(df.cofunctional$gn.a)
  # df.cofunctional$gn.b <- toupper(df.cofunctional$gn.b)
  # df.cofunctional <- subset(df.cofunctional, df.cofunctional$gn.a %in% v.tgs)
  # df.cofunctional <- subset(df.cofunctional, df.cofunctional$gn.b %in% v.tgs)
  # m.cofunctional <- acast(df.cofunctional, gn.a~gn.b, value.var = "value")
  # m.cofunctional[is.na(m.cofunctional)] <- 0
  # class(m.cofunctional) <- "numeric"  
  # m.cofunctional <- as(m.cofunctional, "CsparseMatrix")
  # 
  # # expand aranet matrix to cover all tgs
  # mat.RN <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
  # gns.rows <- intersect(rownames(mat.RN), rownames(m.cofunctional))
  # gns.cols <- intersect(colnames(mat.RN), colnames(m.cofunctional))
  # mat.RN[gns.rows, gns.cols] <- m.cofunctional[gns.rows, gns.cols]
  # 
  # update here to account for symmetry
  mat.RN <- mapping_intervals(mat.RN, min(mat.RN), max(mat.RN), 0, 1) 
  
  #五、STRING-network用于validation（function evidence）
  df.string <- read.table(paste('BEELINE-Networks/Networks',spec,'STRING-network.csv',sep='/'),
                          stringsAsFactors = FALSE,header = TRUE,sep=',',col.names = c("gn.a", "gn.b"))
  df.string["value"] <- 1
  mat.string <- df_link_to_mat2(df.string,v.tgs,v.tgs)
  mat.string[mat.string != 0] <- 1
  
  #六、specific-chip-seq-network用于validation（regulatory evidence）
  chip_dset=as.character(unlist(strsplit(dset, split = "-")))[1]
  df.chip=read.table(file = paste('BEELINE-Networks/Networks',spec,paste(chip_dset,'ChIP-seq-network.csv',sep='-'),sep='/'),
                     stringsAsFactors = FALSE,header = TRUE,sep=',',col.names = c("gn.a", "gn.b","value"))
  df.chip["value"] <- 1
  mat.chip <- df_link_to_mat2(df.chip,v.tfs,v.tgs)
  mat.chip[mat.chip != 0] <- 1
  
  return(list(mat.conserved_motif_plus_expression=m.initial_grn[v.tfs,], 
              v.tfs=v.tfs, v.tgs=v.tgs, mat.GS = m.GS, mat.FunctionalAssociation=mat.RN,
              mat.GO.tgs=mat.tgs.cofunctional_benchmark,mat.string=mat.string,mat.chip=mat.chip))
  
}

#GRACE_load_datasets.R
request_dataset <- function(){
  # browser()
  # set paths to novel datasets 
  lst.data <- vector(mode = "list", length = 8)
  
  message("loading datasets")
  
  # load regulatory evidence (2 column dataframe with experimental regulatory links)
  #reference_TFTF_network.txt for SCODE-datasets
  # df.GS <- read.table(paste(dir,'reference_TFTF_network.txt',sep = '/'), header = FALSE,col.names = c('TF','Target','TFid','Targetid'), sep = "\t", quote = "", stringsAsFactors = FALSE)
  df.GS <- read.table(paste('BEELINE-Networks/Networks',spec,'Non-Specific-Chip-seq-network.csv',sep='/'),
                      stringsAsFactors = FALSE,header = TRUE,sep=',',col.names = c("TF", "Target"))
  m.GS <- acast(df.GS, TF~Target)
  m.GS[is.na(m.GS)] <- 0
  m.GS[m.GS != 0] <- 1
  class(m.GS) <- "numeric"   
  
  ### load DNA binding set 
  
  # the filtered regulatory network based on integrating binding with developmental gene expression ( as computed in GRACE_initieal_GRN_inference.R )
  m.initial_grn <- readRDS(paste(rdir,'m.grn_rf.rds',sep = '/'))
  v.tfs <- rownames(m.initial_grn)
  v.tgs <- colnames(m.initial_grn)
  v.gns <- unique(c(v.tfs, v.tgs))
  
  # 不用co-functional
  # load function evidence  (2 column dataframe with experimental gene ontology annotations)
  df.cofunctional_benchmark <- read.table('co-function/MouseNetV2_GS_symbol.txt', header = FALSE,col.names = c("gn.a", "gn.b"), sep = "\t", stringsAsFactors = FALSE)
  df.cofunctional_benchmark$value <-rep(1,dim(df.cofunctional_benchmark)[1]) 
  mat.tgs.cofunctional_benchmark <- df_link_to_mat2(df.cofunctional_benchmark,v.tfs,v.tgs)
  # names(df.cofunctional_benchmark) <- c("gn.a", "gn.b")
  # df.cofunctional_benchmark$gn.a <- toupper(df.cofunctional_benchmark$gn.a)
  # df.cofunctional_benchmark$gn.b <- toupper(df.cofunctional_benchmark$gn.b)
  # df.cofunctional_benchmark <- subset(df.cofunctional_benchmark, df.cofunctional_benchmark$gn.a %in% v.gns)
  # df.cofunctional_benchmark <- subset(df.cofunctional_benchmark, df.cofunctional_benchmark$gn.b %in% v.gns)
  # df.cofunctional_benchmark["val"] <- 1
  # m.cofunctional_benchmark <- acast(df.cofunctional_benchmark, gn.a~gn.b, value.var = "val")
  # m.cofunctional_benchmark[is.na(m.cofunctional_benchmark)] <- 0
  # class(m.cofunctional_benchmark) <- "numeric"
  # m.cofunctional_benchmark <- as(m.cofunctional_benchmark, "CsparseMatrix")
  # 
  # mat.tgs.cofunctional_benchmark <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
  # gns.rows <- intersect(rownames(mat.tgs.cofunctional_benchmark), rownames(m.cofunctional_benchmark))
  # gns.cols <- intersect(colnames(mat.tgs.cofunctional_benchmark), colnames(m.cofunctional_benchmark))
  # mat.tgs.cofunctional_benchmark[gns.rows, gns.cols] <- m.cofunctional_benchmark[gns.rows, gns.cols]
  # 
  # browser()
  # load cofunctional network (2 column dataframe)
  df.cofunctional <- read.table('co-function/MouseNetV2_symbol.txt', header  = FALSE,col.names = c("gn.a", "gn.b", "value"), sep = "\t", stringsAsFactors = FALSE)
  mat.RN <- df_link_to_mat2(df.cofunctional,v.tfs,v.tgs)
  # names(df.cofunctional) <- c("gn.a", "gn.b", "value")
  # df.cofunctional$gn.a <- toupper(df.cofunctional$gn.a)
  # df.cofunctional$gn.b <- toupper(df.cofunctional$gn.b)
  # df.cofunctional <- subset(df.cofunctional, df.cofunctional$gn.a %in% v.tgs)
  # df.cofunctional <- subset(df.cofunctional, df.cofunctional$gn.b %in% v.tgs)
  # m.cofunctional <- acast(df.cofunctional, gn.a~gn.b, value.var = "value")
  # m.cofunctional[is.na(m.cofunctional)] <- 0
  # class(m.cofunctional) <- "numeric"  
  # m.cofunctional <- as(m.cofunctional, "CsparseMatrix")
  # 
  # # expand aranet matrix to cover all tgs
  # mat.RN <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
  # gns.rows <- intersect(rownames(mat.RN), rownames(m.cofunctional))
  # gns.cols <- intersect(colnames(mat.RN), colnames(m.cofunctional))
  # mat.RN[gns.rows, gns.cols] <- m.cofunctional[gns.rows, gns.cols]
  # 
  # update here to account for symmetry
  mat.RN <- mapping_intervals(mat.RN, min(mat.RN), max(mat.RN), 0, 1) 
  
  #STRING-network用于validation（function evidence）
  df.string <- read.table(paste('BEELINE-Networks/Networks',spec,'STRING-network.csv',sep='/'),
                          stringsAsFactors = FALSE,header = TRUE,sep=',',col.names = c("gn.a", "gn.b"))
  df.string$value <- rep(1,dim(df.string)[1])
  mat.RN <- df_link_to_mat2(df.cofunctional,v.tfs,v.tgs)
  
  #specific-chip-seq-network用于validation（regulatory evidence）
  chip_dset=as.character(unlist(strsplit(dset, split = "-")))[1]
  df.chip=read.table(file = paste('BEELINE-Networks/Networks',spec,paste(chip_dset,'ChIP-seq-network.csv',sep='-'),sep='/'),
                     stringsAsFactors = FALSE,header = TRUE,sep=',',col.names = c("gn.a", "gn.b","value"))
  df.chip$Specific=rep(1,dim(df.chip)[1])
  mat.chip <- df_link_to_mat2(df.chip,v.tfs,v.tgs)
  
  return(list(mat.conserved_motif_plus_expression=m.initial_grn[v.tfs,], 
              v.tfs=v.tfs, v.tgs=v.tgs, mat.GS = m.GS, mat.FunctionalAssociation=mat.RN,
              mat.GO.tgs=mat.tgs.cofunctional_benchmark,mat.string=mat.string,mat.chip=mat.chip))
  
}

mapping_intervals <- function(input, input_min, input_max, output_min, output_max){
  output = output_min + ((output_max - output_min) / (input_max - input_min)) * (input - input_min)
  return(output)
}

df_link_to_mat2<-function(df,rname,cname){
  #筛除非rname、cname中包含的基因
  df$gn.a <- toupper(df$gn.a)
  df$gn.b <- toupper(df$gn.b)
  df <- subset(df, df$gn.a %in% rname)
  df <- subset(df, df$gn.b %in% cname)
  # df<-df[duplicated(df)==FALSE,]  #去重
  # df$value=rep(1,dim(df)[1])
  
  # df转换为矩阵
  library(reshape2)
  mat.temp <- acast(df, gn.a~gn.b, value.var = "value")
  mat.temp[is.na(mat.temp)] <- 0
  class(mat.temp) <- "numeric"  
  library('Matrix')
  # mat.temp <- as(mat.temp, "CsparseMatrix") #稀疏矩阵
  mat.temp <- as(mat.temp, "dgCMatrix") #稀疏矩阵
  
  #扩充为rname*cname的矩阵
  mat <- Matrix(0, nrow = length(rname), ncol = length(cname), dimnames = list(rname, cname))
  gns.rows <- intersect(rownames(mat.temp), rownames(mat))
  gns.cols <- intersect(colnames(mat.temp), colnames(mat))
  mat[gns.rows, gns.cols] <- mat.temp[gns.rows, gns.cols]
  return(mat)
}


# only jaccard and th.min.coreg.tfs = 1
compute_model_assessment_parallel_optimized <- function(mat.p.grace_ensemble=mat.p.grace_ensemble, lst.benchmarks=lst.benchmarks){
  
  mat.p.crf= mat.p.grace_ensemble        
  mat.GS.grn=lst.benchmarks$mat.regulatory_evidence
  mat.bp.tgs=lst.benchmarks$mat.cofunctional_evidence
  
  
  idx.tfs <- which(rowSums(mat.p.crf) > 0)
  idx.tgs <- which(colSums(mat.p.crf) > 0)
  mat.p.crf <- mat.p.crf[idx.tfs,idx.tgs, drop = FALSE]
  
  ## atrm recovery
  tf <- intersect(rownames(mat.p.crf),rownames(mat.GS.grn)) 
  tgs <- intersect(colnames(mat.p.crf), colnames(mat.GS.grn))
  n.gs.grn <- sum(mat.p.crf[tf,tgs] * mat.GS.grn[tf, tgs])
  
  #mat.grn.MR <- jaccard(t((mat.p.crf)))
  mat.grn.MR <- jaccard(t(as.matrix(mat.p.crf)))
  mat.grn.MR <- as.matrix(mat.grn.MR); 
  diag(mat.grn.MR) <- 0;
  rownames(mat.grn.MR) <- colnames(mat.grn.MR) <- colnames(mat.p.crf) 
  mat.grn.MR <- as(mat.grn.MR, "CsparseMatrix")
  
  mat.grn.MR@x[mat.grn.MR@x < 0.5] <- 0
  mat.grn.MR@x[mat.grn.MR@x >= 0.5] <- 1
  
  tgs <- intersect(colnames(mat.p.crf), colnames(mat.bp.tgs))
  
  mat.grn.MR <- mat.grn.MR[tgs,tgs, drop = FALSE]
  mat.bp.tgs <- mat.bp.tgs[tgs,tgs, drop = FALSE]
  
  n.coreg.tg_pairs <- sum(mat.grn.MR) 
  
  # bp pairs (subset)
  idx <- which(mat.grn.MR == 1)
  n.coreg.bp_tg_pairs <- sum(mat.bp.tgs[idx])
  
  idx <- which(colSums(mat.p.crf) > 0)
  v.tgs.set <- colnames(mat.p.crf)[idx]
  n.tg_pairs.bp <- sum(mat.bp.tgs[v.tgs.set, v.tgs.set])
  
  n.tg_pairs.bp.max <- sum(mat.bp.tgs)
  
  n.pred <- sum(mat.p.crf)
  
  rm(mat.p.crf)
  #rm(mat.grn.MR)
  #rm(mat.grn.MR.count)
  rm(mat.GS.grn)
  rm(mat.bp.tgs)
  
  return(list(n.pred=n.pred, n.gs.grn=n.gs.grn, n.coreg.bp_tg_pairs=n.coreg.bp_tg_pairs,n.coreg.tg_pairs=n.coreg.tg_pairs,
              n.tg_pairs.bp=n.tg_pairs.bp, n.tg_pairs.bp.max=n.tg_pairs.bp.max, mat.grn.MR=mat.grn.MR))
  
}

compute_grace_complete_assessment1 <- function(mat.grace=mat.grace, 
                                               lst.benchmarks=lst.benchmarks,
                                               lst.validation=NULL,
                                               code.validation = c("")){
  mat.p.crf= mat.grace        
  mat.GS.grn=lst.benchmarks$mat.regulatory_evidence
  mat.bp.tgs=lst.benchmarks$mat.cofunctional_evidence
  
  idx.tfs <- which(rowSums(mat.p.crf) > 0)
  idx.tgs <- which(colSums(mat.p.crf) > 0)
  mat.p.crf <- mat.p.crf[idx.tfs,idx.tgs, drop = FALSE]
  
  #regulatory-evidence
  ## atrm recovery
  tf <- intersect(rownames(mat.p.crf),rownames(mat.GS.grn)) 
  tgs <- intersect(colnames(mat.p.crf), colnames(mat.GS.grn))
  n.gs.grn <- sum(mat.p.crf[tf,tgs] * mat.GS.grn[tf, tgs])
  
  #co-functional evidence
  mat.grn.MR <- jaccard(t(as.matrix(mat.p.crf)))
  mat.grn.MR <- as.matrix(mat.grn.MR); 
  diag(mat.grn.MR) <- 0;
  rownames(mat.grn.MR) <- colnames(mat.grn.MR) <- colnames(mat.p.crf) 
  mat.grn.MR <- as(mat.grn.MR, "CsparseMatrix")
  
  mat.grn.MR@x[mat.grn.MR@x < 0.5] <- 0
  mat.grn.MR@x[mat.grn.MR@x >= 0.5] <- 1
  tgs <- intersect(colnames(mat.p.crf), colnames(mat.bp.tgs))
  mat.grn.MR <- mat.grn.MR[tgs,tgs, drop = FALSE]
  mat.bp.tgs <- mat.bp.tgs[tgs,tgs, drop = FALSE]
  n.coreg.tg_pairs <- sum(mat.grn.MR) 
  
  # bp pairs (subset)
  idx <- which(mat.grn.MR == 1)
  n.coreg.bp_tg_pairs <- sum(mat.bp.tgs[idx])
  
  idx <- which(colSums(mat.p.crf) > 0)
  v.tgs.set <- colnames(mat.p.crf)[idx]
  n.tg_pairs.bp <- sum(mat.bp.tgs[v.tgs.set, v.tgs.set])
  
  n.tg_pairs.bp.max <- sum(mat.bp.tgs)
  
  n.pred <- sum(mat.p.crf)
  
  
  lst.results <- list(n.pred=n.pred, n.gs.grn=n.gs.grn, n.coreg.bp_tg_pairs=n.coreg.bp_tg_pairs,n.coreg.tg_pairs=n.coreg.tg_pairs,
                      n.tg_pairs.bp=n.tg_pairs.bp, n.tg_pairs.bp.max=n.tg_pairs.bp.max, mat.grn.MR=mat.grn.MR)
  # subset into regulatory and co-functional evidence 
  if(!is.null(lst.validation)){
    
    lst.validation.results <- numeric(length(lst.validation))
    names(lst.validation.results) <- names(lst.validation)
    
    tgs <- intersect(colnames(mat.grace), colnames(mat.grn.MR))
    mat.grn.MR <- mat.grn.MR[tgs,tgs]
    idx <- which(mat.grn.MR == 1)
    
    tfs <- rownames(mat.grace)
    
    for(v in 1:length(lst.validation)){
      
      if(code.validation[v] == "regulatory_evidence"){
        
        lst.validation.results[v] <- sum(as.numeric(lst.validation[[v]][tfs, tgs] * mat.grace[tfs, tgs]))
        
      }else{
        
        mat.MR.validation.tgs <- lst.validation[[v]][rownames(mat.grn.MR), colnames(mat.grn.MR)]
        lst.validation.results[v] <- sum(mat.MR.validation.tgs[idx])
        
        rm(mat.MR.validation.tgs)
      }
      
    }
    
    lst.results <- c(lst.results, lst.validation.results)
  }
  
  return(lst.results)
}

compute_grace_complete_assessment <- function(mat.grace=mat.grace, 
                                              lst.benchmarks=lst.benchmarks,
                                              lst.validation=NULL,
                                              code.validation = c("")){
  
  lst.res <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.grace, lst.benchmarks=lst.benchmarks)
  
  
  lst.results <- (list(  n.pred=lst.res$n.pred,
                         n.gs.grn=lst.res$n.gs.grn,
                         
                         n.coreg.bp_tg_pairs=lst.res$n.coreg.bp_tg_pairs,
                         n.coreg.tg_pairs=lst.res$n.coreg.tg_pairs,
                         
                         n.tg_pairs.bp=lst.res$n.tg_pairs.bp,
                         n.tg_pairs.bp.max=lst.res$n.tg_pairs.bp.max))
  
  # subset into regulatory and co-functional evidence 
  
  if(!is.null(lst.validation)){
    
    lst.validation.results <- numeric(length(lst.validation))
    names(lst.validation.results) <- names(lst.validation)
    
    tgs <- intersect(colnames(mat.grace), colnames(lst.res$mat.grn.MR))
    mat.grn.MR <- lst.res$mat.grn.MR[tgs,tgs]
    idx <- which(mat.grn.MR == 1)
    
    tfs <- rownames(mat.grace)
    
    for(v in 1:length(lst.validation)){
      
      if(code.validation[v] == "regulatory_evidence"){
        
        lst.validation.results[v] <- sum(as.numeric(lst.validation[[v]][tfs, tgs] * mat.grace[tfs, tgs]))
        
      }else{
        
        mat.MR.validation.tgs <- lst.validation[[v]][rownames(mat.grn.MR), colnames(mat.grn.MR)]
        lst.validation.results[v] <- sum(mat.MR.validation.tgs[idx])
        
        rm(mat.MR.validation.tgs)
      }
      
    }
    
    lst.results <- c(lst.results, lst.validation.results)
  }
  
  rm(lst.res)
  
  return(lst.results)
}


#GRACE_helperfunctions.R
compute_fmeasures_regulatory_network <- function(mat.grn=mat.gene_regulatory_network, 
                                                 
                                                 lst.benchmarks=lst.benchmarks,
                                                 lst.validation=NULL,
                                                 code.validation = c(""),
                                                 
                                                 n.samples = 1000,
                                                 n.cpus = 2){
  # browser()
  strt<-Sys.time() 
  cl<-makeCluster(min(n.samples, n.cpus))
  registerDoParallel(cl)
  
  ### running multiple bootstrap samples in parallel
  lst.res <- foreach(i = (n.samples - 1):0, .packages=c("Matrix", "reshape2")) %dopar% {
    
    #source("GRACE_Athaliana_evaluation.R")
    # source("GRACE_helperfunctions.R")
    source("GRACE_function.R")
    # i=1
    cut <- 1/n.samples * i
    # cat("Processing... ", 1 - round(cut, digits = 2) , "%", "\r"); flush.console()    
    
    v.grn <- as.numeric(mat.grn[mat.grn>0])
    th.grn <- as.numeric(quantile(v.grn, cut))
    
    mat.grn.tmp <- as.matrix(mat.grn)
    mat.grn.tmp[mat.grn.tmp < th.grn] <- 0
    mat.grn.tmp[mat.grn.tmp >= th.grn] <- 1
    
    browser()
    lst.results <- compute_grace_complete_assessment(mat.grace=mat.grn.tmp, 
                                                     lst.benchmarks=lst.benchmarks,
                                                     lst.validation=lst.validation,
                                                     code.validation = code.validation)
    
    lst.results <- c(lst.results, th.grn)
    lst.results <- c(lst.results, cut)
    n.idx <- length(lst.results) - 1
    names(lst.results)[n.idx:(n.idx+1)] <- c("th.grn", "cut")
    
    rm(mat.grn.tmp)
    
    return(lst.results)
  }
  
  
  stopCluster(cl)
  print(Sys.time()-strt)
  
  df.rate_density <- as.data.frame(matrix(0, nrow=n.samples, ncol = length(lst.res[[1]])))
  names(df.rate_density) <- names(lst.res[[1]])
  
  for(i in n.samples:1){
    for(j in 1:ncol(df.rate_density)){
      df.rate_density[i,j] <- lst.res[[i]][[j]]
    }
  }
  
  df.rate_density <- df.rate_density[order(df.rate_density$n.pred),]
  df.rate_density <- subset(df.rate_density, df.rate_density$th.grn != 0)
  
  return(df.rate_density)
}


jaccard <- function(m) {
  ## common values:
  A = tcrossprod(m) #交叉积：m %*% t(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  ## only non-zero values of common
  Aim = A[im]
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  ) 
  return( J )
}

#GRACE_helperfunctions.R
plot_fmeasure_scores <- function(df.rate_density){
  
  max.gs.grn <- max(df.rate_density$n.gs.grn)
  max.coreg.tgs.bp <- max(df.rate_density$n.coreg.bp_tg_pairs)
  max(df.rate_density$n.tg_pairs.bp)
  
  # extract ratio of precision x recall
  Precision.gs.grn <- (df.rate_density$n.gs.grn/df.rate_density$n.pred)
  #Recall.coreg.tgs.bp <- df.rate_density$n.coreg.tg_pairs.BP/max(df.rate_density$n.coreg.tg_pairs.BP)
  
  ratio.gs.grn_vs_bp.normalization <- 1 / max(Precision.gs.grn) #/ max(Recall.coreg.tgs.bp)
  max.coreg.bp <- max(df.rate_density$n.coreg.bp_tg_pairs)
  
  ####
  x <- df.rate_density$n.pred / max(df.rate_density$n.pred)
  
  y1 <- (df.rate_density$n.gs.grn/df.rate_density$n.pred)
  y1[is.na(y1)] <- 0
  y1 <- y1/max(y1) # normalization
  #plot(x,y1, type = "l") # training (ATRM)
  
  y2 <- df.rate_density$n.coreg.bp_tg_pairs/max(df.rate_density$n.coreg.bp_tg_pairs)
  y2[is.na(y2)] <- 0
  y2 <- y2#/max(y2)
  #lines(x, y2, col = "red") # training (GO-BP TG<->TG)
  
  beta.precision = 0.5
  fscore_beta <- (1 + beta.precision^2)*(y1*y2/((beta.precision^2 * y1)+y2))
  fscore_beta[is.na(fscore_beta)] <- 0
  fscore_beta_05 <- fscore_beta#/max(fscore_beta)
  
  beta.traditional = 1 #(44/10087)/(175/98596)
  fscore_beta <- (1 + beta.traditional^2)*(y1*y2/((beta.traditional^2 * y1)+y2))
  fscore_beta[is.na(fscore_beta)] <- 0
  fscore_beta_1 <- fscore_beta#/max(fscore_beta)
  
  df.rate_density["fscore_beta_05_normalized"] <- fscore_beta_05/max(fscore_beta_05)
  df.rate_density["fscore_beta_1_normalized"] <- fscore_beta_1/max(fscore_beta_1)
  
  df.rate_density["fscore_beta_05"] <- fscore_beta_05
  df.rate_density["fscore_beta_01"] <- fscore_beta_1
  
  ####
  
  data <- data.frame(Predictions = x,
                     Precision = y1, 
                     Recall = y2, 
                     fscore_05 = df.rate_density$fscore_beta_05_normalized,
                     fscore_1 = df.rate_density$fscore_beta_1_normalized)
  
  idx.max_05 <- max(which(data$fscore_05 == max(data$fscore_05)))
  idx.max_1 <- max(which(data$fscore_1 == max(data$fscore_1)))
  
  
  data <- data.frame(x = c(x,x,x,x),
                     y = c(y1, y2, df.rate_density$fscore_beta_05_normalized, df.rate_density$fscore_beta_1_normalized),
                     z = c(rep("Precision (number of known regulatory links)", length(x)),
                           rep("Recall (number of coregulated gene pairs with co-functional annotation)", length(x)),
                           rep("F-score (\u03B2 = 0.5)", length(x)),
                           rep("F-score (\u03B2 = 1.0)", length(x))))
  # 
  #                          rep(expression("f"[0.5]*"-score"), length(x)),
  #                          rep(expression("f"[1]*"-score"), length(x))))
  
  data$z <- factor(data$z, levels=unique(as.character(data$z)) )
  
  #cairo_pdf("example.pdf", width = 7.5, height = 6, family= "Helvetica")#"DejaVu Sans")
  
  (plot.PR <- ggplot(data, aes(x=x, y=y)) +
      geom_line(mapping=aes(colour=z,  linetype = z)) + 
      theme_bw() + 
      
      theme(legend.title=element_blank()) + 
      theme(legend.key = element_blank()) + 
      
      theme(
        plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
      ) +
      
      scale_color_manual(values = c("red","steelblue","darkgray","black")) +
      scale_linetype_manual(values=c("solid", "solid","dashed", "dashed")) + 
      
      ylab("Normalized scores") +
      xlab("Fraction of predicted links / coregulated gene pairs") +  
      #scale_y_log10() +
      theme(legend.text= element_text(size=12)) + 
      theme(axis.text=element_text(size=12, colour = "black"), axis.title=element_text(size=12)) +
      theme(legend.justification=c(1,0), legend.position=c(1,0)))
  
  #dev.off()
  
  return(list(df.rate_density=df.rate_density , plot.PR=plot.PR, ratio.gs.grn_vs_bp.normalization=ratio.gs.grn_vs_bp.normalization, max.coreg.bp=max.coreg.bp))
  
}


#GRACE.R
contruct_grace_meta_connectivity_network <- function(mat.grn=mat.grn,  mat.FunctionalAssociation=mat.FunctionalAssociation,
                                                     th.preselection = 0.0, n.cores = 2){
  
  mat.RN <- mat.FunctionalAssociation
  v.grn.significant <- th.preselection
  
  df.grn <- as.data.frame(as.table(mat.grn), stringsAsFactors = FALSE) 
  names(df.grn)[1:3] <- c("TF","TG","v.grn")
  df.grn <- subset(df.grn, df.grn$v.grn > 0)
  df.grn <- df.grn[order(-df.grn$v.grn),] 
  
  if(v.grn.significant < df.grn$v.grn[nrow(df.grn)]){
    v.grn.significant <- df.grn$v.grn[nrow(df.grn)]
  }
  df.grn <- subset(df.grn, df.grn$v.grn >= v.grn.significant)
  
  tfs <- unique(df.grn$TF) 
  tgs <- unique(df.grn$TG)
  
  v.grn.node <- df.grn$v.grn
  n.links <- nrow(df.grn)
  df.grn["link.id"] <- seq(1, n.links) 
  rownames(df.grn) <- df.grn$link.id
  
  tf.map <- new.env(hash=T, parent=emptyenv())
  for(j in 1:length(tfs)){ tf.map[[tfs[j]]] <- j}
  df.grn["tf.id"] <- 0
  for(j in 1:length(tfs)){
    #cat("Processing... ", round(j/length(tfs) * 100, digits = 2) , "%", "\r"); flush.console()   
    idx <- which(df.grn$TF == tfs[j]); df.grn$tf.id[idx] <- tf.map[[tfs[j]]]
  }
  
  tg.map <- new.env(hash=T, parent=emptyenv())
  for(j in 1:length(tgs)){ tg.map[[tgs[j]]] <- j}
  df.grn["tg.id"] <- 0
  for(j in 1:length(tgs)){    
    #cat("Processing... ", round(j/length(tgs) * 100, digits = 2) , "%", "\r"); flush.console()   
    idx <- which(df.grn$TG == tgs[j]); df.grn$tg.id[idx] <- tg.map[[tgs[j]]]
  }
  
  mat.link_id <- acast(df.grn, TF~TG, value.var="link.id")
  mat.link_id[is.na(mat.link_id)] <- 0; 
  mat.link_id <- as(mat.link_id, "CsparseMatrix")
  
  ### IMPORTANT - Reorder tfs and tgs in all relevant matrices
  mat.link_id <- mat.link_id[tfs, tgs]
  mat.RN <- mat.RN[tgs, tgs]; 
  
  v.links.grn <- df.grn$link.id
  
  adj.links <- build_link_ranking_connectivity_matrix(df.grn=df.grn, mat.RN=mat.RN, tfs=tfs, tgs=tgs, n.links=n.links, n.cores = n.cores, v.th.rank_cutoff=1)[v.links.grn, v.links.grn]
  rm(mat.RN)
  
  return(list(df.grn=df.grn, adj.links=adj.links, v.grn.significant=v.grn.significant))
}

#GRACE.R
build_link_ranking_connectivity_matrix <- function(df.grn=df.grn, mat.RN=mat.RN, tfs=tfs, tgs=tgs, n.links=n.links, v.th.rank_cutoff=1, n.cores = 2){
  
  mat.RN.sparse <- Matrix(mat.RN)
  mat.ranking <- Matrix(0, nrow = length(tfs), ncol = n.links, dimnames = list(tfs, as.character(seq(1:n.links))))
  for(i in 1:length(tfs)){
    #cat("Processing... ", round(i/length(tfs) * 100, digits = 2) , "%", "\r"); flush.console()    
    df.grn.set <- subset(df.grn, df.grn$TF %in% tfs[i])
    rownames(df.grn.set) <- df.grn.set$link.id
    mat.ranking[i, df.grn.set$link.id] <- df.grn.set$v.grn  #df.grn.set$rank.global ##
  }
  
  v.links_tgs.map <- df.grn$link.id
  names(v.links_tgs.map) <- df.grn$tg.id
  v.links_tgs.map <- v.links_tgs.map[order(v.links_tgs.map)]
  
  v.cutoff <- build_estimate_rank_difference_cutoff(mat.ranking=mat.ranking, tfs=tfs, select = select, v.th.rank_cutoff = v.th.rank_cutoff)
  n.total <- length(tfs)
  n.each <- floor(n.total/n.cores)
  
  cl<-makeCluster(n.cores)
  registerDoParallel(cl)
  strt<-Sys.time()
  lst.adj <- foreach(r = 1:n.cores, .packages=c("Matrix", "reshape2")) %dopar% {    
    v.start <- n.each * (r - 1) + 1
    v.stop <- n.each * r
    tfs.set <- tfs[v.start:v.stop]
    mat.ranks <- Matrix(0, nrow = n.links, ncol = n.links, dimnames = list(as.character(seq(1:n.links)), as.character(seq(1:n.links))))
    for(i in 1:length(tfs.set)){
      cat("Processing... ", round(i/length(tfs.set) * 100, digits = 2) , "%", "\r"); flush.console()    
      v.ranks <- mat.ranking[tfs.set[i],]
      v.ranks <- v.ranks[v.ranks > 0]
      if(length(v.ranks) > 1){  
        df.links <- as.data.frame(t(combn(v.ranks,2)), stringsAsFactors = FALSE)
        v.rank_combs <- abs(df.links$V1 - df.links$V2) 
        idx <- which(v.rank_combs <= v.cutoff) #seq(1:length(v.rank_combs)) #
        if(length(idx) > 0){
          v.rank_combs <- v.rank_combs[idx]  
          df.links <- as.data.frame(t(combn(names(v.ranks),2)), stringsAsFactors = FALSE)
          df.links <- df.links[idx,]
          
          l1 <- as.numeric(df.links$V1)
          l2 <- as.numeric(df.links$V2)
          
          ind.pcc <- cbind(as.numeric(names(v.links_tgs.map)[l1]), as.numeric(names(v.links_tgs.map)[l2]))
          v.pccs <- mat.RN.sparse[ind.pcc]
          
          ind <- cbind(l1,l2)
          mat.ranks[ind] <- v.pccs / (1 + v.rank_combs)
          #mat.ranks[ind] <- 1 / (1 + v.rank_combs) # modified to test influence of Aranet
          #mat.ranks[ind] <- (v.rank_combs) / (1 +v.pccs) # modified to test inverse (negative test)
        }
      }
    }
    mat.ranks
  }
  stopCluster(cl)
  mat.ranks_links <- Reduce('+',lst.adj)
  
  ### residuals
  v.start <- n.each * n.cores + 1
  v.stop <- n.total 
  if(v.start <= v.stop){
    tfs.set <- tfs[v.start:v.stop]
    for(i in 1:length(tfs.set)){
      #cat("Processing... ", round(i/length(tfs.set) * 100, digits = 2) , "%", "\r"); flush.console()    
      v.ranks <- mat.ranking[tfs.set[i],]
      v.ranks <- v.ranks[v.ranks > 0]
      if(length(v.ranks) > 1){  
        df.links <- as.data.frame(t(combn(v.ranks,2)), stringsAsFactors = FALSE)
        v.rank_combs <- abs(df.links$V1 - df.links$V2) 
        idx <- which(v.rank_combs <= v.cutoff) # seq(1:length(v.rank_combs)) #
        if(length(idx) > 0){
          v.rank_combs <- v.rank_combs[idx]
          df.links <- as.data.frame(t(combn(names(v.ranks),2)), stringsAsFactors = FALSE)
          df.links <- df.links[idx,]
          l1 <- as.numeric(df.links$V1)
          l2 <- as.numeric(df.links$V2)
          
          ind.pcc <- cbind(as.numeric(names(v.links_tgs.map)[l1]), as.numeric(names(v.links_tgs.map)[l2]))
          v.pccs <- mat.RN.sparse[ind.pcc]
          
          ind <- cbind(l1,l2)
          mat.ranks_links[ind] <- v.pccs / (1 + v.rank_combs) #对应公式（1）
          #mat.ranks_links[ind] <- 1 / (1 + v.rank_combs) # modified to test influence of Aranet
          #mat.ranks_links[ind] <- (v.rank_combs) / (v.pccs + 1) # modified to test inverse (negative test)
        } 
      }
    }
  }
  mat.ranks_links.target_connectivity <- mat.ranks_links + t(mat.ranks_links)
  
  rm(mat.RN.sparse)
  rm(mat.ranking)
  rm(lst.adj)
  
  mat.ranks_links <- mat.ranks_links.target_connectivity
  
  return(mat.ranks_links)
}

build_estimate_rank_difference_cutoff <- function(mat.ranking=mat.ranking, tfs=tfs, select = "global", v.th.rank_cutoff = 0.05){
  
  v.rank_combs <- numeric()
  for(i in 1:length(tfs)){
    v.ranks <- mat.ranking[tfs[i],]
    v.ranks <- v.ranks[v.ranks > 0]
    if(length(v.ranks) > 1){  
      df.links <- as.data.frame(t(combn(v.ranks,2)), stringsAsFactors = FALSE)
      v.rank_combs <- c(v.rank_combs, abs(df.links$V1 - df.links$V2))
    }
  }
  v.cutoff <- as.numeric(quantile(v.rank_combs, v.th.rank_cutoff))
  return(v.cutoff) 
}

prepare_meta_modules <- function(df.grn=df.grn, adj.links=adj.links){
  
  mat.grn.specs <- as.matrix(df.grn)
  
  v.links.singletons <- rownames(adj.links)[which(rowSums(adj.links) == 0)] 
  
  ## remove singleton nodes
  adj.links.nonsingleton <- adj.links[!rownames(adj.links) %in% v.links.singletons, !colnames(adj.links) %in% v.links.singletons]
  #adj.links.singleton <- adj.links[rownames(adj.links) %in% v.links.singletons, colnames(adj.links) %in% v.links.singletons]
  
  ### STEP 2 - compute Markov clustering     
  #print("extract link modules...")
  g.coreg <- graph.adjacency(adj.links.nonsingleton, weighted = TRUE, mode = "undirected")
  clust.coreg <- clusters(g.coreg) # igraph based pre-MCL extraction of connected components s
  v.coreg_modules <- unique(clust.coreg$membership)
  n.crfs <- length(v.coreg_modules)
  
  lst.crfs <- vector(mode = "list", length = n.crfs)
  lst.crf.prim_based_links <- vector(mode = "list", length = n.crfs)
  lst.v.grn.node <- vector(mode = "list", length = n.crfs)
  lst.tf.ids <- vector(mode = "list", length = n.crfs)
  lst.tg.ids <- vector(mode = "list", length = n.crfs)
  
  for(i in 1:length(v.coreg_modules)){ 
    # cat("Processing... ", round(i/length(v.coreg_modules) * 100, digits = 2) , "%", "\r"); flush.console()  
    
    v.links <- names(which(clust.coreg$membership == v.coreg_modules[i]))
    names(v.links) <- seq(1:length(v.links))
    
    adj.links.set <- adj.links.nonsingleton[v.links,v.links]
    
    ru.adj.links <- adj.links.nonsingleton[v.links,v.links]
    ru.adj.links@x <- 1 - ru.adj.links@x
    g.ru <- graph.adjacency(ru.adj.links, weighted = TRUE, mode = "undirected")
    g.ru.mst <- minimum.spanning.tree(g.ru, algorithm= "prim")  #minimum.spanning.tree 最小生成树算法
    
    df.ru.edges <- get.edges(g.ru.mst, E(g.ru.mst))
    ru.adj.links <- Matrix(0, nrow = length(v.links), ncol = length(v.links), dimnames = list(v.links,v.links))
    
    l1 <- as.numeric(names(v.links[df.ru.edges[,1]]))
    l2 <- as.numeric(names(v.links[df.ru.edges[,2]]))
    ru.adj.links[cbind(l1,l2)] <- adj.links.set[cbind(l1,l2)]
    rm(adj.links.set)
    
    crf <- make.crf(ru.adj.links, n.states = 2)
    
    lst.crfs[[i]] <- crf
    lst.crf.prim_based_links[[i]] <- ru.adj.links[crf$edges]
    
    v.links <- names(which(clust.coreg$membership == v.coreg_modules[i]))
    mat.set <- mat.grn.specs[v.links,]
    v.grn.node <- as.numeric(mat.set[,3]) # instead of rank : 6
    lst.v.grn.node[[i]] <- v.grn.node
    
    tf.ids <- as.numeric(mat.set[,5])
    tg.ids <- as.numeric(mat.set[,6])
    
    lst.tf.ids[[i]] <- tf.ids
    lst.tg.ids[[i]] <- tg.ids
  } 
  
  return(list(v.coreg_modules=v.coreg_modules, lst.crfs=lst.crfs, lst.crf.prim_based_links=lst.crf.prim_based_links, g.coreg=g.coreg, 
              lst.v.grn.node=lst.v.grn.node, lst.tf.ids=lst.tf.ids, lst.tg.ids=lst.tg.ids))
}

#GRACE_optimization_routines.R
evaluate_mode_parameter <- function(lst.grace_support=lst.grace_support, 
                                    lst.modules=lst.modules, 
                                    truereg=truereg, 
                                    mat.cofunctional_evidence=mat.cofunctional_evidence, 
                                    
                                    gamma.gridSearch=gamma.gridSearch,
                                    lambda.gridSearch=lambda.gridSearch,
                                    th.percentage.hyperparameters=th.percentage.hyperparameters,
                                    
                                    ratio.gs.grn_vs_bp.normalization = ratio.gs.grn_vs_bp.normalization,
                                    max.coreg.bp = max.coreg.bp,
                                    
                                    n.sample_size=n.sample_size, 
                                    max.call = 100, 
                                    v.seed = 1234, 
                                    beta = 1.0){
  
  # source("GRACE_optimization_routines.R")
  # source("GRACE.R")
  # source("GRACE_helperfunctions.R")
  
  #### Bootstrapping (based on regulators based meta module sampling (0.632))
  
  # set.seed(v.seed)
  
  tfs.mrf <- unique(lst.grace_support$df.grn$TF) 
  tgs.mrf <- unique(lst.grace_support$df.grn$TG)
  
  lst.bootstrap_sample <- generate_bootstrapping_sample(lst.grace_support=lst.grace_support,lst.modules=lst.modules, truereg=truereg, mat.cofunctional_evidence=mat.cofunctional_evidence, n.sample_size = n.sample_size)
  lst.grace_support.bootstrap_sample <- lst.bootstrap_sample$lst.grace_support.bootstrap_sample
  lst.modules.bootstrap_sample <- lst.bootstrap_sample$lst.modules.bootstrap_sample
  
  lst.gs.grn_bp.train <- list(lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.GS.bootstrap_sample, lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.cofunctional_evidence.bootstrap_sample)    #list(mat.GS.train, mat.BP.coreg.train)
  names(lst.gs.grn_bp.train) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")
  
  ## Step 1: Grid Search (Coarse) optimization to identify finer range for Simulated annealing optimization ##
  strt<-Sys.time() 
  #lst.hyperparameter <- gridSearchOptimization(lst.gs.grn_bp.train=lst.gs.grn_bp.train,gamma.gridSearch, lambda.gridSearch, th.percentage.hyperparameters, ratio.gs.grn_vs_bp.normalization = ratio.gs.grn_vs_bp.normalization, max.coreg.bp = max.coreg.bp)
  
  # add begin grid search
  df.gridSearch <- expand.grid(gamma.gridSearch, lambda.gridSearch)
  names(df.gridSearch) <- c("gamma", "lambda")
  
  df.gridSearch["f_score"] <- 0
  df.gridSearch["n.predictions"] <- 0
  df.gridSearch["n.gs.grn"] <- 0
  df.gridSearch["n.coreg.bp.tgs"] <- 0
  
  for(j in 1:nrow(df.gridSearch)){
    
    cat("Processing... ", round(j/nrow(df.gridSearch) * 100, digits = 2) , "%", "\r"); flush.console()  
    v.alpha = v.gamma = df.gridSearch$gamma[j]; v.lambda= df.gridSearch$lambda[j];
    
    if(v.alpha >= v.gamma){
      
      lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields (lst.grace_support=lst.grace_support.bootstrap_sample, lst.modules=lst.modules.bootstrap_sample,
                                                                 v.alpha = v.alpha, v.gamma = v.gamma, v.lambda = v.lambda, b.prim=TRUE, tfs = tfs.mrf, tgs = tgs.mrf) 
      
      mat.p.crf <- lst.mat.p.crf$mat.p.crf # get if equal to zero ... , could be based on bootstrap sample, no value in there... 
      
      if(!is.null(mat.p.crf) & !is.null(dim(mat.p.crf)) & sum(mat.p.crf) > 0){
        
        lst.eval <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.p.crf,
                                                                lst.benchmarks=lst.gs.grn_bp.train)
        
        
        y1 <- (lst.eval$n.gs.grn/lst.eval$n.pred)
        y1[is.na(y1)] <- 0
        if(b.normalize_precision){
          y1 <- ratio.gs.grn_vs_bp.normalization * y1  
        }
        
        y2 <- lst.eval$n.coreg.bp_tg_pairs/max.coreg.bp
        y2[is.na(y2)] <- 0
        
        fscore_beta <- (1 + beta^2)*(y1*y2/((beta^2 * y1)+y2))
        fscore_beta[is.na(fscore_beta)] <- 0
        
        fscore <- fscore_beta
        n.gs.grn <- lst.eval$n.gs.grn
        n.coreg.bp.tgs <- lst.eval$n.coreg.bp_tg_pairs
        n.pred <- lst.eval$n.pred
        
      }else{
        fscore <- - 100
        n.gs.grn <- - 1
        n.coreg.bp.tgs <- - 1
        n.pred <- - 1
      }
      
    }else{
      fscore <- - 100
      n.gs.grn <- - 1
      n.coreg.bp.tgs <- - 1
      n.pred <- - 1
    }
    
    df.gridSearch$f_score[j] <- fscore
    df.gridSearch$n.predictions[j] <- n.pred
    df.gridSearch$n.gs.grn[j] <- n.gs.grn
    df.gridSearch$n.coreg.bp.tgs[j] <- n.coreg.bp.tgs
    #print(paste("alpha", v.alpha, ", gamma", v.gamma, ", lambda", v.lambda, ", f_score", fscore, "n.predictions", n.pred, "n.gs.grn",n.gs.grn ,"n.coreg.bp.tgs", n.coreg.bp.tgs))
  }
  
  # return Simulated annealing range 
  v.fscore.selection <- quantile(df.gridSearch$f_score, th.percentage.hyperparameters)
  df.gridSearch <- subset(df.gridSearch, df.gridSearch$f_score >= v.fscore.selection)
  
  v.double_fscores <- names(which(table(df.gridSearch$f_score) > 1))
  if(length(v.double_fscores) > 0){
    for(l in 1:length(v.double_fscores)){
      v.lambda.tmp <- subset(df.gridSearch, df.gridSearch$f_score == v.double_fscores[l])$lambda
      v.lambda.tmp <- v.lambda.tmp[!v.lambda.tmp %in% min(v.lambda.tmp)]
      
      idx.remove <- which(df.gridSearch$f_score == v.double_fscores[l] & df.gridSearch$lambda %in% v.lambda.tmp)
      df.gridSearch <- df.gridSearch[-idx.remove,]
    }
  }
  
  v.gamma.selection <- unique(df.gridSearch$gamma)
  v.lambda.selection <- unique(df.gridSearch$lambda)
  
  w  <- numeric(2)
  lower_limits <- numeric(2)
  upper_limits <- numeric(2)
  
  # gamma - single value handling
  if(length(v.gamma.selection) == 1){
    idx.gamma <- which(gamma.gridSearch == v.gamma.selection)
    if(idx.gamma == 1){
      lower_limits[[1]] <- v.gamma.selection
      upper_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma - 1]) # inverse indexing because of v.grn ordering
    }else if(idx.gamma == length(gamma.gridSearch)){
      lower_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma + 1])
      upper_limits[[1]] <- v.gamma.selection
    }else{
      lower_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma + 1])
      upper_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma - 1])
    }
    w[[1]] <- v.gamma.selection
  }else{
    w[[1]] <- min(v.gamma.selection)
    lower_limits[[1]] <- min(v.gamma.selection)
    upper_limits[[1]] <- max(v.gamma.selection)
  }
  
  # lambda - single value handling
  if(length(v.lambda.selection) == 1){
    idx.lambda <- which(lambda.gridSearch == v.lambda.selection)
    if(idx.lambda == 1){
      lower_limits[[2]] <- v.lambda.selection
      upper_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda + 1])
    }else if(idx.lambda == length(lambda.gridSearch)){
      lower_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda - 1])
      upper_limits[[2]] <- v.lambda.selection
    }else{
      lower_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda - 1])
      upper_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda + 1])
    }
    w[[2]] <- v.lambda.selection
  }else{
    w[[2]] <- min(v.lambda.selection)
    lower_limits[[2]] <- min(v.lambda.selection)
    upper_limits[[2]] <- max(v.lambda.selection)
  }
  
  #return(list(w=w,lower_limits=lower_limits,upper_limits=upper_limits))
  
  w.intials     <- w
  lower_limits <- lower_limits
  upper_limits <- upper_limits
  
  
  #   w.intials     <- lst.hyperparameter$w
  #   lower_limits <- lst.hyperparameter$lower_limits
  #   upper_limits <- lst.hyperparameter$upper_limits
  #   
  
  print(Sys.time()-strt)
  
  #   #     w     <- c(gamma.upper, 1.5) #0.5) # start at upper end - therefore, the one with smallest 
  #   #     lower <- c(gamma.lower, 2.0) #0.01) 
  #   #     upper <- c(gamma.upper, 2.5) # 1.0)
  #   
  ## Step2: Simulated annealing ##
  
  
  
  strt<-Sys.time() 
  
  # begin 
  evaluate_SA <- function(w) {
    
    v.alpha = v.gamma = w[1]; v.lambda= w[2]; 
    if(v.alpha >= v.gamma){
      lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields (lst.grace_support=lst.grace_support.bootstrap_sample, lst.modules=lst.modules.bootstrap_sample,
                                                                 v.alpha = v.alpha, v.gamma = v.gamma, v.lambda = v.lambda, b.prim=TRUE, tfs = tfs.mrf, tgs = tgs.mrf) 
      
      mat.p.crf <- lst.mat.p.crf$mat.p.crf
      
      if(!is.null(mat.p.crf) & !is.null(dim(mat.p.crf)) & sum(mat.p.crf) > 0){
        
        
        # use the fast, simpler function
        lst.eval <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.p.crf,
                                                                lst.benchmarks=lst.gs.grn_bp.train)
        
        #     lst.eval <- compute_fscore_model_assessment(mat.p.grace_ensemble=mat.p.crf,
        #                                                 lst.benchmarks=lst.gs.grn_bp.train,
        #                                                 th.min.coreg.tfs = th.min.coreg.tfs, 
        #                                                 b.jaccard = b.jaccard,
        #                                                 beta = beta)
        
        y1 <- (lst.eval$n.gs.grn/lst.eval$n.pred)
        y1[is.na(y1)] <- 0
        if(b.normalize_precision){
          y1 <- ratio.gs.grn_vs_bp.normalization * y1  
        }
        
        y2 <- lst.eval$n.coreg.bp_tg_pairs/max.coreg.bp
        y2[is.na(y2)] <- 0
        
        fscore_beta <- (1 + beta^2)*(y1*y2/((beta^2 * y1)+y2))
        fscore_beta[is.na(fscore_beta)] <- 0
        
        fscore <- fscore_beta
        n.gs.grn <- lst.eval$n.gs.grn
        n.coreg.bp.tgs <- lst.eval$n.coreg.bp_tg_pairs
        n.pred <- lst.eval$n.pred
        
      }else{
        fscore <- - 100
        n.gs.grn <- - 1
        n.coreg.bp.tgs <- - 1
        n.pred <- - 1
      }
      
    }else{
      fscore <- - 100
      n.gs.grn <- - 1
      n.coreg.bp.tgs <- - 1
      n.pred <- - 1
    }
    
    print(paste("alpha", v.alpha, ", gamma", v.gamma, ", lambda", v.lambda, ", f_score", fscore, "n.predictions", n.pred, "n.gs.grn",n.gs.grn ,"n.coreg.bp.tgs", n.coreg.bp.tgs))
    returnVal = fscore * -1 
    returnVal
  }
  ## end define SA
  
  out <- GenSA(w.intials, lower = lower_limits, upper = upper_limits, fn = evaluate_SA , control=list(max.call=max.call)) 
  #out <- GenSA(w, lower = lower, upper = upper, fn = evaluate_SA , control=list(max.call=max.call, smooth = FALSE, simple.function = TRUE))
  results <- out[c("value", "par", "counts")]
  weights <- results$par # go for smallest lambda
  print(Sys.time()-strt)
  
  
  #### construct model using learned weights
  alpha <- gamma <- weights[1]
  lambda <- weights[2]
  
  lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields (lst.grace_support=lst.grace_support.bootstrap_sample, lst.modules=lst.modules.bootstrap_sample,
                                                             v.alpha = alpha, v.gamma = gamma, v.lambda = lambda, b.prim=TRUE, tfs = tfs.mrf, tgs = tgs.mrf) 
  
  mat.p.crf <- lst.mat.p.crf$mat.p.crf
  
  if(!is.null(mat.p.crf) & !is.null(dim(mat.p.crf)) & sum(mat.p.crf) > 0){
    
    lst.gs.grn_bp <- list(truereg, mat.cofunctional_evidence)
    names(lst.gs.grn_bp) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")
    
    lst.eval <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.p.crf,
                                                            lst.benchmarks=lst.gs.grn_bp)
    
    #     lst.eval <- compute_fscore_model_assessment(mat.p.grace_ensemble=mat.p.crf,
    #                                                 lst.benchmarks=lst.gs.grn_bp,  
    #                                                 th.min.coreg.tfs = th.min.coreg.tfs, 
    #                                                 b.jaccard = b.jaccard,
    #                                                 beta = beta)
    
    
    y1 <- (lst.eval$n.gs.grn/lst.eval$n.pred)
    y1[is.na(y1)] <- 0
    if(b.normalize_precision){
      y1 <- ratio.gs.grn_vs_bp.normalization * y1  
    }
    
    y2 <- lst.eval$n.coreg.bp_tg_pairs/max.coreg.bp
    y2[is.na(y2)] <- 0
    
    fscore_beta <- (1 + beta^2)*(y1*y2/((beta^2 * y1)+y2))
    fscore_beta[is.na(fscore_beta)] <- 0
    
    fscore <- fscore_beta
    n.gs.grn <- lst.eval$n.gs.grn
    n.coreg.bp.tgs <- lst.eval$n.coreg.bp_tg_pairs
    n.pred <- lst.eval$n.pred
    
    print(paste("alpha", alpha, ", gamma", gamma, ", lambda", lambda, ", f_score", fscore, "n.predictions", n.pred, "n.gs.grn",n.gs.grn ,"n.coreg.bp.tgs", n.coreg.bp.tgs))
    
    
    lst.grace_models <- list(weights=c(alpha, gamma, lambda), fscore=fscore, n.pred=n.pred,n.gs.grn=n.gs.grn,n.coreg.bp.tgs=n.coreg.bp.tgs)
  }else{
    lst.grace_models <- list(weights=c(-1, -1, -1), fscore=-100, n.pred=-1,n.gs.grn=-1,n.coreg.bp.tgs=-1)  
  }
  
  return(lst.grace_models)
}

#修改删除模拟退火SA过程
evaluate_mode_parameter1 <- function(lst.grace_support=lst.grace_support, 
                                     lst.modules=lst.modules, 
                                     truereg=truereg, 
                                     mat.cofunctional_evidence=mat.cofunctional_evidence, 
                                     
                                     gamma.gridSearch=gamma.gridSearch,
                                     lambda.gridSearch=lambda.gridSearch,
                                     th.percentage.hyperparameters=th.percentage.hyperparameters,
                                     
                                     ratio.gs.grn_vs_bp.normalization = ratio.gs.grn_vs_bp.normalization,
                                     max.coreg.bp = max.coreg.bp,
                                     
                                     n.sample_size=n.sample_size, 
                                     max.call = 100, 
                                     v.seed = 1234, 
                                     beta = 1.0){
  
  # source("GRACE_optimization_routines.R")
  # source("GRACE.R")
  # source("GRACE_helperfunctions.R")
  
  #### Bootstrapping (based on regulators based meta module sampling (0.632))
  
  # set.seed(v.seed)
  
  tfs.mrf <- unique(lst.grace_support$df.grn$TF) 
  tgs.mrf <- unique(lst.grace_support$df.grn$TG)
  
  lst.bootstrap_sample <- generate_bootstrapping_sample(lst.grace_support=lst.grace_support,lst.modules=lst.modules, truereg=truereg, mat.cofunctional_evidence=mat.cofunctional_evidence, n.sample_size = n.sample_size)
  lst.grace_support.bootstrap_sample <- lst.bootstrap_sample$lst.grace_support.bootstrap_sample
  lst.modules.bootstrap_sample <- lst.bootstrap_sample$lst.modules.bootstrap_sample
  
  lst.gs.grn_bp.train <- list(lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.GS.bootstrap_sample, lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.cofunctional_evidence.bootstrap_sample)    #list(mat.GS.train, mat.BP.coreg.train)
  names(lst.gs.grn_bp.train) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")
  
  ## Step 1: Grid Search (Coarse) optimization to identify finer range for Simulated annealing optimization ##
  strt<-Sys.time() 
  #lst.hyperparameter <- gridSearchOptimization(lst.gs.grn_bp.train=lst.gs.grn_bp.train,gamma.gridSearch, lambda.gridSearch, th.percentage.hyperparameters, ratio.gs.grn_vs_bp.normalization = ratio.gs.grn_vs_bp.normalization, max.coreg.bp = max.coreg.bp)
  
  # add begin grid search
  df.gridSearch <- expand.grid(gamma.gridSearch, lambda.gridSearch)
  names(df.gridSearch) <- c("gamma", "lambda")
  
  df.gridSearch["f_score"] <- 0
  df.gridSearch["n.predictions"] <- 0
  df.gridSearch["n.gs.grn"] <- 0
  df.gridSearch["n.coreg.bp.tgs"] <- 0
  
  for(j in 1:nrow(df.gridSearch)){
    
    cat("Processing... ", round(j/nrow(df.gridSearch) * 100, digits = 2) , "%", "\r"); flush.console()  
    v.alpha = v.gamma = df.gridSearch$gamma[j]; v.lambda= df.gridSearch$lambda[j];
    
    if(v.alpha >= v.gamma){
      
      lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields (lst.grace_support=lst.grace_support.bootstrap_sample, lst.modules=lst.modules.bootstrap_sample,
                                                                 v.alpha = v.alpha, v.gamma = v.gamma, v.lambda = v.lambda, b.prim=TRUE, tfs = tfs.mrf, tgs = tgs.mrf) 
      
      mat.p.crf <- lst.mat.p.crf$mat.p.crf # get if equal to zero ... , could be based on bootstrap sample, no value in there... 
      
      if(!is.null(mat.p.crf) & !is.null(dim(mat.p.crf)) & sum(mat.p.crf) > 0){
        
        lst.eval <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.p.crf,
                                                                lst.benchmarks=lst.gs.grn_bp.train)
        
        
        y1 <- (lst.eval$n.gs.grn/lst.eval$n.pred)
        y1[is.na(y1)] <- 0
        if(b.normalize_precision){
          y1 <- ratio.gs.grn_vs_bp.normalization * y1  
        }
        
        y2 <- lst.eval$n.coreg.bp_tg_pairs/max.coreg.bp
        y2[is.na(y2)] <- 0
        
        fscore_beta <- (1 + beta^2)*(y1*y2/((beta^2 * y1)+y2))
        fscore_beta[is.na(fscore_beta)] <- 0
        
        fscore <- fscore_beta
        n.gs.grn <- lst.eval$n.gs.grn
        n.coreg.bp.tgs <- lst.eval$n.coreg.bp_tg_pairs
        n.pred <- lst.eval$n.pred
        
      }else{
        fscore <- - 100
        n.gs.grn <- - 1
        n.coreg.bp.tgs <- - 1
        n.pred <- - 1
      }
      
    }else{
      fscore <- - 100
      n.gs.grn <- - 1
      n.coreg.bp.tgs <- - 1
      n.pred <- - 1
    }
    
    df.gridSearch$f_score[j] <- fscore
    df.gridSearch$n.predictions[j] <- n.pred
    df.gridSearch$n.gs.grn[j] <- n.gs.grn
    df.gridSearch$n.coreg.bp.tgs[j] <- n.coreg.bp.tgs
    #print(paste("alpha", v.alpha, ", gamma", v.gamma, ", lambda", v.lambda, ", f_score", fscore, "n.predictions", n.pred, "n.gs.grn",n.gs.grn ,"n.coreg.bp.tgs", n.coreg.bp.tgs))
  }
  
  # return Simulated annealing range 
  v.fscore.selection <- quantile(df.gridSearch$f_score, th.percentage.hyperparameters)
  df.gridSearch <- subset(df.gridSearch, df.gridSearch$f_score >= v.fscore.selection)
  
  v.double_fscores <- names(which(table(df.gridSearch$f_score) > 1))
  if(length(v.double_fscores) > 0){
    for(l in 1:length(v.double_fscores)){
      v.lambda.tmp <- subset(df.gridSearch, df.gridSearch$f_score == v.double_fscores[l])$lambda
      v.lambda.tmp <- v.lambda.tmp[!v.lambda.tmp %in% min(v.lambda.tmp)]
      
      idx.remove <- which(df.gridSearch$f_score == v.double_fscores[l] & df.gridSearch$lambda %in% v.lambda.tmp)
      df.gridSearch <- df.gridSearch[-idx.remove,]
    }
  }
  
  v.gamma.selection <- unique(df.gridSearch$gamma)
  v.lambda.selection <- unique(df.gridSearch$lambda)
  
  w  <- numeric(2)
  lower_limits <- numeric(2)
  upper_limits <- numeric(2)
  
  # gamma - single value handling
  if(length(v.gamma.selection) == 1){
    idx.gamma <- which(gamma.gridSearch == v.gamma.selection)
    if(idx.gamma == 1){
      lower_limits[[1]] <- v.gamma.selection
      upper_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma - 1]) # inverse indexing because of v.grn ordering
    }else if(idx.gamma == length(gamma.gridSearch)){
      lower_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma + 1])
      upper_limits[[1]] <- v.gamma.selection
    }else{
      lower_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma + 1])
      upper_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma - 1])
    }
    w[[1]] <- v.gamma.selection
  }else{
    w[[1]] <- min(v.gamma.selection)
    lower_limits[[1]] <- min(v.gamma.selection)
    upper_limits[[1]] <- max(v.gamma.selection)
  }
  
  # lambda - single value handling
  if(length(v.lambda.selection) == 1){
    idx.lambda <- which(lambda.gridSearch == v.lambda.selection)
    if(idx.lambda == 1){
      lower_limits[[2]] <- v.lambda.selection
      upper_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda + 1])
    }else if(idx.lambda == length(lambda.gridSearch)){
      lower_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda - 1])
      upper_limits[[2]] <- v.lambda.selection
    }else{
      lower_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda - 1])
      upper_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda + 1])
    }
    w[[2]] <- v.lambda.selection
  }else{
    w[[2]] <- min(v.lambda.selection)
    lower_limits[[2]] <- min(v.lambda.selection)
    upper_limits[[2]] <- max(v.lambda.selection)
  }
  alpha <- gamma <- w[1]
  lambda <- w[2]
  #return(list(w=w,lower_limits=lower_limits,upper_limits=upper_limits))
  
  w.intials     <- w
  lower_limits <- lower_limits
  upper_limits <- upper_limits
  
  
  #   w.intials     <- lst.hyperparameter$w
  #   lower_limits <- lst.hyperparameter$lower_limits
  #   upper_limits <- lst.hyperparameter$upper_limits
  #   
  
  print(Sys.time()-strt)
  
  #   #     w     <- c(gamma.upper, 1.5) #0.5) # start at upper end - therefore, the one with smallest 
  #   #     lower <- c(gamma.lower, 2.0) #0.01) 
  #   #     upper <- c(gamma.upper, 2.5) # 1.0)
  #   
  ## Step2: Simulated annealing ##
  
  
  
  
  lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields (lst.grace_support=lst.grace_support.bootstrap_sample, lst.modules=lst.modules.bootstrap_sample,
                                                             v.alpha = alpha, v.gamma = gamma, v.lambda = lambda, b.prim=TRUE, tfs = tfs.mrf, tgs = tgs.mrf) 
  
  mat.p.crf <- lst.mat.p.crf$mat.p.crf
  
  if(!is.null(mat.p.crf) & !is.null(dim(mat.p.crf)) & sum(mat.p.crf) > 0){
    
    lst.gs.grn_bp <- list(truereg, mat.cofunctional_evidence)
    names(lst.gs.grn_bp) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")
    
    lst.eval <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.p.crf,
                                                            lst.benchmarks=lst.gs.grn_bp)
    
    #     lst.eval <- compute_fscore_model_assessment(mat.p.grace_ensemble=mat.p.crf,
    #                                                 lst.benchmarks=lst.gs.grn_bp,  
    #                                                 th.min.coreg.tfs = th.min.coreg.tfs, 
    #                                                 b.jaccard = b.jaccard,
    #                                                 beta = beta)
    
    
    y1 <- (lst.eval$n.gs.grn/lst.eval$n.pred)
    y1[is.na(y1)] <- 0
    if(b.normalize_precision){
      y1 <- ratio.gs.grn_vs_bp.normalization * y1  
    }
    
    y2 <- lst.eval$n.coreg.bp_tg_pairs/max.coreg.bp
    y2[is.na(y2)] <- 0
    
    fscore_beta <- (1 + beta^2)*(y1*y2/((beta^2 * y1)+y2))
    fscore_beta[is.na(fscore_beta)] <- 0
    
    fscore <- fscore_beta
    n.gs.grn <- lst.eval$n.gs.grn
    n.coreg.bp.tgs <- lst.eval$n.coreg.bp_tg_pairs
    n.pred <- lst.eval$n.pred
    
    print(paste("alpha", alpha, ", gamma", gamma, ", lambda", lambda, ", f_score", fscore, "n.predictions", n.pred, "n.gs.grn",n.gs.grn ,"n.coreg.bp.tgs", n.coreg.bp.tgs))
    
    
    lst.grace_models <- list(weights=c(alpha, gamma, lambda), fscore=fscore, n.pred=n.pred,n.gs.grn=n.gs.grn,n.coreg.bp.tgs=n.coreg.bp.tgs)
  }else{
    lst.grace_models <- list(weights=c(-1, -1, -1), fscore=-100, n.pred=-1,n.gs.grn=-1,n.coreg.bp.tgs=-1)  
  }
  
  return(lst.grace_models)
}

#GRACE
generate_bootstrapping_sample <- function(lst.grace_support=lst.grace_support,lst.modules=lst.modules, truereg=truereg, mat.cofunctional_evidence=mat.cofunctional_evidence, n.sample_size = 0.632, b.validation = FALSE){
  # browser()
  mat.GS.bootstrap_sample <- truereg
  df.GS <- as.data.frame(which(mat.GS.bootstrap_sample == 1, arr.ind = TRUE))
  names(df.GS) <- c("TF", "TG")
  idx.tfs.gs <- unique(df.GS$TF)
  v.tfs.gs <- rownames(mat.GS.bootstrap_sample)[idx.tfs.gs]
  
  # sample from gold standard
  v.filter <- sample(nrow(df.GS), round(nrow(df.GS) * (n.sample_size)))
  
  if(b.validation){
    # select the complement for validation 
    v.filter.train <- v.filter
    v.filter <- seq(1:nrow(df.GS))[!seq(1:nrow(df.GS)) %in% v.filter.train]
  }
  
  # sample of the gold standard transcription factors (needs to be non-overlapping separated into training and testing)
  v.tfs.gs.bootstrap_sample <- unique(df.GS$TF[v.filter])
  v.tfs.gs.bootstrap_sample <- rownames(mat.GS.bootstrap_sample)[v.tfs.gs.bootstrap_sample]
  
  # sample from data (then replace ) - this can be identical to the training, as this is a semisupervised setting
  v.tfs <- (unique(lst.grace_support$df.grn$TF))
  names(v.tfs) <- unique(lst.grace_support$df.grn$tf.id)
  n.tfs <- length(v.tfs)
  v.tfs.ids <- unique(lst.grace_support$df.grn$tf.id)
  names(v.tfs.ids) <- (unique(lst.grace_support$df.grn$TF))
  
  n.tfs.bootstrap_sample <- round(n.tfs * n.sample_size, digits = 0)
  n.tfs.bootstrap_sample <- n.tfs.bootstrap_sample - length(v.tfs.gs.bootstrap_sample)
  
  v.tfs.bootstrap_sample.gs_plus <- sample(v.tfs[!v.tfs %in% v.tfs.gs], n.tfs.bootstrap_sample, replace = FALSE)
  n.tfs.bootstrap_sample.gs_plus <- length(v.tfs.bootstrap_sample.gs_plus)
  
  # >= 0.632 of features (regulators) - final selection
  v.tfs.bootstrap_sample <- unique(c(v.tfs.bootstrap_sample.gs_plus, v.tfs.gs.bootstrap_sample))
  v.tfs.filter <- v.tfs[!v.tfs %in% v.tfs.bootstrap_sample]
  
  mat.GS.bootstrap_sample[v.tfs.filter,] <- 0 # this might be re-designed
  
  # sample datasets using regulator bootstrap samples (adjust targets for max coreg bp pairs computation)
  lst.grace_support.bootstrap_sample <- lst.grace_support
  lst.grace_support.bootstrap_sample$df.grn <- subset(lst.grace_support.bootstrap_sample$df.grn, lst.grace_support.bootstrap_sample$df.grn$TF %in% v.tfs.bootstrap_sample)
  
  v.tgs.bootstrap_sample <- unique(lst.grace_support.bootstrap_sample$df.grn$TG)
  v.tgs.cofunctional_evidence <- colnames(mat.cofunctional_evidence)
  v.tgs.filter <- v.tgs.cofunctional_evidence[!v.tgs.cofunctional_evidence %in% v.tgs.bootstrap_sample]
  
  mat.cofunctional_evidence.bootstrap_sample <- mat.cofunctional_evidence
  if(length(v.tgs.filter)>0){
    mat.cofunctional_evidence.bootstrap_sample[v.tgs.filter,v.tgs.filter] <- 0
  }
  ###
  
  v.link.bootstrap_sample <- lst.grace_support.bootstrap_sample$df.grn$link.id
  lst.grace_support.bootstrap_sample$adj.links <- lst.grace_support.bootstrap_sample$adj.links[v.link.bootstrap_sample,v.link.bootstrap_sample]
  
  ## more regulators than regulators?
  lst.modules.bootstrap_sample <- lst.modules
  idx.tfs_ids.bootstrap_sample <- unlist(lapply(lst.modules.bootstrap_sample$lst.tf.ids, function(m) unlist(m[[1]])))
  
  
  idx.tfs_ids.bootstrap_sample <- which(idx.tfs_ids.bootstrap_sample %in% as.numeric(v.tfs.ids[v.tfs.bootstrap_sample]))
  
  lst.modules.bootstrap_sample$v.coreg_modules <- lst.modules$v.coreg_modules[idx.tfs_ids.bootstrap_sample]
  lst.modules.bootstrap_sample$lst.crfs <- lst.modules$lst.crfs[idx.tfs_ids.bootstrap_sample]
  lst.modules.bootstrap_sample$lst.crf.prim_based_links <- lst.modules$lst.crf.prim_based_links[idx.tfs_ids.bootstrap_sample]
  lst.modules.bootstrap_sample$lst.v.grn.node <- lst.modules$lst.v.grn.node[idx.tfs_ids.bootstrap_sample]
  lst.modules.bootstrap_sample$lst.tf.ids <- lst.modules$lst.tf.ids[idx.tfs_ids.bootstrap_sample]
  lst.modules.bootstrap_sample$lst.tg.ids <- lst.modules$lst.tg.ids[idx.tfs_ids.bootstrap_sample]
  
  
  return(list(lst.grace_support.bootstrap_sample=lst.grace_support.bootstrap_sample, lst.modules.bootstrap_sample=lst.modules.bootstrap_sample, 
              lst.benchmark.bootstrap_sample=list(mat.GS.bootstrap_sample=mat.GS.bootstrap_sample, mat.cofunctional_evidence.bootstrap_sample=mat.cofunctional_evidence.bootstrap_sample)))
}

#GRACE.R
reevaluate_links_with_MarkovRandomFields <- function(lst.grace_support=lst.grace_support, lst.modules=lst.modules,
                                                     v.alpha=v.alpha, v.gamma=v.gamma, v.lambda=v.lambda,  b.prim = TRUE,
                                                     tfs = unique(lst.grace_support$df.grn$TF), tgs = unique(lst.grace_support$df.grn$TG)){
  
  th.gamma <- v.gamma
  th.alpha <- v.alpha
  
  df.grn <- lst.grace_support$df.grn
  adj.links <- lst.grace_support$adj.links
  th.preselection <- lst.grace_support$v.grn.significant 
  
  g.coreg <- lst.modules$g.coreg
  v.coreg_modules <- lst.modules$v.coreg_modules
  n.crfs <- length(v.coreg_modules)
  lst.crfs.prim <- lst.modules$lst.crfs
  lst.crf.prim_based_links <- lst.modules$lst.crf.prim_based_links
  lst.v.grn.node <- lst.modules$lst.v.grn.node
  
  lst.tf.ids <- lst.modules$lst.tf.ids
  lst.tg.ids <- lst.modules$lst.tg.ids
  
  #   mat.link_id <- acast(df.grn, TF~TG, value.var="link.id")
  #   mat.link_id[is.na(mat.link_id)] <- 0;  
  #   mat.link_id <- as(mat.link_id, "CsparseMatrix")
  #   
  #   ### IMPORTANT - Reorder tfs and tgs in all relevant matrices
  #   mat.link_id <- mat.link_id[tfs, tgs]
  #   
  # v.links.grn <- df.grn$link.id
  
  ## singleton nodes                         
  v.links.singletons <- rownames(adj.links)[which(Matrix::rowSums(adj.links) == 0)] 
  
  ## remove singleton nodes
  adj.links.nonsingleton <- adj.links[!rownames(adj.links) %in% v.links.singletons, !colnames(adj.links) %in% v.links.singletons]
  
  clust.coreg <- clusters(g.coreg) # igraph based extraction of connected components s
  
  lst.crfs <- vector(mode = "list", length = n.crfs)
  for(i in 1:length(v.coreg_modules)){ 
    
    v.links <- names(which(clust.coreg$membership == v.coreg_modules[i]))
    adj.links.set <- adj.links.nonsingleton[v.links,v.links]
    
    v.links %in% rownames(adj.links.nonsingleton)
    v.links %in% rownames(adj.links)
    
    if(b.prim){
      crf <- lst.crfs.prim[[i]]
      v.edge_pots <- exp(v.lambda * lst.crf.prim_based_links[[i]]) # ru.adj.links[crf$edges])
    } else{
      crf <- make.crf(adj.links.set, n.states = 2)
      v.edge_pots <- exp(v.lambda * adj.links.set[crf$edges])
    }
    ## edge potentials
    for(j in 1:crf$n.edges){ # only lambda
      crf$edge.pot[[j]] <- matrix(c(v.edge_pots[j], 1, 1, v.edge_pots[j]), ncol = 2, byrow = FALSE)
    } 
    lst.crfs[[i]] <- crf
  }
  
  ### Compute MRF
  if(length(v.links.singletons) > 0){
    df.grn.singletons <- subset(df.grn, df.grn$link.id %in% v.links.singletons)
    #mat.p.crf_precomputed <- acast(df.grn.singletons, TF~TG, value.var = "rank.global")
    mat.p.crf_singleton <- acast(df.grn.singletons, TF~TG, value.var = "v.grn")
    mat.p.crf_singleton[is.na(mat.p.crf_singleton)] <- 0
    mat.p.crf_singleton <- exp(mat.p.crf_singleton - th.gamma)
    mat.p.crf_singleton[mat.p.crf_singleton <= 1] <- 0
    mat.p.crf_singleton[mat.p.crf_singleton > 1] <- 1
    #mat.p.crf_singleton <- as(mat.p.crf_singleton, "CsparseMatrix")
  }
  
  
  mat.p.crf <- matrix(0, nrow = length(tfs), ncol = length(tgs), dimnames = list(tfs, tgs)) 
  v.counter.crf <- numeric(length(v.coreg_modules))
  
  for(i in 1:length(v.coreg_modules)){ 
    v.grn.node <- lst.v.grn.node[[i]]
    crf <- lst.crfs[[i]]
    crf$node.pot[,1] <- exp(v.grn.node - th.gamma)
    crf$node.pot[,2] <- 1 # mean normalized to 0
    
    idx.gamma <- which(v.grn.node > th.alpha)
    #idx.gamma <- which(crf$node.pot[,1] > 1) # everything > 0.5 is automatically included (p = 1)
    clamped <- rep(0, length(crf$node.pot[,1]))
    clamped[idx.gamma] <- 1  
    
    if(!all(clamped == 1)){
      #infer.conditional配分函数和边际概率的计算
      p.crf <- infer.conditional(crf, clamped, infer.tree)$node.bel[,1] # infer marginals / node beliefs 
    }else{
      v.counter.crf[i] <- 1
      p.crf <- rep(1.0, length(clamped))
    }
    
    p.crf <- ifelse(p.crf > 0.5, 1, 0) 
    tf.ids <- lst.tf.ids[[i]]
    tg.ids <- lst.tg.ids[[i]]
    ind = (tg.ids-1)*nrow(mat.p.crf) + tf.ids
    mat.p.crf[ind] <- p.crf
    
  }
  
  if(length(v.links.singletons) > 0){
    tf.pre <- rownames(mat.p.crf_singleton)
    tg.pre <- colnames(mat.p.crf_singleton)
    mat.p.crf[tf.pre, tg.pre] <- mat.p.crf[tf.pre, tg.pre] + mat.p.crf_singleton
  }
  
  if(all(v.counter.crf == 1)){
    print("Note: all links to be re-evaluated based on functional association data already in pre-selection, no changes")
  }
  
  return(list(mat.p.crf=mat.p.crf, th.gamma=th.gamma))
}
