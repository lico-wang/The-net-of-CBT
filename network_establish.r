library(igraph)
library(cellchat)
library(CellCall)
chat_TF1=read.table('./new_ligand_receptor_TFs.txt',header = T)
chat_TF2=read.table('./new_ligand_receptor_TFs_extended.txt',header = T)
chat_TF3=read.table('./new_ligand_receptor_TFs_homology.txt',header = T)
chat_TF4=read.table('./new_ligand_receptor_TFs_homology_extended.txt',header = T)
chat_TFall=rbind(chat_TF1,chat_TF2)
chat_TFall=rbind(chat_TFall,chat_TF3)
chat_TFall=rbind(chat_TFall,chat_TF4)
TF_TG1=read.table('./tf_target.txt',header = T)
TF_TG2=read.table('./tf_target_homology.txt',header = T)
TF_TGall=rbind(TF_TG1,TF_TG2)

ligand_target_matrix <- readRDS("./ligand_target_matrix.rds")
ligand_target_matrix=ligand_target_matrix[rownames(ligand_target_matrix) %in% unique(TF_TGall$TF_Symbol),]
load('net.Rdata')
net_target=ligand_target_matrix[,colnames(ligand_target_matrix) %in% unique(net$ligand)]

#SCENIC算出来的转录因子活跃度与net_target中存在的转录因子做对应
regu=row.names(regulonAUC)
regul=unlist(lapply(regu, function(x){substr(x,1,nchar(x)-3)}))
row.names(regulonAUC)=regul
celltype='interested cell'
net_cell=row.names(meta0M)[meta0M$Celltype %in% celltype]
net_reguAUC=regulonAUC[,colnames(regulonAUC) %in% net_cell]
net_reguAUC=net_reguAUC[row.names(net_reguAUC) %in% row.names(net_target),]
net_target=net_target[row.names(net_target) %in% row.names(net_reguAUC),]
#基线转录因子活跃度
baseAUC0=rowMeans(regulonAUC@assays@data$AUC)
netAUC0=rowMeans(net_reguAUC@assays@data$AUC)
netAUC0=as.data.frame(netAUC0);baseAUC0=as.data.frame(baseAUC0)
#根据SCENIC确定转录因子活性
baseAUC0$TF=row.names(baseAUC0)
baseAUC0=baseAUC0[match(row.names(netAUC0),baseAUC0$TF),]
netAUC0$exist=FALSE
for (i in 1:nrow(netAUC0)){
  if (netAUC0$netAUC0[i]>=baseAUC0$baseAUC0[i]){
    netAUC0$exist[i]=TRUE
  }
}
netAUC0=netAUC0[netAUC0$exist,]
remain_TF=unlist(lapply(row.names(netAUC0), function(x){substr(x,1,nchar(x)-3)}))
net_target=net_target[row.names(net_target) %in% remain_TF,]
#根据每条通讯机制下的转录因子的大小确定转录因子的归属
threshold=as.data.frame(apply(net, 1,function(x){quantile(x,0.9)}))
net_target_remain=nettarget
nettarget_remain[,]=FALSE
for (i in 1:nrow(net_target)){
  for (j in 1:ncol(net_target)){
      if(net_target[i,j]>=threshold[i,1]){
        net_target_remain[i,j]=TRUE
      }
  }}
#将得到的ligand - TF 找到对应的TG 再根据regulons_incidMat0中的TG的基线活跃度
#确定是否被调控 阳性的是1
net_TG=TF_TGall[TF_TGall$TF_Symbol %in%  row.names(net_target_remain),]
TGnames=unlist(lapply(row.names(regulons_incidMat0), function(x){substr(x,1,nchar(x)-3)}))
row.names(regulons_incidMat0)=TGnames
net_incide=regulons_incidMat0[row.names(regulons_incidMat0) %in% row.names(net_target_remain),]
net_incide=as.data.frame(net_incide)
net_TG=net_TG[net_TG$Target_Symbol %in% colnames(regulons_incidMat0),]
net_TG$exist=FALSE
for (i in 1:nrow(net_TG)){
  tmptf=net_TG$TF_Symbol[i]
  tmptg=net_TG$Target_Symbol[i]
  if (net_incide[row.names(net_incide) %in% tmptf,colnames(net_incide) %in% tmptg]==1){
    net_TG$exist[i]=TRUE
  }
}
net_TG=net_TG[net_TG$exist,]
net0M$s_node=paste(net0M$ligand,net0M$source,sep='_')
net0M$t_node=paste(net0M$receptor,net0M$ligand,sep='_')
edge1=net0M[,colnames(net0M) %in% c('s_node','t_node','prob')]
edge1=edge1[,c(2,3,1)]
sourcenode=edge1$s_node
sourcenode=data.frame(node=sourcenode,Z=rep(4,length(sourcenode)))
targetnode=edge1$t_node
targetnode=data.frame(node=targetnode,Z=rep(3,length(targetnode)))
node1=rbind(sourcenode,targetnode)
node1=node1[!duplicated(node1),]
cc_TF=net0M[,c('ligand','t_node')]
cc_TF=cc_TF[!duplicated(cc_TF),]
C_TF=data.frame(t_node=rep(cc_TF$t_node[1],length(net0M_target_remain[,cc_TF$ligand[1]][net0M_target_remain[,cc_TF$ligand[1]]>0])),
                    TF=names(net0M_target_remain[,cc_TF$ligand[1]][net0M_target_remain[,cc_TF$ligand[1]]>0]))
for(i in 2:nrow(cc_TF)){
  tmp_C_TF=data.frame(t_node=rep(cc_TF$t_node[i],length(net0M_target_remain[,cc_TF$ligand[i]][net0M_target_remain[,cc_TF$ligand[i]]>0])),
  TF=names(net0M_target_remain[,cc_TF$ligand[i]][net0M_target_remain[,cc_TF$ligand[i]]>0]))
  C_TF=rbind(C_TF,tmp_C_TF)
}
edge2=C_TF
node2=data.frame(node=edge2$TF,Z=rep(2,nrow(edge2)))
node2=node2[!duplicated(node2),]
edge3=net0M_TG[,c(1,2)]
node3=data.frame(node=edge3$Target_Symbol,Z=rep(1,nrow(edge3)))
node3=node3[!duplicated(node3),]
edge1_2=edge1[,c(1,2)]
colnames(edge1_2)=c('source','target')
edge1_2$weight=rep(4,nrow(edge1_2))
edge2_2=edge2[,c(1,2)]
colnames(edge2_2)=c('source','target')
edge2_2$weight=rep(3,nrow(edge2_2))
colnames(edge3)=c('source','target')
edge3$weight=rep(2,nrow(edge3))
TG_edge=read.table('./string_interactions_short.tsv')
TG_edge=TG_edge[,c(1,2)]
colnames(TG_edge)=c('source','target')
TG_edge$weight=rep(1,nrow(TG_edge))
TF_edge=read.table('./string_interactions_short2.tsv')
TF_edge=TF_edge[,c(1,2)]
colnames(TF_edge)=c('source','target')
TF_edge$weight=rep(3,nrow(TF_edge))
edgenetT0M=rbind(edge1_2,edge2_2)
edgenetT0M=rbind(edgenetT0M,edge3)
edgenetT0M=rbind(edgenetT0M,TF_edge)
edgenetT0M=rbind(edgenetT0M,TG_edge)
edgenetT0M=edgenetT0M[!duplicated(edgenetT0M),]
edgenetT0M$ID=c(1:nrow(edgenetT0M))
edgenetT0M=edgenetT0M[,c(4,1,2,3)]

nodenetT0M=rbind(node1,node2)
nodenetT0M$label=nodenetT0M$node
node3$label=rep('',nrow(node3))
nodenetT0M=rbind(nodenetT0M,node3)
nodenetT0M$ID=nodenetT0M$node
nodenetT0M=nodenetT0M[,c(4,1,2,3)]
nodenetT0M=nodenetT0M[!duplicated(nodenetT0M$node),]

diffnode=unique(c(edgenetT0M$source,edgenetT0M$target))[! unique(c(edgenetT0M$source,edgenetT0M$target)) %in% nodenetT0M$node]
edgenetT0M=edgenetT0M[!edgenetT0M$source %in% diffnode,]
edgenetT0M=edgeCnetT0M[!edgenetT0M$target %in% diffnode,]
edges=edgenetT0M[,c(2,3)]
colnames(edges)=c('from','to')
nodes=nodenetT0M[,c(2,3)]
net=graph_from_data_frame(d=edges,
              vertices=nodes,directed = T)

