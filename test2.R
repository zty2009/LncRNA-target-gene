

lncrna.feature=read.table('lncrna.feature.txt',stringsAsFactor=F)
for (i in 1:max(lncrna.feature[,1])){
	index=which(lncrna.feature[,1]==i)
	lncrna.feature[index,2]=lncrna.feature[index,2]/max(lncrna.feature[index,2])
	}
lncrna.feature[,1]=lncrna.feature[,1]/max(lncrna.feature[,1])

gene.feature=read.table('gene.feature.txt',stringsAsFactor=F)
gene.net=read.table('gene-net1.txt',stringsAsFactor=F)
for (i in 1:max(gene.feature[,1])){
	index=which(gene.feature[,1]==i)
	gene.feature[index,2]=gene.feature[index,2]/max(gene.feature[index,2])
	}
gene.feature[,1]=gene.feature[,1]/max(gene.feature[,1])
##############GCN layer############
alllncrna=read.table('all.lncrna.txt',stringsAsFactor=F)
allgene=read.table('allgene.txt',stringsAsFactor=F)
lncrna.gene.net=read.table('lncrna-gene.net.txt',stringsAsFactor=F)
gene.net=gene.net+diag(1,nrow(gene.net),nrow(gene.net))
D=colSums(gene.net)
D_hat=diag(D^(-0.5))
gene.gcn=D_hat%*%as.matrix(gene.net)%*%D_hat%*%as.matrix(gene.feature)


pos=read.table('lncrna_target.txt',fill=T,stringsAsFactor=F)
posdata=matrix(0,nrow(pos),ncol(gene.gcn)*2+2)
title=matrix(0,nrow(pos),2)
for (i in 1:nrow(pos)){
	index1=which(pos[i,3]==alllncrna)
	index2=which(pos[i,4]==allgene)
	if (length(c(index1,index2))==2){
	posdata[i,]=unlist(c(lncrna.feature[index1,],lncrna.gene.net[index1,index2],gene.gcn[index2,],lncrna.gene.net[index1,index2]))
	title[i,]=unlist(c(pos[i,3],pos[i,4]))
	}
	}
index=which(posdata[,1]!=0)
posdata=posdata[index,]
title=title[index,]
write.table(file='posdata2.txt',posdata,append=T,row.names=F,col.names=F,quote=F)
write.table(file='postitle2.txt',title,append=T,row.names=F,col.names=F,quote=F)

for (j in 1:nrow(posdata)){
	set.seed(j+9)
	lnsample1=sample.int(nrow(alllncrna), size = 1)
	gene=pos[which(pos[,1]==alllncrna[lnsample1,1]),2]
	index=apply(as.matrix(gene),1,function(x){which(allgene==x)})
	index=unlist(index)
	s=1:nrow(allgene)
	if (length(index)!=0){
		s=s[-index]
		}
	set.seed(j+9)
	lnsample2=sample(x=s, size = 1)
	neg=unlist(c(alllncrna[lnsample1,1],allgene[lnsample2,1],lncrna.feature[lnsample1,],lncrna.gene.net[lnsample1,lnsample2],gene.gcn[lnsample2,],lncrna.gene.net[lnsample1,lnsample2]))
	dim(neg)=c(1,36)
	write.table(file='neg.2.10.txt',neg,append=T,row.names=F,col.names=F,quote=F)
	}


library(keras)
negdata=read.table('neg.2.10.txt',stringsAsFactor=F)
negdata=as.matrix(negdata)

title.neg=negdata[,1:2]
negdata=negdata[,3:ncol(negdata)]

negindex=cut(sample(1:nrow(negdata)),10, labels=F)
posindex=cut(sample(1:nrow(posdata)),10, labels=F)



for (kk in 1:10){

	test.index1=which(posindex==kk)
	test.index2=which(negindex==kk)
	train.index1=which(posindex!=kk)
	train.index2=which(negindex!=kk)

	x_train=rbind(posdata[train.index1,],negdata[train.index2,])
	pp=nrow(x_train)
	dim(x_train)=c(pp,2,ncol(lncrna.feature)+1,1)
	y_train=to_categorical(c(rep(1,length(train.index1)),rep(0,length(train.index2))))
	
	x_test=rbind(posdata[test.index1,],negdata[test.index2,])
	y_test=to_categorical(c(rep(1,length(test.index1)),rep(0,length(test.index2))))
	ti=rbind(title[test.index1,],title.neg[test.index2,])
	label=c(rep(1,length(test.index1)),rep(0,length(test.index2)))
	ti=cbind(ti,label)
	pp=nrow(x_test)
	dim(x_test)=c(pp,2,ncol(lncrna.feature)+1,1)

	
	
model <- keras_model_sequential() 
model %>% 

  layer_conv_2d(filters = 64, kernel_size = c(1, 4), activation = "tanh", input_shape = c(2,ncol(lncrna.feature)+1,1)) %>% 
  layer_max_pooling_2d(pool_size = c(1, 2)) %>% 
  layer_dropout(rate = 0.25) %>% 

  layer_conv_2d(filters = 128, kernel_size = c(2, 2), activation = "tanh") %>% 
  layer_max_pooling_2d(pool_size = c(1, 2)) %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_conv_2d(filters = 256, kernel_size = c(1, 2), activation = "tanh") %>% 
  layer_max_pooling_2d(pool_size = c(1, 1)) %>% 
  layer_dropout(rate = 0.25) %>% 

  layer_flatten() %>% 
  layer_dense(units = 512, activation = 'tanh') %>% 
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 2, activation = 'sigmoid')


model %>% compile(
  loss = 'binary_crossentropy',
  optimizer = optimizer_adam(lr = 0.001),
  metrics = c('accuracy')
)

history <- model %>% fit(
  x_train, y_train, 
  epochs = 30, batch_size = 512, 
  validation_split = 0.2,
  verbose = 2
)
result=predict(model,x_test)[,2]
re=cbind(ti,result)
write.table(file='result2.10.txt',re,append=T,row.names=F,col.names=F,quote=F)

}


library(pROC)
library(PRROC)
a=read.table('result2.10.txt',stringsAsFactor=F)
roc(a[,3],a[,4])
pr.curve(a[a[,3]==1,4],a[a[,3]==0,4])



