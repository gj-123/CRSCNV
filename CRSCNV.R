library(PythonInR)
library(readr)
RD_filename <- read_table2("RD_filename.txt",col_names = FALSE)
data <- read_table2(RD_filename$X1,col_names = FALSE)
#depth:readcount  index:xia biao 
index<-vector(mode = "numeric",length = 0)
depth<-vector(mode = "numeric",length = 0)
index<-data$X1 
depth<-data$X2 
#binstart  binend binsize:1000 
binstart <- vector(mode =  "numeric" , length = length(index))
binend <- vector(mode =  "numeric" , length = length(index))
j<-0
for (i in index) {
  j<-j+1
  binstart[j]<-(i*1000+1)
  binend[j]<-(i*1000+1000)
}
write(depth, file = "average_gc_1.txt",ncolumns =1,sep = " ")
pyConnect(RD_filename$X2)
pyExecfile("CRSCNV_TV.py")
#denoise  
data2 <- read_table2("average_tv.txt",col_names = FALSE)
average_2_new<-data2$X1
#get mode
average_2_sort<-sort(average_2_new)
average_2_min<-mean(average_2_sort[1:1000]) 
average_2_max<-mean(average_2_sort[(length(average_2_new)-1000+1):length(average_2_new)]) 
average_2_mode<-as.numeric(names(table(average_2_new))[which.max(table(average_2_new))])
average_2_sub<-average_2_new-rep(average_2_mode,length(average_2_new))
average_2_new1<-abs(average_2_max/average_2_min)*average_2_sub
#qiu 10 ge segment de jun zhi index  matrix_data1
matrix_data1<-matrix(data = NA,nrow = floor(length(index)/10),ncol =10)
matrix_data3<-matrix(data = NA,nrow = floor(length(index)/10),ncol =10)
binstart_1_matrix<-matrix(data = NA,nrow = floor(length(index)/10),ncol =10)
binend_1_matrix<-matrix(data = NA,nrow = floor(length(index)/10),ncol =10)
k<-0
for(i in seq(1,length(index)-length(index)%%10,floor(length(index)/10)))
{
    l<-0
    k <- k+1
    for(j in 0:(floor(length(index)/10)-1))
    {
      l<-l+1
      #matrix_data1:store each segment element
      matrix_data1[l,k]<-average_2_new1[i+j]
      #matrix_data3:store each segment index
      matrix_data3[l,k]<-i+j
      binstart_1_matrix[l,k]<-binstart[i+j]
      binend_1_matrix[l,k]<-binend[i+j]
    }
}
#matrix_data4:each col data matrix_data4
matrix_data4<-matrix(data = NA,nrow = length(index)-length(index)%%10-floor(length(index)/10),ncol =10)
for (i in 1:10) {
    c<-vector(mode="numeric",length=0)
    c<-setdiff(1:(length(index)-length(index)%%10),as.vector(matrix_data3[,i]))
    for (j in 1:length(c)) {
      matrix_data4[j,i]<-average_2_new1[c[j]]
    }
}
average1 <- vector(mode = "numeric",length = 10)
standard_deviation <- vector(mode = "numeric",length = 10)
for (i in 1:10) {
    average1[i]<-mean(matrix_data4[,i])
    standard_deviation[i]<-sqrt(var(matrix_data4[,i]))
}
#erf funtion
erf <- function(z){
  t = 1.0 / (1.0 + 0.5 * abs(z));
  
  
  ans = 1 - t * exp( -z*z   -   1.26551223 +
                       t * ( 1.00002368 +
                               t * ( 0.37409196 + 
                                       t * ( 0.09678418 + 
                                               t * (-0.18628806 + 
                                                      t * ( 0.27886807 + 
                                                              t * (-1.13520398 + 
                                                                     t * ( 1.48851587 + 
                                                                             t * (-0.82215223 + 
                                                                                    t * ( 0.17087277))))))))));
  if (z >= 0) return(ans)
  else        return(-ans)
}
#pvalue
t<-matrix(data = NA,nrow = floor(length(index)/10),ncol =10)
for(i in 1:10){
  for (j in 1:floor(length(index)/10)) {
    t[j,i] = (matrix_data1[j,i]-average1[i])/standard_deviation[i]
  }
}
#pvaluelist:matrix_data2 store each bin pvalue
matrix_data2<-matrix(data = NA,nrow = floor(length(index)/10),ncol =10)
matrix_data2l<-matrix(data = NA,nrow = floor(length(index)/10),ncol =10)
matrix_data2r<-matrix(data = NA,nrow = floor(length(index)/10),ncol =10)
for(i in 1:10){
  for (j in 1:floor(length(index)/10)) {
    #left
    matrix_data2l[j,i] = 0.5*(1-erf(-t[j,i]*sqrt(0.5)))
    #right
    matrix_data2r[j,i] = 1-0.5*(1-erf(-t[j,i]*sqrt(0.5)))
  }
}
matrix_data2_indexlist<-list()
matrix_data2_list<-list()
binstart_1_list<-list()
binend_1_list<-list()
for (i in 1:10) {
  matrix_data2_indexlist[[i]]<-which(matrix_data2r[,i]<0.05)
}
m<-10
for (i in 1:10) {
  m<-m+1
  matrix_data2_indexlist[[m]]<-which(matrix_data2l[,i]<0.05)
}
#right
for (i in 1:10) {
    k<-0
    matrix_data2_vector <- vector(mode = "numeric",length = 0)
    binstart_1_vector <- vector(mode = "numeric",length = 0)
    binend_1_vector <- vector(mode = "numeric",length = 0)
    for (j in matrix_data2_indexlist[[i]]) {
      k<-k+1
      matrix_data2_vector[k]<-matrix_data2[j,i]
      binstart_1_vector[k]<-binstart_1_matrix[j,i]
      binend_1_vector[k]<-binend_1_matrix[j,i]
    }
    matrix_data2_list[[i]]<-matrix_data2_vector
    binstart_1_list[[i]]<-binstart_1_vector
    binend_1_list[[i]]<-binend_1_vector
}
#left
for (i in 1:10) {
  k<-0
  matrix_data2_vector <- vector(mode = "numeric",length = 0)
  binstart_1_vector <- vector(mode = "numeric",length = 0)
  binend_1_vector <- vector(mode = "numeric",length = 0)
  for (j in matrix_data2_indexlist[[i+10]]) {
    k<-k+1
    matrix_data2_vector[k]<-matrix_data2[j,i]
    binstart_1_vector[k]<-binstart_1_matrix[j,i]
    binend_1_vector[k]<-binend_1_matrix[j,i]
  }
  matrix_data2_list[[i+10]]<-matrix_data2_vector
  binstart_1_list[[i+10]]<-binstart_1_vector
  binend_1_list[[i+10]]<-binend_1_vector
}
#matrix_data2_list1 store pvalue
matrix_data2_list1 <- vector(mode = "numeric",length = 0)
binstart_1_list1 <- vector(mode = "numeric",length = 0)
binend_1_list1 <- vector(mode = "numeric",length = 0)
matrix_data2_list1<-unlist(matrix_data2_list)
binstart_1_list1<-unlist(binstart_1_list)
binend_1_list1<-unlist(binend_1_list)
matrix_data5<-matrix(data = NA,nrow = length(matrix_data2_list1),ncol =2)
matrix_data5<-cbind(binend_1_list1,matrix_data2_list1)
matrix_data5_1<-matrix(data = NA,nrow = length(matrix_data2_list1),ncol =3)
matrix_data5_1<-cbind(binstart_1_list1,matrix_data5)
#write(t(matrix_data5_1), file = "matrix_data5_1.txt",ncolumns =3,sep = " ")
#get variation result
if(length(matrix_data2_list1)!=0){
    #Define variation region1
    variation_type<-vector(mode="character",length=nrow(matrix_data5_1))
    charname<-vector(mode="character",length=nrow(matrix_data5_1))
    #get boundary
    for (i in 1:nrow(matrix_data5_1)) {
      binstartindex<-0
        binstartindex<-which(binstart==matrix_data5_1[i,1])
        if(average_2_new1[binstartindex]>mean(average_2_new1)){
          variation_type[i]<-"gain"
        }else{
          variation_type[i]<-"loss"
        }
        charname[i]<-"chr21"
    }
    #total
    variation<-cbind(matrix_data5_1[,2],variation_type)
    variation1<-cbind(matrix_data5_1[,1],variation)
    variation2<-cbind(charname,variation1)
    #write(t(variation2), file = "variation2.txt",ncolumns =4,sep = " ")
    if(nrow(matrix_data5_1)>=2){
      for (i in 1:(nrow(matrix_data5_1)-1)){
        if((as.numeric(variation2[i,3])+1==as.numeric(variation2[i+1,2]))&&(variation2[i,4]==variation2[i+1,4])&&(variation2[i,1]==variation2[i+1,1])){
          variation2[i+1,2]<-variation2[i,2]
          variation2[i,]=""
        }
      }
      variation2_new<-variation2[which(variation2[,1]!=""),]
      write(t(variation2_new), file = "result.txt",ncolumns =4,sep = " ")
    }else{
      write(t(variation2), file = "result.txt",ncolumns =4,sep = " ")
    }
    
}else{
  variation2_new<-""
  write(t(variation2_new), file = "result.txt",ncolumns =4,sep = " ")
}


