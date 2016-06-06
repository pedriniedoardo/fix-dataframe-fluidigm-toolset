#read the file
setwd("C:/Users/edoardo/Desktop/seattle experiments/single cells analysis/160420 single cell HL60 experiment/single cell analysis vincristine")
df<-read.csv("data mock vs 48h.csv",skip = 11)
#is better to be sure to remove the factor and transform the variable in character
df$ID<-as.character(df$ID)
df$Name<-as.character(df$Name)
df$Name.1<-as.character(df$Name.1)
##################################
##################################
#prepare a better label for the samples
l<-rep(101:196,each=96)
l2<-rep(1:96,each=96)
#part from the df
t(as.data.frame(strsplit(df$Name,split = "_")))
first<-t(as.data.frame(strsplit(df$Name,split = "_")))[,1]
second<-t(as.data.frame(strsplit(df$Name,split = "_")))[,2]
#set up the complete new label
label<-paste0(first,second,"_",l,"-",l2)
#change the old label with the new one
df_2<-df
df_2$treatment<-df_2$Name
df_2$Name<-as.character(label)
##################################
##################################
#rearrange the df
#sort out all the 100 cell and the all human rna control
df_3<-df_2[df_2$rConc==1,]
#collapse the tecnical replicate
px<-{}
#1
for(i in unique(df_3$Name)){
  print(i)
  p<-df_3[df_3$Name==i,]
  #2
  for(l in unique(df_3$Name.1)){
    print(l)
    p2<-p[p$Name.1==l,]
    if(p2$Call[1]=="Fail" & p2$Call[2]=="Fail"){
      px<-rbind.data.frame(px,p2[1,])
    }
    else if(p2$Call[1]=="Fail" & p2$Call[2]!="Fail"){
      px<-rbind.data.frame(px,p2[2,])
    }
    else if(p2$Call[1]!="Fail" & p2$Call[2]=="Fail"){
      px<-rbind.data.frame(px,p2[1,])
    }
    else if(p2$Call[1]!="Fail" & p2$Call[2]!="Fail"){
      p3<-p2
      p3[1,7]<-mean(p2$Value)
      px<-rbind.data.frame(px,p3[1,])
    }
  }
}
#remove the column calibrated conc
px<-px[-8]
#change the label of the assays and sample in the ID column
S<-rep("S",nrow(px))
samp<-c(rep(paste0(0,1:9),each=48),rep(10:length(unique(px$Name)),each=48))
A<-rep("A",nrow(px))
assay<-rep(c(paste0(0,1:9),10:48),length(unique(px$Name)))
ID_new<-paste0(S,samp,"-",A,assay)
px$ID<-ID_new
#add the column Peak Call
px$`Peak Call`[px$Comments=="ok"]<-"Pass"
px$`Peak Call`[px$Comments!="ok"]<-"Fail"
##################################
##################################
#produce the final df and write them
#df treated
px_treat<-px[px$treatment=="48h_low"|px$treatment=="48h_high",]
#df mock
px_mock<-px[px$treatment=="mock_low"|px$treatment=="mock_high",]
#very important remove the dots from the column names!
colnames(px)<-c("ID","Name","Type","rConc","Name","Type","Value","Quality","Call","Threshold","In Range","Out Range","Peak Ratio","Comments","treatment","Peak Call")
colnames(px_mock)<-c("ID","Name","Type","rConc","Name","Type","Value","Quality","Call","Threshold","In Range","Out Range","Peak Ratio","Comments","treatment","Peak Call")
colnames(px_treat)<-c("ID","Name","Type","rConc","Name","Type","Value","Quality","Call","Threshold","In Range","Out Range","Peak Ratio","Comments","treatment","Peak Call")
#write the df
write.csv(px,row.names = FALSE,"C:/Users/edoardo/Desktop/seattle experiments/single cells analysis/160420 single cell HL60 experiment/single cell analysis vincristine/dataframe fluidigm.csv")
write.csv(px_treat,row.names = FALSE,"C:/Users/edoardo/Desktop/seattle experiments/single cells analysis/160420 single cell HL60 experiment/single cell analysis vincristine/dataframe fluidigm treat.csv")
write.csv(px_mock,row.names = FALSE,"C:/Users/edoardo/Desktop/seattle experiments/single cells analysis/160420 single cell HL60 experiment/single cell analysis vincristine/dataframe fluidigm mock.csv")
#prepare a file with the unique name of the sample and the label of the relative control
sampleMock<-unique(px_mock$Name)
sampleTreat<-unique(px_treat$Name)
labelMock<-rep("mock",length(sampleMock))
labelTreat<-rep("48h vinc 3 nM",length(sampleTreat))
name_sample<-cbind("SampleID"=c(sampleMock,sampleTreat),"GroupID"=c(labelMock,labelTreat))
write.table(name_sample,quote = FALSE,sep="\t",row.names = FALSE,"C:/Users/edoardo/Desktop/seattle experiments/single cells analysis/160420 single cell HL60 experiment/single cell analysis vincristine/name_sample.txt")
