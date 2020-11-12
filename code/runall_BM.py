
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time

#/lustre2/home/tangqianzi/software/mirnylab-hiclib-8c146707b6f6/examples/pipeline_2017/

class generate:

	def __init__ (self,path='/Lustre02/data/hic/forJinplot/crossspecies/exprs_sele/'):

		self.path = path
		self.jobfolder = 'runall_BM.nochickenjobs'

	def generate (self):

                ## handle=open(self.path+'/tissuenames.txt','r')
                ## handle2=open('/Lustre02/data/hic/forJinOrth/mRNAOrth/chicken/finalforcornew.txt','r')
                #handle3=open(self.path+'/datafiles/myconvert.txt','r')

                samples=['Subcutaneous_adipose']

                tissue2species={'Subcutaneous_adipose':['cat','dog','human','pig','rabbit','rat','sheep']}

	        os.system('mkdir -p '+self.path+'/'+self.jobfolder)
	        os.system('mkdir -p '+self.path+'/output/')
	        os.system('mkdir -p '+self.path+'/Rscripts/')

	        myoutput=open(self.path+'/'+self.jobfolder+'/runthis_pipeline.sh','w')

                #======= cat ===============================

	        myoutput1=open(self.path+'/'+self.jobfolder+'/runthis_cat.sh','w')
	        myoutput2=open(self.path+'/'+self.jobfolder+'/alljobs_cat.sh','w')

	        for tissue in samples:
	               for species in tissue2species[tissue]:
	                       handle=open(self.path+'/Rscripts/'+tissue+'_'+species+'.R','w')

	                       print>>handle,'''setwd("'''+self.path+'''/")
library(phytools)
input_trees<-c("/Lustre02/data/hic/forJinplot/crossspecies/exprs_sele/'''+species+'''.tre")

input_csvs<-c("forselection_'''+tissue+'''_mean.csv")
input_ses<-c("forselection_'''+tissue+'''_sd.csv")

myfitBM<-function(input_tree,input_csv,input_se){

tree <- read.simmap(input_tree, format = "phylip")
data_ori <- read.csv(file=input_csv,header=T,row.names=1)
#diff from previous codes here!!!
data_ori<-t(data_ori)

se_ori <- read.csv(file=input_se,header=T,row.names=1)
#diff from previous codes here!!!
se_ori<-t(se_ori)

delete_col<-c()

myfilter<-function(inputfilter){


!isTRUE(all.equal.numeric(as.vector(inputfilter),rep(0,'''+str(len(tissue2species[tissue]))+''')))

}

myfilter_result <- sapply( 1:ncol(data_ori), function(j) myfilter( data_ori[,j] ) )

data<-data_ori[,myfilter_result]

#mytestse<-matrix(rep(c(0,0.01,0.01,0,0.01),10),nrow=5,ncol=10)
mytestsenew<-se_ori[,myfilter_result]

BMstep<-function(inputstep,inputname,inputse,count){

res <- try(

      brownie.lite(tree, inputstep)
      #brownie.lite(tree, inputstep, se=inputse)

  )

  if(!inherits(res, "try-error"))
{
print (count)

fitBMi <- brownie.lite(tree, inputstep)
#fitBMi <- brownie.lite(tree, inputstep, se=inputse)

parameters<-c(fitBMi$logL1,fitBMi$k1,fitBMi$sig2.multiple[1],fitBMi$sig2.multiple[2],fitBMi$logL.multiple,fitBMi$k2,fitBMi$P.chisq)
para_addname<-c(inputname,parameters)
return(para_addname)

}
}

results <- sapply( 1:ncol(data), function(j) BMstep( data[,j], colnames(data)[j], mytestsenew[,j], j ) )

results<-Filter(Negate(is.null), results)

results <- matrix(unlist(results), ncol = 8, byrow = TRUE)

colnames(results)<-c('genename','logL1','k1','BMrate_background','BMrate_specificspecies','logL.multiple','k2','P.chisq')

write.table(results,file="'''+self.path+'/output/'+tissue+'-'+species+'''.BMresults.txt",sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

}

myfitBM(input_trees[1],input_csvs[1],input_ses[1])
'''


	                       handle.close()

	                       print>>myoutput2,'''export PATH=/Lustre01/share/software/R-3.5.0/bin/:$PATH
export R_LIBS=/Lustre01/tangqianzi/software/RlibsR35new/:$R_LIBS
export LD_LIBRARY_PATH=/Lustre01/tangqianzi/software/anaconda2new/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/Lustre01/tangqianzi/software/:$LD_LIBRARY_PATH'''
	                       print>>myoutput2,'Rscript '+self.path+'/Rscripts/'+tissue+'_'+species+'.R'

	        print>>myoutput1,'/usr/bin/perl /Lustre01/tangqianzi/software/scripts/qsub-sgenew.pl --maxjob 30 --lines 5 --jobprefix runBM --convert no --resource nodes=1:ppn=1,mem=20g '+self.path+'/'+self.jobfolder+'/alljobs_cat.sh'

	        myoutput1.close()
	        myoutput2.close()

	        print>>myoutput,'sh '+self.path+'/'+self.jobfolder+'/runthis_cat.sh'

		myoutput.close()


def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	generate().generate()


if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
