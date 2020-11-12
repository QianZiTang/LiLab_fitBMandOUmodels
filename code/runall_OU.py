
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
		self.jobfolder = 'runall_OU.nochickenjobs'

	def generate (self):

                ## handle=open(self.path+'/tissuenames.txt','r')
                ## handle2=open('/Lustre02/data/hic/forJinOrth/mRNAOrth/chicken/finalforcornew.txt','r')
                #handle3=open(self.path+'/datafiles/myconvert.txt','r')

                samples=['Subcutaneous_adipose']

                tissue2species={'Subcutaneous_adipose':['cat']}

	        os.system('mkdir -p '+self.path+'/'+self.jobfolder)
	        os.system('mkdir -p '+self.path+'/Rscripts/')
	        os.system('mkdir -p '+self.path+'/output/')

	        myoutput=open(self.path+'/'+self.jobfolder+'/runthis_pipeline.sh','w')

                #======= cat ===============================

	        myoutput1=open(self.path+'/'+self.jobfolder+'/runthis_cat.sh','w')
	        myoutput2=open(self.path+'/'+self.jobfolder+'/alljobs_cat.sh','w')

	        for tissue in samples:
	               for species in tissue2species[tissue]:
	                       for m in range(1,10):

	                               handle=open(self.path+'/Rscripts/'+tissue+'_'+species+'.OU.'+str(m)+'.R','w')

	                               print>>handle,'''setwd("'''+self.path+'''/cut3parts/")
library(geiger)
library(phytools)

input_trees<-c("/Lustre02/data/hic/forJinplot/crossspecies/exprs_sele/'''+species+'''.tre")

input_csvs<-c("forselection_'''+tissue+'''_mean_'''+str(m)+'''.csv")
input_ses<-c("forselection_'''+tissue+'''_sd_'''+str(m)+'''.csv")

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

#mytestse<-matrix(rep(c(NA,1,1,NA,1),10),nrow=5,ncol=10)
mytestsenew<-se_ori[,myfilter_result]

BMstep<-function(inputstep,inputname,inputse,count){


res <- try(

      fitContinuous(tree, inputstep, ncores=2, model="OU")
      #fitContinuous(tree, inputstep, SE=inputse, ncores=2, model="OU")

  )

  if(!inherits(res, "try-error"))
{
print (count)

fitBMi <- fitContinuous(tree, inputstep, ncores=2, model="OU")

#fitBMi <- fitContinuous(tree, inputstep, SE=inputse, ncores=2, model="OU")

parameters<-c(fitBMi$opt$aic,fitBMi$opt$lnL,fitBMi$opt$k,fitBMi$opt$aicc)
para_addname<-c(inputname,parameters)

return(para_addname)


}
}


results <- sapply( 1:ncol(data), function(j) BMstep( data[,j], colnames(data)[j], mytestsenew[,j], j ) )

results<-Filter(Negate(is.null), results)

results <- matrix(unlist(results), ncol = 5, byrow = TRUE)

colnames(results)<-c('genename','AIC','lnL','k','AICC')

write.table(results,file="'''+self.path+'/output/'+tissue+'-'+species+'''-'''+str(m)+'''.OUresults.txt",sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

}

myfitBM(input_trees[1],input_csvs[1],input_ses[1])
'''


	                               handle.close()

	                               print>>myoutput2,'''export PATH=/Lustre01/share/software/R-3.5.0/bin/:$PATH
export R_LIBS=/Lustre01/tangqianzi/software/RlibsR35new/:$R_LIBS
export LD_LIBRARY_PATH=/Lustre01/tangqianzi/software/anaconda2new/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/Lustre01/tangqianzi/software/:$LD_LIBRARY_PATH'''
	                               print>>myoutput2,'Rscript '+self.path+'/Rscripts/'+tissue+'_'+species+'.OU.'+str(m)+'.R'

	        print>>myoutput1,'/usr/bin/perl /Lustre01/tangqianzi/software/scripts/qsub-sgenew.pl --maxjob 70 --lines 5 --jobprefix runOU --convert no --resource nodes=1:ppn=1,mem=20g '+self.path+'/'+self.jobfolder+'/alljobs_cat.sh'

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
