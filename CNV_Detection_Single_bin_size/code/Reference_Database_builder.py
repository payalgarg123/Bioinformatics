import argparse
import numpy as np
import glob
import pandas as pd
import os
import time as tm


************************************************************************************************Required Arguments***************************************************************************************************************************************************

parser=argparse.ArgumentParser(description="refernce data of gc value for all normalization files")
parser.add_argument('-N','--normalization_file_location',type=str,metavar="",required=True,help='location of the normalization files to be in this program')
parser.add_argument('-M','--normalization_method',type=str,metavar=["gc","reads"],required=True,help=' normalization method')
parser.add_argument('-R','--refernce_data_base',type=str,metavar='',required=True,help='location of refernce data base')
args=parser.parse_args()

startTime=tm.time()
appended_data=[]
data_raw=[]

#####################################################################################################1.Function#####################################################################################################################################################

'''' For normalization  used gc method. In this method we can used gc values to normalized count
data_1_22= data[data['space'].str.contains("chrX|chrY")!=True] data_1_22 is not contain chrX and ChrY
d= data_1_22['reads'].mean()
dgc=data_1_22.loc[(data_1_22['gc'] >= df.gc[i]-0.1) & (data_1_22['gc'] <= df.gc[i]+0.1), 'reads'].mean() 
ARC_gc_count(df.reads[i]*(d/dgc))------ method'''


if args.normalization_method=='gc':
    def normalization_reads(data,new_s):
        data_1_22= data[data['space'].str.contains("chrX|chrY")!=True]
        d= data_1_22['reads'].mean()
        ARC=[]
        for i in range (len(data)):
            if(data.gc[i]==-1):
                dgc=d
            else:
                dgc=data_1_22.loc[(data_1_22['gc'] >= data.gc[i]-0.1) & (data_1_22['gc'] <= data.gc[i]+0.1), 'reads'].mean()
            ARC.append(data.reads[i]*(d/dgc))
        data['ARC'+new_s]= ARC
        data1=data[['ARC'+new_s]]
        return (data1)

 ''' for normalization used read method.In this method we can use the reads value to create a new column normalization read count
  df['data_r']=(df['reads']/data_1_22['reads'].sum())*1000-------Method'''
      
elif args.normalization_method=='reads':
        def normalization_reads(data,new_s):
            data_1_22= data[data['space'].str.contains("chrX|chrY")!=True]
            data['nor_reads'+new_s]=(data['reads']/(data_1_22['reads'].sum()))*1000
            data1=data[['nor_reads'+new_s]]
            return (data1)
else:
    print(error)

        
*********************************************************************************************************************************************************************************************************************************************************************


         
for infile in glob.glob(args.normalization_file_location+"/*tsv"): ####used glob to read multiple files
        data = pd.read_csv(infile, sep='\t',usecols=['reads','gc','space'])
        new_s=os.path.basename(infile) #create basename path
        new_s=new_s.replace('.Normalization.txt','')
        arc_data=normalization_reads(data,new_s)
        appended_data.append(arc_data)
        
df2=pd.read_csv(infile, sep='\t', usecols = ['space','start','end','width'])
# concat the data
appended_data = pd.concat(appended_data, axis=1)
in_file=df2.join(appended_data)
in_file["nor_mean"]=appended_data.mean(axis=1) ####calculate mean
in_file["nor_std"]=appended_data.std(axis=1)   ####calculate standard deviation
new_database=in_file[['space','start','end','nor_mean','nor_std']]
            
new_database.to_csv(args.refernce_data_base,sep="\t",index=False,header=True,index_label=None)  ####save reference database #output
totalTime=tm.time() - startTime
print(totalTime/60)
