import numpy as np
import glob
import pandas as pd
import os
from multiprocessing import Pool
import time as tm
import argparse


##############################################################################################Required arguments#####################################################################################################################################################

parser=argparse.ArgumentParser(description="refernce database creator of different bin size 500kb to 10000kb")
parser.add_argument('-N','--normalization_file_location',type=str,required=True,help='Required: location of the normalization files whose  bin size is 500kb')
parser.add_argument('-R','--output_dir',type=str,required=True,help='Required: path of folder name where all refernce database file of different size should be saved/created')
args=parser.parse_args()
startTime=tm.time()




#############################################################################################1.FUNCTION##############################################################################################################################################################


''''normalization function used to formed a new column normalization read count by using the method read count 
    #####data['nor_reads'+new_s]=(data['reads']/(data_1_22['reads'].sum()))*1000###############
    *******************data_1_22 is that which not contain chrx and chry************************************'''

def normalization_reads(data,new_s):
    data_1_22= data[data['space'].str.contains("chrX|chrY")!=True]
    data['nor_reads'+new_s]=(data['reads']/(data_1_22['reads'].sum()))*1000
    data1=data[['nor_reads'+new_s]]
    return (data1)
	
###################################################################################################2.FUNCTION########################################################################################################################################################

'''' reference datase function is used to create the database of different bin size this is the output which we used in the cnv detection programme
    *********************************************function return new database of different bin size*************************************************'''

def reference_database(nums):
    appended_data=[]
    for infile in glob.glob('args.normalization_file_location+"/*tsv'):
            data = pd.read_csv(infile, sep='\t',usecols=['space','start','end','reads'])
            group = data.groupby('space')
            chromosome_list=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
            data_final=pd.DataFrame()
            for  v in chromosome_list:
                df=group.get_group(v).reset_index()
                #df is a dataframe of group by data
                data1= pd.DataFrame()
                lists=[]
                lists1=[]
                reads1=[]
                spaces=[]
                i=0
                while i<(len(df)-nums):
                    start_bin =df.start[i]
                    lists.append(start_bin)
                    space=df.space[i]
                    spaces.append(space)
                    end_bin=df.end[i+nums]
                    lists1.append(end_bin)
                    reads=df.loc[range(i,i+nums+1),'reads'].sum()
                    reads1.append(reads)
                    i+=nums+1
                remain_line=len(df) % (nums+1)
                if(remain_line>0):
                    start_bin =df.start[i]
                    lists.append(start_bin)
                    space=df.space[i]
                    spaces.append(space)
                    end_bin=df.end[i+remain_line-1]
                    lists1.append(end_bin)
                    reads=reads=df.loc[range(i,i+remain_line),'reads'].sum()
                    reads1.append(reads)
                    width=width1[1]
                    width1.append(width)
                data1['space']=spaces
                data1['start']=lists
                data1['end']=lists1
                data1['width']=data1.end-data1.start+1
                data1['reads']=reads1
                data_final=pd.concat([data_final,data1],axis=0)
            data_final=data_final.reset_index()
            data_final= data_final.drop(['index'],axis=1)
            #data_final is data of different size of bin 
            new_s=os.path.basename(infile)
            new_s=new_s.replace('.Normalization.txt','')
            arc_data=normalization_reads(data_final,new_s)
            appended_data.append(arc_data)
    appended_data = pd.concat(appended_data, axis=1)
    data_final=data_final[['space','start','end','width']]
    in_file=data_final.join(appended_data)
    in_file["nor_mean"]=appended_data.mean(axis=1)
    in_file["nor_std"]=appended_data.std(axis=1)
    new_database=in_file[['space','start','end','width','nor_mean','nor_std']]
    return(new_database)
	
	
*********************************************************************************************************************************************************************************************************************************************************************


values=[0,1,3,5,7,9,11,13,15,17,19]
###Use values to create a data of different bin size taking normal files of size 500kb
for j in values:
    path= args.output_dir
    df_1000 =reference_database(j)
    output_file = os.path.join(path,'file_1000'+str(j)+'.txt')
    df_1000.to_csv(output_file,index=False,sep="\t")
    ###df_1000 we have also 500kb and 1 mb to 10 mb bin size files
totalTime=tm.time() - startTime
print(totalTime/60)


