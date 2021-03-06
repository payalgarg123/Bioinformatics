
# coding: utf-8

# In[1]:


import numpy as np
import glob
import pandas as pd
import os
import seaborn as sns; sns.set(style='whitegrid')
import matplotlib.pyplot as plt
import argparse
import time as tm

****************************************************************Required Arguments************************************************************************************

parser=argparse.ArgumentParser(description="CNV deTection and HeatMap of different bin size")
parser.add_argument('-N','--sample_file',type=str,required=True,help='Required: location of the patients read count file at 500kb bin size')
parser.add_argument('-R','--reference_database_folder',type=str,required=True,help='reference database of folder where multiple bins reference file were generated by Reference_Database_builder.py programe ') 
parser.add_argument('-Z','--z_score_table',type=str,required=True,help='name of output file containing z-score  values for various bin size')
parser.add_argument('-W','--window_size',metavar="",type=int,help='rolling size')
parser.add_argument('-S','--min_cnv_size',type=int,help='minimum continous cnv periods')
parser.add_argument('-O','--output_file',type=str,required=True,help="final report of continous deletion and duplication")
parser.add_argument('-P','--plotting_file',type=str,required=True,help="Heatmap of multiple bin z-scores for human Genome")
args=parser.parse_args()
startTime=tm.time()


######################################################################1.function#######################################################################################

'''' This function is used for normalization reads count
data_1_22 is the data which not contain chrX and chrY  
data['nor_reads']=(data['reads']/(data_1_22['reads'].sum()))*1000 -----Method'''

def normalization_reads1(data):
    data_1_22= data[data['space'].str.contains("chrX|chrY")!=True]
    data['nor_reads']=(data['reads']/(data_1_22['reads'].sum()))*1000
    return (data)

#####################################################################2.function######################################################################################

''''This function is used to create a data variable different size 1mb to 10 mb include 500kb=0.5mb
1000000=1mb=1000kb'''


def cnv_data(nums):
    data = pd.read_csv(args.sample_file, sep='\t',usecols=['space','start','end','reads'])
    group = data.groupby('space')
    data_final1=pd.DataFrame()
    chromosome_list=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    for  v in chromosome_list:
        df=group.get_group(v).reset_index()
        data1= pd.DataFrame()
        lists=[]
        lists1=[]
        reads1=[]
        spaces=[]
        width1=[]
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
            width=end_bin-start_bin+1
            width1.append(width)
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
        data1['width']=width1
        data1['reads']=reads1
        
        data_final1=pd.concat([data_final1,data1],axis=0)
    data_final1=data_final1.reset_index()
    data_final1= data_final1.drop(['index'],axis=1)
    normalization_reads1(data_final1)
    return (data_final1)

#############################################################3.function###############################################################################################

'''' this function is used to detect the cnv(copy number variation) using the continous bin with z-score +/-1.2 and lists the result in a txt'''  


def continous_bin(df,data,condition):
    i=0
    count=0
    data_list=[]
    cont_data=pd.DataFrame()
    un_cont_data=pd.DataFrame()
    while i<(len(data)-1):
        start_bin=data.start[i]
        end_bin=data.start[i+1]
        i+=1
        
        
        
        if end_bin-start_bin==df.width[i]:
            count+=1
             
        elif count>=args.min_cnv_size:
            real_start=data.start[i-count-1]
            real_end=data.end[i-1]
            con_bin_mean=data[(data['start'] >= real_start) & (data['start'] <= data.start[i-1])]
            cont_data=pd.concat([cont_data,con_bin_mean])
            total_c_mean=con_bin_mean['z'].mean()
            data_list.append([condition,data.space[i],real_start,real_end,total_c_mean])
            count=0
        
        else:
            real_start=data.start[i-count-1]
            real_end=data.end[i-1]
            un_con=data[(data['start'] >=real_start ) & (data['start'] <= data.start[i-1])]
            un_cont_data=pd.concat([un_cont_data,un_con])
            count=0
        
        
        if i== len(data)-1:
            if count>=args.min_cnv_size:
                real_start=data.start[i-count]
                real_end=data.end[i]
                con_bin_mean=data[(data['start'] >= real_start) & (data['start'] <= data.start[i])]
                cont_data=pd.concat([cont_data,con_bin_mean])
                un_cont_data=un_cont_data
                total_c_mean=con_bin_mean['z'].mean()
                data_list.append([condition,data.space[i],real_start,real_end,total_c_mean])
				Myfile=open(args.output_file,'a')
                for lines in data_list:
                    Myfile.write(str(lines))
                    Myfile.write(" ")
                Myfile.close()
                count=0
                #return(cont_data)
        
            else:
                real_start=data.start[i-count]
                real_end=data.end[i]
                un_con=data[(data['start'] >=real_start ) & (data['start'] <= data.start[i])]
                un_cont_data=pd.concat([un_cont_data,un_con])
                cont_data=cont_data
                data_list=data_list
				Myfile=open(args.output_file,'a')
                for lines in data_list:
                    Myfile.write(str(lines))
                    Myfile.write(" ")
                Myfile.close()
                count=0
                
                
#####################################################################4.function######################################################################################## 

'''create a new column index_plot for plotting '''

               
def index_p(df2):
    index_plot=[]
    if df2.width[1]==500000:
        for k in (df2.index):
            index_plot1=((11-df2.width[1]*0/1000000)/(len(df2)))*k
            index_plot.append(index_plot1)
    else:
        for k in (df2.index):
            index_plot1=(((11-df2.width[1]/1000000)/(len(df2)))*k)
            index_plot.append(index_plot1)
    df2['index_plot']=index_plot
    return(df2)

***********************************************************************************************************************************************************************



values=[0,1,3,5,7,9,11,13,15,17,19]
df_c=pd.DataFrame()
for j in  values :
    finals=cnv_data(j)
    main_con=pd.DataFrame()
    for infiles in glob.glob('C:/Users/EDGC_IN_01/reference database/file_1000'+str(j)+'.txt'):
        r_data = pd.read_csv(infiles,sep='\t',usecols=['space','nor_mean','nor_std'])
        finals['zscore']=((finals.nor_reads - r_data['nor_mean'])/r_data['nor_std']).fillna(0)
        groups = finals.groupby('space')
        df_h=pd.DataFrame()
        chromosome_list=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
        #chromosome_list=['chr1']
        for  v in chromosome_list:
            df2=groups.get_group(v).reset_index()
            df2= df2.drop(['index'],axis=1)
            df2.loc[:,'z']=df2['zscore'].rolling(window=args.window_size,min_periods=1,center=True).mean()
            df2=index_p(df2)
            
            df_duplication=df2[df2.z>=1.2].reset_index()
            df_deletion=df2[df2.z<=-1.2].reset_index()
            continous_bin(df2,df_duplication,'duplication')
            continous_bin(df2,df_deletion,'Deletion')
            df_h=pd.concat([df_h,df2],axis=0)
    df_c=pd.concat([df_c,df_h],axis=0)
df_c.to_csv(args.z_score_table) #database of z-score 
group_whole=df_c.groupby('space')
rowlength = group_whole.ngroups//2 #for subplot divided into  parts (2,12)
fig,axes = plt.subplots(figsize=(100,20),nrows=2, ncols=rowlength,gridspec_kw=dict(hspace=1))
chromosome_list=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
for g,ax in zip(chromosome_list,axes.flatten()):
    group=group_whole.get_group(g)
    f = group.pivot('width','index_plot',"z") #make a database in pivot wich we used in plotting
    f.sort_index(ascending=False, inplace=True)
    g1=sns.heatmap(f.bfill(axis=1), cmap='RdBu',vmax=3,vmin=-3,square=True,xticklabels=False, cbar=False,ax=ax) #seaborn library
    g1.set_aspect(1./ax.get_data_ratio())
    g1.set_title(g,fontdict=30)
    g1.set_xlabel('#bins')
    g1.set_ylabel('size of bins')
    fig = g1.get_figure()
    fig.savefig(args.plotting_file)
                


# In[ ]:





# In[120]:


values=[0,1,3,5,7,9,11,13,15,17,19]
df_c=pd.DataFrame()
for j in  values :
    finals=cnv_data(j)
    main_con=pd.DataFrame()
    for infiles in glob.glob('C:/Users/EDGC_IN_01/reference database/file_1000'+str(j)+'.txt'):
        r_data = pd.read_csv(infiles,sep='\t',usecols=['space','nor_mean','nor_std'])
        finals['zscore']=((finals.nor_reads - r_data['nor_mean'])/r_data['nor_std']).fillna(0)
        groups = finals.groupby('space')
        df_h=pd.DataFrame()
        chromosome_list=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
        #chromosome_list=['chr1']
        for  v in chromosome_list:
            df2=groups.get_group(v).reset_index()
            df2= df2.drop(['index'],axis=1)
            df2.loc[:,'z']=df2['zscore'].rolling(window=10,min_periods=1,center=True).mean()
            df2=index_p(df2)
            #print(df2)
            #print(len(df2),v,df2.width[1])
            
            df_duplication=df2[df2.z>=1.2].reset_index()
            df_deletion=df2[df2.z<=-1.2].reset_index()
            continous_bin(df2,df_duplication,'duplication')
            continous_bin(df2,df_deletion,'Deletion')
            df_h=pd.concat([df_h,df2],axis=0)
    df_c=pd.concat([df_c,df_h],axis=0)
group_whole=df_c.groupby('space')
rowlength = group_whole.ngroups//2
fig,axes = plt.subplots(figsize=(100,20),nrows=2, ncols=rowlength,gridspec_kw=dict(hspace=1))
chromosome_list=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
for g,ax in zip(chromosome_list,axes.flatten()):
    group=group_whole.get_group(g)
    flights = group.pivot('width','index_plot',"z")
    flights.sort_index(ascending=False, inplace=True)
    g1=sns.heatmap(flights.bfill(axis=1), cmap='RdBu',vmax=3,vmin=-3,square=True,xticklabels=False, cbar=False,ax=ax)
    g1.set_aspect(1./ax.get_data_ratio())
    g1.set_title(g,fontdict=30)
    g1.set_xlabel('#bins')
    g1.set_ylabel('size of bins')
    fig = g1.get_figure()
    fig.savefig('different_bin_size_plot.png')
        

