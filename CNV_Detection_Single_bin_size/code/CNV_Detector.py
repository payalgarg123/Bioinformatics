
# coding: utf-8

# In[19]:


import numpy as np
import glob
import pandas as pd
import os
from matplotlib import pyplot as plt
import seaborn as sns
from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.models.widgets import Tabs, Panel
from bokeh.layouts import column
from bokeh.models.glyphs import Text
from bokeh.plotting import reset_output
from bokeh.models import HoverTool
from bokeh.models.widgets import Select
import argparse
import time as tm


******************************************************************************************Required Arguments********************************************************************************************************************************************************


parser=argparse.ArgumentParser(description="CNV deTection and plotting")
parser.add_argument('-N','--normal_file_location',type=str,required=True,help='location of the normal file to be in this program')
parser.add_argument('-M','--normalization_method',type=str,metavar=['gc','reads'],required=True,help='normalization method')
parser.add_argument('-R','--reference_database',type=str,required=True,help='reference database of normalization files')
parser.add_argument('-Z','--z_score_data',type=str,required=True,help='z score database of normalization file')
parser.add_argument('-E','--width_size',type=int,help='difference between start bin and end bin in kb')
parser.add_argument('-W','--window_size',metavar="",type=int,help='rolling size')
parser.add_argument('-S','--min_cnv_size',type=int,help='minimum cnv periods')
parser.add_argument('-O','--output_file',type=str,required=True,help="final report of continous deletion and duplication")
parser.add_argument('-P','--plotting_file',type=str,required=True,help="plots")

args=parser.parse_args()


startTime=tm.time()

output_file(args.plotting_file) #saved html file of plots


#########################################################################################################1.Function#################################################################################################################################################


'''' For normalization  used gc method. In this method we can used gc values to normalized count
data_1_22= data[data['space'].str.contains("chrX|chrY")!=True] data_1_22 is not contain chrX and ChrY
d= data_1_22['reads'].mean()
dgc=data_1_22.loc[(data_1_22['gc'] >= df.gc[i]-0.1) & (data_1_22['gc'] <= df.gc[i]+0.1), 'reads'].mean() 
ARC_gc_count(df.reads[i]*(d/dgc))------ method'''

if args.normalization_method=='gc':
    def normalization_reads(data,df):
        data_1_22= data[data['space'].str.contains("chrX|chrY")!=True]
        d= data_1_22['reads'].mean()
        ARC=[]
        for i in range (len(df)):
           if(df.gc[i]==-1):
               dgc=d
           else:
              dgc=data_1_22.loc[(data_1_22['gc'] >= df.gc[i]-0.1) & (data_1_22['gc'] <= df.gc[i]+0.1), 'reads'].mean()
           ARC.append(df.reads[i]*(d/dgc))
        df['data_r']= ARC
        return(df)

  ''' for normalization used read method.In this method we can use the reads value to create a new column normalization read count
  df['data_r']=(df['reads']/data_1_22['reads'].sum())*1000-------Method'''
  
elif args.normalization_method=='reads':
        def normalization_reads(data,df):
            data_1_22= data[data['space'].str.contains("chrX|chrY")!=True]
            df['data_r']=(df['reads']/data_1_22['reads'].sum())*1000
            return (df)
else:
    print(error)

###############################################################################################################2.Function############################################################################################################################################

'''continous function in this function we can check copy number variations.
it is divided in two parts deletion and duplication'''


def continous_bin(df,data,condition,condition1,color,cont_color):
    i=0
    count=0
    data_list=[]
    cont_data=pd.DataFrame()
    un_cont_data=pd.DataFrame()
    while i<(len(data)-1):
        #print("I am in continous_bin loop")
        start_bin=data.start[i]
        end_bin=data.start[i+1]
        i+=1
        
        
        
        if end_bin-start_bin==args.width_size:
            #print(start_bin,end_bin)
            count+=1
             
        elif count>=args.min_cnv_size:
            #print(count)
            real_start=data.start[i-count-1]
            real_end=data.end[i-1]
            con_bin_mean=data[(data['start'] >= real_start) & (data['start'] <= data.start[i-1])]
            cont_data=pd.concat([cont_data,con_bin_mean])
            graph_plotting(cont_data,cont_data.z,cont_color,condition1)
            total_c_mean=con_bin_mean['z'].mean()
            data_list.append([condition,data.space[i],real_start,real_end,total_c_mean])
            count=0
        
        else:
            real_start=data.start[i-count-1]
            real_end=data.end[i-1]
            un_con=data[(data['start'] >=real_start ) & (data['start'] <= data.start[i-1])]
            un_cont_data=pd.concat([un_cont_data,un_con])
            graph_plotting(un_cont_data,un_cont_data.z,color,condition)
            count=0
        
        
        if i== len(data)-1:
            if count>=args.min_cnv_size:
                real_start=data.start[i-count]
                real_end=data.end[i]
                con_bin_mean=data[(data['start'] >= real_start) & (data['start'] <= data.start[i])]
                cont_data=pd.concat([cont_data,con_bin_mean])
                un_cont_data=un_cont_data
                graph_plotting(cont_data,cont_data.z,cont_color,condition1)
                #print(un_cont_data)
                total_c_mean=con_bin_mean['z'].mean()
                data_list.append([condition,data.space[i],real_start,real_end,total_c_mean])
                Myfile=open(args.output_file,'a')
                for lines in data_list:
                    Myfile.write(str(lines))
                    Myfile.write(" ")
                Myfile.close()
                count=0
                
        
            else:
                #print("I was here too")
                real_start=data.start[i-count]
                real_end=data.end[i]
                un_con=data[(data['start'] >=real_start ) & (data['start'] <= data.start[i])]
                un_cont_data=pd.concat([un_cont_data,un_con])
                cont_data=cont_data
                data_list=data_list
                graph_plotting(un_cont_data,un_cont_data.z,color,condition)
                Myfile=open(args.output_file,'a')
                for lines in data_list:
                    Myfile.write(str(lines))
                    Myfile.write(" ")
                Myfile.close()
				
                count=0
				
###############################################################################################3.Function############################################################################################################################################################

'''Bokeh library used in this function '''
            
           
        
def graph_plotting(data_graph,y_axis,color,condition):
    p.title.text_font_size = '30pt'
    p.title.align = 'center'
    p.outline_line_color = "navy"
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_alpha = 0.5
    p.circle(data_graph.size1,y_axis,size=4,line_color=color, fill_color=color, fill_alpha=0.9,legend=condition)
    p.legend.location = 'top_right'
    p.legend.click_policy = 'hide'             
        
            

*********************************************************************************************************************************************************************************************************************************************************************


data = pd.read_csv(args.normal_file_location, sep='\t', usecols=['space','start','end','gc','reads'])
data1 = pd.read_csv(args.reference_database, sep='\t', usecols=['space','nor_mean','nor_std']) #reference_data_base 
group = data.groupby('space')
groups=data1.groupby('space')
data['size1']=data.end/1000000
chromosome_list=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
z_score_data =pd.DataFrame()
plots=[]
for  v in chromosome_list:
    #print(v)
    df=group.get_group(v).reset_index()
    df1=groups.get_group(v).reset_index()
    df= df.drop(['index'],axis=1)
    normalization_reads(data,df)
	#create figure for plots
    p=figure(plot_width=800, plot_height=700, x_axis_label='size', y_axis_label='z',x_range=([df.size1.min(),df.size1.max()]),y_range=(-4,5),title='relation between z and size of bin',tools=['box_select', 'reset', 'save','wheel_zoom','zoom_in','pan'])
    p.add_tools(HoverTool()) #hovertool 
    df['zscore']=((df['data_r'] - df1['nor_mean'])/df1['nor_std']).fillna(0)
    df.loc[:,'z']=df['zscore'].rolling(window=args.window_size,min_periods=1,center=True).mean()
    z_score_data=pd.concat([z_score_data,df])
    df_normal=df[(df.z<=1.5) & (df.z>=-1.5)]
    graph_plotting(df_normal,df_normal.z,'black','normal') #plotting function
    df_duplication=df[df.z>=1.5].reset_index()
    df_deletion=df[df.z<=-1.5].reset_index()
    continous_bin(df,df_duplication,'duplication','con_duplication','blue','purple')
    continous_bin(df,df_deletion,'deletion','con_Deletion','pink','red')
    p_panels = Panel(child=p, title=v) #create panels
    plots.append(p_panels)
tabs1 = Tabs(tabs=plots,width =820) #create tabs
show(tabs1) #show plots
z_score_data.to_csv(args.z_score_data) #saved zscore database
totalTime=tm.time() - startTime
print(totalTime/60)
reset_output()

