#import xlrd  #read the excel data
import numpy as np
from itertools import islice
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
matplotlib.use('Agg')

#fname = "*.xls"
#bk = xlrd.open_workbook(fname)
#shxrange = range(bk.nsheets)
#sh = bk.sheet_by_name("Sheet1")
#nrows = sh.nrows
#ncols = sh.ncols
# extract the data

#row_list = []
#for i in range(1,nrows)
#    row_data = sh.row_values(i)
#    row_list.append(row_data)

def import_trio(trio_name_file):
    trio_name = open(trio_name_file,'r')
    lines = trio_name.readlines()
    trio_name_list=[]
    for line in lines:
        line = line.strip()
        trio_name_list.append(line)
    trio_name.close()
    return trio_name_list
#print trio_name_list

def import_sample(sample_name_file):
    sample_name = open(sample_name_file,'r')
    lines = sample_name.readlines()
    sample_name_list=[]
    for line in lines:
       line = line.strip()
       sample_name_list.append(line)
    sample_name.close()
    return sample_name_list
#print len(sample_name_list)

def trio_3_person(trio_name_list,sample_name_list):
    trio_3_person=[]
    for i in trio_name_list:
       m=0
       x=i+"-1"
       y=i+"-2"
       z=i+"-3"
       if x in sample_name_list:
            m=m+1
       if y in sample_name_list:
            m=m+1
       if z in sample_name_list:
            m=m+1
       if m==3 :
           trio_3_person.append(i)
    return trio_3_person
#print new_trio
#print len(new_trio)

def import_all_data_trio_index(all_data_file,trio_3_person,normal_trio_list,disease_fa_list,disease_ma_list):
    all_data = open(all_data_file,'r')
    lines = all_data.readlines()
    trio_3_person_nor_index=[]
    trio_3_person_fa_index=[]
    trio_3_person_ma_index=[]
    head = lines[0].split("\t")
    for i in trio_3_person:
       x=i+"-1"
       y=i+"-2"
       z=i+"-3"
       if i in normal_trio_list:
          trio_3_person_nor_index.append(head.index(x)+1)
          trio_3_person_nor_index.append(head.index(y)+1)   
          trio_3_person_nor_index.append(head.index(z)+1)    
       elif i in disease_fa_list:
          trio_3_person_fa_index.append(head.index(x)+1)
          trio_3_person_fa_index.append(head.index(y)+1)
          trio_3_person_fa_index.append(head.index(z)+1)
       elif i in disease_ma_list:
          trio_3_person_ma_index.append(head.index(x)+1)
          trio_3_person_ma_index.append(head.index(y)+1)
          trio_3_person_ma_index.append(head.index(z)+1)
       else:
          print("ERROR:Trio %s is not  normal or disease"%i)
    print(len(trio_3_person_nor_index),len(trio_3_person_fa_index),len(trio_3_person_ma_index))
    all_data.close()
    return trio_3_person_nor_index,trio_3_person_fa_index,trio_3_person_ma_index

def met_nor_result(all_data,trio_3_person_norm_index,outfile1,trio_3_person_fa_index,outfile2,trio_3_person_ma_index,outfile3,outfile4):

    outfile_fa = open(outfile2,'w')
    all_modes=['111', '110', '11-1', '101', '100', '10-1', '1-11', '1-10', '1-1-1', '011', '010', '01-1', '001', '000', '00-1', '0-11', '0-10', '0-1-1',  '-111', '-110', '-11-1', '-101', '-100', '-10-1', '-1-11', '-1-10', '-1-1-1']
    fa_matrix = [[0 for col in range(27)]for row in range(252699)]
    fa_matrix[0][0]='111'
    fa_matrix[0][1]='110'
    fa_matrix[0][2]='11-1'
    fa_matrix[0][3]='101'
    fa_matrix[0][4]='100'
    fa_matrix[0][5]='10-1'
    fa_matrix[0][6]='1-11'
    fa_matrix[0][7]='1-10'
    fa_matrix[0][8]='1-1-1'
    fa_matrix[0][9]='011'
    fa_matrix[0][10]='010'
    fa_matrix[0][11]='01-1'
    fa_matrix[0][12]='001'
    fa_matrix[0][13]='000'
    fa_matrix[0][14]='00-1'
    fa_matrix[0][15]='0-11'
    fa_matrix[0][16]='0-10'
    fa_matrix[0][17]='0-1-1'
    fa_matrix[0][18]='-111'
    fa_matrix[0][19]='-110'
    fa_matrix[0][20]='-11-1'
    fa_matrix[0][21]='-101'
    fa_matrix[0][22]='-100'
    fa_matrix[0][23]='-10-1'
    fa_matrix[0][24]='-1-11'
    fa_matrix[0][25]='-1-10'
    fa_matrix[0][26]='-1-1-1'

    all_data = open(all_data,'r')
    n = 1
    for line in islice(all_data,2,None):
        line = line.strip().split("\t")
        for i in range(0,len(trio_3_person_fa_index),3):
          #mode=lines[l][trio_3_person_norm_index[i]]+lines[l][trio_3_person_norm_index[i+1]]+lines[l][trio_3_person_norm_index[i+2]]
          mode=line[trio_3_person_fa_index[i]]+line[trio_3_person_fa_index[i+1]]+line[trio_3_person_fa_index[i+2]]
          #111 110 11-1 101 100 10-1 1-11 1-10 1-1-1 011 010 01-1 001 000 00-1 0-11 0-10 0-1-1  -111 -110 -11-1 -101 -100 -10-1 -1-11 -1-10 -1-1-1 
          t=all_modes.index(mode)
          fa_matrix[n][t] += 1
        n = n+1
    np.savetxt(outfile_fa,fa_matrix,fmt='%s')
    #np.savetxt(outfile_in,inhert_matrix,fmt='%s')
    all_data.close()
    outfile_fa.close()
    return fa_matrix,all_modes

def figure(met_result_matrix,all_modes):

    data = np.array(met_result_matrix[1:])
    data_sum = np.sum(data,axis=0)
    fig = plt.figure(figsize=(10.0,5.0))
    plt.xticks(fontsize=6)
    plt.bar(all_modes,data_sum,0.4,color="green")
    #sns_plot = sns.heatmap(data)
    plt.xlabel("Inheritance Mode")
    plt.ylabel("Count")
    plt.title("Trio (Father)")

    plt.show()
    plt.savefig("Father_all_modes_bar.jpg",dpi=300)
    plt.savefig("Father_all_modes_bar.pdf", bbox_inches='tight')

if __name__ == "__main__":
    
    trio_name = import_trio('trio_name.txt')
    sample_name = import_sample('sample_name.txt')
    normal_trio_name = import_trio('normal_trio.txt')
    disease_fa_name = import_trio('disease_father_trio.txt')
    disease_ma_name = import_trio('disease_mather_trio.txt')
    trio_3_person = trio_3_person(trio_name,sample_name)
    trio_3_person_nor_index,trio_3_person_fa_index,trio_3_person_ma_index =  import_all_data_trio_index('all_data_new_header.txt',trio_3_person,normal_trio_name,disease_fa_name,disease_ma_name)
    result_nor_file = 'met_single_probe_norm_all_modes.txt'
    result_fa_file = 'met_single_probe_father_all_modes.txt'
    result_ma_file = 'met_single_probe_mather_all_modes.txt'
    result_in_file = 'father_mather_normal_candidate_3fold_12_1_2_all_modes.txt'
    met_result_matrix,all_modes = met_nor_result('all_data_new.txt',trio_3_person_nor_index,result_nor_file,trio_3_person_fa_index,result_fa_file,trio_3_person_ma_index,result_ma_file,result_in_file)
    #met_nor_result('tmp_all_data_new.txt',trio_3_person_nor_index,result_nor_file,trio_3_person_fa_index,result_fa_file,trio_3_person_ma_index,result_ma_file,result_in_file)
    figure(met_result_matrix,all_modes)
