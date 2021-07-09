#import xlrd  #read the excel data
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import plot,savefig
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

def import_all_data_trio_index(all_data,trio_3_person,normal_trio_list,disease_fa_list,disease_ma_list):
    all_data = open(all_data,'r')
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
          print ("ERROR:Trio %s is not  normal or disease"%i)
    #print len(trio_3_person_nor_index),len(trio_3_person_fa_index),len(trio_3_person_ma_index)
    all_data.close()
    return trio_3_person_nor_index,trio_3_person_fa_index,trio_3_person_ma_index

def met_nor_result(all_data,trio_3_person_norm_index,outfile1,trio_3_person_fa_index,outfile2,trio_3_person_ma_index,outfile3,outfile4):
    all_data = open(all_data,'r')
    #outfile_norm = open(outfile1,'w')
    #outfile_fa = open(outfile2,'w')
    #outfile_ma = open(outfile3,'w')
    outfile_in = open(outfile4,'w')
    lines = all_data.readlines()
    #normal_matrix = [['11-1',0],['-1-11',0]] father_matrix = [['1-11',0],['-11-1',0]] mather_matrix = [['-111',0],['1-1-1',0]]
    normal_matrix = [[0 for col in range(2)]for row in range(len(lines)-1)]
    father_matrix = [[0 for col in range(2)]for row in range(len(lines)-1)]
    mather_matrix = [[0 for col in range(2)]for row in range(len(lines)-1)]
    inhert_matrix = [[0 for col in range(8)]for row in range(len(lines)-1)]
    normal_matrix[0][0] = '-1-11'; normal_matrix[0][1] = '11-1';
    father_matrix[0][0] = '1-11'; father_matrix[0][1] = '-11-1';
    mather_matrix[0][0] = '-111'; mather_matrix[0][1] = '1-1-1';
    inhert_matrix[0][0] = 'high'; inhert_matrix[0][1] = 'low';
    inhert_matrix[0][2] = 'normal_-1-11'; inhert_matrix[0][3] = 'normal_11-1';
    inhert_matrix[0][4] = 'father_1-11'; inhert_matrix[0][5] = 'father_-11-1';
    inhert_matrix[0][6] = 'mather_-111'; inhert_matrix[0][7] = 'mather_1-1-1';

    n=1;
    for l in range(2,len(lines)):  #line of all_data_new 
       lines[l] = lines[l].split("\t")
       for i in range(0,len(trio_3_person_norm_index),3):
          lines[l][trio_3_person_norm_index[i]]=float(lines[l][trio_3_person_norm_index[i]])
          lines[l][trio_3_person_norm_index[i+1]]=float(lines[l][trio_3_person_norm_index[i+1]])
          lines[l][trio_3_person_norm_index[i+2]]=float(lines[l][trio_3_person_norm_index[i+2]])
          if (lines[l][trio_3_person_norm_index[i]]==-1 and lines[l][trio_3_person_norm_index[i+1]]==-1 and lines[l][trio_3_person_norm_index[i+2]]==1):
              normal_matrix[n][0]+=1
          elif (lines[l][trio_3_person_norm_index[i]]== 1 and lines[l][trio_3_person_norm_index[i+1]]==1 and lines[l][trio_3_person_norm_index[i+2]]==-1):
              normal_matrix[n][1]+=1
          else:
              pass
       for i in range(0,len(trio_3_person_fa_index),3):
          lines[l][trio_3_person_fa_index[i]]=float(lines[l][trio_3_person_fa_index[i]])
          lines[l][trio_3_person_fa_index[i+1]]=float(lines[l][trio_3_person_fa_index[i+1]])
          lines[l][trio_3_person_fa_index[i+2]]=float(lines[l][trio_3_person_fa_index[i+2]])
          if (lines[l][trio_3_person_fa_index[i]]==1 and lines[l][trio_3_person_fa_index[i+1]]==-1 and lines[l][trio_3_person_fa_index[i+2]]==1):
              father_matrix[n][0]+=1
          elif (lines[l][trio_3_person_fa_index[i]]==-1 and lines[l][trio_3_person_fa_index[i+1]]==1 and lines[l][trio_3_person_fa_index[i+2]]==-1):
              father_matrix[n][1]+=1
          else:
              pass
       for i in range(0,len(trio_3_person_ma_index),3):
          lines[l][trio_3_person_ma_index[i]]=float(lines[l][trio_3_person_ma_index[i]])
          lines[l][trio_3_person_ma_index[i+1]]=float(lines[l][trio_3_person_ma_index[i+1]])
          lines[l][trio_3_person_ma_index[i+2]]=float(lines[l][trio_3_person_ma_index[i+2]])    
          if (lines[l][trio_3_person_ma_index[i]]==-1 and lines[l][trio_3_person_ma_index[i+1]]==1 and lines[l][trio_3_person_ma_index[i+2]]==1):
              mather_matrix[n][0]+=1
          elif (lines[l][trio_3_person_ma_index[i]]==1 and lines[l][trio_3_person_ma_index[i+1]]==-1 and lines[l][trio_3_person_ma_index[i+2]]==-1):
              mather_matrix[n][1]+=1
          else:
              pass
       if(((normal_matrix[n][0]>=12 and normal_matrix[n][1]==0) or (normal_matrix[n][1]>=3 and normal_matrix[n][0]/normal_matrix[n][1]>=3)) and ((father_matrix[n][0]>=1 and father_matrix[n][1]==0) or (father_matrix[n][1]!=0 and father_matrix[n][0]/father_matrix[n][1]>=3)) and ((mather_matrix[n][0]>=1 and mather_matrix[n][1]==0) or (mather_matrix[n][1]!=0 and mather_matrix[n][0]/mather_matrix[n][1]>=3))):
          inhert_matrix[n][0]=lines[l][6];
          inhert_matrix[n][2]=normal_matrix[n][0];inhert_matrix[n][3]=normal_matrix[n][1];
          inhert_matrix[n][4]=father_matrix[n][0];inhert_matrix[n][5]=father_matrix[n][1];
          inhert_matrix[n][6]=mather_matrix[n][0];inhert_matrix[n][7]=mather_matrix[n][1];
       elif (((normal_matrix[n][1]>=12 and normal_matrix[n][0]==0) or (normal_matrix[n][0]>=3 and normal_matrix[n][1]/normal_matrix[n][0]>=3)) and ((father_matrix[n][1]>=1 and father_matrix[n][0]==0) or (father_matrix[n][0]!=0 and father_matrix[n][1]/father_matrix[n][0]>=3)) and ((mather_matrix[n][1]>=1 and mather_matrix[n][0]==0) or (mather_matrix[n][0]!=0 and mather_matrix[n][1]/mather_matrix[n][0]>=3))):
          inhert_matrix[n][1]=lines[l][6]
          inhert_matrix[n][2]=normal_matrix[n][0];inhert_matrix[n][3]=normal_matrix[n][1];
          inhert_matrix[n][4]=father_matrix[n][0];inhert_matrix[n][5]=father_matrix[n][1];
          inhert_matrix[n][6]=mather_matrix[n][0];inhert_matrix[n][7]=mather_matrix[n][1];
       else:
          pass
       n+=1;
    #np.savetxt(outfile_norm,normal_matrix,fmt='%s')
    #np.savetxt(outfile_fa,father_matrix,fmt='%s')
    #np.savetxt(outfile_ma,mather_matrix,fmt='%s')
    np.savetxt(outfile_in,inhert_matrix,fmt='%s')
    all_data.close()
    #outfile_norm.close()
    #outfile_fa.close()
    #outfile_ma.close()
    #return met_result_matrix

def figure(met_result):
    x=range(1,28)
    col_sum = []
    sum=0
    for i in range(0,27):
        
        col_sum.append(sum(met_result[:,i]))
    plot(x,col_sum)
    savefig('single_probe.jpg')

if __name__ == "__main__":
    
    trio_name = import_trio('trio_name.txt')
    sample_name = import_sample('sample_name.txt')
    normal_trio_name = import_trio('normal_trio.txt')
    disease_fa_name = import_trio('disease_father_trio.txt')
    disease_ma_name = import_trio('disease_mather_trio.txt')
    trio_3_person = trio_3_person(trio_name,sample_name)
    trio_3_person_nor_index,trio_3_person_fa_index,trio_3_person_ma_index =  import_all_data_trio_index('all_data_new.txt',trio_3_person,normal_trio_name,disease_fa_name,disease_ma_name)
    result_nor_file = 'met_single_probe_norm_multi_output.txt'
    result_fa_file = 'met_single_probe_father_multi_output.txt'
    result_ma_file = 'met_single_probe_mather_multi_output.txt'
    result_in_file = 'father_mather_normal_candidate_3fold_12_1_2_multi_output.txt'
    met_nor_result('all_data_new.txt',trio_3_person_nor_index,result_nor_file,trio_3_person_fa_index,result_fa_file,trio_3_person_ma_index,result_ma_file,result_in_file)
    
