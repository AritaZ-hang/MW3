from multiprocess import Pool
import pybedtools
import pandas as pd
import os
import sys
import pysam
#set CPU num
cpu = 24

dir=sys.argv[1]
os.chdir(dir)

#read filenames
filenames = os.listdir()

#subscript
loc = list(range(0, len(filenames)))

#total num
len_file = len(filenames)

def Merge(dict1, dict2): 
    res = {**dict1, **dict2} 
    return res 
    
def correlation(set_a,set_b):
    unions = len(set_a.union(set_b))
    intersections = len(set_a.intersection(set_b))
    return 1. * intersections / unions


def get_dic(file):
    bf = pysam.AlignmentFile(file, 'rb', check_sq=False)
    bc_list = []
    for line in bf:
        tn = line.get_tag('KB')
        bc_list.append(tn)
    bc=set(bc_list)
    dic = {file[13:31]:bc}
    return(dic)
    
pool_file = Pool(cpu)
dic = pool_file.map(func = get_dic, iterable = filenames)
pool_file.close()
pool_file.join()

dic_all={}  
for i in dic:
    dic_all=Merge(dic_all,i)    

def calcuateJaccard_Pool(i):
    global dic_all
    global len_file
    all_jaccard = pd.DataFrame(columns=['match_a', 'match_b','jaccard'])
    if i< len_file:
        for j in range(i+1,len_file):
            a=list(dic_all.keys())[i]
            b=list(dic_all.keys())[j]
            set_a=dic_all[a]
            set_b=dic_all[b]
            result = correlation(set_a,set_b)
            tmp={'match_a': a,'match_b':b,'jaccard':result}
            tmp=pd.DataFrame.from_dict(tmp,orient='index').T
            all_jaccard=pd.concat([all_jaccard, tmp])
        return(all_jaccard)

pool_jaccard = Pool(cpu)
temp_res = pool_jaccard.map(func = calcuateJaccard_Pool, iterable = loc)
pool_jaccard.close()
pool_jaccard.join()
all_jaccard_res = pd.DataFrame(columns=['match_a', 'match_b','jaccard'])
for i in temp_res:
    all_jaccard_res = pd.concat([all_jaccard_res,i])
all_jaccard_res.to_csv('../2w_jaccard.csv', index = 0)
