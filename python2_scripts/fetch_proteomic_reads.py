import shelve
import collections

def merge_dicts(dict1,dict2):
    for readlen in dict2:
        if readlen not in dict1:
            dict1[readlen] = dict2[readlen]
        else:
            for pos in dict2[readlen]:
                if pos in dict1[readlen]:
                    dict1[readlen][pos] += dict2[readlen][pos]
                else:
                    dict1[readlen][pos] = dict2[readlen][pos]
    return dict1
            


def get_reads(min_read, max_read, tran, user_files, organism,tranlen,minscore):

    master_dict = {}
    for i in range(0,tranlen):
        master_dict[i] = 0



    for filename in user_files:
        openshelve = shelve.open("/home/DATA/GWIPS_viz/proteomics_shelves/{}/{}.shelf".format(organism, filename))
        if tran in openshelve:
            for readlen in openshelve[tran]:
                if readlen > min_read:
                    for pos in openshelve[tran][readlen]:
                        for score in openshelve[tran][readlen][pos]:
                            if score >= minscore:
                                count = openshelve[tran][readlen][pos][score]
                                for i in range((pos),(pos)+readlen,3):
                                    if i not in master_dict:
                                        master_dict[i] = 0
                                    master_dict[i] += count

            
   

    sorted_master_dict = collections.OrderedDict()
    for key in sorted(master_dict.keys()):
        sorted_master_dict[key] = master_dict[key]

    
    return sorted_master_dict
        

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


