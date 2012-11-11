import math
from StringIO import StringIO
import numpy
import csv
import operator
import random
from scipy.stats import stats
import itertools
import gzip  #needed for gzipped files
import collections
import re   # needed for splitting strings
import os   # needed for file names

# CONSTANTS would go here, if any - eg, CONSTANT={'value1': num, 'value2':num2, etc.}
FNAME_MED_TO_CID=os.path.dirname(__file__)+"/static/label_mapping.tsv.gz"
FNAME_CID_TO_AE=os.path.dirname(__file__)+"/static/meddra_freq_parsed.tsv.gz"

def check_medlist(variables): 
    """

    Check a list of meds and sum up burden

    Variables is a dictionary passed from server.py
 
    """

# take csv list passed of meds
    #complist=[x.strip() for x in variables['Druglist'].replace('\n',',').split(',')]  
    complist=[x for x in variables['Druglist'].replace('\n',',').split(',')]  
    complist=filter(None,complist)
    complist=[y.lstrip(" ").split(" ")[0] for y in complist]
    print("complist",complist)
# map to CID for lookup

    matcher_dict={}
    backmatch_dict={}
    matchedcid=[]
    matcheddrugs=[]
    matcheddrug_to_cid={}

    with gzip.open(FNAME_MED_TO_CID) as gzfile:
        medpairs=csv.reader(gzfile,delimiter='\t')    
        for row in medpairs:

            if ((row[0] in complist) or (row[1] in complist)):
                if (not row[1] in matcher_dict) and (not row[1] in matcheddrugs and not row[0] in matcheddrugs) :
                    matcher_dict[row[1]]= row[3]
                    backmatch_dict[row[3]]=row[1]
                    matcheddrugs.append(row[1])
                    matchedcid.append(row[3])
    print("matchedlist:",matcher_dict)
# make aelist from comparator
    if variables['Comparator']=="Psychiatry":
        aelist= load_aefilelist("CNS_psychiatric.txt")  
    elif variables['Comparator']=="Cognition":
        aelist=load_aefilelist("CNS_cognition.txt")
    elif variables['Comparator']=="Other Neurologic":
        aelist=load_aefilelist("CNS_other_neurologic.txt")
    elif variables['Comparator']=="All CNS":
        aelist=load_aefilelist("CNS_full.txt")
    elif variables['Comparator']=="Fall Risk":
        aelist=load_aefilelist("CNS_fall_risk.txt")     
    else : aelist=[variables['Comparator']]    
            
# read in AE megafile - it's gzipped...
    property_dict={}
    list_by_ae={}
    list_by_drug={}
    
    with gzip.open(FNAME_CID_TO_AE) as gzfile:
        drug_properties=csv.reader(gzfile,delimiter='\t')
        for row in drug_properties:
            if (row[9]=="PT" and not(row[5]=="placebo")) and not row[0]=='' :
                # need to reassign frequencies (exact) to categories
                   #postmarketing, rare, infrequent, frequent, or exact #
                freqnum=parse_frequency(row[6])
                if variables['Option_1']==0 and freqnum>0: freqnum=.1   # treats any ae as freq=.1
                if freqnum>0 and row[0] in property_dict:
                    if not(row[11] in property_dict[row[0]]):
                        property_dict[row[0]][row[11]]=freqnum
                elif freqnum>0 :
                    property_dict[row[0]]={}
                    property_dict[row[0]][row[11]]=freqnum
                    
    # now remove drugs which are not in dictionary
    drug_not_in_dictionary=[]
    for cid in matchedcid:
        if not property_dict.has_key(cid): 
            drug_not_in_dictionary.append(backmatch_dict[cid])
            matchedcid.remove(cid)
            matcheddrugs.remove(backmatch_dict[cid])
                
    #now figure out p450 interactions!
    modifiers_p450={}
    subs_p450={}
    multiplier={}
    
    p450_results=map_p450(matcheddrugs)
    modifiers_p450=p450_results['modifiers']
    subs_p450=p450_results['substrates']
    multiplier=p450_results['multiplier']
    
    print("mods_p450",modifiers_p450)
    
    # now calculate burden score
    list_by_ae={}
    list_by_drug={}

    # loop over all AE's in list to query
    for cid in matchedcid:
        for ae in aelist:
            if not property_dict.has_key(cid): drug_not_in_dictionary.append(backmatch_dict[cid])
            elif ae in property_dict[cid] :
                if ae in list_by_ae:
                    list_by_ae[ae][backmatch_dict[cid]]=property_dict[cid][ae]*multiplier[backmatch_dict[cid]]
                else :
                    list_by_ae[ae]={}
                    list_by_ae[ae][backmatch_dict[cid]]=property_dict[cid][ae]*multiplier[backmatch_dict[cid]] 
                    
                if backmatch_dict[cid] in list_by_drug:
                    list_by_drug[backmatch_dict[cid]][ae]=property_dict[cid][ae]*multiplier[backmatch_dict[cid]] 
                else:
                    list_by_drug[backmatch_dict[cid]]={}
                    list_by_drug[backmatch_dict[cid]][ae]=property_dict[cid][ae]*multiplier[backmatch_dict[cid]] 
    print("not_in_dict",drug_not_in_dictionary)
    
    #if we want to add a warning for high placebo rate, add it here.

    
    # now sum up freq burden or risk, by AE
    print("show list_by_ae",list_by_ae)
    ae_score={}
    for ae in list_by_ae:
        aeburden=0
        aeburden=sum(list_by_ae[ae].itervalues())
        ae_score[ae]=aeburden
        
    drug_score={}    
    for drug in matcher_dict:
        drugburden=0
        if drug in list_by_drug:
            drugburden=sum(list_by_drug[drug].itervalues()) 
        drug_score[drug]=drugburden
    print(drug_score)
    # now sum up overall burden (all AE's)
    ae_total=sum(ae_score.itervalues())     
               
    # here's where we can add custom annotation by drug
    #FNAME_DRUG_ANNOTATION="none"
    annotation_by_drug={}
    #for drug in matched_drugs:
    #    annotation_by_drug[drug]=annotation[drug]
    
# now return results    
    return {
        'matched_drugs': matcheddrugs,
        'modifiers_p450':modifiers_p450,
        'subs_p450':subs_p450,
        'listed_CID': "temp",
        'list_by_drug':list_by_drug,
        'list_by_ae':list_by_ae,
        'annotation_by_drug':annotation_by_drug,         
        'ae_score':ae_score,
        'drug_score':drug_score,
        'ae_total':ae_total,
    }

"""
this function categorizes AE's and writes a new gz file
"""

def load_aefilelist (aelistfile):
# take a text file of AE's, remap to a category-list, save a new gzipped file

      # read in list of AE's to match - just a text file
      # files from Victor - column 1 (0-index) is name, col 0 is #, col 3 is high-level term
    aelist2=[]   
    ae_file=csv.reader(open(os.path.dirname(__file__)+"/static/"+aelistfile),delimiter='\t')
    for row in ae_file:
       if not row[1] in aelist2:
           aelist2.append(row[1])
    return aelist2

"""
this function parses frequencies from screwy file
"""

def parse_frequency(freqtemp):
    if freqtemp=="postmarketing": freqnum=.01
    elif freqtemp=="rare": freqnum=.01
    elif freqtemp=="infrequent": freqnum=.05
    elif freqtemp=="frequent": freqnum=.25
    elif freqtemp=="potential": freqnum=.001     # potential should be nearly zero!
    elif "-" in freqtemp or "to" in freqtemp: 
        if freqtemp[-1:]=="%" : freqtemp=freqtemp[:-1]
        freqnums=re.split('-|to',freqtemp.replace(" ",""))
        freqnum=(float(freqnums[0])+float(freqnums[1]))/200
    else : freqnum=float(freqtemp[:-1])/100
    
    return freqnum

"""
this function takes a list of meds and figures out p450 status
returns a dictionary of drug->AE multiplier (ie, 2=double AE's, 0.5=half-normal AE's)
"""

def map_p450(list_of_meds):
    CYP450_MODIFIERS="cyp450_mods.txt"
    CYP450_SUBSTRATES="cyp450_substrates.txt"
    p450_substrates={}
    p450_modifiers={}
    
    multiplier={}
    for med in list_of_meds:
        multiplier[med]=1
        
    p450_panel={'1A2':1, '2B6':1, '2C8':1, '2C9': 1, '2C19':1, '2D6': 1, '2E1':1, '3A457': 1}

    # now read in inhibitors/inducers: format is DRUGNAME P450 MULTIPLIER
    cyp450_mods=csv.reader(open(os.path.dirname(__file__)+"/static/"+CYP450_MODIFIERS),delimiter='\t')
    for row in cyp450_mods:
       if row[0] in list_of_meds:
           p450_panel[row[1]]=p450_panel[row[1]]*float(row[2])                     #modify p450 status panel
           p450_modifiers[row[0]]=row[1]                                    #add to list of modifiers
    
    # now generate multipliers for med side effects
    cyp450_subs=csv.reader(open(os.path.dirname(__file__)+"/static/"+CYP450_SUBSTRATES),delimiter='\t')
    for row in cyp450_subs:
       if row[0] in list_of_meds:
           multiplier[row[0]]=multiplier[row[0]]*p450_panel[row[1]]         #lookup p450 key, multiply  
           p450_substrates[row[0]]=row[1]                                   #add to list of substrates

    return {
        'modifiers': p450_modifiers,
        'substrates':p450_substrates,
        'multiplier': multiplier
    }