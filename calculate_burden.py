import math
from StringIO import StringIO
import numpy
import csv
import operator
import random
from scipy.stats import ks_2samp
from scipy.stats import stats
import itertools
import gzip
import collections
import re

# CONSTANTS would go here, if any - eg, CONSTANT={'value1': num, 'value2':num2, etc.}
FNAME_MED_TO_CID="static/label_mapping.tsv.gz"
FNAME_CID_TO_AE="static/meddra_freq_parsed.tsv.gz"

def check_medlist(variables): 
    """

    Check a list of meds and sum up burden

    Variables is a dictionary passed from server.py
 
    """

# take csv list passed of meds
    complist=[x.strip() for x in variables['Druglist'].replace('\n',',').split(',')]  

# map to CID for lookup

    matcher_dict={}
    backmatch_dict={}
    matchedcid=[]
    matcheddrugs=[]
    matcheddrug_to_cid={}

    with gzip.open(FNAME_MED_TO_CID) as gzfile:
        medpairs=csv.reader(gzfile,delimiter='\t')    
        for row in medpairs:

            if (row[0] in complist) or (row[1] in complist):
                if not row[1] in matcher_dict :
                    matcher_dict[row[1]]= row[3]
                    backmatch_dict[row[3]]=row[1]
                    matcheddrugs.append(row[1])
                    matchedcid.append(row[3])
 
# make aelist from comparator
    if variables['Comparator']=="All CNS":
        aelist=["Confusion","Delirium","Agitation"]
    else : aelist=[variables['Comparator']]    
            
# read in AE megafile - it's gzipped...
    property_dict={}
    list_by_ae={}
    list_by_drug={}
    
    with gzip.open(FNAME_CID_TO_AE) as gzfile:
        drug_properties=csv.reader(gzfile,delimiter='\t')
        for row in drug_properties:
            if row[9]=="PT" and not(row[5]=="placebo") :
                # need to reassign frequencies (exact) to categories
                freqtemp=row[6]
                   #postmarketing, rare, infrequent, frequent, or exact #
                if freqtemp=="postmarketing": freqnum=.01
                elif freqtemp=="rare": freqnum=.01
                elif freqtemp=="infrequent": freqnum=.05
                elif freqtemp=="frequent": freqnum=.25
                elif freqtemp=="potential": freqnum=.001
                elif "-" in freqtemp or "to" in freqtemp: freqnum=.01  #fix this!!
                else : freqnum=float(freqtemp[:-1])/100
                if row[0] in property_dict:
                    if not(row[11] in property_dict[row[0]]):
                        property_dict[row[0]][row[11]]=freqnum
                else :
                    property_dict[row[0]]={}
                    property_dict[row[0]][row[11]]=freqnum
                    
    # now calculate burden score
    list_by_ae={}
    list_by_drug={}
    # loop over all AE's in list to query
    for cid in matchedcid:
        for ae in aelist:
            if ae in property_dict[cid] :
                if ae in list_by_ae:
                    list_by_ae[ae][backmatch_dict[cid]]=property_dict[cid][ae]
                else :
                    list_by_ae[ae]={}
                    list_by_ae[ae][backmatch_dict[cid]]=property_dict[cid][ae] 
                    
                if backmatch_dict[cid] in list_by_drug:
                    list_by_drug[backmatch_dict[cid]][ae]=property_dict[cid][ae] 
                else:
                    list_by_drug[backmatch_dict[cid]]={}
                    list_by_drug[backmatch_dict[cid]][ae]=property_dict[cid][ae] 

    
    #if we want to add a warning for high placebo rate, add it here.
    
    # now sum up freq burden or risk, by AE
    ae_score={}
    for ae in aelist:
        aeburden=0
        if ae in list_by_ae:
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
        'listed_CID': "temp",
        'list_by_drug':list_by_drug,
        'list_by_ae':list_by_ae,
        'annotation_by_drug':annotation_by_drug,         
        'ae_score':ae_score,
        'drug_score':drug_score,
        'ae_total':ae_total,
    }


