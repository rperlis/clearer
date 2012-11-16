from StringIO import StringIO

import csv


import itertools
import gzip  #needed for gzipped files

import re   # needed for splitting strings
import os   # needed for file names

# TO-DO: 
    # confirm numbers corresponding to frequency
    # fix problem w generic/brand mixups
    
# CONSTANTS would go here, if any 
FNAME_MED_TO_CID=os.path.dirname(__file__)+"/static/label_mapping.tsv.gz" 
FNAME_CID_TO_AE=os.path.dirname(__file__)+"/static/meddra_freq_parsed_mini.tsv.gz"
CYP450_MODIFIERS="cyp450_mods.txt"
CYP450_SUBSTRATES="cyp450_substrates.txt"
    
    
"""
this function parses frequencies from screwy file
"""

def parse_frequency(freqtemp):
    if freqtemp=="postmarketing": freqnum=.0005
    elif freqtemp=="rare": freqnum=.0005
    elif freqtemp=="infrequent": freqnum=.005
    elif freqtemp=="frequent": freqnum=.05
    elif freqtemp=="potential": freqnum=.000     # potential should be nearly zero!
    elif "-" in freqtemp or "to" in freqtemp: 
        if freqtemp[-1:]=="%" : freqtemp=freqtemp[:-1]
        freqnums=re.split('-|to',freqtemp.replace(" ",""))
        freqnum=(float(freqnums[0])+float(freqnums[1]))/200
    else : freqnum=float(freqtemp[:-1])/100
    
    return freqnum

property_dict={}
with gzip.open(FNAME_CID_TO_AE) as gzfile:
    drug_properties=csv.reader(gzfile,delimiter='\t')
    for row in drug_properties:
        if (row[9]=="PT" and not(row[5]=="placebo")) and not row[0]=='' :
                # need to reassign frequencies (exact) to categories
                   #postmarketing, rare, infrequent, frequent, or exact #
            freqnum=parse_frequency(row[6])
            if freqnum>0 and row[0] in property_dict:
                if not(row[11] in property_dict[row[0]]):
                    property_dict[row[0]][row[11]]=freqnum
            elif freqnum>0 :
                property_dict[row[0]]={}
                property_dict[row[0]][row[11]]=freqnum

def check_medlist(variables): 
    """

    Check a list of meds and sum up burden

    Variables is a dictionary passed from server.py
 
    """

# take csv list passed of meds
    #complist=[x.strip() for x in variables['Druglist'].replace('\n',',').split(',')]  
    complist=[x for x in variables['Druglist'].replace('\n',',').replace('\r',',').split(',')]  
    complist=filter(None,complist)
    complist=[y.lstrip(" ").split(" ")[0] for y in complist]
    print("complist",complist)
# map to CID for lookup

    matcher_dict={}
    backmatch_dict={}
    matchedcid=[]
    matcheddrugs=[]
    matched_othername=[]

    with gzip.open(FNAME_MED_TO_CID) as gzfile:
        medpairs=csv.reader(gzfile,delimiter='\t')    
        for row in medpairs:

            gname=row[1].upper().split(" ")[0]
            bname=row[0].upper().split(" ")[0]
            if ((gname in complist) or (bname in complist)) and not gname=='':
                print("in complist: gname",gname,"bname",bname)
                if (not gname in matcher_dict) and (not gname in matcheddrugs) and (not bname in matcheddrugs) :
                    matcher_dict[gname]= row[3]
                    backmatch_dict[row[3]]=gname
                    matcheddrugs.append(gname)
                    matched_othername.append(bname)             # hack to address bname and gname switch
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
    elif variables['Comparator']=="Bleeding":
        aelist=load_aefilelist("Heme_bleeding.txt")
    elif variables['Comparator']=="Fall Risk":
        aelist=load_aefilelist("CNS_fall_risk.txt")     
    else : aelist=[variables['Comparator']]    
            
# read in AE megafile - it's gzipped...
    list_by_ae={}
    list_by_drug={}
    
# moved this reading in of dictionary to be compiled with server.
                    
    # now remove drugs which are not in dictionary
    drug_not_in_dictionary=[]
    for cid in matchedcid:
        if not property_dict.has_key(cid): 
            drug_not_in_dictionary.append(backmatch_dict[cid])
            matchedcid.remove(cid)
            matcheddrugs.remove(backmatch_dict[cid])
            del matcher_dict[backmatch_dict[cid]]
    #now figure out p450 interactions!
    modifiers_p450={}
    substrates_p450={}
    multiplier={}
    
    inhibitors_p450,inducers_p450,substrates_p450,multiplier=map_p450(matcheddrugs,matched_othername,variables['Option_2'])
    
    print("mods",modifiers_p450)
    
    # now calculate burden score
    list_by_ae={}
    list_by_drug={}

    # loop over all AE's in list to query
    for cid in matchedcid:
        for ae in aelist:
            if not property_dict.has_key(cid): drug_not_in_dictionary.append(backmatch_dict[cid])
            elif ae in property_dict[cid] :
                freqnumtemp=property_dict[cid][ae]
                if variables['Option_1']==0: freqnumtemp=.01
                if ae in list_by_ae:
                    list_by_ae[ae][backmatch_dict[cid]]=freqnumtemp*multiplier[backmatch_dict[cid]]
                else :
                    list_by_ae[ae]={}
                    list_by_ae[ae][backmatch_dict[cid]]=freqnumtemp*multiplier[backmatch_dict[cid]] 
                    
                if backmatch_dict[cid] in list_by_drug:
                    list_by_drug[backmatch_dict[cid]][ae]=freqnumtemp*multiplier[backmatch_dict[cid]] 
                else:
                    list_by_drug[backmatch_dict[cid]]={}
                    list_by_drug[backmatch_dict[cid]][ae]=freqnumtemp*multiplier[backmatch_dict[cid]] 
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
    print(make_table(list_by_drug,'drug','adverse effect'))  
    return {
        'matched_drugs': matcheddrugs,
        'mods_p450':make_table_list(inhibitors_p450,'Inhibitor','Enzyme') + make_table_list(inducers_p450,'Inducer','Enzyme'),
        'subs_p450':make_table_list(substrates_p450,'Substrate','Enzyme'),
        'list_by_drug':make_table(list_by_drug,'Drug','Adverse Effect',variables['Option_1']),
        'list_by_ae':make_table(list_by_ae,'Adverse effect','Drug',variables['Option_1']),
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

def map_p450(list_of_meds,altnames_meds,use_p450):
    CYP450_MODIFIERS="cyp450_mods.txt"
    CYP450_SUBSTRATES="cyp450_substrates.txt"
    p450_substrates={}
    p450_inhibitors={}
    p450_inducers={}
    
    p450_multiplier={}
    for med in list_of_meds:
        p450_multiplier[med]=1
        
    p450_panel={'1A2':1, '2B6':1, '2C8':1, '2C9': 1, '2C19':1, '2D6': 1, '2E1':1, '3A457': 1}

    # now read in inhibitors/inducers: format is DRUGNAME P450 MULTIPLIER
    cyp450_mods=csv.reader(open(os.path.dirname(__file__)+"/static/"+CYP450_MODIFIERS),delimiter='\t')
    for row in cyp450_mods:
       drugname=row[0].upper().split(" ")[0]
       if drugname in list_of_meds or drugname in altnames_meds:
           p450_panel[row[1]]=p450_panel[row[1]]*float(row[2])                     #modify p450 status panel
           if float(row[2])<1: p450_inducers[drugname]=row[1]                                    #add to list of modifiers
           if float(row[2])>1: p450_inhibitors[drugname]=row[1]
    # now generate multipliers for med side effects
    cyp450_subs=csv.reader(open(os.path.dirname(__file__)+"/static/"+CYP450_SUBSTRATES),delimiter='\t')
    for row in cyp450_subs:
       drugname=row[0].upper().split(" ")[0] 
       if drugname in list_of_meds:
           if use_p450==1: p450_multiplier[drugname]=p450_multiplier[drugname]*p450_panel[row[1]]         #lookup p450 key, multiply  
           p450_substrates[drugname]=row[1]
           if not p450_panel[row[1]]==1 : p450_substrates[drugname]=row[1]+"*"             
     #add to list of substrates
    return p450_inhibitors,p450_inducers,p450_substrates,p450_multiplier
    
def make_table(any_dictionary,col1,col2,showfreq=1):
    htmltext='<table class="table table-condensed table-bordered table-striped"><thead><th style="width:20%">'+col1+'</th><th>'+col2+'</th></thead>'
    htmltext=htmltext+'<tbody>'
    if any_dictionary:
        for key1 in any_dictionary:
            htmltext=htmltext+'<tr><td>'+key1+'</td><td>'
            for key2 in any_dictionary[key1]:
                if showfreq==1:
                    htmltext=htmltext+key2+' ('+str(any_dictionary[key1][key2]*100)+'%), '
                else:
                    htmltext=htmltext+key2+', '
            htmltext=htmltext[:-2]+'</td></tr>'    
    else: htmltext=htmltext+'<tr><td>None</td><td> </td></tr>'
    htmltext=htmltext+'</tbody></table>'
    return htmltext

def make_table_list(any_dictionary,col1,col2):
    htmltext='<table class="table table-condensed table-bordered table-striped"><thead><th style="width:20%">'+col1+'</th><th>'+col2+'</th></thead>'
    htmltext=htmltext+'<tbody>'
    if any_dictionary:
        for key1 in any_dictionary:
            htmltext=htmltext+'<tr><td>'+key1+'</td><td>'+any_dictionary[key1]+'</td></tr>'
    else: htmltext=htmltext+'<tr><td>None</td><td> </td></tr>'
    htmltext=htmltext+'</tbody></table>'
    return htmltext
