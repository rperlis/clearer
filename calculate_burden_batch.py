from StringIO import StringIO
import csv
import itertools
import gzip  #needed for gzipped files

import re   # needed for splitting strings
import os   # needed for file names

# Refactored 11/12 to allow batch input of files
# TO-DO: 
    # confirm numbers corresponding to frequency
    # fix med-naming (fluoxetine <> fluoxetine hydrochloride) --> USE fuzzy wuzzy!
    
# CONSTANTS would go here, if any 
FNAME_MED_TO_CID=os.path.dirname(__file__)+"/static/label_mapping.tsv.gz" 
FNAME_CID_TO_AE=os.path.dirname(__file__)+"/static/meddra_freq_parsed.tsv.gz"
CYP450_MODIFIERS="cyp450_mods.txt"
CYP450_SUBSTRATES="cyp450_substrates.txt"

# note patient data-in is really csv not tsv!!!

PATIENT_FILE=os.path.dirname(__file__)+"/static/patient_data.tsv"
PATIENT_FILE_OUTPUT=os.path.dirname(__file__)+"/static/patient_data_out"

def check_medlist_batch(variables): 
    """

    Check a list of meds and sum up burden

    Variables is a dictionary passed from server.py
 
    """
    
# 1. read in matching dictionary (drug->CID) - we keep all drugs in dictionary, then batch
# 2. read in AE profile lists (ie, panels of AE's we care about)
# 3. read in CID AE list- do not need to keep in memory
# 4. map all pt meds to CID's - let complist be a dictionary
# 5. parse all patients
# 6. spit out file

# 1. read in drug->CID dict

# map to CID for lookup

    matcher_dict={}
    backmatch_dict={}
    matcheddruglist=[]
    
    # format: row0=brand, row1=generic, row3=CID
    # note we currently only match brand!
    with gzip.open(FNAME_MED_TO_CID) as gzfile:
        medpairs=csv.reader(gzfile,delimiter='\t')    
        for row in medpairs:

# first match generic name
            gname=row[1].upper().split(" ")[0]
            if (not gname in matcher_dict) and (not gname in matcheddruglist) :
                matcher_dict[gname]= row[3]
                if not backmatch_dict.has_key(row[3]):    #necessary to prevent duplicate backmatches
                    backmatch_dict[row[3]]=gname
                matcheddruglist.append(gname) #this is here to prevent duplicate entries
# then match brand name                
            bname=row[0].upper().split(" ")[0]
            if (not bname in matcher_dict) and (not bname in matcheddruglist) :
                matcher_dict[bname]= row[3]
                # backmatch_dict[row[3]]=row[1].upper() ----> backmatch only for generics
                matcheddruglist.append(bname) #this is here to prevent duplicate entries

    print("loaded drug match list")
# 2. read in panels of AE's

# aelist will be a dictionary, rather than a plain list as in single-patient version
    #key[drug]->list of AE concept ID's
    aedict={}
    aedict['Psychiatric']= load_aefilelist("CNS_psychiatric.txt")  
    aedict['Cognition']=load_aefilelist("CNS_cognition.txt")
    aedict['Other_neurologic']=load_aefilelist("CNS_other_neurologic.txt")
    aedict['Fall_risk']=load_aefilelist("CNS_fall_risk.txt")     

    print("loaded AE panels")
# 3. read in CID AE list
# read in AE megafile - it's gzipped...
    property_dict={}
    parseddict={}   # this prevents duplicate AEs per drug
    with gzip.open(FNAME_CID_TO_AE) as gzfile:
        drug_properties=csv.reader(gzfile,delimiter='\t')
        for row in drug_properties:
            drugname=row[0].upper().split(" ")[0]
            aename=row[11]
            if not(parseddict.has_key(drugname)): parseddict[drugname]={}
            if not(parseddict[drugname].has_key(aename)):
                if (row[9]=="PT" and not(row[5]=="placebo")) and not drugname=='' :
                    parseddict[drugname][aename]=1
                # need to reassign frequencies (exact) to categories
                   #postmarketing, rare, infrequent, frequent, or exact #
                    freqnum=parse_frequency(row[6])
                    if variables['Option_1']==0 and freqnum>0: freqnum=.1   # treats any ae as freq=.1
                
                # next section modified from single-patient version:
                # loop through aelist and see if the ae here is in our list. 
                # if it is, increment the drug-ae-burden-item for that aelist
                
                    for key in aedict:                  # loop through AE panels
                        if aename in aedict[key]:       #if ae is in a panel,
                            if property_dict.has_key(drugname):
                                 if property_dict[drugname].has_key(key):
                                     property_dict[drugname][key]=property_dict[drugname][key]+freqnum
                                 else : property_dict[drugname][key]=freqnum
                            else : 
                                 property_dict[drugname]={}
                                 property_dict[drugname][key]=freqnum
    print("loaded AE list")                
#4. now read in list of patients and meds. format: pt_id drug_list (tab-separated)
# read by making everything into a list, pop first item from list and call it pt_id, assign remaining list to it.

# files from Victor - column 1 (0-index) is name, col 0 is #, col 3 is high-level term
    pt_drugdict={}   

    pt_drugfile=csv.reader(open(PATIENT_FILE),delimiter=',')
    for row in pt_drugfile:
        rowlist=row
        pt_id=rowlist.pop(0)
        rowlist=[element.upper() for element in rowlist]
        pt_drugdict[pt_id]=rowlist

# now loop through patients - for each patient, make p450 list, use it to parse AE's
    pt_output={}
    pt_output_nop450={}
    for ptname in pt_drugdict:
        p450_inhibitors,p450_inducers,p450_substrates,multiplier=map_p450(pt_drugdict[ptname],1)
        # match meds to cid's for AE lookup
        matchedcid=[]
        for med in pt_drugdict[ptname]:
            if matcher_dict.has_key(med.upper()):
                matchedcid.append(matcher_dict[med.upper()])  
                
        # loop by drug, then by ae category - make one dict with p450, another without
        for cid in matchedcid:
            for aecat in aedict:
                if property_dict.has_key(cid):
                    if property_dict[cid].has_key(aecat):
                        if not pt_output.has_key(ptname): 
                            pt_output[ptname]={}
                            pt_output_nop450[ptname]={}
                        if pt_output[ptname].has_key(aecat):
                            pt_output[ptname][aecat]=pt_output[ptname][aecat]+(property_dict[cid][aecat]*multiplier[backmatch_dict[cid]])
                            pt_output_nop450[ptname][aecat]=pt_output[ptname][aecat]+(property_dict[cid][aecat])
                        else : 
                            pt_output[ptname][aecat]=property_dict[cid][aecat]*multiplier[backmatch_dict[cid]]
                            pt_output_nop450[ptname][aecat]=property_dict[cid][aecat]
    print("matched patients")     
    
# now return results  
    filename=PATIENT_FILE_OUTPUT
    if variables['Option_1']==0: filename=filename+"_nofreq"
    plink_file=file(filename+".tsv","w")
    plink_file.write("Pt_id")
    for aecats in aedict:
        plink_file.write("\t"+aecats)
    for pt_name in pt_output:
        plink_file.write("\n"+pt_name)
        for aecats in aedict:
            if pt_output[pt_name].has_key(aecats) : plink_file.write("\t"+str(pt_output[pt_name][aecats]))
            else : plink_file.write("\t"+'0')
    plink_file.close()

# and again with nop450    
    plink_file=file(filename+"_nop450"+".tsv","w")
    plink_file.write("Pt_id")
    for aecats in aedict:
        plink_file.write("\t"+aecats)
    for pt_name in pt_output_nop450:
        plink_file.write("\n"+pt_name)
        for aecats in aedict:
            if pt_output_nop450[pt_name].has_key(aecats) : plink_file.write("\t"+str(pt_output_nop450[pt_name][aecats]))
            else : plink_file.write("\t"+'0')
    plink_file.close()
    
    return {
        'pts_written': len(pt_output),
        'ae_categories': len(aedict)
    }

"""
this function categorizes AE's and writes a new gz file
"""

def load_aefilelist (aelistfile):

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

def map_p450(list_of_meds,use_p450):
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
       if row[0] in list_of_meds:
           p450_panel[row[1]]=p450_panel[row[1]]*float(row[2])                     #modify p450 status panel
           if float(row[2])<1: p450_inducers[row[0]]=row[1]                                    #add to list of modifiers
           if float(row[2])>1: p450_inhibitors[row[0]]=row[1]
    # now generate multipliers for med side effects
    cyp450_subs=csv.reader(open(os.path.dirname(__file__)+"/static/"+CYP450_SUBSTRATES),delimiter='\t')
    for row in cyp450_subs:
       drugname=row[0].upper().split(" ")[0] 
       if drugname in list_of_meds:
           if use_p450==1: p450_multiplier[drugname]=p450_multiplier[drugname]*p450_panel[row[1]]         #lookup p450 key, multiply  
           p450_substrates[drugname]=row[1]                                   #add to list of substrates
    return p450_inhibitors,p450_inducers,p450_substrates,p450_multiplier
    
def make_table(any_dictionary,col1,col2,showfreq=1):
    htmltext='<table class="table table-condensed table-bordered table-striped"><thead><th>'+col1+'</th><th>'+col2+'</th></thead>'
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
    htmltext='<table class="table table-condensed table-bordered table-striped"><thead><th>'+col1+'</th><th>'+col2+'</th></thead>'
    htmltext=htmltext+'<tbody>'
    if any_dictionary:
        for key1 in any_dictionary:
            htmltext=htmltext+'<tr><td>'+key1+'</td><td>'+any_dictionary[key1]+'</td></tr>'
    else: htmltext=htmltext+'<tr><td>None</td><td> </td></tr>'
    htmltext=htmltext+'</tbody></table>'
    return htmltext
