# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 16:21:28 2020

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code loads in and sorts through the NMDID participant metadata to
    identify appropriate ID numbers for inclusion in this study.
    
    TODO: add inclusion criteria
    
    Potential inclusion:
        - No broken bones for leg/arm (depending on study)
        - History of certain sports (e.g. running, team sports [check file for options])
        - No chromosal abnormalities that might impact bone (e.g. cancer)
    
"""

# %% Import packages

import pandas as pd

# %% Load in metadata

#Load in the main decedents info
df_decedents = pd.read_csv('decedents.csv')

#Load in the supplementary data required
df_abnormalities = pd.read_csv('chromosomal-abnormalities.csv')
df_envHistory = pd.read_csv('env-history.csv')
df_medHistory = pd.read_csv('medical-history.csv')

# %% Extract activities data

#Loop through participant ID's in activities database and separately extract

#Get the unique list of id's
decedentIDs = df_decedents['id'].unique()

#Create a dictionary to sotre activity data in
activityDict = {'id': [], 'activities': []}

for ii in range(len(decedentIDs)):
    
    #Get the current participants activities set
    df_currID = df_envHistory.loc[(df_envHistory['id'] == decedentIDs[ii]) &
                                  (df_envHistory['environmental_factor_type'] == 'Activities'),
                                  ['id','environmental_factor_item']]
    df_currID.reset_index(inplace = True, drop = True)
    
    #Loop through activities and set in a list
    activities = list()
    for aa in range(len(df_currID)):
        activities.append(df_currID['environmental_factor_item'][aa])
        
    #Append id and activities to dictionary
    activityDict['id'].append(decedentIDs[ii])
    activityDict['activities'].append(activities)

#Convert the dictionary to dataframe
df_activities = pd.DataFrame.from_dict(activityDict)

#Merge the activities dataframe with the main decedents dataframe
df_decedents = pd.merge(left = df_decedents, right = df_activities,
                        how = 'left', left_on = 'id', right_on='id')
    
# %% Extract medical data

#Use medical data to extract surgeries and broken bones

#Create a dictionary to sotre medical data in
medDict = {'id': [], 'surgery': [], 'fracture': []}

for ii in range(len(decedentIDs)):
    
    #Get the current participants surgery and fracture history
    df_currID = df_medHistory.loc[(df_medHistory['id'] == decedentIDs[ii]) &
                                  (df_medHistory['medical_history_type'].isin(['Surgery','Broken Bone'])),
                                  ['id','medical_history_type','medical_history_detail']]
    df_currID.reset_index(inplace = True, drop = True)
    
    #Check if current id is empty before progressing
    if len(df_currID) != 0:
        
        #Extract any surgeries
        surgeries = list(df_currID.loc[(df_currID['medical_history_type'] == 'Surgery'),
                                       ['medical_history_detail']].reset_index(drop = True)['medical_history_detail'])
        
        ### TODO: remove if 'None' is the only surgery
        
        #Extract any broken bones
        fractures = list(df_currID.loc[(df_currID['medical_history_type'] == 'Broken Bone'),
                                       ['medical_history_detail']].reset_index(drop = True)['medical_history_detail'])
    
        ### TODO: remove if 'None' is the only broken bone
        
        ### TODO: append to the dictionary...
        
### TODO: left joint surgery and broken bones to decedents dataframe...


