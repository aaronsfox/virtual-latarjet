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
    
    The inclusion criteria used are:
        - Normal BMI
        - No history of arm/leg (or ankle/foot) fracture (depending on project)
        - No history of arm/leg surgery (depending on project)
        - Physically active, based on activities of:
            > american flag or touch football
            > american tackle football
            > baseball
            > basketball
            > cross country skiing
            > dancing
            > exercise machines primarily for cardiorespiratory conditioning
            > exercise machines primarily for muscle strengthening
            > free weights
            > golf
            > ice hockey
            > martial arts
            > mountain climbing, rock climbing and wall climbing
            > other involving cardiorespiratory exercise
            > other involving dancing and rhythmic movements
            > other involving muscle strengthening exercises
            > other involving other sports and athletics played as a team or group
            > other involving other sports and athletics played individually
            > other specified sports and athletics
            > running
            > soccer
            > swimming
            > volleyball (beach) (court)
            > walking, marching and hiking
        - Specific chromosomal abnormality affecting bone:
            > Bone Cancer\xa0(includes Ewing Sarcoma and Osteosarcoma and Malignant Fibrous Histiocytoma)
            > Congenital deformities of feet (for certain studies)
            > Congenital malformations of musculoskeletal system, not elsewhere classified
            > Other congenital musculoskeletal deformities
            > Spina bifida
    
"""

# %% Import packages

import pandas as pd
import numpy as np
from itertools import chain
import random

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
surgeryDict = {'id': [], 'surgery': []}
fractureDict = {'id': [], 'fracture': []}

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
        
        #Check if the list just contains 'None'
        if surgeries and all(elem == 'None' for elem in surgeries):
            #Set the list back to empty
            surgeries = None
        
        #Extract any broken bones
        fractures = list(df_currID.loc[(df_currID['medical_history_type'] == 'Broken Bone'),
                                       ['medical_history_detail']].reset_index(drop = True)['medical_history_detail'])
    
        #Check if the list just contains 'None'
        if fractures and all(elem == 'None' for elem in fractures):
            #Set the list back to empty
            fractures = None
            
        #Append id and medical history to dictionaries
        if surgeries is not None:
            surgeryDict['id'].append(decedentIDs[ii])
            surgeryDict['surgery'].append(surgeries)
        if fractures is not None:
            fractureDict['id'].append(decedentIDs[ii])
            fractureDict['fracture'].append(fractures)

#Convert the dictionary to dataframe
df_surgeries = pd.DataFrame.from_dict(surgeryDict)
df_fractures = pd.DataFrame.from_dict(fractureDict)

#Merge the medical history dataframes with the main decedents dataframe
df_decedents = pd.merge(left = df_decedents, right = df_surgeries,
                        how = 'left', left_on = 'id', right_on='id')
df_decedents = pd.merge(left = df_decedents, right = df_fractures,
                        how = 'left', left_on = 'id', right_on='id')

# %% Extract chromosomal abnormalities

#Create a dictionary to sotre activity data in
abnormalitiesDict = {'id': [], 'abnormalities': []}

for ii in range(len(decedentIDs)):
    
    #Get the current participants activities set
    df_currID = df_abnormalities.loc[(df_abnormalities['id'] == decedentIDs[ii]),
                                     ['id','abnormality_factor']]
    df_currID.dropna(inplace = True)
    df_currID.reset_index(inplace = True, drop = True)
    
    #Remove nones from list too
    abnormalities = list(df_currID['abnormality_factor'])
    abnormalities = [x for x in abnormalities if x != 'None']    
    abnormalities = [x for x in abnormalities if x != 'Unknown']
    
    #First check if dataframe is empty
    if len(abnormalities) != 0:
            
        #Append id and activities to dictionary
        abnormalitiesDict['id'].append(decedentIDs[ii])
        abnormalitiesDict['abnormalities'].append(abnormalities)

#Convert the dictionary to dataframe
df_abnormalities = pd.DataFrame.from_dict(abnormalitiesDict)

#Merge the activities dataframe with the main decedents dataframe
df_decedents = pd.merge(left = df_decedents, right = df_abnormalities,
                        how = 'left', left_on = 'id', right_on='id')

# %% Calculate BMI

#Calculate BMI using living weight and height data
bmi = np.array(df_decedents['living_weight']) / \
    ((np.array(df_decedents['living_height']) / 100)**2)
    
#Append to decedents dataframe
df_decedents['bmi'] = bmi

#Identify BMI category
bmiCat = list()
for bb in range(len(bmi)):
    #Check for nan value
    if np.isnan(bmi[bb]):
        bmiCat.append(np.nan)
    else:
        #Calculate according to category
        if bmi[bb] < 18.5:
            bmiCat.append('underweight')
        elif bmi[bb] >= 18.5 and bmi[bb] < 25:
            bmiCat.append('normal')
        elif bmi[bb] >= 25 and bmi[bb] < 30:
            bmiCat.append('overweight')
        elif bmi[bb] >= 30:
            bmiCat.append('obese')
            
#Append to decedents dataframe
df_decedents['bmi_category'] = bmiCat

# %% Extract participants

#Set the inclusion lists for the different categories

#Fracture criteria
#Unpack all unique values in the list
fractureList = list(set(chain(*list(df_fractures['fracture']))))
#Remove the unwanted fractures in the list (separate arm and leg)
fractureArmInc = [x for x in fractureList if x != 'Arm']
fractureLegInc = [x for x in fractureList if x != 'Leg']
fractureLegInc = [x for x in fractureLegInc if x != 'Ankle/foot']
#Create exclusion fracture list
fractureArmExc = ['Arm']
fractureLegExc = ['Leg','Ankle/foot']

#Surgery criteria
#Unpack all unique values in the list
surgeryList = list(set(chain(*list(df_surgeries['surgery']))))
#Remove the unwanted surgeries in the list (separate arm and leg)
surgeryArmInc = [x for x in surgeryList if x != 'Arms']
surgeryLegInc = [x for x in surgeryList if x != 'Legs']
#Create exclusion for surgery list
surgeryArmExc = ['Arms']
surgeryLegExc = ['Legs']

#Activities criteria
#Unpack all unique values in the list
activityInc = list(set(chain(*list(df_activities['activities']))))
activityExc = list(set(chain(*list(df_activities['activities']))))
desiredActivities = ['american flag or touch football','american tackle football','baseball',
                     'basketball','cross country skiing','dancing','exercise machines primarily for cardiorespiratory conditioning',
                     'exercise machines primarily for muscle strengthening','free weights',
                     'golf','ice hockey','martial arts','mountain climbing, rock climbing and wall climbing',
                     'other involving cardiorespiratory exercise','other involving dancing and other rhythmic movements',
                     'other involving muscle strengthening exercises','other involving other sports and athletics played as a team or group',
                     'other involving other sports and athletics played individually',
                     'other specified sports and athletics','running','soccer','swimming',
                     'volleyball (beach) (court)','walking, marching and hiking']
#Get a list of unwanted activities
for aa in range(len(desiredActivities)):
    activityExc = [x for x in activityExc if x != desiredActivities[aa]]
#Remove unwanted activities from the list
for aa in range(len(activityExc)):
    activityInc = [x for x in activityInc if x != activityExc[aa]]

#Abnormalities criteria
#Unpack all unique values in the list
abnormalityInc = list(set(chain(*list(df_abnormalities['abnormalities']))))
abnormalityExc = list(set(chain(*list(df_abnormalities['abnormalities']))))
undesiredAbnormalities = ['Bone Cancer\xa0(includes Ewing Sarcoma and Osteosarcoma and Malignant Fibrous Histiocytoma)',
                          'Congenital deformities of feet (for certain studies)',
                          'Congenital malformations of musculoskeletal system, not elsewhere classified',
                          'Other congenital musculoskeletal deformities',
                          'Spina bifida']
#Get a list of OK abnormalities
for aa in range(len(undesiredAbnormalities)):
    abnormalityInc = [x for x in abnormalityInc if x != undesiredAbnormalities[aa]]
#Get a list of excluded abnormalities
for aa in range(len(abnormalityInc)):
    abnormalityExc = [x for x in abnormalityExc if x != abnormalityInc[aa]]
    
# %% Upper limb participants    

#Extract the participants who meet the desired criteria
#With the variable methods used here it works best going step-by-step

#Extract normal BMI participants  
df_incParticipantsArm = df_decedents.loc[(df_decedents['bmi_category'] == 'normal'),]    

#Convert empty list to nan for later application
df_incParticipantsArm.mask(df_incParticipantsArm.applymap(str).eq('[]'),
                           inplace = True)

#With desired activities
#Unpack the individual activities into their own dataframe
df_splitActivities = pd.DataFrame(df_incParticipantsArm['activities'].values.tolist(),
                                  index = df_incParticipantsArm.index)
#Set up a boolena mask for the participants who have the criteria
maskActivities = df_splitActivities.astype(str).isin(activityInc).any(axis=1)
#Extract the participants using the mask
df_incParticipantsArm = df_incParticipantsArm[maskActivities]

#Remove undesired fractures
#Convert nans to 'None'
for row in df_incParticipantsArm.loc[df_incParticipantsArm['fracture'].isna(),'fracture'].index:
    df_incParticipantsArm.at[row, 'fracture'] = ['None']
#Unpack the individual fractures into their own dataframe
df_splitFractures = pd.DataFrame(df_incParticipantsArm['fracture'].values.tolist(),
                                 index = df_incParticipantsArm.index)
#Set up a boolean mask for the participants who have the criteria
#Inverting the boolean mask (i.e. ~) is required given using exclusion criteria
maskFractures = ~df_splitFractures.astype(str).isin(fractureArmExc).any(axis=1)
#Extract the participants using the mask
df_incParticipantsArm = df_incParticipantsArm[maskFractures]

#Remove undesired surgeries
#Convert nans to 'None'
for row in df_incParticipantsArm.loc[df_incParticipantsArm['surgery'].isna(),'surgery'].index:
    df_incParticipantsArm.at[row, 'surgery'] = ['None']
#Unpack the individual surgeries into their own dataframe
df_splitSurgeries = pd.DataFrame(df_incParticipantsArm['surgery'].values.tolist(),
                                 index = df_incParticipantsArm.index)
#Set up a boolean mask for the participants who have the criteria
#Inverting the boolean mask (i.e. ~) is required given using exclusion criteria
maskSurgeries = ~df_splitSurgeries.astype(str).isin(surgeryArmExc).any(axis=1)
#Extract the participants using the mask
df_incParticipantsArm = df_incParticipantsArm[maskSurgeries]

#Remove undesired abnormalities
#Convert nans to 'None'
for row in df_incParticipantsArm.loc[df_incParticipantsArm['abnormalities'].isna(),'abnormalities'].index:
    df_incParticipantsArm.at[row, 'abnormalities'] = ['None']
#Unpack the individual activities into their own dataframe
df_splitAbnormalities = pd.DataFrame(df_incParticipantsArm['abnormalities'].values.tolist(),
                                     index = df_incParticipantsArm.index)
#Set up a boolean mask for the participants who have the criteria
#Inverting the boolean mask (i.e. ~) is required given using exclusion criteria
maskAbnormalities = ~df_splitAbnormalities.astype(str).isin(abnormalityExc).any(axis=1)
#Extract the participants using the mask
df_incParticipantsArm = df_incParticipantsArm[maskAbnormalities]

# %% Lower limb participants    

#Extract the participants who meet the desired criteria
#With the variable methods used here it works best going step-by-step

#Extract normal BMI participants  
df_incParticipantsLeg = df_decedents.loc[(df_decedents['bmi_category'] == 'normal'),]    

#Convert empty list to nan for later application
df_incParticipantsLeg.mask(df_incParticipantsLeg.applymap(str).eq('[]'),
                           inplace = True)

#With desired activities
#Unpack the individual activities into their own dataframe
df_splitActivities = pd.DataFrame(df_incParticipantsLeg['activities'].values.tolist(),
                                  index = df_incParticipantsLeg.index)
#Set up a boolena mask for the participants who have the criteria
maskActivities = df_splitActivities.astype(str).isin(activityInc).any(axis=1)
#Extract the participants using the mask
df_incParticipantsLeg = df_incParticipantsLeg[maskActivities]

#Remove undesired fractures
#Convert nans to 'None'
for row in df_incParticipantsLeg.loc[df_incParticipantsLeg['fracture'].isna(),'fracture'].index:
    df_incParticipantsLeg.at[row, 'fracture'] = ['None']
#Unpack the individual fractures into their own dataframe
df_splitFractures = pd.DataFrame(df_incParticipantsLeg['fracture'].values.tolist(),
                                 index = df_incParticipantsLeg.index)
#Set up a boolean mask for the participants who have the criteria
#Inverting the boolean mask (i.e. ~) is required given using exclusion criteria
maskFractures = ~df_splitFractures.astype(str).isin(fractureLegExc).any(axis=1)
#Extract the participants using the mask
df_incParticipantsLeg = df_incParticipantsLeg[maskFractures]

#Remove undesired surgeries
#Convert nans to 'None'
for row in df_incParticipantsLeg.loc[df_incParticipantsLeg['surgery'].isna(),'surgery'].index:
    df_incParticipantsLeg.at[row, 'surgery'] = ['None']
#Unpack the individual surgeries into their own dataframe
df_splitSurgeries = pd.DataFrame(df_incParticipantsLeg['surgery'].values.tolist(),
                                 index = df_incParticipantsLeg.index)
#Set up a boolean mask for the participants who have the criteria
#Inverting the boolean mask (i.e. ~) is required given using exclusion criteria
maskSurgeries = ~df_splitSurgeries.astype(str).isin(surgeryLegExc).any(axis=1)
#Extract the participants using the mask
df_incParticipantsLeg = df_incParticipantsLeg[maskSurgeries]

#Remove undesired abnormalities
#Convert nans to 'None'
for row in df_incParticipantsLeg.loc[df_incParticipantsLeg['abnormalities'].isna(),'abnormalities'].index:
    df_incParticipantsLeg.at[row, 'abnormalities'] = ['None']
#Unpack the individual activities into their own dataframe
df_splitAbnormalities = pd.DataFrame(df_incParticipantsLeg['abnormalities'].values.tolist(),
                                     index = df_incParticipantsLeg.index)
#Set up a boolean mask for the participants who have the criteria
#Inverting the boolean mask (i.e. ~) is required given using exclusion criteria
maskAbnormalities = ~df_splitAbnormalities.astype(str).isin(abnormalityExc).any(axis=1)
#Extract the participants using the mask
df_incParticipantsLeg = df_incParticipantsLeg[maskAbnormalities]

# %% Export participants

#Upper limb
df_incParticipantsArm.to_csv('extractedParticipants_UpperLimb.csv', index = False)


#Lower limb

#All
df_incParticipantsLeg.to_csv('extractedParticipants_LowerLimb.csv',index = False)

#Extract a set of 30 male and female participants
#Note there is only 11 female, so we will take 11 and 19
#Get the females
df_legFemale = df_incParticipantsLeg.loc[df_incParticipantsLeg['sex_code'] == 'Female',]
#Get male id's
legMale_id = list(df_incParticipantsLeg.loc[df_incParticipantsLeg['sex_code'] == 'Male',]['id'])
#Randomly select 19 male participants from ID list
#Set seed and select 19 participants
random.seed(1234)
legMale_selected = random.sample(legMale_id, 19)
#Extract male dataframe
df_legMale = df_incParticipantsLeg.loc[(df_incParticipantsLeg['sex_code'] == 'Male') &
                                       (df_incParticipantsLeg['id'].isin(legMale_selected)),]
#Concatenate male and female databases
df_selectParticipantsLeg = pd.concat([df_legFemale,df_legMale])
#Print to file
df_selectParticipantsLeg.to_csv('selectedParticipants_LowerLimb.csv',index = False)


#Extract the unique ID numbers across the two dataframes for NMDID searching
uniqueParticipants = pd.concat([df_incParticipantsArm,df_incParticipantsLeg])['deidentified_record_number'].unique()

#Export participant identifier numbers
np.savetxt('extractedParticipants_recordNumbers.txt',
           uniqueParticipants, delimiter=',')

# %% ----- End of participantFinder.py ----- %% #