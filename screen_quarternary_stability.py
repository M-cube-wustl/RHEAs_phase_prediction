#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import math
from pymatgen import Composition
from pymatgen import MPRester
import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, Element
from matminer.data_retrieval.retrieve_MP import MPDataRetrieval
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
mpdr = MPDataRetrieval(api_key='sTzfI3LdPo2O5J8h') # or MPDataRetrieval(api_key=YOUR_API_KEY here)


# In[ ]:


#import matrixs
df_sqs_2=pd.read_excel('./enthalpy_data_and_predictions/pairwise_mixing_enthalpy.xlsx', sheet_name='our work')
df_sqs_2.set_index('Unnamed: 0', inplace=True)


# In[ ]:


import itertools

def binary():
    first_el = el
    second_el = el
    return sorted(["{}-{}".format(*sorted(pair))
                   for pair in itertools.product(first_el, second_el)])

def ternary():
    """
    Return a sorted list of chemical systems
        of the form [...,"Li-Be-S",...,"Al-Zr-Se",...]
    """
    first_el = el
    second_el = el
    third_el = el
    return sorted(["{}-{}-{}".format(*sorted(triple))
                   for triple in itertools.product(first_el, second_el, third_el)])


# In[ ]:


df_3=pd.read_excel('./enthalpy_data_and_predictions/mixing_enthalpy_predictions.xlsx', sheet_name ='ternary')
df_4=pd.read_excel('./enthalpy_data_and_predictions/mixing_enthalpy_predictions.xlsx', sheet_name ='quarternary')


# In[ ]:


for index, row in df_4.iterrows():
    e1=row['e1']
    e2=row['e2']
    e3=row['e3']
    e4=row['e4']
    df_pd = pd.DataFrame()
    #quarternary
    df_pd.at[0,'comp']= row['comp']
    df_pd.at[0,'Hf']= row['enthalpy(prediction)']
    df_pd.at[0,'S']= row['entropy']
    #ternary
    df_a=df_3[df_3['comp'].str.match(e1+e2+e3)]
    df_a.reset_index(drop=True, inplace=True)
    df_b=df_3[df_3['comp'].str.match(e1+e2+e4)]
    df_b.reset_index(drop=True, inplace=True)
    df_c=df_3[df_3['comp'].str.match(e1+e3+e4)]
    df_c.reset_index(drop=True, inplace=True)
    df_d=df_3[df_3['comp'].str.match(e2+e3+e4)]
    df_d.reset_index(drop=True, inplace=True)
    df_pd.at[1,'comp']= df_a.loc[0,'comp']
    df_pd.at[1,'Hf']= df_a.loc[0,'enthalpy(prediction)']
    df_pd.at[1,'S']= df_a.loc[0,'entropy']
    df_pd.at[2,'comp']= df_b.loc[0,'comp']
    df_pd.at[2,'Hf']= df_b.loc[0,'enthalpy(prediction)']
    df_pd.at[2,'S']= df_b.loc[0,'entropy']
    df_pd.at[3,'comp']= df_c.loc[0,'comp']
    df_pd.at[3,'Hf']= df_c.loc[0,'enthalpy(prediction)']
    df_pd.at[3,'S']= df_c.loc[0,'entropy']
    df_pd.at[4,'comp']= df_d.loc[0,'comp']
    df_pd.at[4,'Hf']= df_d.loc[0,'enthalpy(prediction)']
    df_pd.at[4,'S']= df_d.loc[0,'entropy']
    #binary
    df_pd.at[5,'comp']= e1+e2
    df_pd.at[5,'Hf']= df_sqs_2[e1][e2]
    df_pd.at[5,'S']=5.9730802545007364e-05
    df_pd.at[6,'comp']= e1+e3
    df_pd.at[6,'Hf']= df_sqs_2[e1][e3]
    df_pd.at[6,'S']=5.9730802545007364e-05
    df_pd.at[7,'comp']= e1+e4
    df_pd.at[7,'Hf']= df_sqs_2[e1][e4]
    df_pd.at[7,'S']=5.9730802545007364e-05
    df_pd.at[8,'comp']= e2+e3
    df_pd.at[8,'Hf']= df_sqs_2[e2][e3]
    df_pd.at[8,'S']=5.9730802545007364e-05
    df_pd.at[9,'comp']= e2+e4
    df_pd.at[9,'Hf']= df_sqs_2[e2][e4]
    df_pd.at[9,'S']=5.9730802545007364e-05
    df_pd.at[10,'comp']= e3+e4
    df_pd.at[10,'Hf']= df_sqs_2[e3][e4]
    df_pd.at[10,'S']=5.9730802545007364e-05
    #elements
    df_pd.at[11,'comp']= e1
    df_pd.at[11,'Hf']= 0
    df_pd.at[11,'S']= 0
    df_pd.at[12,'comp']= e2
    df_pd.at[12,'Hf']= 0
    df_pd.at[12,'S']= 0
    df_pd.at[13,'comp']= e3
    df_pd.at[13,'Hf']= 0
    df_pd.at[13,'S']= 0
    df_pd.at[14,'comp']= e4
    df_pd.at[14,'Hf']= 0
    df_pd.at[14,'S']= 0
    #query binary intermetallic from the materials project
    el={e1,e2,e3,e4}
    df_1 = mpdr.get_dataframe({'chemsys': {'$in': binary()}}, 
                            properties=['pretty_formula', "formation_energy_per_atom", "e_above_hull"])
    df_2 = mpdr.get_dataframe({'chemsys': {'$in': ternary()}}, 
                            properties=['pretty_formula', "formation_energy_per_atom", "e_above_hull"])
    df_1 = df_1.append(df_2, ignore_index=True)
    df_in=pd.DataFrame()
    df_in['comp']=df_1['pretty_formula']
    df_in['Hf']=df_1['formation_energy_per_atom']
    df_in['S']=0
    df_final=df_pd.append(df_in, ignore_index=True)
    #build phase diagram
    comps = df_final['comp']
    Ef=df_final['Hf']-2000*df_final['S']
    mg_comp=[None]*len(comps)
    for i in range(len(comps)):
        mg_comp[i]=Composition(comps[i])
        entries3=[None]*len(mg_comp)
    for i in range(len(mg_comp)):
        entries3[i]=PDEntry(composition=mg_comp[i], energy=Ef[i])
    phase = PhaseDiagram(entries3)
    #get decomposition and e_aobve_hull
    test = PDEntry(composition=mg_comp[0], energy=Ef[0])
    print(phase.get_decomp_and_e_above_hull(test))
    #df_4.loc[index,'phase']=phase.get_decomp_and_e_above_hull(test)


# In[7]:


df_4.to_excel("4-element.xlsx",sheet_name='2000K') 

