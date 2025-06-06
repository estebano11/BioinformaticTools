#%%
# Script to get reactions from EC or KO
import pandas as  pd
import numpy as np
import os
import Bio
#%%
# import the table 
file = 'metabolism.parsed.tsv'
vmh_file = 'vmh_reactions.tsv'

vmh_df = pd.read_csv(vmh_file,sep='\t', index_col=0)
df = pd.read_csv(file,sep='\t')
df.head()

# %%
# I take only the aminoacids
aas = df[df['MAIN'] == 'aminoacids']
# I verify the KOs present in the predefined tables
kosInfile = aas.KO.unique()
kosVMH = vmh_df.keggorthology.unique()

presentes = set(kosInfile) & set(kosVMH)
diferencia = len(presentes)
porcentaje = np.round((diferencia/len(kosInfile)) * 100,2)

print (f'Hay una representatión de {diferencia} reacciones en un total \n de {len(kosInfile)} ({porcentaje}%)')

# %%
# Take out the rows that we have
for ko in kosInfile:
    idxs = list(aas[aas['KO'] == ko].index)
    if ko in presentes:
        rxns = list(vmh_df[vmh_df['keggorthology'] == ko]['seed'].values)
        for idx in idxs:
            aas.loc[idx,'RXN'] = rxns
    else:
        pass
#%%
Bio.KEGG.
