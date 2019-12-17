# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 10:36:29 2019

@author: fkurtog
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from itertools import compress, islice
import urllib.parse
import urllib.request


Phi_A=[]
Phi_B=[]

with open('PHISTO_HerpesVirus_1_p_p_interaction.txt') as f:
    while True:
        A=f.readline()
        if A == '':
            break
        A=A.replace('\n','')
        B=A.split(' ')
        Phi_A.append(B[0])
        Phi_B.append(B[1])

G= nx.Graph()
for p in range(len(Phi_A)):
    G.add_edge(Phi_A[p],Phi_B[p])

#centrality = nx.eigenvector_centrality(G)
centrality = nx.eigenvector_centrality_numpy(G)
#nx.draw(G)

cntr_list=np.asarray(list(centrality.values()))*100
_, bins = np.histogram(np.log10(cntr_list + 1), bins='auto')
plt.figure(1)
plt.hist(cntr_list, bins=10**bins);
plt.yscale('log', nonposy='clip')
Phi_title = 'PHI Data - ' + str(G.number_of_edges()) + ' interactions, & ' + str(G.number_of_nodes()) + ' proteins'
plt.title(Phi_title)

thresh= 5

thresh_cntr_list = [cntr_list>thresh]
thresh_cntr_list = thresh_cntr_list[0]

Thresh_hold_prots = list(G.nodes())
Thresh_hold_prots=list(compress(Thresh_hold_prots, thresh_cntr_list))

Possible_Human_Targets=[]
Pathogen_Protens=[]
for a in Thresh_hold_prots:
    if a in Phi_A:
        Pathogen_Protens.append(a)
        print('Pathogen')
    elif a in Phi_B:
        Possible_Human_Targets.append(a)
        print('Human')
    else:
        print('None')


#%%
query=''
for k in Possible_Human_Targets:
    query = query + k +  " " 

query=query[:-1]
#%%
url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ACC',
'to': '	ENSEMBL_PRO_ID',
'format': 'tab',
'query': query
}
data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)

with urllib.request.urlopen(req) as f:
   response = f.read()
   
#%%
Enesmbl_Respnse= response.decode('utf-8')
Enesmbl_Respnse=Enesmbl_Respnse.split('\n')
Enesmbl_H_prot=[]
for line in Enesmbl_Respnse:
    if line == '':
        break
    Enesmbl_H_prot.append(line.split('\t')[1])
    
Enesmbl_H_prot=Enesmbl_H_prot[1:]




#%%
Hhi_A=[]
Hhi_B=[]
with open('STRING_p_p_interaction_999_thr.txt') as f:
    while True:
        A=f.readline()
        if A == '':
            break
        A=A.replace('\n','')
        B=A.split(' ')
        Hhi_A.append(B[0])
        Hhi_B.append(B[1])
#%%
G_HHI= nx.Graph()
for j in range(len(Hhi_A)):
    G_HHI.add_edge(Hhi_A[j],Hhi_B[j])

#centrality = nx.eigenvector_centrality(G)
HHI_centrality = nx.eigenvector_centrality_numpy(G_HHI)
#%%
HHI_cntr_list=np.asarray(list(HHI_centrality.values()))*100
_, bins = np.histogram(np.log10(HHI_cntr_list + 1),bins=50)
plt.figure(2)
plt.hist(HHI_cntr_list, bins=10**bins);
plt.yscale('log', nonposy='clip')



Hhi_title = 'HHI Data - ' + str(G_HHI.number_of_edges()) + ' interactions, & ' + str(G_HHI.number_of_nodes()) + ' proteins'
plt.title(Hhi_title)


#%%
HHI_thresh= 0.01

HHI_thresh_cntr_list = [HHI_cntr_list<HHI_thresh]
HHI_thresh_cntr_list = HHI_thresh_cntr_list[0]

HHI_Thresh_hold_prots = list(G_HHI.nodes())
HHI_Thresh_hold_prots=list(compress(HHI_Thresh_hold_prots, HHI_thresh_cntr_list))

Mutual_Proteins =[]

for l in Enesmbl_H_prot:
    if l in HHI_Thresh_hold_prots:
        print('Found one')
        Mutual_Proteins.append(l)
