# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 10:36:29 2019

@author: fkurtog
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
Phi_A=[]
Phi_B=[]

with open('PHISTO_p_p_interaction.txt') as f:
    while True:
        A=f.readline()
        if A == '':
            break
        A=A.replace('\n','')
        B=A.split(' ')
        Phi_A.append(B[0])
        Phi_B.append(B[1])
#%%
G= nx.Graph()
for p in range(len(Phi_A)):
    G.add_edge(Phi_A[p],Phi_B[p])

#centrality = nx.eigenvector_centrality(G)
centrality = nx.eigenvector_centrality_numpy(G)
#nx.draw(G)

cntr_list=np.asarray(list(centrality.values()))
plt.figure(1)
plt.hist(cntr_list)
Phi_title = 'PHI Data - ' + str(len(Phi_B)) + ' interactions, & 2882 proteins'
plt.title(Phi_title)

thresh= 0.10

thresh_cntr_list = cntr_list[cntr_list>thresh]


#%%
Hhi_A=[]
Hhi_B=[]
with open('STRING_p_p_interaction_990_thr.txt') as f:
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
for p in range(len(Hhi_A)):
    G_HHI.add_edge(Hhi_A[p],Hhi_B[p])

#centrality = nx.eigenvector_centrality(G)
HHI_centrality = nx.eigenvector_centrality_numpy(G_HHI)
#%%
HHI_cntr_list=np.asarray(list(HHI_centrality.values()))
plt.figure(2)
plt.hist(HHI_cntr_list,alpha=1)

HHI_thresh= 0.000000000000001

HHI_thresh_cntr_list = HHI_cntr_list[HHI_cntr_list>HHI_thresh]
Hhi_title = 'HHI Data - ' + str(len(Hhi_B)) + ' interactions, & 3783 proteins'
plt.title(Hhi_title)
