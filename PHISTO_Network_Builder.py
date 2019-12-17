# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 10:54:09 2019

@author: fkurtog & jkrzesni
"""


file_name='PHISTO_HerpesVirus_1.csv'
line_counter = 0
link_counter =0
Links_1=[]
protein_A=[]
protein_B=[]
with open(file_name) as f:
    A=f.readline()
    while True:
        A=f.readline()
        A=A.replace('\n','')
        B=A.split(',')
        line_counter += 1
        if (A  == ''):
            break
        if (B[2][1:-1]+'_'+B[4][1:-1] not in Links_1) and (B[4][1:-1]+'_'+B[2][1:-1] not in Links_1):
            Links_1.append(B[2][1:-1]+'_'+B[4][1:-1])
            link_counter += 1
            protein_A.append(B[2][1:-1])
            protein_B.append(B[4][1:-1])
        else:
            print('Redundant Interaction found')
            
with open("PHISTO_HerpesVirus_1_p_p_interaction.txt", 'w') as f:
    for i in range(len(protein_A)):
        f.write("%s %s\n" % (protein_A[i], protein_B[i]))
    f.close()
