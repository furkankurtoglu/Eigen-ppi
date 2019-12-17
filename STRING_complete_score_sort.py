#import sys 
from itertools import islice
import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'


#file_name = sys.argv[1]
#c_s_thresh =  int(sys.argv[2]) #combined interaction score threshold

file_name='STRING_9606.protein.links.full.v11.0.txt'
c_s_thresh = 990

n = 1 #slice size
protein_dict = {}
ppi_counter = 0
ppi_counter_thresh = 0
with open(file_name) as f:
    while True:
        next_n_lines = list(islice(f, n))
        ppi_counter += 1
        if not next_n_lines:
            break
        for line in next_n_lines:
            line_split = line.split()
            if line_split[0] == "protein1":
                break
            elif (int(line_split[len(line_split)-1]) >= c_s_thresh):
                protein_dict[line_split[0][5:]] = line_split[1][5:]

f.close()


#%%
query1=''
query2=''
keys = protein_dict.keys()
for geneA in keys:
    geneB = protein_dict[geneA]
    query1 = query1 + geneA + " "
    query2 = query2 + geneB + " "
#%%

params = {
'from': 'ENSEMBL_PRO_ID',
'to': '	ACC',
'format': 'tab',
'query': query1
}
data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)

with urllib.request.urlopen(req) as f:
   response = f.read()
Protein_A= response.decode('utf-8')


#%%
params = {
'from': 'ENSEMBL_PRO_ID',
'to': '	ACC',
'format': 'tab',
'query': query2
}
data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)

with urllib.request.urlopen(req) as f:
   response = f.read()
Protein_B= response.decode('utf-8')

#%%
Prot_A=[]
Prot_A_lines = Protein_A.split('\n')
for line in Prot_A_lines:
    if line == '':
        break
    Prot_A.append(line.split('\t')[1])

Prot_B=[]
Prot_B_lines = Protein_B.split('\n')
for line in Prot_B_lines:
    if line == '':
        break
    Prot_B.append(line.split('\t')[1])

Prot_A.remove('To')
Prot_B.remove('To')

#%%

f = open("STRING_p_p_interaction_"+str(c_s_thresh)+"_thr.txt", 'w')
for protein1 in protein_dict.keys():
    f.write("%s %s\n" % (protein1, protein_dict[protein1]))
    ppi_counter_thresh += 1 

print(str(ppi_counter_thresh),'ppi out of',str(ppi_counter),'ppi is above threshold')
f.close()
