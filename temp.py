import shelve
import json

with shelve.open('/Users/hui/Scripts/plojo/plojo-nior/hplc_data','rb') as f:
    # data = f['index']
    data=f['ams41']


print(json.dumps(data,indent=4))


print(1)

a="""CGCCCTCGTCCCATCTC TTT GCCAGTGAGTCGGATCTC CGCATATCTGC GAACACCAACCGAGAACG"""

swap='NNNNNN'

mut=[]
for i in range(0,19-len(swap),3):
    new = a[0:22+i]+swap+a[22+i+len(swap):]
    mut.append(new)
print('\n'.join(mut))
