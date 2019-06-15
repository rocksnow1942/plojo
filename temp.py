import shelve
import json

with shelve.open('/Users/hui/Scripts/plojo/plojo-nior/hplc_data','rb') as f:
    # data = f['index']
    data=f['ams41']


print(json.dumps(data,indent=4))


print(1)
