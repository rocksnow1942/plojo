import shelve
from os import path

new = 'new_data'
file_name = 'hplc_data'
file_path = '/Users/hui/Documents/PycharmProjects/HPLCviewer/data_storage'
upload_data_folder = '/Users/hui/Documents/PycharmProjects/HPLCviewer/data_upload'
temp_position = '/Users/hui/Documents/PycharmProjects/HPLCviewer/data_storage/temp'


to_check = 'ams2-run2'

run = to_check.split('-')[0]
load = {}
with shelve.open(path.join(file_path,file_name)) as hd:
    for key,item in hd.items():
        load.update({key:item})

with shelve.open(path.join(file_path,new)) as f:
    for key,item in load.items():
        f[key] = item








print(data_exp)
print(a)
print(b)


import json
import shelve

with open('data.json','rt') as f:
    data = json.load(f)
ams_2 = {'ams2':{'name':'Old HPLC Data','date':'20190218','author':'jones'}}
rd_exp = {}
rd_exp_raw = {}
for key,item in data.items():
    if item.get('fit_method','none') == 'hplc':
        run_name = item['name']
        run_date = item['date']
        run_time = item['concentration']
        run_start = run_time[0]
        run_end = run_time[-1]
        run_len = len(run_time)
        delta_t = run_time[1]-run_time[0]
        if delta_t > 0.006:
            curve ='A'
        else:
            curve = 'B'
        run_sig = item['signal']
        if run_name not in rd_exp.keys():
            rd_exp_raw.update({run_name:{curve:{'time':(run_start,run_end,run_len),'signal':run_sig}}})
            if curve == 'A':
                rd_exp.update({run_name:{'date':run_date,'name':run_name,'speed':0.5,curve:{'extcoef':40}}})
            else:
                rd_exp.update({run_name:{'date':run_date,'name':run_name,'speed':0.5,curve:{}}})
        else:
            rd_exp_raw[run_name].update({curve:{'time':(run_start,run_end,run_len),'signal':run_sig}})
            if curve == 'A':
                rd_exp[run_name].update({curve:{'extcoef':40}})
            else:
                rd_exp[run_name].update({curve:{}})

exp = {}
raw = {}

for i,j in enumerate(rd_exp.keys()):
    exp.update({'ams2-run'+str(i):rd_exp[j]})
    raw.update({'ams2-run'+str(i):rd_exp_raw[j]})


with shelve.open('hplc_data',writeback=True) as hd:
    hd['index'].update(ams_2)
    hd['ams2'] = exp
    hd['ams2raw'] = raw
    hd.sync()
