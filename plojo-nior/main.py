from bokeh.plotting import figure, ColumnDataSource
from bokeh.io import curdoc
from bokeh.models import HoverTool, Plot,CustomJS,LinearAxis, Range1d,LabelSet,Arrow,VeeHead
from bokeh.models.glyphs import Text
from bokeh.models.widgets import Panel, Tabs, Button, TextInput, Select, MultiSelect, RadioButtonGroup, PasswordInput, PreText, DataTable, TableColumn, TextAreaInput,Dropdown,Div,CheckboxGroup
from bokeh.layouts import widgetbox, row, column, layout
from bokeh.palettes import Category10
import shelve
import datetime,time
import os
import glob
from io import BytesIO
import base64
import pandas as pd
import numpy as np
import copy


file_name = 'hplc_data'
file_path = os.path.dirname(__file__)
upload_data_folder = '/Users/hui/Documents/PycharmProjects/HPLCviewer/data_upload'
temp_position = file_path


global data_index,info_deque_holder,current_time,temp_data_to_save,raw_data,axis_label_dict,user_pwd,copyed_runs
user_pwd = {'hui':['h']}
axis_label_dict = {'A': 'DAD signal / mA.U.', 'B':'Solvent B Gradient %' ,'P': 'Pressure / bar', 'default': 'Arbitrary Unit' }
current_time = datetime.datetime.now()
info_deque_holder = ['Welcome!']*3
upload_file_source = ColumnDataSource({'file_contents':[],'file_name':[]})
curve_type_options = ['A','B','P','A1','A2','A3']
copyed_runs = []

view_data_plot = figure(plot_width=1200,plot_height=300)
view_data_plot.annulus(x=[1, 2, 3,4,5], y=[1, 2, 3,4,5], color="hotpink",inner_radius=0.2, outer_radius=0.5)
updatenote = ColumnDataSource(dict(x=[1],y=[3.7],text=['Update:\n 1.Add Repair tool to repair previous broken data due to integration bug. \n2. Add tools to copy and paste run.\n 3. Integration bugs fixed.']))
update_ = Text(x='x', y='y', text='text',text_font_size='12pt',text_font='helvetica')
view_data_plot.add_glyph(updatenote,update_)


class Data():
    """
    class to interact with shelve data storage.
    data stored in shelvs as follows:
    index : a dict contain experiment information; { ams5:{'name': 'experiment name','date': "20190101",'author':'jones'}, ams23:{}}
    meta information and raw data of each experiment is under: "amsXX" and "amsXXraw"
    meta information example:
            'ams39':{
            "ams39-run0": {
                "speed": 1.0,
                "date": "20190409",
                "name": "PHENORP YPEGBS3 2-5V1RXN 0-5MMBS3 35C 1UG 5-60 28M 1MM",
                "A": {
                    "y_label": "DAD signal / mA.U.",
                    "extcoef": 33.0
                },
                "B": {
                    "y_label": "Solvent B Gradient %"
                }
            },
            "ams39-run1": { }, "ams39-run2": { }}
    raw data example:
        ams39raw:{'ams39-run0':{'A':{'time':[0.0,39.99,6000]},'signal':[0,0.1,...]},'B':{'time':},'ams39-run1'}
    when load in to Data Class,
    index is loaded to data.index
    meta information is loaded to data.experiment = {'ams39':{}}
    raw data is loaded to data.experiment_raw = {'ams39':{}}  !!! caution: 'raw' tag is discarded during loading.
    """
    def __init__(self,data_index):
        self.index = data_index
        self.experiment = {} # {ams0:{ams0-run1:{date: ,time: , A:{}}}}
        self.experiment_to_save = {}
        self.experiment_raw = {}
        self.experiment_load_hist = []
        self.max_load = 200
        self.cds = {}
    def next_index(self,n):
        entry = list(self.index.keys())
        if not entry:
            entry = ['ams0']
        entry = sorted(entry, key=lambda x: int(x.split('-')[0][3:]), reverse=True)[0]
        entry_start = int(entry.split('-')[0][3:])+1
        new_entry_list=['ams'+str(i) for i in range(entry_start, entry_start+n)]
        return new_entry_list
    def next_run(self,index,n):
        entry = list(self.experiment[index].keys())
        if not entry:
            entry_start = 0
        else:
            entry = sorted(entry, key=lambda x: int(x.split('-')[1][3:]), reverse=True)[0]
            entry_start = int(entry.split('-')[1][3:])+1
        new_entry_list=[index+'-run'+str(i) for i in range(entry_start, entry_start+n)]
        return new_entry_list
    def load_experiment(self,new):
        if self.max_load < len(self.experiment.keys()):
            to_delete =[]
            for i in self.experiment_load_hist[0:100]:
                if i not in self.experiment_to_save.keys():
                    del self.experiment[i]
                    del self.experiment_raw[i]
                    to_delete.append(i)
            self.experiment_load_hist = [i for i in self.experiment_load_hist if i not in to_delete ]
        new_load = list(set(new)-raw_data.experiment.keys())
        if new_load:
            with shelve.open(os.path.join(file_path,file_name)) as hd:
                for i in new_load:
                    raw_data.experiment[i] = hd.get(i,{})
                    raw_data.experiment_raw[i] = hd.get(i+'raw',{})
                    self.experiment_load_hist.append(i)


with shelve.open(os.path.join(file_path,file_name),writeback=False) as hd:
    data_index = hd['index']
    raw_data = Data(data_index)


# functions

def info_deque(text):
    global info_deque_holder
    j = len(info_deque_holder)
    info_deque_holder.append(str(j)+'>>'+text)
    result = '\n'.join(info_deque_holder[-3:])
    return result

def experiment_menu_generator(items):
    menu = []
    dates = []
    for i in items:
        name = raw_data.index[i].get('name', 'No_Name')
        name = (name[:25]+'...') if len(name)>30 else name
        menu.append(' '.join([raw_data.index[i].get('date', '000000')[-4:],name]))
        dates.append(int(raw_data.index[i].get('date', '00000000')))
    result = list(zip(items, menu))
    result = [i for _,i in sorted(zip(dates,result),reverse=True)]
    return result


def runs_menu_generator(items):
    menu = []
    value_list = []
    for i in items:
        for j in raw_data.experiment[i].keys():
            value_list.append(j)
            curve = list(raw_data.experiment[i][j].keys())
            curve = [i for i in curve if i in curve_type_options]
            name = raw_data.experiment[i][j].get('name','No_name')
            name = (name[:35]+'...') if len(name)>35 else name
            menu.append(' '.join([j,raw_data.experiment[i][j].get('date','No_date')[-4:],name,str(curve)]))
    result = list(zip(value_list, menu))
    result = sorted(result, key=lambda x: (int(x[0].split('-')[0][3:]),int(x[0].split('-')[1][3:])), reverse=True)
    return result



# widgets
tools_menu = [('Copy Analysis','copy'),('Paste Analysis','paste'),None,('Integrate','integrate'),('Annotate','annotate'),None,('Cut Runs','copy_run'),('Paste Runs','paste_run'),('Repair Broken Data','check')]
mode_selection = RadioButtonGroup(labels=['Upload', 'View', 'Analyze'], button_type='success', width=250)
edit_dropdown = Dropdown(label='Tool Kit', button_type='success',value='integrate',menu = tools_menu,width=150)
info_box = PreText(text='Welcome!',width=400)
top_div_2 = Div(text='',width=30)
top_div_4 = Div(text='',width=30)
load_button = Button (label = 'Load Data',button_type='success',width = 150)
save_button = Button(label = 'Save Data',button_type='danger',width = 150)
top_row_= row(mode_selection,edit_dropdown,top_div_2,info_box,load_button,top_div_4,save_button)


# upload widgets
sd_file_upload = Button(label='Upload Data from File', button_type='success')
sd_folder_upload = Button(label='Upload Data from folder', button_type='success')
sd_experiment = TextInput(title='New Experiment Name:',value='None')
sd_author = TextInput(title='Author', value='none')
sd_date = TextInput(title='Experiment Date',
                    value=current_time.strftime('%Y%m%d'))
sd_run_speed = TextInput(title='Run Speed ml/min',value='0.5')
sd_ext_coef = TextInput(title='Extinction Coefficient ug/ml/OD',value='33')
sd_create_new_exp = Button(label='Create New Experiment', button_type='success')
sd_experiment_list = MultiSelect(title='List of Experiments',options=experiment_menu_generator(raw_data.index.keys()),size=15,width=300)
sd_run_list = MultiSelect(title='List of runs',options=runs_menu_generator([]),size=15,width=400)
sd_save_data = Button(label='Save Input Data', button_type='danger')

sd_row_1 = row(column(sd_author,sd_run_speed),column(sd_date,sd_ext_coef),column(sd_experiment,sd_create_new_exp))
sd_bottom_row = row(sd_file_upload,sd_folder_upload,sd_save_data)



# view data widgets




vd_axis_option_menu = curve_type_options.copy()
vd_axis_option_menu.append('none')
vd_primary_axis_option = Select(title='Primary Axis',value='A',options=vd_axis_option_menu,width=170)
vd_secondary_axis_option = Select(title='Secondary Axis',value='B',options=vd_axis_option_menu,width=170)
vd_anotation_option = MultiSelect(title='Primary',value=['integrate-area'],options=[('none','None'),('integrate-line','Integrate Line'),('integrate-area','Integrate Area'),('integrate-percent','Integrate Percent'),('integrate-mass','Integrate Mass'),('integrate-label','Integrate Label'),('annotate','Annotation')],size=7,width=150)
vd_anotation_option_s = MultiSelect(title='Secondary',value=['none'],options=[('none','None'),('integrate-line','Integrate Line'),('integrate-area','Integrate Area'),('integrate-percent','Integrate Percent'),('integrate-mass','Integrate Mass'),('integrate-label','Integrate Label'),('annotate','Annotation')],size=7,width=150)
vd_secondary_axis_range = TextInput(title = 'Secondary Axis Range',value='0-100',width=150)
vd_div_0 = Div(text='',width=50)
vd_plot_backend = Select(title='Plot Format', value='PNG',options = ['PNG','SVG'],width = 170)
vd_offset_option = TextInput(title = 'Secondary Axis Offset',value='0',width=150)
vd_plot_options = row(widgetbox(vd_primary_axis_option,vd_secondary_axis_option,vd_secondary_axis_range,vd_offset_option,width=200),column(row(vd_anotation_option,vd_anotation_option_s),vd_plot_backend))
# vd_plot_options = row(column(vd_primary_axis_option,vd_secondary_axis_option,vd_secondary_axis_range),vd_div_0,column(row(vd_anotation_option,vd_anotation_option_s),vd_offset_option))
vd_selection_tab = Tabs(active=0, width=800, height=280, tabs=[Panel(child=row(sd_experiment_list,sd_run_list), title="Experiment"),Panel(child=vd_plot_options,title='Plot Options')])
vd_div_1 = Div(text='',width=50)
vd_exp_name = TextInput(title='Experiment Name : ')
vd_exp_tag = TextAreaInput(title='Experiment Note :', rows=5, cols=35, max_length=50000)
vd_run_name = TextInput(title='Run name :')
vd_run_ext_coef = TextInput(title='Primary Extinction Coefficient ug/ml/OD : ')
vd_run_speed = TextInput(title='Run Speed (ml/min):')
vd_run_date = TextInput(title='Run Date :')
vd_run_note = TextInput(title='Run Tag :')
vd_info_widgets = widgetbox(vd_exp_name,vd_exp_tag,vd_run_name,vd_run_note,vd_run_ext_coef,vd_run_speed,vd_run_date)
search_field_options = [('exp_name','Experiment Name'),('exp_tag','Experiment Note'),('exp_date','Experiment Date')]
vd_search_field = MultiSelect(title='Search field:', value=['exp_name','exp_tag'], options=search_field_options, size=3,width = 300)
vd_save_info_button = Button(label='Save Info Changes', button_type='danger',width=300)
vd_delete_data_button = Button(label='Delete Selected Data', button_type='success',width=300)
vd_search_keyword = TextInput(title='Keyword filter',width=300)
vd_search_botton = Button(label = 'Search Experiment',button_type='success',width=300)
vd_div_2 = Div(text='',width=105)
vd_div_3 = Div(text='',width=25)
vd_div_4 = Div(text='',width=50)
vd_div_5 = Div(text='',width=150,height=50)
vd_button_widgets = widgetbox(vd_search_botton,vd_save_info_button,vd_delete_data_button)
vd_search_box = row(vd_div_4,column(vd_search_keyword,vd_search_field),vd_div_3,column(vd_div_2,vd_button_widgets))

#analyze layout widgets

# integrate tools

it_start_x = TextInput(title='Integration Start X',width=70,value='none')
it_start_y = TextInput(title='Integration Start Y',width=70,value='0')
it_end_x = TextInput(title='Integration End X',width=70,value='none')
it_end_y = TextInput(title='Integration End Y',width=70,value='0')
it_integration_list = MultiSelect(title='Integration List',value=[],options=[],size=15,width=350)
it_add_button = Button(label='New Integration', button_type='success',width=190)
it_delete_button = Button(label='Delete Integration', button_type='success',width=190)
it_update_button = Button(label='Update Integration', button_type='success',width=190)
it_integration_name = TextInput(title='Integration Name',width=150,value='none')
it_ext_coef= TextInput(title='Extinction Coef. (DO)',width=150,value='none')
it_run_speed = TextInput(title='Run Speed ml/min (DO)',width=150,value='none')
it_div_1 = Div(text='',width=150)
it_div_2=Div(text='',width=150)
it_div_3 = Div(text='',width=100)
it_div_4 = Div(text='',width=70)
it_div_5 = Div(text='',width=30)
it_tool_box = row(it_integration_list,it_div_4,column(row(it_start_x,it_div_1,it_end_x),row(it_start_y,it_div_2,it_end_y),row(it_add_button,it_div_5,it_update_button),it_delete_button),it_div_3,column(it_integration_name,it_ext_coef,it_run_speed))


# annotate tools

an_label = TextInput(title='New Annotation Label',value='none')
an_x = TextInput(title='Annotation Position',value='none')
an_height = TextInput(title='Annotation Height',value='none')
an_list = MultiSelect(title='Annotation List',value=[],options=[],size=15,width=350)
an_div_4 = Div(text='',width=70)
an_div_1 = Div(text='',width=60)
an_div_0 = Div(text='',width=70)
an_add_button = Button(label='New Annotation', button_type='success')
an_delete_button = Button(label='Delete Annotation', button_type='success')
an_update_button = Button(label='Update Annotation', button_type='success')
annotate_tool_box = row(an_list,an_div_4,column(an_label,an_x,an_height),an_div_1,column(an_div_0,an_add_button,an_update_button,an_delete_button))



global analysis_temp_keep
analysis_temp_keep = {}


def copy_analysis():
    """
    copy analysis will copy the selected annotation options for primary_axis:
    i.e. if both integration and annotation are selected, both will be copyed to clipboard.
    """
    if len(sd_run_list.value)!=1:
        info_box.text = info_deque('Select Only 1 run to copy analysis.')
        raise ValueError ('sele 1')
    run_index = sd_run_list.value[0]
    index= run_index.split('-')[0]
    curve = vd_primary_axis_option.value
    target = raw_data.experiment[index][run_index]
    analysis = set([i.split('-')[0] for i in vd_anotation_option.value if i!='none'])
    if len(analysis) == 0:
        info_box.text = info_deque('Select primary_axis analysis in plot options.')
        raise ValueError ('sele 1')
    info_box.text = info_deque('Copying...')
    temp={}
    for key in analysis:
        if key == 'integrate':
            inte_para = target.get(curve,{}).get(key,{}).get('inte_para',[]).copy()
            label = target.get(curve,{}).get(key,{}).get('label',[]).copy()
            temp.update({'integrate':{'inte_para':inte_para,'label':label,'integrate_gap_x':[],'integrate_gap_y':[],'label_cordinate_x':[],'label_cordinate_y':[],'area':[]}})
        elif key == 'annotate':
            height = target.get(curve,{}).get(key,{}).get('height',[]).copy()
            x = target.get(curve,{}).get(key,{}).get('x',[]).copy()
            label = target.get(curve,{}).get(key,{}).get('label',[]).copy()
            temp.update({'annotate':{'height':height,'x':x,'label':label,'y':[]}})
        info_box.text = info_deque('{} copyed.'.format(key))
        global analysis_temp_keep
        analysis_temp_keep = copy.deepcopy(temp)

def check_run_analysis(run_index):
    """
    validate integration and annotation parameter
    generate result.
    """
    index= run_index.split('-')[0]
    target = raw_data.experiment[index][run_index]
    repair=False
    for k,i in target.items():
        if k in curve_type_options:
            time_max = raw_data.experiment_raw[index][run_index][k]['time'][1]
            for key,item in i.items():
                if key == 'integrate':
                    to_keep= [i for i in zip(item['inte_para'],item['label']) if i[0][1]<=time_max]
                    if len(to_keep) != len(item['inte_para']):
                        repair = True
                    item['inte_para']=[i[0] for i in to_keep]
                    item['label']=[i[1] for i in to_keep]
                    item.update({'integrate_gap_x':[],'integrate_gap_y':[],'label_cordinate_x':[],'label_cordinate_y':[],'area':[]})
                    generate_integration_from_para(run_index,k)
                elif key == 'annotate':
                    to_keep= [i for i in zip(item['height'],item['x'],item['label']) if i[1]<=time_max]
                    if len(to_keep) != len(item['x']):
                        repair = True
                    item.update({'height':[i[0] for i in to_keep],'x':[i[1] for i in to_keep],'label':[i[2] for i in to_keep],'y':[]})
                    annotate_generator(run_index,k)
    return repair

def check_selected_runs():
    """
    check the analysis of selected runs.
    """
    entry = sd_run_list.value
    if mode_selection.active != 1:
        info_box.text = info_deque('Go to View tab to Select Experiments/Runs.')
        raise ValueError('go to view tab.')
    if not entry:
        to_repair = [i[0] for i in sd_run_list.options]
        info_box.text = info_deque('No runs selected, Repair {} runs in all selected experiments...'.format(len(to_repair)))
    else:
        to_repair = entry
        info_box.text = info_deque('Start repair {} runs...'.format(len(to_repair)))
    count=0
    for i in to_repair:
        _=check_run_analysis(i)
        if _:
            count+=1
            raw_data.experiment_to_save.update({i.split('-')[0]:'sync'})
            info_box.text = info_deque('Run {} repaired.'.format(i))
    info_box.text = info_deque('Done! {} runs repaired.'.format(count))
    if count:
        info_box.text = info_deque('!!!DON\'T forget to SAVE!!!')

def paste_analysis():
    """
    paste analysis to selected runs.
    will only paste down integration within allowed run time.
    """
    global analysis_temp_keep
    run_index = sd_run_list.value
    if not run_index:
        info_box.text = info_deque('Select runs to paste to.')
        raise ValueError('select run to paste.')
    index= [i.split('-')[0] for i in run_index]
    curve = vd_primary_axis_option.value
    for i,ri in zip(index,run_index):
        target = raw_data.experiment[i][ri]
        if target.get(curve,'none') == 'none':
            info_box.text = info_deque('{} doesn\'t have curve {}'.format(ri,curve))
        else:
            target[curve].update(copy.deepcopy(analysis_temp_keep))
            _=check_run_analysis(ri)
        raw_data.experiment_to_save.update({i:'sync'})
        info_box.text = info_deque('Pasted to {}.'.format(ri))
    sync_plot(run_index,**vd_plot_options_fetcher())



def copy_selected_run():
    global copyed_runs
    entry = sd_run_list.value
    if mode_selection.active != 1:
        info_box.text = info_deque('Go to View tab to cut.')
        raise ValueError('go to view tab to cut.')
    copyed_runs = entry.copy()
    if not entry:
        info_box.text = info_deque('No runs selected, clipboard cleared.')
    else:
        info_box.text = info_deque('{} runs in clipboard.'.format(len(copyed_runs)))

def paste_selected_run():
    global copyed_runs
    if mode_selection.active != 1:
        info_box.text = info_deque('Go to View tab to Paste.')
        raise ValueError('go to view tab to Paste.')
    target = sd_experiment_list.value
    if len(target) != 1:
        info_box.text = info_deque('Selecte only 1 Experiment to Paste to.')
        raise ValueError('SELECT 1 exp to paste.')
    if len(copyed_runs)==0:
        info_box.text = info_deque('Clipboard is empty, Cut runs first.')
        raise ValueError('SELECT before paste.')
    target = target[0]
    new_run = raw_data.next_run(target,len(copyed_runs))
    try:
        for i,j in zip(copyed_runs,new_run):
            exp = i.split('-')[0]
            raw_data.experiment[target].update({j:raw_data.experiment[exp].pop(i)})
            raw_data.experiment_raw[target].update({j:raw_data.experiment_raw[exp].pop(i)})
            raw_data.experiment_to_save.update({exp:'upload'})
        raw_data.experiment_to_save.update({target:'upload'})
        sd_run_list.options = runs_menu_generator([target])
        info_box.text = info_deque('{} runs moved to {}.'.format(len(copyed_runs),target))
        copyed_runs = []
    except:
        load_button_cb()
        info_box.text = info_deque('Error occured during pasting. No changes were made.')
        copyed_runs=[]

def edit_dropdown_cb(attr,old,new):
    try:
        if new == 'integrate':
             analyze_layout.children[2].children = [it_tool_box]
        elif new == 'annotate':
             analyze_layout.children[2].children = [annotate_tool_box]
        elif new == 'copy':
            copy_analysis()
        elif new == 'paste':
            paste_analysis()
        elif new == 'copy_run':
            copy_selected_run()
        elif new == 'paste_run':
            paste_selected_run()
        elif new == 'check':
            check_selected_runs()
        else:
            pass
    except:
        pass
    if new != 'none':
        edit_dropdown.value = 'none'


edit_dropdown.on_change('value',edit_dropdown_cb)

def it_para_check(xe,run_index,curve,xs=0):
    """
    check if the x value of integration interval is within the time of current run.
    """
    index= run_index.split('-')[0]
    time_ = raw_data.experiment_raw[index][run_index][curve]['time']
    if xs>=xe or xe>time_[1]:
        info_box.text = info_deque('Time parameter error, action aborted.')
        raise ValueError('time para error.')



def it_add_button_cb():
    """
    add integration
    """
    run_index = sd_run_list.value[0]
    index= run_index.split('-')[0]
    curve = vd_primary_axis_option.value
    label = it_integration_name.value
    try:
        # raw_data.experiment[index][run_index].update(extcoef=float(it_ext_coef.value))
        xs= float(it_start_x.value)
        xe= float(it_end_x.value)
        ys = float(it_start_y.value)
        ye = float(it_end_y.value)
    except:
        info_box.text = info_deque('enter valid numbers')
        raise ValueError('wrong input')
    it_para_check(xe,run_index,curve,xs)
    inte_para= (xs,xe,ys,ye)
    target = raw_data.experiment[index][run_index][curve]
    if not target.get('integrate',{}):
        target.update(integrate={'inte_para':[inte_para],'integrate_gap_x':[],'integrate_gap_y':[],'label_cordinate_x':[],'label_cordinate_y':[],'area':[],'label':[label]})
    else:
        target['integrate']['inte_para'].append(inte_para)
        target['integrate']['label'].append(label)
    generate_integration_from_para(run_index,curve)
    plot_para = vd_plot_options_fetcher()
    plot_para.update(analysis=['integrate-area','integrate-mass','integrate-label','integrate-percent'])
    sync_plot([run_index],**plot_para)
    it_integration_list.options = it_integration_list_menu_generator(run_index,curve)
    it_integration_list.value = [it_integration_list.options[-1][0]]
    raw_data.experiment_to_save.update({index:'sync'})



it_add_button.on_click(it_add_button_cb)


def it_update_button_cb():
    it_delete_button_cb()
    it_add_button_cb()

it_update_button.on_click(it_update_button_cb)


def generate_integration_from_para(run_index,curve):
    index= run_index.split('-')[0]
    target = raw_data.experiment[index][run_index][curve]['integrate']
    time_ = raw_data.experiment_raw[index][run_index][curve]['time']
    delta_t = (time_[1]-time_[0])/(time_[2]-1)
    signal=raw_data.experiment_raw[index][run_index][curve]['signal']
    update_start = len(target['integrate_gap_x'])
    for xs,xe,ys,ye in target['inte_para'][update_start:]:
        xs_i = int(round((xs-time_[0])/delta_t))
        xe_i = int(round((xe-time_[0])/delta_t))
        ym = signal[int((xs_i+xe_i)/2)]
        area = sum(signal[xs_i:xe_i])*(delta_t)-(ys+ye)*(xe-xs)/2
        target['integrate_gap_x'].append([ xs,xs ,xe ,xe ])
        target['integrate_gap_y'].append([ signal[xs_i],ys ,ye ,signal[xe_i] ])
        target['label_cordinate_x'].append((xs+xe)/2)
        target['label_cordinate_y'].append(ym)
        target['area'].append(area)

# call backs


def it_integration_list_menu_generator(run_index,curve):
    index= run_index.split('-')[0]
    inte_para = raw_data.experiment[index][run_index].get(curve,{}).get('integrate',{}).get('inte_para',[])
    inte_label = raw_data.experiment[index][run_index].get(curve,{}).get('integrate',{}).get('label',[])
    menu = []
    for (i,j),k in zip(enumerate(inte_para),inte_label):
        xs = str(round(j[0],3))
        xe =str(round(j[1],3))
        ys =str(round(j[2],3))
        ye =str(round(j[3],3))
        value = k + ' ('+xs+','+ys+') -> ('+xe+','+ye+')'
        menu.append((str(i),value))
    return menu

def an_list_menu_generator(run_index,curve):
    index= run_index.split('-')[0]
    x_position = raw_data.experiment[index][run_index].get(curve,{}).get('annotate',{}).get('x',[])
    height = raw_data.experiment[index][run_index].get(curve,{}).get('annotate',{}).get('height',[])
    label = raw_data.experiment[index][run_index].get(curve,{}).get('annotate',{}).get('label',[])
    menu = []
    for (i,j),k,l in zip(enumerate(x_position),height,label):

        value = l + ' -> ' + str(round(j,3)) +' - ' + str(int(k))
        menu.append((str(i),value))
    return menu


def an_list_cb(attr,old,new):
    if len(new)>0:
        if len(new)>1:
            info_box.text = info_deque('Only first sele used.')
        run_index = sd_run_list.value[0]
        index= run_index.split('-')[0]
        curve = vd_primary_axis_option.value
        an_index = int(new[0])
        label = raw_data.experiment[index][run_index].get(curve,{}).get('annotate',{}).get('label',[])[an_index]
        x_position = raw_data.experiment[index][run_index].get(curve,{}).get('annotate',{}).get('x',[])[an_index]
        height = raw_data.experiment[index][run_index].get(curve,{}).get('annotate',{}).get('height',[])[an_index]
        an_label.value = str(label)
        an_x.value = str(x_position)
        an_height.value = str(height)
    else:
        pass

an_list.on_change('value',an_list_cb)


def annotate_generator(run_index,curve):
    index= run_index.split('-')[0]
    target =  raw_data.experiment[index][run_index][curve]['annotate']
    signal=raw_data.experiment_raw[index][run_index][curve]['signal']
    time_ = raw_data.experiment_raw[index][run_index][curve]['time']
    delta_t = (time_[1]-time_[0])/(time_[2]-1)
    x_position = target['x']
    update_start = len(target['y'])
    for i in x_position[update_start:]:
        y_position = signal[int(round((i-time_[0])/delta_t))]
        target['y'].append(y_position)



def an_add_button_cb():
    run_index = sd_run_list.value[0]
    index= run_index.split('-')[0]
    curve = vd_primary_axis_option.value
    label = it_integration_name.value
    try:
        x_position = float(an_x.value)
        height = float(an_height.value)
        label = an_label.value
        signal=raw_data.experiment_raw[index][run_index][curve]['signal']
        time_ = raw_data.experiment_raw[index][run_index][curve]['time']
        delta_t = (time_[1]-time_[0])/(time_[2]-1)
        y_position = signal[int(round((x_position-time_[0])/delta_t))]
    except:
        info_box.text = info_deque('enter valid numbers')
        raise ValueError('wrong input')
    target =  raw_data.experiment[index][run_index][curve]
    it_para_check(x_position,run_index,curve)
    if not target.get('annotate',{}):
        target.update(annotate={'height':[height],'label':[label],'x':[x_position],'y':[y_position]})
    else:
        target['annotate']['height'].append(height)
        target['annotate']['label'].append(label)
        target['annotate']['x'].append(x_position)
        target['annotate']['y'].append(y_position)
    plot_para = vd_plot_options_fetcher()
    plot_para.update(analysis=['annotate'])
    sync_plot([run_index],**plot_para)
    an_list.options = an_list_menu_generator(run_index,curve)
    an_list.value = [an_list.options[-1][0]]
    raw_data.experiment_to_save.update({index:'sync'})

an_add_button.on_click(an_add_button_cb)

def an_delete_button_cb():
    to_delete=an_list.value
    if len(to_delete)!=1:
        info_box.text = info_deque('select 1 to delete')
        raise ValueError('sele 1')
    run_index = sd_run_list.value[0]
    index= run_index.split('-')[0]
    to_delete = int(to_delete[0])
    curve = vd_primary_axis_option.value
    target = raw_data.experiment[index][run_index].get(curve,{}).get('annotate',{})
    for key in target.keys():
        try:
            del target[key][to_delete]
        except:
            pass
    an_list.options = an_list_menu_generator(run_index,curve)
    an_list.value = []
    raw_data.experiment_to_save.update({index:'sync'})
    plot_para = vd_plot_options_fetcher()
    plot_para.update(analysis=['annotate'])
    sync_plot([run_index],**plot_para)

def an_update_button_cb():
    an_delete_button_cb()
    an_add_button_cb()



an_delete_button.on_click(an_delete_button_cb)
an_update_button.on_click(an_update_button_cb)

def it_integration_list_cb(attr,old,new):
    if len(new)>0:
        if len(new)>1:
            info_box.text = info_deque('Only first sele used.')
        run_index = sd_run_list.value[0]
        index= run_index.split('-')[0]
        curve = vd_primary_axis_option.value
        int_index = int(new[0])
        inte_para = raw_data.experiment[index][run_index].get(curve,{}).get('integrate',{}).get('inte_para',[])[int_index]
        inte_label = raw_data.experiment[index][run_index].get(curve,{}).get('integrate',{}).get('label',[])[int_index]
        it_start_x.value = str(inte_para[0])
        it_end_x.value = str(inte_para[1])
        it_start_y.value = str(inte_para[2])
        it_end_y.value = str(inte_para[3])
        it_integration_name.value = inte_label
    else:
        pass

def it_delete_button_cb():
    to_delete=it_integration_list.value
    if len(to_delete)!=1:
        info_box.text = info_deque('select 1 to delete')
        raise ValueError('sele 1')
    run_index = sd_run_list.value[0]
    index= run_index.split('-')[0]
    to_delete = int(to_delete[0])
    curve = vd_primary_axis_option.value
    target = raw_data.experiment[index][run_index].get(curve,{}).get('integrate',{})
    for key in target.keys():
        try:
            del target[key][to_delete]
        except:
            pass
    it_integration_list.options = it_integration_list_menu_generator(run_index,curve)
    it_integration_list.value = []
    raw_data.experiment_to_save.update({index:'sync'})
    plot_para = vd_plot_options_fetcher()
    plot_para.update(analysis=['integrate-area','integrate-mass','integrate-label'])
    sync_plot([run_index],**plot_para)





it_integration_list.on_change('value',it_integration_list_cb)
it_delete_button.on_click(it_delete_button_cb)


def sync_vd_exp_info_widgets(sele):
    global info_change_keep
    vd_exp_name.value = str([raw_data.index[i].get('name','No_Name') for i in sele]).strip('\'[]')
    vd_exp_tag.value = str([raw_data.index[i].get('tag','No_note') for i in sele]).strip('\'[]')
    info_change_keep.update(exp_name=False,exp_tag=False)



def sync_vd_run_info_widgets(sele):
    global info_change_keep
    exp = [i.split('-')[0] for i in sele]
    curve = vd_primary_axis_option.value
    vd_run_ext_coef.value = str([raw_data.experiment[i][j].get(curve,{'extcoef':'No_Curve'}).get('extcoef','No_Ext') for i,j in zip(exp,sele)]).strip('\'[]')
    vd_run_name.value = str([raw_data.experiment[i][j].get('name','No_Name') for i,j in zip(exp,sele)]).strip('\'[]')
    vd_run_speed.value = str([raw_data.experiment[i][j].get('speed','No_Speed') for i,j in zip(exp,sele)]).strip('\'[]')
    vd_run_date.value = str([raw_data.experiment[i][j].get('date','No_Date') for i,j in zip(exp,sele)]).strip('\'[]')
    vd_run_note.value = str([raw_data.experiment[i][j].get('note','No_Tag') for i,j in zip(exp,sele)]).strip('\'[]')
    info_change_keep.update(run_name=False,run_date=False,run_note=False,run_speed=False,run_extcoef=False)


global info_change_keep
info_change_keep = dict.fromkeys(['exp_name','run_extcoef','exp_tag','run_name','run_speed','run_date','run_note'],False)

def vd_save_info_button_cb():
    global info_change_keep
    selected_exp = sd_experiment_list.value
    selected_run = sd_run_list.value
    curve = vd_primary_axis_option.value
    try:
        new_input = dict(zip(['exp_name','exp_tag','run_name','run_extcoef','run_speed','run_date','run_note'],[vd_exp_name.value,vd_exp_tag.value,vd_run_name.value,vd_run_ext_coef.value,vd_run_speed.value,vd_run_date.value,vd_run_note.value]))
        for key,item in info_change_keep.items():
            if item:
                if key.split('_')[0]=='exp':
                    for i in selected_exp:
                        raw_data.index[i].update({key.split('_')[1]:new_input[key]})
                else:
                    for j in selected_run:
                        if key == 'run_speed':
                            raw_data.experiment[j.split('-')[0]][j].update({key.split('_')[1]:float(new_input[key])})
                        elif key == 'run_extcoef':
                            raw_data.experiment[j.split('-')[0]][j][curve].update({key.split('_')[1]:float(new_input[key])})
                        else:
                            raw_data.experiment[j.split('-')[0]][j].update({key.split('_')[1]:new_input[key]})
                        raw_data.experiment_to_save.update({j.split('-')[0]:'sync'})
        sd_save_data_cb()
        info_change_keep = dict.fromkeys(info_change_keep.keys(),False)
    except:
        info_box.text = info_deque('Wrong input format.')
    sd_experiment_list.options = experiment_menu_generator(raw_data.index.keys())
    sd_run_list.options = runs_menu_generator(selected_exp)



def vd_exp_name_cb(attr,old,new):
    info_change_keep.update(exp_name=True)
def vd_exp_tag_cb(attr,old,new):
    info_change_keep.update(exp_tag=True)
def vd_run_name_cb(attr,old,new):
    info_change_keep.update(run_name=True)
def vd_run_speed_cb(attr,old,new):
    info_change_keep.update(run_speed=True)
def vd_run_date_cb(attr,old,new):
    info_change_keep.update(run_date=True)
def vd_run_note_cb(attr,old,new):
    info_change_keep.update(run_note=True)
def vd_run_ext_coef_cb(attr,old,new):
    info_change_keep.update(run_extcoef=True)

vd_save_info_button.on_click(vd_save_info_button_cb)
vd_exp_name.on_change('value',vd_exp_name_cb)
vd_run_ext_coef.on_change('value',vd_run_ext_coef_cb)
vd_exp_tag.on_change('value',vd_exp_tag_cb)
vd_run_name.on_change('value',vd_run_name_cb)
vd_run_speed.on_change('value',vd_run_speed_cb)
vd_run_date.on_change('value',vd_run_date_cb)
vd_run_note.on_change('value',vd_run_note_cb)


def generate_cds(run_index,**kwargs):
    offset_curve,offset = kwargs.get('offset',('A',0))
    index = run_index.split('-')[0]
    run_speed = raw_data.experiment[index][run_index].get('speed',0)
    run = raw_data.experiment[index][run_index]
    run_raw = raw_data.experiment_raw[index][run_index]
    run_dict = dict.fromkeys(set(curve_type_options))
    integration_dict = dict.fromkeys(set(curve_type_options))
    annotate = dict.fromkeys(set(curve_type_options))
    for curve in curve_type_options:
        curve_dict = run.get(curve,{})
        extinction_coeff = curve_dict.get('extcoef',0)
        curve_dict_raw = run_raw.get(curve,{})
        curve_time = np.linspace(*curve_dict_raw.get('time',(0,0,1)))+offset*int(curve==offset_curve)
        curve_signal = curve_dict_raw.get('signal',[0])
        curve_integrate = curve_dict.get('integrate',{})
        integrate_gap_x=curve_integrate.get('integrate_gap_x',[])
        integrate_gap_y=curve_integrate.get('integrate_gap_y',[])
        label_cordinate_x = curve_integrate.get('label_cordinate_x',[])
        label_cordinate_y = curve_integrate.get('label_cordinate_y',[])
        label_ = curve_integrate.get('area',[])
        label = curve_integrate.get('label',[])
        label_area = ['{:.4g}'.format(i) for i in label_]
        integrate_percent = ['{:.2f}%'.format(i*100/(sum(label_))) for i in label_]
        label_mass = ['{:.4g}ug'.format(i*extinction_coeff*run_speed/1000) for i in label_]
        temp_cds = ColumnDataSource( {'time':curve_time,'signal':curve_signal})
        integration_cds = ColumnDataSource({'integrate_percent':integrate_percent,'integrate_gap_x':integrate_gap_x,'integrate_gap_y':integrate_gap_y,'label_cordinate_x':label_cordinate_x,'label_cordinate_y':label_cordinate_y,'label_area':label_area,'label_mass':label_mass,'label':label})
        run_dict.update({curve:temp_cds})
        integration_dict.update({curve:integration_cds})
        # generate annotate
        annotate_dict = curve_dict.get('annotate',{})
        anno_height,anno_x,anno_y,anno_label = [annotate_dict.get(i,[]) for i in ['height','x','y','label']]
        label_data_source = ColumnDataSource({'label_x':anno_x,'label_y':[i+j for i,j in zip(anno_y,anno_height)],'label':anno_label})
        arrow_position = [i for i in zip(anno_x,anno_y,anno_height) ]
        annotate.update({curve:{'label':label_data_source,'arrow':arrow_position}})
    return {'run':run_dict,'integrate':integration_dict,'annotate':annotate}



def plot_generator(plot_list=[],**kwargs):
    plot_backend = kwargs.get('plot_backend','canvas')
    primary_axis = kwargs.get('primary_axis','A')
    secondary_axis=kwargs.get('secondary_axis','B')
    secondary_axis_range = kwargs.get('secondary_axis_range',(0,100))
    color_={'A':'magenta','B':'green','P':'royalblue','A1':'teal','A2':'deepskyblue','A3':'red','none':'red'}
    tools_list = "pan,ywheel_zoom,xwheel_zoom,box_zoom,save,crosshair,reset"
    if len(plot_list) > 1:
        alpha = 0.75
        use_colorp = True
    else:
        use_colorp = False
        alpha = 0.9
    p=figure(plot_width=1200,plot_height=300,tools=tools_list,toolbar_location='above')
    p.toolbar.active_inspect = None
    p.output_backend=plot_backend
    p.xaxis.axis_label = 'Run Time /min'
    p.yaxis.axis_label = axis_label_dict.get(primary_axis[0],'AU')
    p.extra_y_ranges = {'secondary':Range1d(*secondary_axis_range)}
    p.add_layout(LinearAxis(y_range_name='secondary',axis_label=axis_label_dict.get(secondary_axis[0],'AU')),'right')
    analysis = kwargs.get('analysis',[])
    analysis_s = kwargs.get('analysis_s',[])
    analysis = set([ k for j in [i.split('-') for i in analysis] for k in j])
    analysis_s = set([ k for j in [i.split('-') for i in analysis_s] for k in j])
    angle = 0.6
    x_offset = -10
    for i,run_index in enumerate(plot_list):
        cds_dict = generate_cds(run_index,**kwargs)
        tag = raw_data.experiment[run_index.split('-')[0]][run_index].get('note','')
        if use_colorp:
            p_color_touse = Category10[10][i % 10]
            s_color_touse = Category10[10][(i+3) % 10]
        else:
            p_color_touse = color_[primary_axis]
            s_color_touse = color_[secondary_axis]

        if secondary_axis != 'none':
            p_render_2=p.line('time','signal',source=cds_dict['run'][secondary_axis],alpha=alpha,line_width=1.5,color=s_color_touse,legend=run_index+' '+secondary_axis+' '+tag,y_range_name='secondary')
            p_data_hover_2 = HoverTool(tooltips=[('Time: ', '@time{0.000}'), ('Value: ', '@signal{0,0.00}')],renderers=[p_render_2],mode='vline')
            p.add_tools(p_data_hover_2)

        if primary_axis != 'none':
            p_render = p.line('time','signal',source=cds_dict['run'][primary_axis],alpha=alpha,color=p_color_touse,line_width=1.5,legend=run_index+' '+primary_axis+' '+tag)
            p_data_hover = HoverTool(tooltips=[('Time: ', '@time{0.000}'), ('Value: ', '@signal{0,0.00}')],renderers=[p_render],mode='vline')
            p.add_tools(p_data_hover)

        if secondary_axis != 'none':
            integrate_offset = 15
            if 'integrate' in analysis_s:
                p.multi_line(y_range_name='secondary',xs='integrate_gap_x',ys='integrate_gap_y',source=cds_dict['integrate'][secondary_axis],line_width=1,color='black')
            if 'area' in analysis_s:
                int_labels = LabelSet(y_range_name='secondary',x='label_cordinate_x',y='label_cordinate_y',angle=angle, text='label_area', level='glyph',text_font_size='9pt', x_offset=x_offset, y_offset=integrate_offset, source=cds_dict['integrate'][secondary_axis], render_mode='canvas')
                p.add_layout(int_labels)
                integrate_offset +=22
            if 'percent' in analysis_s:
                int_labels = LabelSet(y_range_name='secondary',x='label_cordinate_x',y='label_cordinate_y',angle=angle, text='integrate_percent', level='glyph',text_font_size='9pt', x_offset=x_offset, y_offset=integrate_offset, source=cds_dict['integrate'][secondary_axis], render_mode='canvas')
                p.add_layout(int_labels)
                integrate_offset +=22
            if 'mass' in analysis_s:
                int_labels = LabelSet(y_range_name='secondary',x='label_cordinate_x',y='label_cordinate_y',angle=angle, text='label_mass', level='glyph',text_font_size='9pt', x_offset=x_offset, y_offset=integrate_offset, source=cds_dict['integrate'][secondary_axis], render_mode='canvas')
                p.add_layout(int_labels)
                integrate_offset +=22
            if 'label' in analysis_s:
                int_labels = LabelSet(y_range_name='secondary',x='label_cordinate_x',y='label_cordinate_y', angle=angle,text='label', level='glyph',text_font_size='9pt', x_offset=x_offset, y_offset=integrate_offset, source=cds_dict['integrate'][secondary_axis], render_mode='canvas')
                p.add_layout(int_labels)
                integrate_offset +=22
            if 'annotate' in analysis_s:
                for j in cds_dict['annotate'][secondary_axis]['arrow']:
                    p.add_layout(Arrow(y_range_name='secondary',end=VeeHead(size=8), line_color=s_color_touse,x_start=j[0], y_start=j[1]+j[2], x_end=j[0], y_end=j[1]))
                anno_label = LabelSet(y_range_name='secondary',x='label_x',y='label_y',angle=angle, text='label', level='glyph',text_font_size='9pt', x_offset=0, y_offset=5, source=cds_dict['annotate'][secondary_axis]['label'], render_mode='canvas')
                p.add_layout(anno_label)

        if primary_axis != 'none':
            integrate_offset = 15
            if 'integrate' in analysis:
                p.multi_line(xs='integrate_gap_x',ys='integrate_gap_y',source=cds_dict['integrate'][primary_axis],line_width=1,color='black')
            if 'area' in analysis:
                int_labels = LabelSet(x='label_cordinate_x',y='label_cordinate_y',angle=angle, text='label_area', level='glyph',text_font_size='9pt', x_offset=x_offset, y_offset=integrate_offset, source=cds_dict['integrate'][primary_axis], render_mode='canvas')
                p.add_layout(int_labels)
                integrate_offset +=22
            if 'percent' in analysis:
                int_labels = LabelSet(x='label_cordinate_x',y='label_cordinate_y', angle=angle,text='integrate_percent', level='glyph',text_font_size='9pt', x_offset=x_offset, y_offset=integrate_offset, source=cds_dict['integrate'][primary_axis], render_mode='canvas')
                p.add_layout(int_labels)
                integrate_offset +=22
            if 'mass' in analysis:
                int_labels = LabelSet(x='label_cordinate_x',y='label_cordinate_y',angle=angle, text='label_mass', level='glyph',text_font_size='9pt', x_offset=x_offset, y_offset=integrate_offset, source=cds_dict['integrate'][primary_axis], render_mode='canvas')
                p.add_layout(int_labels)
                integrate_offset +=22
            if 'label' in analysis:
                int_labels = LabelSet(x='label_cordinate_x',y='label_cordinate_y', angle=angle,text='label', level='glyph',text_font_size='9pt', x_offset=x_offset, y_offset=integrate_offset, source=cds_dict['integrate'][primary_axis], render_mode='canvas')
                p.add_layout(int_labels)
                integrate_offset +=22
            if 'annotate' in analysis:
                for j in cds_dict['annotate'][primary_axis]['arrow']:
                    p.add_layout(Arrow(end=VeeHead(size=8), line_color=p_color_touse,x_start=j[0], y_start=j[1]+j[2], x_end=j[0], y_end=j[1]))
                anno_label = LabelSet(x='label_x',y='label_y',angle=angle, text='label', level='glyph',text_font_size='9pt', x_offset=0, y_offset=5, source=cds_dict['annotate'][primary_axis]['label'], render_mode='canvas')
                p.add_layout(anno_label)
    p.legend.click_policy = 'hide'
    p.legend.border_line_alpha = 0
    p.legend.background_fill_alpha = 0.1
    return p


def sync_plot(plot_list,**kwargs):
    p = plot_generator(plot_list,**kwargs)
    plot_row.children = [p]


def vd_plot_options_fetcher():
    pa = vd_primary_axis_option.value
    sa = vd_secondary_axis_option.value
    annotation = vd_anotation_option.value
    annotation_s = vd_anotation_option_s.value
    plot_format = 'canvas' if vd_plot_backend.value=='PNG' else 'svg'
    try:
        a=vd_secondary_axis_range.value.split('-')
        if len(a)>2:
            a = ('-'+a[1],a[2])
        sa_range = tuple(map(float,a))
        assert len(sa_range)==2
        offset = (sa,float(vd_offset_option.value))
    except:
        info_box.text = info_deque('enter valid plot options')
        sa_range = (0,100)
        offset = ('A',0)
    return dict(offset=offset,plot_backend=plot_format,primary_axis=pa,secondary_axis=sa,analysis=annotation,analysis_s=annotation_s,secondary_axis_range=sa_range)

def sd_run_list_cb(attr,old,new):
    if mode_selection.active != 0 and new:
        sync_plot(new,**vd_plot_options_fetcher())
        sync_vd_run_info_widgets(new)
    else:
        pass

def vd_plot_option_cb(attr,old,new):
    selection = sd_run_list.value
    sd_run_list_cb(1,2,selection)



def vd_search_botton_cb():
    keyword = vd_search_keyword.value.split()
    keyword = keyword if keyword else [""]
    field = vd_search_field.value
    search_hit = []
    for key, item in raw_data.index.items():
        for field_ in field:
            hit_found = False
            for i in keyword:
                if i.lower() in item.get(field_.split('_')[1],'').lower():
                    search_hit.append(key)
                    hit_found = True
                    break
            if hit_found:
                break
    sd_experiment_list.options = experiment_menu_generator(search_hit)
    sd_experiment_list.value = []



sd_file_upload.callback = CustomJS(args=dict(file_source=upload_file_source), code = """
function read_file(filename) {
    var reader = new FileReader();
    reader.onload = load_handler;
    reader.onerror = error_handler;
    // readAsDataURL represents the file's data as a base64 encoded string
    reader.readAsDataURL(filename);
}

function load_handler(event) {
    var b64string = event.target.result;
    file_source.data = {'file_contents' : [b64string], 'file_name':[input.files[0].name]};
    file_source.trigger("change");
}

function error_handler(evt) {
    if(evt.target.error.name == "NotReadableError") {
        alert("Can't read file!");
    }
}

var input = document.createElement('input');
input.setAttribute('type', 'file');
input.onchange = function(){
    if (window.FileReader) {
        read_file(input.files[0]);
    } else {
        alert('FileReader is not supported in this browser');
    }
}
input.click();
""")


def mode_selection_cb(attr,old,new):
    if new == 0:
        display_layout.children = upload_layout.children
    elif new == 1:
        display_layout.children = vd_layout.children
    elif new == 2:
        display_layout.children = analyze_layout.children
        if len(sd_run_list.value) ==0:
            info_box.text = info_deque ('Select a run to start!')
            raise ValueError('sele 1')
        else:
            info_box.text = info_deque ('Only first sele used!')
            run_index=sd_run_list.value[0]
            index= run_index.split('-')[0]
            sd_run_list.value = [run_index]
            curve = vd_primary_axis_option.value
        load_analysis_para(index,run_index,curve)
    else:
        pass

def load_analysis_para(index,run_index,curve):
    it_integration_list.options = it_integration_list_menu_generator(run_index,curve)
    an_list.options = an_list_menu_generator(run_index,curve)
    it_ext_coef.value =str(raw_data.experiment[index][run_index].get(curve,{}).get('extcoef','None'))
    it_run_speed.value = str(raw_data.experiment[index][run_index].get('speed','None'))



def sd_data_generator():
    if upload_file_source.data['file_name'][0] == 'nofile':
        result_dict = 'none'
    else:
        try:
            raw_contents = upload_file_source.data['file_contents'][0]
            prefix, b64_contents = raw_contents.split(",", 1)
            file_contents = base64.b64decode(b64_contents)
            data = file_contents.decode("utf-16")
            file_name = upload_file_source.data['file_name'][0].split('.')[0].split('_')
            data = data.strip('\n')
            data = [list(map(float,i.strip('\r').split('\t'))) for i in data.split('\n')]
            data = list(map(list, zip(*data)))
            data_range = (data[0][0],data[0][-1],len(data[0]))
            data_signal = data[1]
            data_date=file_name[0]
            data_name=file_name[1]
            data_key = file_name[2].upper()
            run_speed = float(sd_run_speed.value)
            ext_coef = float(sd_ext_coef.value)
            if data_key not in curve_type_options:
                raise ValueError('not correct curve type')
            meta = {'speed':run_speed,'date':data_date,'name':data_name,data_key:{'y_label':axis_label_dict.get(data_key[0],'default')}}
            if data_key[0]=='A':
                meta[data_key].update(extcoef=ext_coef)
            raw = {data_key:{'time':data_range,'signal':data_signal}}
            result_dict = {'meta':meta,'raw':raw}
        except:
            info_box.text = info_deque('Wrong uploaded.')
            result_dict = 'none'
    return result_dict



def data_dict_to_exp(data,index,run_index):
    if run_index == 'new':
        run_index = raw_data.next_run(index,1)[0]
    else:
        if raw_data.experiment[index][run_index]['name'] != data['meta']['name'] or raw_data.experiment[index][run_index]['date'] != data['meta']['date']:
            info_box.text = info_deque('Input file name didn\'t match selected run name')
            raise ValueError('name mismatch')
    curve = list(data['raw'].keys())
    if run_index not in raw_data.experiment[index].keys():
        raw_data.experiment[index].update({run_index:data['meta']})
        raw_data.experiment_raw[index].update({run_index:data['raw']})
    else:
        raw_data.experiment[index][run_index].update(data['meta'])
        raw_data.experiment_raw[index][run_index].update(data['raw'])
    raw_data.experiment_to_save.update({index:'upload'})
    info_box.text = info_deque('Curve {} uploaded to run {}.'.format(curve[0],run_index))
    sd_run_list.options = runs_menu_generator(sd_experiment_list.value)


def upload_file_source_cb(attr,old,new):
    data = sd_data_generator()
    if data != 'none':
        index = sd_experiment_list.value
        if len(index)!=1:
            info_box.text = info_deque('Select 1 experiment to upload.')
            raise ValueError('sele 1')
        index = index[0]
        try:
            run_index = sd_run_list.value[0]
        except:
            run_index = 'new'
        data_dict_to_exp(data,index,run_index)
        upload_file_source.data={'file_contents':['none'],'file_name':['nofile']}

def sd_save_data_cb():
    cur = time.time()
    file_list = sorted(glob.glob(os.path.join(temp_position,file_name+'*')))
    if len(file_list)>90:
        os.remove(file_list[0])
        os.remove(file_list[1])
        os.remove(file_list[2])
    if len(file_list)>9:
        last_save = os.path.getmtime(file_list[-1])
    else:
        last_save = 0.0
    if cur-last_save> 10800:
        source_1 = os.path.join(file_path,file_name)
        dest = os.path.join(temp_position,file_name+'_'+datetime.datetime.now().strftime('%Y%m%d%H%M'))
        with shelve.open(source_1) as old:
            with shelve.open(dest) as new:
                for key,item in old.items():
                    new[key] = old[key]
    with shelve.open(os.path.join(file_path,file_name),writeback=False) as hd:
        hd['index'] = raw_data.index
        for key,item in raw_data.experiment_to_save.items():
            if item == 'sync':
                hd[key]=raw_data.experiment[key]
            elif item == 'del':
                del hd[key]
                del hd[key+'raw']
            elif item == 'upload':
                hd[key]=raw_data.experiment[key]
                hd[key+'raw']=raw_data.experiment_raw[key]
    info_box.text = info_deque('{} experiments have been saved'.format(len(raw_data.experiment_to_save.keys())))
    raw_data.experiment_to_save = {}


def sd_experiment_list_cb(attr,old,new):
    raw_data.load_experiment(new)
    new_option = runs_menu_generator(new)
    if mode_selection.active ==0:
        a=[('new','new')]
        a.extend(new_option)
        sd_run_list.options = a
    else:
        sd_run_list.options = new_option
    sd_run_list.value = []
    sync_vd_exp_info_widgets(new)


def sd_create_new_exp_cb():
    if sd_experiment.value == 'None':
        info_box.text = info_deque('Enter new exp name.')
        raise ValueError('enter new')
    new_exp = raw_data.next_index(1)
    index_update = {new_exp[0]:{'name':sd_experiment.value,'date':sd_date.value,'author':sd_author.value}}
    raw_data.index.update(index_update)
    raw_data.experiment.update({new_exp[0]:{}})
    raw_data.experiment_raw.update({new_exp[0]:{}})
    raw_data.experiment_to_save.update({new_exp[0]:'upload'})
    new_exp_opt = experiment_menu_generator(raw_data.index.keys())
    sd_experiment_list.options=new_exp_opt
    sd_experiment_list.value = [new_exp_opt[0][0]]
    # raw_data.index_to_save.add(new_exp[0])



def load_csv_to_experiment(run_list,index):
    temp_result_dict={}
    temp_raw_dict={}
    run_speed = float(sd_run_speed.value)
    ext_coef = float(sd_ext_coef.value)
    for i in run_list:
        file_name = os.path.basename(i)
        data_date,data_name,data_key = file_name.split('.')[0].split('_')
        with open(i,'rb') as f:
            data = pd.read_csv(BytesIO(f.read().decode('utf-16').encode('utf-8')),names=['time','signal'],sep='\t')
        data_range = (float(data.iloc[0,0]),float(data.iloc[-1,0]),len(data['time']))
        data_signal = list(data['signal'])
        temp_key = data_date+data_name
        if temp_key in temp_result_dict.keys():
            temp_result_dict[temp_key].update({data_key:{}})
            if data_key[0]=='A':
                temp_result_dict[temp_key][data_key].update(extcoef=ext_coef)
            temp_raw_dict[temp_key].update({data_key:{'time':data_range,'signal':data_signal}})
        else:
            temp_result_dict.update({temp_key:{'date':data_date,'name':data_name,'speed':run_speed,data_key:{}}})
            temp_raw_dict.update({temp_key:{data_key:{'time':data_range,'signal':data_signal}}})
    run_number = len(temp_result_dict.keys())
    run_entry = raw_data.next_run(index,run_number)
    for i,j in zip(run_entry,temp_result_dict.keys()):
        raw_data.experiment[index].update({i:temp_result_dict[j]})
        raw_data.experiment_raw[index].update({i:temp_raw_dict[j]})
        info_box.text = info_deque('Run {} uploaded to experiment {}.'.format(run_entry,index))
    raw_data.experiment_to_save.update({index:'upload'})
    sd_run_list.options = runs_menu_generator([index])


def load_folder_to_experiment(exp_list):
    new_exp_entry = raw_data.next_index(len(exp_list))
    for i,j in zip(exp_list,new_exp_entry):
        exp_date,exp_name = os.path.basename(i).split('_')
        raw_data.index.update({j:{'date':exp_date,'name':exp_name,'author':sd_author.value}})
        raw_data.experiment.update({j:{}})
        raw_data.experiment_raw.update({j:{}})
        run_list = glob.glob(os.path.join(i,'*.[cC][sS][vV]'))
        load_csv_to_experiment(run_list,j)
    sd_experiment_list.options=experiment_menu_generator(raw_data.index.keys())

def sd_folder_upload_cb():
    """
    call back function for upload data from folder button
    """
    search=os.path.join(upload_data_folder,'[1-2][0][0-9]*')
    file_list = glob.glob(search)
    run_list = []
    exp_list = []
    for i in file_list:
        if os.path.isdir(i):
            exp_list.append(i)
        else:
            run_list.append(i)
    current_index = sd_experiment_list.value
    if run_list:
        if len(current_index)!=1:
            info_box.text = info_deque('Select 1 experiment to upload.')
            raise ValueError('sele 1')
        load_csv_to_experiment(run_list,current_index[0])
    else:
        info_box.text = info_deque('No CSV file found.')
    if exp_list:
        load_folder_to_experiment(exp_list)
    else:
        info_box.text = info_deque('No experiment folder found.')


def vd_delete_data_button_cb():
    if vd_selection_tab.active!=0:
        info_box.text = info_deque('Go to Experiment Tab')
        raise ValueError('sele')
    exp_to_del = sd_experiment_list.value
    run_to_del = sd_run_list.value
    if run_to_del:
        for i in run_to_del:
            exp = i.split('-')[0]
            _ = raw_data.experiment[exp].pop(i,None)
            _ = raw_data.experiment_raw[exp].pop(i,None)
            raw_data.experiment_to_save.update({exp:'upload'})
            sd_run_list.options = runs_menu_generator(exp_to_del)
            sd_run_list.value = []
    else:
        for i in exp_to_del:
            _= raw_data.index.pop(i,None)
            _= raw_data.experiment.pop(i,None)
            _= raw_data.experiment_raw.pop(i,None)
            raw_data.experiment_to_save.update({i:'del'})
            sd_experiment_list.options = experiment_menu_generator(raw_data.index.keys())
            sd_experiment_list.value = []

def load_button_cb():
    global raw_data
    with shelve.open(os.path.join(file_path,file_name),writeback=False) as hd:
        data_index = hd['index']
    raw_data = Data(data_index)
    info_box.text = info_deque('Data re-loaded.')
    sd_experiment_list.options = experiment_menu_generator(raw_data.index.keys())
    sd_experiment_list.value = []


# add callbacks
mode_selection.on_change('active',mode_selection_cb)
upload_file_source.on_change('data',upload_file_source_cb)
sd_experiment_list.on_change('value',sd_experiment_list_cb)
sd_create_new_exp.on_click(sd_create_new_exp_cb)
sd_save_data.on_click(sd_save_data_cb)
sd_folder_upload.on_click(sd_folder_upload_cb)
vd_delete_data_button.on_click(vd_delete_data_button_cb)
load_button.on_click(load_button_cb)
save_button.on_click(sd_save_data_cb)
sd_run_list.on_change('value',sd_run_list_cb)
vd_search_botton.on_click(vd_search_botton_cb)
vd_anotation_option.on_change('value',vd_plot_option_cb)
vd_anotation_option_s.on_change('value',vd_plot_option_cb)
vd_plot_backend.on_change('value',vd_plot_option_cb)
vd_offset_option.on_change('value',vd_plot_option_cb)
vd_primary_axis_option.on_change('value',vd_plot_option_cb)
vd_secondary_axis_option.on_change('value',vd_plot_option_cb)
vd_secondary_axis_range.on_change('value',vd_plot_option_cb)





############login layout
plot_login = Plot(title=None, plot_width=600, plot_height=201,
                  min_border=0, toolbar_location=None)
glyph_0_ = ColumnDataSource(dict(x=[0], y=[0], text=['']))
glyph_0 = Text(x='x', y='y', text='text')
glyph_00_ = ColumnDataSource(dict(x=[6], y=[3.1], text=['']))
glyph_00 = Text(x=6, y=4, text='text')
glyph_1_ = ColumnDataSource(
    dict(x=[0], y=[2.5], text=['Aptitude Medical Systems']))
glyph_1 = Text(x='x', y='y', text='text', text_color="orangered",
               text_font_size='33pt', text_font_style='bold', text_font='helvetica')
glyph_2_ = ColumnDataSource(dict(x=[0.84,2.245,3], y=[0,0.99,0.99], text=['PL@j@-nior','~','~']))
glyph_2 = Text(x='x', y='y', text='text', text_color="darkturquoise",
               text_font_size='60pt', text_font='cursive',text_font_style='bold')
plot_login.add_glyph(glyph_0_, glyph_0)
plot_login.add_glyph(glyph_00_, glyph_00)
plot_login.add_glyph(glyph_1_, glyph_1)
plot_login.add_glyph(glyph_2_, glyph_2)
login_text = PreText(
    text="Aptitude Medical Systems, Inc.\n\n      LOG IN TO USE\n")
login_user = TextInput(placeholder="username", title="User Name:")
login_pwd = PasswordInput(placeholder="password", title="Password:")
login_btn = Button(label="LOGIN", width=150, button_type='success')
time_text="""<div style="text-align:center;padding:1em 0;"><h3><span style="color:gray;">Current local time in</span>
<br />Santa Barbara, United States</a></h3>
<iframe src="http://free.timeanddate.com/clock/i6nvx52g/n1050/szw110/szh110/hocbbb/hbw6/
cf100/hgr0/fas16/fdi64/mqc000/mqs4/mql20/mqw2/mqd94/mhc000/mhs3/mhl20/mhw2/mhd94/mmc000/mml10/mmw1/mmd94/hmr7/hsc000/hss1/hsl90"
frameborder="0" width="110" height="110"></iframe>
</div>"""
login_info = Div(text=time_text)

botton_spacer = PreText(text=('ho'), width=1200, height=100)


def login_btn_callback():
    global user_pwd
    user = login_user.value
    password = login_pwd.value
    if password in user_pwd.get(user, ['aptitude','ams']):
        mode_selection.active = 1
        sd_author.value = user
    else:
        login_text.text = 'Wrong Username/Password \n'

def login_pwd_cb(attr,old,new):
    login_btn_callback()

def refresh_time_cb():
    ltime = datetime.datetime.now()
    botton_spacer.text = ('\n'*2+' '*83+'Today is : '+ltime.strftime('%Y-%m-%d, %a')+','
                          + '\tCurrent time: '+ltime.strftime('%I:%M:%S %p')+"\n\n"+' '*107+'Aptitude Medical Systems, Inc.')



login_btn.on_click(login_btn_callback)
login_pwd.on_change('value',login_pwd_cb)


#layouts

plot_row = row(view_data_plot)



analyze_layout = layout([top_row_],[plot_row],[it_tool_box], [botton_spacer])#[analysis_too_box]



vd_layout = layout([top_row_],[plot_row],[column(vd_selection_tab,vd_div_5,vd_search_box),vd_div_1,vd_info_widgets],[botton_spacer])

upload_layout = layout([top_row_],[sd_row_1],[sd_experiment_list,sd_run_list],[sd_bottom_row], [botton_spacer])

display_layout = layout([plot_login], [login_info, column(
    login_text, login_user, login_pwd, login_btn)])


curdoc().add_periodic_callback(refresh_time_cb, 1000)
curdoc().add_root(display_layout)
