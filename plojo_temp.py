from bokeh.plotting import figure, ColumnDataSource
from bokeh.io import curdoc
from bokeh.models import HoverTool, Plot
from bokeh.models.glyphs import Text
from bokeh.models.widgets import Panel, Tabs, Button, TextInput, Select, MultiSelect, RadioButtonGroup, PasswordInput, PreText, DataTable, TableColumn, TextAreaInput
from bokeh.layouts import widgetbox, row, column, layout
from bokeh.palettes import Category10
import json
import numpy as np
import fit_binding as fb
import datetime
from os import path




# common variables:
global info_change_keep, raw_data, file_save_location, file_name, current_time, info_deque_holder
info_deque_holder = ['Welcome!']*3
info_change_keep = dict.fromkeys(
    ['name','author', 'date','tag', 'note', 'fit_method', 'assay_type'], False)

current_time = datetime.datetime.now()


file_save_location = '/Users/hui/Documents/PycharmProjects/plojo/data_storage'
file_name = 'data.json'


with open(path.join(file_save_location, file_name), 'rt') as f:
    raw_data = json.load(f)


assay_type_options = ['beads_kd', 'beads_ic50', 'beads_ric50',
                      'elisa_kd', 'elisa_ic50', 'elisa_ric50', 'N/A']
search_field_options = [('name','Experiment Name'),('date', 'Experiment Date'), ('tag', 'Tag'), ('note', 'Note'), ('author', 'Author'),
                        ('fit_method', 'Fit Method')]


plot_ = figure(plot_width=600, plot_height=400)
plot_.annulus(x=[1, 2, 3], y=[1, 2, 3], color="hotpink",
              inner_radius=0.2, outer_radius=0.5)
norm_plot = figure(plot_width=600, plot_height=400)
norm_plot.annulus(x=[1, 2, 3], y=[1, 2, 3], color="orangered",
              inner_radius=0.2, outer_radius=0.5)


# functions

def info_deque(text):
    global info_deque_holder
    j = len(info_deque_holder)
    info_deque_holder.append(str(j)+'>>'+text)
    result = '\n'.join(info_deque_holder[-3:])
    return result

def plot_generator(source=[], **kwargs):
    global raw_data
    color_ = Category10[10]
    tools_list = "pan,wheel_zoom,save,crosshair,reset"
    p = figure(x_axis_type='log', tools=tools_list)
    p.xaxis.axis_label = 'Concentration /nM'
    p.yaxis.axis_label = 'Raw Signal A.U.'
    p.title.text = 'Raw Data'
    p.title_location = 'above'
    pn = figure(x_axis_type='log', x_range=p.x_range, tools=tools_list)
    pn.xaxis.axis_label = 'Concentration /nM'
    pn.yaxis.axis_label = 'Normalized Signal A.U.'
    pn.title.text = 'Normalized Data '
    pn.title_location = 'above'
    fit_method_count=[0,0]
    for i, j in zip(source, color_):
        data = raw_data[i]
        name=raw_data[i].get('name','No_Name')
        fit_result = fb.generate_cds(data, **kwargs)
        if fit_result['fit_method']=='kd':
            fit_method_count[0] +=1
        else:
            fit_method_count[1] +=1
        raw_data[i].update({'fit_para': fit_result['fit_para']})
        p.circle(
            'x', 'y', source=fit_result['raw_data'], color=j, line_width=2, legend=name)
        p_data_hover = HoverTool(
            tooltips=[('x: ', '@x{0.00}'), ('y: ', '@y{0,0}')])
        p.add_tools(p_data_hover)
        p_temp = p.line(
            'x', 'y', source=fit_result['fit_data'], color=j, line_width=2, legend=name)
        p_hover = HoverTool(renderers=[p_temp], tooltips=fit_result['tooltip'])
        p.add_tools(p_hover)
        pn.circle(
            'x', 'yn', source=fit_result['raw_data'], color=j, line_width=2, legend=name)
        pn_temp = pn.line(
            'x', 'yn', source=fit_result['fit_data'], color=j, line_width=2, legend=name)
        pn_hover = HoverTool(
            renderers=[pn_temp], tooltips=fit_result['tooltip'])
        pn.add_tools(pn_hover)
        p.cross('concentration', 'signal',
                source=fit_result['outlier'], angle=80, size=20, line_width=2, color='red')
    if kwargs.get('mark_outlier', False):
        p.cross('x', 'y', source=kwargs['mark_outlier'],
                angle=80, size=20, line_width=2, color='red')
    if fit_method_count[0]==0:
        p.legend.location = 'top_right'
        pn.legend.location = 'top_right'
    elif fit_method_count[1]==0:
        p.legend.location = 'top_left'
        pn.legend.location = 'top_left'
    elif fit_method_count[0]>fit_method_count[1]:
        p.legend.location = 'center_left'
        pn.legend.location = 'center_left'
    else:
        p.legend.location = 'center_right'
        pn.legend.location = 'center_right'
    p.legend.click_policy = 'hide'
    p.plot_height = 400
    p.plot_width = 600
    pn.legend.click_policy = 'hide'
    pn.plot_height = 400
    pn.plot_width = 600
    return p, pn


def fit_method_para_reader(fit_method, selected_items=[]):
    if fit_method == 'kd':
        cf_fit_bound.children = [cf_fit_method,
                                 cf_kd_bound, cf_Fmax_bound, cf_Fmin_bound]
        cf_fit_method.value = 'kd'
        cf_kd_bound.value, cf_Fmax_bound.value, cf_Fmin_bound.value = bounds_formatter(
            selected_items, fit_method)
    elif fit_method == 'ic_50':
        cf_fit_bound.children = [cf_fit_method,
                                 cf_ic_50_bound, cf_Fmax_bound, cf_Fmin_bound]
        cf_fit_method.value = 'ic_50'
        cf_ic_50_bound.value, cf_Fmax_bound.value, cf_Fmin_bound.value = bounds_formatter(
            selected_items, fit_method)
    elif fit_method == 'ric_50':
        cf_fit_bound.children = [cf_fit_method, cf_r0_bound, cf_v0_bound,
                                 cf_kd_r_bound, cf_kd_a_bound, cf_Fmax_bound, cf_Fmin_bound]
        cf_fit_method.value = 'ric_50'
        cf_r0_bound.value, cf_v0_bound.value, cf_kd_r_bound.value, cf_kd_a_bound.value, cf_Fmax_bound.value, cf_Fmin_bound.value = bounds_formatter(
            selected_items, fit_method)
    else:
        cf_fit_bound.children = [cf_fit_method]
        cf_fit_method.value = 'various'


def bounds_formatter(sele, fit_method):
    """
    format bounds of a selection of raw data to string foramt.
    """

    if len(set([raw_data[i].get('fit_method', 'kd') for i in sele])) > 1:
        # fit_method = raw_data[sele[0]].get('fit_method','kd')
        info_box.text = info_deque('Selection not unique.')
    bounds = [raw_data[i].get('bounds', {}) for i in sele]

    def list_formatter(bounds, tag):
        r = ','.join([str(i.get(tag)[0])+'-'+str(i.get(tag)[1]) if i.get(tag,
                                                                         'default') != 'default' else 'default' for i in bounds])
        return r
    kd = list_formatter(bounds, 'kd')
    Fmax = list_formatter(bounds, 'Fmax')
    Fmin = list_formatter(bounds, 'Fmin')
    ic_50 = list_formatter(bounds, 'ic_50')
    r_0 = list_formatter(bounds, 'r_0')
    v_0 = list_formatter(bounds, 'v_0')
    kd_r = list_formatter(bounds, 'kd_r')
    kd_a = list_formatter(bounds, 'kd_a')
    if fit_method == 'kd':
        result = (kd, Fmax, Fmin)
    elif fit_method == 'ic_50':
        result = (ic_50, Fmax, Fmin)
    elif fit_method == 'ric_50':
        result = (r_0, v_0, kd_r, kd_a, Fmax, Fmin)
    else:
        pass
    return result




def bounds_generator():
    """
    read bounds input and return list of bounds.
    {'ams0':{'kd':[1,10],'Fmin':[1,10]}}
    """
    def parser(string):
        list_string = string.split(',')
        for idx, text in enumerate(list_string):
            if '-' in text:
                try:
                    list_string[idx] = [float(i) if i != '' else [
                                              0, np.inf][j] for j, i in enumerate(text.split('-'))]
                    assert len(list_string[idx]) == 2
                except:
                    info_box.text = info_deque(
                        'Bound is not a range, use default instead.')
                    raise NameError('no valid input')
            elif text == 'default' or text == '':
                list_string[idx] = 'default'
            else:
                info_box.text = info_deque('Bound is not a range,try again.')
                raise NameError('not valid input')
        return list_string

    def dict_getter(sele_list, valuelist, namelist):
        sele_length = len(sele_list)
        for idx,i in enumerate(valuelist):
            if len(i)==sele_length:
                pass
            elif len(i)<sele_length:
                info_box.text = info_deque('bounds too short, fill with last element.')
                valuelist[idx].extend([i[-1] for j in range(sele_length-len(i))])
            else:
                info_box.text = info_deque('bounds too long,omit rest.')
                valuelist[idx]=valuelist[idx][0:sele_length]
        a = [dict(zip(namelist, [j[i_]for j in valuelist]))
             for i_ in range(sele_length)]
        return dict(zip(sele_list, a))
    selected_items = cf_focused_select_data.value
    fit_method = cf_fit_method.value
    kd = parser(cf_kd_bound.value)
    ic_50 = parser(cf_ic_50_bound.value)
    r0 = parser(cf_r0_bound.value)
    v0 = parser(cf_v0_bound.value)
    kdr = parser(cf_kd_r_bound.value)
    kda = parser(cf_kd_a_bound.value)
    Fmax = parser(cf_Fmax_bound.value)
    Fmin = parser(cf_Fmin_bound.value)
    if fit_method == 'kd':
        dict_of_bounds = dict_getter(
            selected_items, [kd, Fmax, Fmin], ['kd', 'Fmax', 'Fmin'])
    elif fit_method == 'ic_50':
        dict_of_bounds = dict_getter(selected_items, [ic_50, Fmax, Fmin], [
                                     'ic_50', 'Fmax', 'Fmin'])
    elif fit_method == 'ric_50':
        dict_of_bounds = dict_getter(selected_items, [r0, v0, kdr, kda, Fmax, Fmin], [
                                     'r_0', 'v_0', 'kd_r', 'kd_a', 'Fmax', 'Fmin'])
    else:
        pass
    return dict_of_bounds


def outlier_generator():
    """
    read the value of selected points on plot, return dict.
    """
    outlier = cf_outlier.value
    conc = [float(i.split('->')[0]) for i in outlier]
    sig = [float(i.split('->')[1]) for i in outlier]
    outlier_dict = {}
    outlier_dict[cf_focused_select_data.value[0]] = {
        'concentration': conc, 'signal': sig}
    return outlier_dict


def menu_generator(items):
    menu = []
    for i in items:
        menu.append(' '.join([i, raw_data[i].get(
            'name', 'No_Name'), raw_data[i].get('tag', 'No_Tag'), raw_data[i].get('date', 'N/A'),raw_data[i].get('assay_type', 'N/A')]))
    result = list(zip(items, menu))
    result = sorted(result, key=lambda x: int(
        x[0].split('-')[0][3:]), reverse=True)
    return result

def save_data():
    global raw_data
    with open(path.join(file_save_location, file_name), 'rt') as f:
        old_data = json.load(f)
    if raw_data == old_data:
        info_box.text = info_deque('No change was made.')
        pass
    elif raw_data:
        now = datetime.datetime.now()
        with open(path.join(file_save_location, file_name), 'wt') as f:
            json.dump(raw_data, f, indent=4)
            print('read')
        with open(path.join(file_save_location, 'temp', (file_name[:-5]+'_'+now.strftime('%Y%m%d%H%M%S')+'.json')), 'wt') as f:
            json.dump(old_data, f, indent=4)
            print('save')
        info_box.text = info_deque('Data saved.')
    else:
        info_box.text = info_deque('Load Data Before Saving!')



def load_data():
    print('-1')
    info_box.text = info_deque('Start Loading data...')
    print('0')
    global raw_data
    with open(path.join(file_save_location, file_name), 'rt') as f:
        raw_data = json.load(f)
    print('1')
    cf_select_data.options = menu_generator(list(raw_data.keys()))
    print('2')
    cf_select_data.value = []
    print('2')
    info_box.text = info_deque('Data Loaded.')

def sync_vd_info_widgets(selected_items):
    vd_author.value = str([raw_data[i].get('author', 'No_Author')
                           for i in selected_items]).strip('[]')
    vd_date.value = str([raw_data[i].get('date', 'No_Date')
                         for i in selected_items]).strip('[]')
    vd_name.value = str([raw_data[i].get('name', 'No_Name')
                           for i in selected_items]).strip('[]')
    vd_tag.value = str([raw_data[i].get('tag', 'No_Tag') for i in selected_items]).strip('[]')
    vd_note.value = str([raw_data[i].get('note', 'No_Note')
                         for i in selected_items]).strip('[]')
    vd_fit_method.value = str(
        ([raw_data[i].get('fit_method', 'N/A') for i in selected_items])).strip('[]')

    if len(set([raw_data[i].get('assay_type', 'N/A') for i in selected_items])) == 1:
        assay_type = raw_data[selected_items[0]].get('assay_type', 'N/A')
    else:
        assay_type = 'N/A'
    vd_assay_type.value = assay_type
    if len(set([raw_data[i].get('fit_method', 'kd') for i in selected_items])) == 1:
        fit_method_para_reader(raw_data[selected_items[0]].get(
            'fit_method', 'kd'), selected_items)
    else:
        fit_method_para_reader('various')
    global info_change_keep
    info_change_keep = dict.fromkeys(info_change_keep, False)
    if len(selected_items) == 1:
        conc = raw_data[selected_items[0]]['concentration']
        sig = raw_data[selected_items[0]]['signal']
        cf_outlier.options = [str(i)+' -> '+str(j) for i, j in zip(conc, sig)]
        cf_outlier.value = []
    else:
        cf_outlier.options = []
        cf_outlier.value = []



def syn_plot_display(selected_items, **kwargs):
    raw_plot, norm_plot = plot_generator(selected_items, **kwargs)
    read_layout.children[1].children[0].children[0] = raw_plot
    read_layout.children[1].children[0].children[1] = norm_plot
    curve_fit_layout.children[1].children[0].children[0] = raw_plot
    curve_fit_layout.children[1].children[0].children[1] = norm_plot
    te = []
    for i in selected_items:
        value = raw_data[i].get('fit_para', {})
        tag = raw_data[i].get('tag', 'N/A')
        t = i+'-'+tag+' -> '
        for key, item in value.items():
            if key in ['kd', 'ic_50', 'r_0', 'v_0', 'kd_a', 'kd_r']:
                temp = '{}:{:.3f}nM '.format(key, item)
            else:
                temp = '{}:{:.1f} '.format(key, item)
            t = t+temp
        te.append(t)
    fit_para_result = '\n'.join(te)
    cf_parameter_viewer.text = fit_para_result


def input_to_rawdata():
    entry = list(raw_data.keys())
    entry = sorted(entry, key=lambda x: int(
        x.split('-')[0][3:]), reverse=True)[0]
    entry_start = int(entry.split('-')[0][3:])
    author = sd_author.value
    date = sd_date.value
    assay_type = sd_assay_type.value

    if assay_type == 'N/A':
        fit_method = 'kd'
    else:
        fit_method = assay_type.split('_')[1]
    data = sd_data.value
    data = data.strip('\n')
    print(repr(data))
    data = [i.split('\t') for i in data.split('\n')]

    save_entry_list = [
        'ams'+str(i) for i in range(entry_start+1, entry_start+len(data[0]))]

    data_re = list(map(list, zip(*data)))
    if data_re[0][0:3] != ['tag', 'note', 'nM_Name']:
        info_box.text = info_deque(
            'Wrong Order: tag, note, nM_Name.')
        raise ValueError('input format')
    concentration = data_re[0][3:]
    result_dict = {}
    for key, value in zip(save_entry_list, data_re[1:]):
        dict_to_save = {}
        dict_to_save.update(tag=value[0], note=value[1],name=value[2],author=author, date=date, assay_type=assay_type, fit_method=fit_method)
        conc_ = []
        signal = []
        for i, j in zip(concentration, value[3:]):
            if i == '' or j == '':
                pass
            else:
                try:
                    conc_.append(float(i))
                    signal.append(float(j))
                except:
                    info_box.text = info_deque('Error: concentration or signal contain non-number values.')
        if len(conc_) == len(signal) and len(signal) > 3:
            dict_to_save.update(concentration=conc_, signal=signal)
        else:
            info_box.text = info_deque('Error: value format not correct or too few data points')
            raise ValueError('not length')
        result_dict[key] = dict_to_save
    print(result_dict)
    info_box.text = info_deque('Format check passed.')
    return result_dict, save_entry_list


# widgets


button_mode = RadioButtonGroup(labels=[
                               'Save Changes', 'View Data', 'Curve Fitting'], active=0, button_type='warning', width=300)
info_box = PreText(text='Welcome!')
button_load = Button(label='Load Data', button_type='success')
button_save = Button(label='Save Data', button_type='danger')
botton_spacer = PreText(text=('ho'), width=1200, height=100)
cf_assay_type = MultiSelect(
    title='Assay Type', value=assay_type_options, options=assay_type_options, size=5)
cf_filter = TextInput(title='Keyword filter')
cf_search_field = MultiSelect(title='Search field:', value=[
                              'tag', 'note', 'fit_method', 'fix_comp'], options=search_field_options, size=7)
cf_select_data = MultiSelect(title='Select Experiment Data To View',
                             options=menu_generator(raw_data.keys()), size=10, width=600)
cf_outlier = MultiSelect(
    title='Mark Data Points As Outlier', options=[], size=10, width=600)
cf_parameter_viewer = PreText(text='Fitting Result:\n', width=600, height=160)
cf_focused_select_data = MultiSelect(
    title='Data of interest', options=[], size=10, width=600)
data_outlier_tab = Tabs(active=0, width=600, height=230, tabs=[Panel(child=cf_select_data, title="Experiment"), Panel(
    child=cf_focused_select_data, title='To Plot'), Panel(child=cf_outlier, title="Raw Data"), Panel(child=cf_parameter_viewer, title='Fitting Result')])
vd_save_info = Button(label='Save Info Changes', button_type='success')
vd_delete_data = Button(label='Delete Selected Data', button_type='danger')
options = widgetbox(cf_assay_type, cf_filter, cf_search_field,
                    vd_save_info, vd_delete_data, sizing_mode='stretch_both')
vd_name = TextInput(title='Caution! Edits to all selected items.\nExperiment Name')
vd_author = TextInput(title='Author')
vd_date = TextInput(title='Experiment Date')
vd_tag = TextAreaInput(title='Tag:', rows=3, cols=35, max_length=50000)
vd_note = TextAreaInput(title='Note', rows=5, cols=35, max_length=50000)
vd_fit_method = TextInput(title='Fit Method, edit here won\'t have effect.')
vd_assay_type = Select(title='Assay Type', value='N/A', options=[
                       'beads_kd', 'beads_ic50', 'beads_ric50', 'elisa_kd', 'elisa_ic50', 'elisa_ric50', 'N/A'])
vd_data_info = widgetbox(vd_name,vd_assay_type,vd_tag, vd_note,vd_author, vd_date, vd_fit_method)


######## curve_fitting layout

cf_alias_name = TextInput(title='Make alias for selected data:', value='copy')
cf_make_alias = Button(label='Creat Alias', button_type='primary')
cf_mark_outlier = Button(label='Mark Outlier', button_type='success')
cf_select_data_widgets = widgetbox(
    cf_assay_type, cf_filter, cf_search_field, cf_alias_name, cf_make_alias, cf_mark_outlier)
cf_fit_method = Select(title='Fit method:', value='kd',
                       options=fb.supported_fit_method)
cf_kd_bound = TextInput(title='Kd range in nM:', value='default')
cf_ic_50_bound = TextInput(title='IC50 range in nM:', value='default')
cf_Fmax_bound = TextInput(title='Max Signal AU:', value='default')
cf_Fmin_bound = TextInput(title='Min Signal AU:', value='default')
cf_r0_bound = TextInput(title='Receptor Conc. in nM:', value='default')
cf_v0_bound = TextInput(title='VEGF Conc. in nM:', value='default')
cf_kd_r_bound = TextInput(title='VEGF-Receptor Kd in nM:', value='default')
cf_kd_a_bound = TextInput(title='VEGF-Aptamer Kd in nM:', value='default')
cf_fit_bound = widgetbox(cf_fit_method, cf_kd_bound,
                         cf_Fmax_bound, cf_Fmin_bound)
cf_fit_plot_data = Button(label='Fit and Plot', button_type='success')
cf_restore = Button(label='Reset Fitting Parameter', button_type='warning')

#########save layout

sd_author = TextInput(title='Author', value='none')
sd_date = TextInput(title='Experiment Date',
                    value=current_time.strftime('%Y%m%d'))
sd_assay_type = Select(title='Assay Type', value='beads_kd', options=[
                       'beads_kd', 'beads_ic50', 'beads_ric50', 'elisa_kd', 'elisa_ic50', 'elisa_ric50', 'N/A'])
sd_row_1 = row(sd_author, sd_date, sd_assay_type)
sd_data = TextAreaInput(title='Experiment Data:',
                        rows=30, cols=140, max_length=500000)
sd_type_check = Button(label='Proof Read Data', button_type='success')
sd_save_data = Button(label='Save Input Data', button_type='danger')
source = ColumnDataSource({})
sd_data_table = DataTable(width=600, height=280, source=source)


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
glyph_2_ = ColumnDataSource(dict(x=[0.94], y=[0.1], text=['p l o j o']))
glyph_2 = Text(x='x', y='y', text='text', text_color="darkturquoise",
               text_font_size='70pt', text_font='chalkboard',text_font_style='bold')
plot_login.add_glyph(glyph_0_, glyph_0)
plot_login.add_glyph(glyph_00_, glyph_00)
plot_login.add_glyph(glyph_1_, glyph_1)
plot_login.add_glyph(glyph_2_, glyph_2)
login_text = PreText(
    text="Aptitude Medical Systems, Inc.\n\n      LOG IN TO USE\n")
login_user = TextInput(placeholder="username", title="User Name:")
login_pwd = PasswordInput(placeholder="password", title="Password:")
login_btn = Button(label="LOGIN", width=150, button_type='success')
login_info = PreText(text=('Today is : '+current_time.strftime('%Y-%m-%d, %a')
                           + '\nCurrent time: '+current_time.strftime('%I:%M:%S %p')))


# call back functions


def mode_callback(attr, old, new):
    value = button_mode.active
    if value == 1:
        display_layout.children = read_layout.children
    elif value == 0:
        display_layout.children = save_layout.children
    else:
        display_layout.children = curve_fit_layout.children


def cf_select_data_cb(attr, old, new):
    selected_items = cf_select_data.value
    # if len(selected_items) == 1:
    #     conc = raw_data[selected_items[0]]['concentration']
    #     sig = raw_data[selected_items[0]]['signal']
    #     cf_outlier.options = [str(i)+' -> '+str(j) for i, j in zip(conc, sig)]
    #     cf_outlier.value = []
    # else:
    #     cf_outlier.options = []
    #     cf_outlier.value = []
    cf_focused_select_data.options = menu_generator(selected_items)
    # cf_focused_select_data.value = []
    sync_vd_info_widgets(selected_items)


def update_select_menu(attr, old, new):
    assay_type = cf_assay_type.value
    search_field = cf_search_field.value
    keyword = cf_filter.value
    search_hit = []
    for key, item in raw_data.items():
        item_type = item.get('assay_type', 'N/A')
        if item_type in assay_type:
            for field_ in search_field:
                if keyword.lower() in item.get(field_, '').lower():
                    search_hit.append(key)
                    break
    menu = menu_generator(search_hit)
    cf_select_data.options = menu


def cf_focused_select_data_cb(attr, old, new):
    info_box.text = info_deque('Drawing...')
    selected_items = cf_focused_select_data.value
    sync_vd_info_widgets(selected_items)
    syn_plot_display(selected_items)
    info_box.text = info_deque('New plot generated.')


def vd_delete_data_cb():
    global current_time
    current = datetime.datetime.now()
    if (current-current_time).total_seconds() > 2:
        current_time = datetime.datetime.now()
        info_box.text = info_deque(
            'Caution! Click again within 2s to delete data.')
        raise ValueError('click again')
    to_delete = cf_focused_select_data.value
    if len(to_delete) == 0:
        info_box.text = info_deque('Select Data in RoI tab to delete')
        raise ValueError('no select')
    for i in to_delete:
        del raw_data[i]
    info_box.text = info_deque('Selected data was deleted.')
    cf_select_data_menu = cf_select_data.options
    selected_data = cf_select_data.value
    cf_select_data.options = [
        i for i in cf_select_data_menu if i[0] not in to_delete]
    cf_select_data.value = [i for i in selected_data if i not in to_delete]


def vd_author_callback(attr, old, new):
    info_change_keep.update(author=True)

def vd_name_callback(attr, old, new):
    info_change_keep.update(name=True)


def vd_date_callback(attr, old, new):
    info_change_keep.update(date=True)


def vd_tag_cb(attr, old, new):
    info_change_keep.update(tag=True)


def vd_note_cb(attr, old, new):
    info_change_keep.update(note=True)

#
# def vd_fit_method_cb(attr, old, new):
#     info_change_keep.update(fit_method=True)


def vd_assay_type_cb(attr, old, new):
    info_change_keep.update(assay_type=True)


def vd_save_info_callback():
    if data_outlier_tab.active ==0:
        selected_items = cf_select_data.value
    elif data_outlier_tab.active ==1:
        selected_items = cf_focused_select_data.value
    else:
        info_box.text = info_deque('Go to Experiment/DOI tab to change info.')
        raise ValueError(' go to experiment tab')
    global info_change_keep
    info_box.text = info_deque('saving info...')
    new_input = dict(zip(['name','author', 'date', 'tag', 'note', 'fit_method', 'assay_type'], [vd_name.value,
                     vd_author.value, vd_date.value, vd_tag.value, vd_note.value, vd_fit_method.value, vd_assay_type.value]))
    update_dict = {i: new_input[i] for i, j in info_change_keep.items() if j}
    for i in selected_items:
        raw_data[i].update(update_dict)
    info_change_keep = dict.fromkeys(info_change_keep, False)
    saved_field = list(update_dict.keys())
    cf_select_data.options = menu_generator([i[0] for i in cf_select_data.options])
    cf_focused_select_data.options = menu_generator([i[0] for i in cf_focused_select_data.options])
    info_box.text = info_deque(str(saved_field).strip('[]')+' saved.')

def cf_restore_callback():
    selected_items = cf_focused_select_data.value
    if len(selected_items) != 1:
        info_box.text = info_deque('Select single experiment to restore.')
        raise ValueError('select 1 item')
    raw_data[selected_items[0]].update(fit_para={}, bounds={}, outlier={
                                       'concentration': [], 'signal': []})
    cf_outlier.value = []
    info_box.text = info_deque('Fitting parameters have been reset.')


def cf_fit_method_callback(attr, old, new):
    choice = cf_fit_method.value
    selected_items = cf_focused_select_data.value
    fit_method_para_reader(choice, selected_items)


def cf_fit_plot_callback():
    info_box.text = info_deque('start fitting...')
    selected_items = cf_focused_select_data.value
    fit_method = cf_fit_method.value
    bounds = bounds_generator()
    info_box.text = info_deque('bounds data loaded...')
    outlier = outlier_generator()
    info_box.text = info_deque('outlier data loaded...')
    for i in selected_items:
        raw_data[i].update(fit_method=fit_method,
                           bounds=bounds.get(i, {}), fit_para={}, outlier=outlier.get(i, {'concentration': [], 'signal': []}))
    syn_plot_display(selected_items, new_fit=True)
    info_box.text = info_deque('Fitting complted.')


def make_alias_callback():
    info_box.text = info_deque('creating copy...')
    alias = cf_focused_select_data.value[0] + '-' + cf_alias_name.value
    if raw_data.get(alias, None) != None:
        info_box.text = info_deque('This alias already exists.')
        raise KeyError('key exist for this alias.')
    raw_data[alias] = raw_data[cf_focused_select_data.value[0]]
    cf_select_data.options.append(alias)
    cf_select_data.value.append(alias)
    cf_focused_select_data.value = [alias]
    info_box.text = info_deque('Data copy created.')



def cf_mark_outlier_callback():
    """
    mark outlier out
    """
    info_box.text = info_deque('Marking outlier ...')
    sele = cf_focused_select_data.value
    outlier = outlier_generator()
    x = outlier[sele[0]].get('concentration', [])
    y = outlier[sele[0]].get('signal', [])
    source = ColumnDataSource(dict(x=x, y=y))
    syn_plot_display(sele, mark_outlier=source)
    info_box.text = info_deque('Outlier marked.')


def sd_type_check_cb():
    data = sd_data.value
    data = data.strip('\n')
    entry = list(raw_data.keys())
    entry = sorted(entry, key=lambda x: int(
        x.split('-')[0][3:]), reverse=True)[0]
    entry_start = int(entry.split('-')[0][3:])

    data = [i.split('\t') for i in data.split('\n')]
    save_entry_list = [
        'ams'+str(i) for i in range(entry_start+1, entry_start+len(data[0]))]
    data_re = list(map(list, zip(*data)))
    data_field = ['Field']
    data_field.extend(save_entry_list)
    datasource = ColumnDataSource(dict(zip(data_field, data_re)))
    columns = [TableColumn(field=i, title=i, width=90) for i in data_field]
    display_layout.children[3].children[0] = DataTable(
        source=datasource, columns=columns, width=600, height=280, fit_columns=False)
    _, _a = input_to_rawdata()


def sd_save_data_cb():
    info_box.text = info_deque('Start saving...')
    result, save_entry_list = input_to_rawdata()
    result_conc_list = []
    for key, item in result.items():
        result_conc_list.append(item['signal'])
        for i, j in raw_data.items():
            if item['signal'] == j['signal']:
                info_box.text = info_deque(
                    'Library record(s) collision found in input data.')
                raise KeyError('same data as in raw data')
    for idx, conc in enumerate(result_conc_list[:-1]):
        for con_2 in result_conc_list[(idx+1):]:
            if conc == con_2:
                info_box.text = info_deque(
                    'Duplicate data was found in input.')
                raise ValueError('duplicate found')
    raw_data.update(result)
    cf_select_data.options = menu_generator(save_entry_list)
    save_data()
    info_box.text = info_deque('Data has successfully saved.')


def login_btn_callback():
    user_pwd = {'aptitude': 'aptitude', 'newuser': 'aptitude','hui':'h'}
    user = login_user.value
    password = login_pwd.value
    if user_pwd.get(user, 'aptitude') == password:
        button_mode.active = 1
        sd_author.value = user
    else:
        login_text.text = 'Wrong Username/Password \n'


def refresh_time_cb():
    ltime = datetime.datetime.now()
    botton_spacer.text = ('\n'*2+' '*83+'Today is : '+ltime.strftime('%Y-%m-%d, %a')+','
                          + '\tCurrent time: '+ltime.strftime('%I:%M:%S %p')+"\n\n"+' '*107+'Aptitude Medical Systems, Inc.')
    login_info.text = ('Today is : '+ltime.strftime('%Y-%m-%d, %a')
                       + '\nCurrent time: '+ltime.strftime('%I:%M:%S %p'))


# assign callback to elements
button_load.on_click(load_data)
button_save.on_click(save_data)
button_mode.on_change('active', mode_callback)
cf_focused_select_data.on_change('value', cf_focused_select_data_cb)
cf_select_data.on_change('value', cf_select_data_cb)
cf_assay_type.on_change('value', update_select_menu)
cf_filter.on_change('value', update_select_menu)
cf_search_field.on_change('value', update_select_menu)
vd_delete_data.on_click(vd_delete_data_cb)
vd_author.on_change('value', vd_author_callback)
vd_date.on_change('value', vd_date_callback)
vd_name.on_change('value',vd_name_callback)
vd_tag.on_change('value', vd_tag_cb)
vd_note.on_change('value', vd_note_cb)
vd_assay_type.on_change('value', vd_assay_type_cb)
vd_save_info.on_click(vd_save_info_callback)
cf_restore.on_click(cf_restore_callback)
cf_fit_method.on_change('value', cf_fit_method_callback)
cf_fit_plot_data.on_click(cf_fit_plot_callback)
cf_make_alias.on_click(make_alias_callback)
cf_mark_outlier.on_click(cf_mark_outlier_callback)
sd_type_check.on_click(sd_type_check_cb)
sd_save_data.on_click(sd_save_data_cb)
login_btn.on_click(login_btn_callback)

# layouts

save_layout = layout([button_mode, info_box, button_load, button_save], [sd_row_1], [
                     sd_data], [sd_data_table, column(sd_type_check, sd_save_data)], [botton_spacer])

read_layout = layout([button_mode, info_box, button_load, button_save], [column(
    plot_, norm_plot), column(data_outlier_tab, row(options, vd_data_info))], [botton_spacer],)

curve_fit_layout = layout([button_mode, info_box, button_load, button_save], [column(plot_, norm_plot), column(
    data_outlier_tab, row(cf_select_data_widgets, column(cf_fit_bound, cf_fit_plot_data, cf_restore)))], [botton_spacer])

display_layout = layout([plot_login], [login_info, column(
    login_text, login_user, login_pwd, login_btn)])


# run

curdoc().add_periodic_callback(refresh_time_cb, 1000)
curdoc().add_root(display_layout)
