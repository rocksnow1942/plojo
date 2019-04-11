from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import HoverTool,Slider,RangeSlider
from bokeh.models.widgets import Button, TextInput,PreText
from bokeh.layouts import widgetbox
import numpy as np
from plojo_module import plojo_data_init

class Kd_simu():
    def __init__(self):
        para_text = PreText(text='Binding Parameters')
        self.slider_kd = Slider(title='log(Kd) (nM)', start=-3, end=3, step=0.01, value=0)
        self.slider_conc = Slider(title='log(S_0) (nM)', start=-3, end=3, step=0.01, value=0)
        self.slider_ns = Slider(title='Non Specific Binding (AU/nM)', start=0, end=20, step=0.01, value=0)
        self.slider_Fminmax = RangeSlider(title = 'Singal Range AU', start =0, end =100000,step =1000,value =(1000,10000))
        para_text_ = PreText(text='Binding Parameters')
        self.input_kd = TextInput(title='Binding Kd (nM)',value='1')
        self.input_conc = TextInput(title='Fixed Component Conc. (nM)',value='1')
        self.input_ns = TextInput(title='Non Specific Binding (AU/nM)',value='1')
        rand_text = PreText(text='Randomization Parameters')
        self.slider_conc_range =RangeSlider(title= 'Experiment Conc. (log nM) : ',  start = -3, end=3, step=0.01,value=(-2,2))
        self.slider_points = Slider(title='Number of concentration points: ', start = 3, end=20, step=1,value=9)
        self.slider_sets = Slider(title='Number of replicates:',start=1,end=9,step=1,value=1)
        self.slider_linear = Slider(title='Linear randomization', start = 0, end=25, step=0.1, value=0.2)
        self.slider_proportional = Slider(title='Proportional randomization', start = 0, end=25, step=0.1,value=0.2)
        self.refresh_button = Button(label='Refresh random data',button_type='success')
        self.add_button = Button(label='Add data to plojo',button_type='success')
        self.name_input = TextInput(title='Create Name for the data',value='Simu_Kd')
        initial_para = self.fetch_para(False)
        self.fit_data = ColumnDataSource(data=self.generate_data(**self.fetch_para(False)))
        self.raw_data = ColumnDataSource(data=self.generate_data(**self.fetch_para(True)))
        tools_list = "pan,ywheel_zoom,xwheel_zoom,save,reset"
        self.p = figure(x_axis_label='Concentration (nM)', y_axis_label='Signal A.U.', x_axis_type='log',tools=tools_list)
        self.p.title.text = 'Binding Kd = {:.3g} nM; Fixed S_0 = {:.3g} nM'.format(initial_para['kd'],initial_para['s_0'])
        self.p.title_location='above'
        self.p.line('x','y', source=self.fit_data,line_color='green',line_width=2) #legend='Binding Curve')
        self.p.circle('x','y',source=self.raw_data,color='red',line_width=3)
        hover_tool_1 = HoverTool(tooltips=[('Aptamer/nM', '@x{0.00}'),( 'Signal', '@y{0}' )],mode='vline')
        self.p.add_tools(hover_tool_1)
        self.p.plot_height = 400
        self.p.plot_width = 700
        opt_input = widgetbox(para_text_,self.input_kd,self.input_conc,self.input_ns)
        options = widgetbox(para_text,self.slider_kd,self.slider_conc,self.slider_ns,self.slider_Fminmax,self.name_input)
        rand_opt = widgetbox(rand_text,self.slider_conc_range,self.slider_points,self.slider_sets,self.slider_linear,self.slider_proportional,self.refresh_button,self.add_button)
        self.slider_kd.on_change('value',self.callback)
        self.slider_conc.on_change('value',self.callback)
        self.slider_ns.on_change('value',self.callback)
        self.slider_Fminmax.on_change('value',self.callback)
        self.slider_conc_range.on_change('value',self.callback)
        self.slider_points.on_change('value',self.callback)
        self.slider_linear.on_change('value',self.callback)
        self.slider_proportional.on_change('value',self.callback)
        self.slider_sets.on_change('value',self.callback)
        self.refresh_button.on_click(self.refresh_button_cb)
        self.add_button.on_click(self.add_button_cb)
        self.layout = ([self.p],[opt_input,options,rand_opt])

    def callback(self,attr,old,new):
        self.fit_data.data = self.generate_data(**self.fetch_para(False))
        self.raw_data.data = self.generate_data(**self.fetch_para(True))
        para = self.fetch_para(True)
        Fmax = self.slider_Fminmax.value[1]
        kd = 10 ** self.slider_kd.value
        self.input_conc.value = str(10 ** self.slider_conc.value)
        self.input_kd.value = str(kd)
        self.input_ns.value = str(self.slider_ns.value * Fmax/(kd*5000))
        self.p.title.text = self.title_generator(para)

    def generate_data(self,kd=1,s_0=1,ns=0,Fmax=10000,Fmin=1000,start=None,end=None,point=100,randomize =False,linear=0.001,proportional=0.001,sets=1,**kwargs):
        if start == None and end == None:
            i_0 = np.geomspace(kd / 1000, kd * 1000, 100)
        else:
            i_0= np.geomspace(start,end,point)
        if randomize:
            dict_ = {'x':np.repeat(i_0,sets), 'y':self.randomizer(self.solve_binding(np.repeat(i_0,sets), kd, s_0,ns,Fmax,Fmin),linear=linear,proportional=proportional,**kwargs)}
        else:
            dict_ = {'x':i_0, 'y':self.solve_binding(i_0, kd, s_0,ns,Fmax,Fmin)}
        return dict_

    def fetch_para(self,randomize):
        kd = 10 ** self.slider_kd.value
        s_0 = 10 ** self.slider_conc.value
        sets=self.slider_sets.value
        Fmax = self.slider_Fminmax.value[1]
        ns = self.slider_ns.value * Fmax/(kd*5000)
        Fmin = self.slider_Fminmax.value[0]
        start = 10**self.slider_conc_range.value[0]
        end = 10**self.slider_conc_range.value[1]
        points = self.slider_points.value
        linear = self.slider_linear.value/100
        proportional = self.slider_proportional.value/100
        if randomize:
            result = dict(kd=kd,s_0=s_0,ns=ns,Fmax=Fmax,Fmin=Fmin,start=start,end=end,point=points,randomize=True,linear=linear,proportional=proportional,sets=sets)
        else:
            result = dict(kd=kd,s_0=s_0,ns=ns,Fmax=Fmax,Fmin=Fmin)
        return result

    def solve_binding(self,i_0, kd_func, s_0func,ns=0,Fmax=10000,Fmin=1000): # i_0 is np array,
        s_free =((s_0func-kd_func-i_0)+np.sqrt((kd_func+i_0)**2+s_0func**2+2*kd_func*s_0func-2*i_0*s_0func))*0.5
        si = s_0func-s_free
        signal = (Fmax-Fmin)*si/s_0func + Fmin + ns*i_0
        return signal

    def randomizer(self,signal,linear=0.001, proportional=0.001,seed=42):
        np.random.seed(seed)
        size = len(signal)
        max_sig = max(signal)
        return signal * np.random.normal(loc=1,scale=proportional,size=size) + max_sig*np.random.normal(loc=0,scale=linear,size=size)

    def title_generator(self,para_dict):
        kd = para_dict['kd']
        s_0 = para_dict['s_0']
        ns = para_dict['ns']
        Fmax = para_dict['Fmax']
        Fmin = para_dict['Fmin']
        title = 'Kd:{:.3g}nM; S_0:{:.3g}nM; NS:{:.3g}AU/nM; Fmax:{:.3g}, Fmin:{:.3g}'.format(kd, s_0,ns,Fmax,Fmin)
        return title

    def refresh_button_cb(self):
        seed = np.random.randint(90000)
        self.raw_data.data = self.generate_data(**self.fetch_para(True),seed=seed)

    def add_button_cb(self):
        plojo_data = plojo_data_init()
        new_entry=plojo_data.next_exp(1)[0]
        new_index = False
        for key in plojo_data.index.keys():
            if 'Simuojo' in key:
                new_index = key
        if not new_index:
            new_index = plojo_data.new_index('Simuojo Data')
        concentration=list(self.raw_data.data['x'])
        signal=list(self.raw_data.data['y'])
        name=self.name_input.value
        para = self.fetch_para(True)
        tags = self.title_generator(para)
        notes = self.title_generator(para) + '; R_l:{:.3g}; R_p:{:.3g}'.format(self.fetch_para(True)['linear'],self.fetch_para(True)['proportional'])
        dict_tosave = dict(date='SimuDate',concentration=concentration,signal=signal,tag=tags,note=notes,fit_method='kd',author='simu',name=name)
        plojo_data.experiment[new_entry]=dict_tosave
        plojo_data.index[new_index].add(new_entry)
        plojo_data.experiment_to_save.update({new_entry:'sync'})
        plojo_data.save_experiment()
        self.p.title.text='Data Saved to plojo.'
