from bokeh.plotting import figure, ColumnDataSource
from bokeh.io import curdoc
from bokeh.models import HoverTool,Slider,RangeSlider
from bokeh.models.widgets import Button, TextInput,PreText
from bokeh.layouts import widgetbox,layout
import numpy as np
from plojo_module import plojo_data_init


class DR_5PL():
    def __init__(self):
        para_text = PreText(text='Model Parameters')
        para_text_ = PreText(text='Model Parameters')
        self.slider_ec_50 = Slider(title='Model EC50 (nM)', start=-3, end=3, step=0.01, value=0)
        self.slider_hill = Slider(title='Hill Coefficient', start=-1, end=1, step=0.01, value=0)
        self.slider_s = Slider(title='Symmetry Factor S', start=-2, end=2, step=0.01, value=0)
        self.slider_Fmax = Slider(title='Fmax value AU',start = 1000, end = 100000,step = 1000, value = 50000)
        self.slider_Fmin = Slider(title='Fmin value AU',start = 0, end = 10000, step = 100, value = 2000)
        self.name_input = TextInput(title='Create Name for the data',value='Simu_DR_5PL')
        rand_text = PreText(text='Randomization Parameters')
        self.slider_conc_range =RangeSlider(title= 'Experiment Conc. (log nM) : ',  start = -3, end=3, step=0.01,value=(-2,2))
        self.slider_points = Slider(title='Number of concentration points: ', start = 3, end=20, step=1,value=9)
        self.slider_sets = Slider(title='Number of replicates:',start=1,end=9,step=1,value=1)
        self.slider_linear = Slider(title='Linear randomization', start = 0, end=25, step=0.1, value=0.2)
        self.slider_proportional = Slider(title='Proportional randomization', start = 0, end=25, step=0.1,value=0.2)
        self.input_ec_50 = TextInput(value="1", title="Model EC50 (nM):")
        self.input_hill = TextInput(value="1", title="Hill Coefficient:")
        self.input_s = TextInput(value="1", title="Symmetry Factor S:")
        refresh_button = Button(label='Refresh Random',button_type='success')
        add_button = Button(label='Add data to plojo',button_type='success')
        para_slider = widgetbox(para_text,self.slider_ec_50,self.slider_hill,self.slider_s,self.slider_Fmax,self.slider_Fmin,self.name_input)
        para_input = widgetbox(para_text_,self.input_ec_50,self.input_hill,self.input_s)
        rand_opt = widgetbox(rand_text,self.slider_conc_range,self.slider_points,self.slider_sets,self.slider_linear,self.slider_proportional,refresh_button,add_button)
        hover_tool_1 = HoverTool(tooltips=[('Aptamer/nM', '@a_0{0.00}'),( 'Signal', '@signal{0.00}' )],mode='vline')
        self.fit_data = ColumnDataSource(data=self.generate_plotdata(**self.fetch_para(False)))
        self.raw_data = ColumnDataSource(data=self.generate_plotdata(**self.fetch_para(True)))
        self.p = figure(x_axis_label='Aptamer (nM)', y_axis_label='Signal / A.U.', x_axis_type='log')
        self.p.add_tools(hover_tool_1)
        self.p.title.text = 'Dose Response 5 Parameter Logistic Model'
        self.p.title_location='above'
        self.p.line('a_0','signal', source=self.fit_data,line_color='green',line_width=2)
        self.p.circle('a_0','signal',source=self.raw_data,color='red',line_width=3)
        self.p.plot_height = 400
        self.p.plot_width = 700
        self.slider_ec_50.on_change('value',self.callback)
        self.slider_hill.on_change('value',self.callback)
        self.slider_s.on_change('value',self.callback)
        self.slider_Fmax.on_change('value',self.callback)
        self.slider_Fmin.on_change('value',self.callback)
        self.slider_conc_range.on_change('value',self.callback)
        self.slider_points.on_change('value',self.callback)
        self.slider_linear.on_change('value',self.callback)
        self.slider_proportional.on_change('value',self.callback)
        self.slider_sets.on_change('value',self.callback)
        refresh_button.on_click(self.refresh_button_cb)
        add_button.on_click(self.add_button_cb)
        self.layout =([self.p],[para_input,para_slider,rand_opt])


    def randomizer(self,signal,linear=0.001, proportional=0.001,seed=42):
        np.random.seed(seed)
        size = len(signal)
        max_sig = max(signal)
        return signal * np.random.normal(loc=1,scale=proportional,size=size) + max_sig*np.random.normal(loc=0,scale=linear,size=size)

    def refresh_button_cb(self):
        seed = np.random.randint(90000)
        self.raw_data.data = self.generate_plotdata(**self.fetch_para(True),seed=seed)
        self.p.title.text = self.title_generator(self.fetch_para(False))

    def DR_5PL(self,x,ec_50,Fmax,Fmin,Hill,S):
        """
        Five parameters logistic regression
        signal = Fmin + (Fmax-Fmin)*(X**(Hill*S))/(X**Hill + EC50**Hill*(2**(1/S)-1))**S
        """
        denominator = (x**Hill+(ec_50)**Hill*(2**(1/S)-1))**S
        signal = Fmax - (Fmax-Fmin)*(x**(Hill*S))/denominator
        return signal


    def generate_plotdata(self,ec_50=1,Fmax=10000,Fmin=1000,hill=1,s=1,start=None,end=None,point=100,randomize=False,proportional=0.001,linear=0.001,sets=1,**kwargs):
        if start == None and end == None:
            a_0 = np.geomspace(ec_50 / 1000, ec_50 * 1000, 100)
        else:
            a_0= np.geomspace(start,end,point)
        signal = self.DR_5PL(a_0,ec_50,Fmax,Fmin,hill,s)
        if randomize:
            signal = self.randomizer(np.repeat(signal,sets),linear=linear,proportional=proportional,**kwargs)
            result_={'a_0':np.repeat(a_0,sets), 'signal':signal}
        else:
            result_ = {'a_0':a_0, 'signal':signal}
        return result_

    def fetch_para(self,randomize):
        ec_50 = 10 ** self.slider_ec_50.value
        hill = 10 ** self.slider_hill.value
        s = 10 ** self.slider_s.value
        Fmax = self.slider_Fmax.value
        Fmin = self.slider_Fmin.value
        sets=self.slider_sets.value
        start = 10**self.slider_conc_range.value[0]
        end = 10**self.slider_conc_range.value[1]
        points = self.slider_points.value
        linear = self.slider_linear.value/100
        proportional = self.slider_proportional.value/100
        if randomize:
            result = dict(ec_50=ec_50,hill=hill,s=s,Fmax=Fmax,Fmin=Fmin,start=start,end=end,point=points,randomize=True,linear=linear,proportional=proportional,sets=sets)
        else:
            result = dict(ec_50=ec_50,hill=hill,s=s,Fmax=Fmax,Fmin=Fmin,)
        return result


    def callback(self,attr,old,new):
        ec_50 = 10 ** self.slider_ec_50.value
        hill = 10 ** self.slider_hill.value
        s = 10 ** self.slider_s.value
        self.input_ec_50.value = str(ec_50)
        self.input_hill.value = str(hill)
        self.input_s.value = str(s)
        self.fit_data.data = self.generate_plotdata(**self.fetch_para(False))
        self.raw_data.data = self.generate_plotdata(**self.fetch_para(True))
        self.p.title.text = self.title_generator(self.fetch_para(False))

    def title_generator(self,para_dict):
        ec_50 = para_dict['ec_50']
        hill = para_dict['hill']
        s = para_dict['s']
        Fmax = para_dict['Fmax']
        Fmin = para_dict['Fmin']
        title = 'EC50:{:.3g}nM; Hill:{:.3g}; S:{:.3g}; Fmax:{:.3g}, Fmin:{:.3g}'.format(ec_50,hill,s,Fmax,Fmin)
        return title

    def add_button_cb(self):
        plojo_data = plojo_data_init()
        new_entry=plojo_data.next_exp(1)[0]
        new_index = False
        for key in plojo_data.index.keys():
            if 'Simuojo' in key:
                new_index = key
        if not new_index:
            new_index = plojo_data.new_index('Simuojo Data')
        concentration=list(self.raw_data.data['a_0'])
        signal=list(self.raw_data.data['signal'])
        name=self.name_input.value
        para = self.fetch_para(True)
        tags = self.title_generator(para)
        notes = self.title_generator(para) + '; R_l:{:.3g}; R_p:{:.3g}'.format(self.fetch_para(True)['linear'],self.fetch_para(True)['proportional'])
        dict_tosave = dict(date='SimuDate',concentration=concentration,signal=signal,tag=tags,note=notes,fit_method='DR_4PL',author='simu',name=name)
        plojo_data.experiment[new_entry]=dict_tosave
        plojo_data.index[new_index].add(new_entry)
        plojo_data.experiment_to_save.update({new_entry:'sync'})
        plojo_data.save_experiment()
        self.p.title.text='Data Saved to plojo.'
