from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import HoverTool,Slider,RangeSlider
from bokeh.models.widgets import Button, TextInput,PreText
from bokeh.layouts import widgetbox
import numpy as np
from plojo_module import plojo_data_init



class ic50_simu():
    def __init__(self):
        para_text = PreText(text='Binding Parameters')
        para_text_ = PreText(text='Binding Parameters')
        self.slider_kd_1 = Slider(title='log(Receptor-VEGF Kdr) (nM)', start=-3, end=3, step=0.01, value=0)
        self.slider_kd_2 = Slider(title='log(Aptamer-VEGF Kda) (nM)', start=-3, end=3, step=0.01, value=0)
        self.slider_r_0 = Slider(title='log(Receptor conc.) (nM)', start=-3, end=3, step=0.01, value=0)
        self.slider_v_0 = Slider(title='log(VEGF conc.) (nM)', start=-3, end=3, step=0.01, value=0)
        self.slider_Fminmax = RangeSlider(title = 'Singal Range AU', start =0, end =100000,step =1000,value =(1000,10000))
        self.name_input = TextInput(title='Create Name for the data',value='Simu_IC50')
        rand_text = PreText(text='Randomization Parameters')
        self.slider_conc_range =RangeSlider(title= 'Experiment Conc. (log nM) : ',  start = -5, end=5, step=0.01,value=(-2,2))
        self.slider_points = Slider(title='Number of concentration points: ', start = 3, end=20, step=1,value=9)
        self.slider_sets = Slider(title='Number of replicates:',start=1,end=9,step=1,value=1)
        self.slider_linear = Slider(title='Linear randomization', start = 0, end=25, step=0.1, value=0.2)
        self.slider_proportional = Slider(title='Proportional randomization', start = 0, end=25, step=0.1,value=0.2)
        self.input_kd_1 = TextInput(value="1", title="Receptor-VEGF Kdr (nM):")
        self.input_kd_2 = TextInput(value="1", title="Aptamer-VEGF Kda (nM):")
        self.input_r_0 = TextInput(value="1", title="Receptor Conc. (nM):")
        self.input_v_0 = TextInput(value="1", title="VEGF Conc. (nM):")
        self.input_ic50 = TextInput(value='2',title='Caculated IC50 (nM):')
        refresh_button = Button(label='Refresh Random',button_type='success')
        add_button = Button(label='Add data to plojo',button_type='success')
        para_slider = widgetbox(para_text,self.slider_kd_1,self.slider_r_0,self.slider_kd_2,self.slider_v_0,self.slider_Fminmax,self.name_input)
        para_input = widgetbox(para_text_,self.input_kd_1,self.input_r_0,self.input_kd_2,self.input_v_0,self.input_ic50)
        rand_opt = widgetbox(rand_text,self.slider_conc_range,self.slider_points,self.slider_sets,self.slider_linear,self.slider_proportional,refresh_button,add_button)
        self.fit_data = ColumnDataSource(data=self.generate_plotdata(**self.fetch_para(False)))
        self.raw_data = ColumnDataSource(data=self.generate_plotdata(**self.fetch_para(True)))
        self.p = figure(x_axis_label='Aptamer Concentration (nM)', y_axis_label='Receptor and Aptamer Signal A.U.', x_axis_type='log')
        self.p.title.text = 'Receptor-VEGF complex / Receptor percentage'
        self.p.title_location='above'
        p_1=self.p.line('a_0','rv_per', source=self.fit_data,line_color='green',line_width=2,alpha=0.75,legend='Receptor')
        hover_tool_1 = HoverTool(renderers=[p_1],tooltips=[('Aptamer/nM', '@a_0{0.000}'),( 'Rcpt Signal', '@rv_per{0}' )],mode='vline')
        self.p.add_tools(hover_tool_1)
        p_2=self.p.line('a_0','av_per', source=self.fit_data,line_color='red',line_width=2,alpha=0.5,legend='Aptamer')
        hover_tool_2 = HoverTool(renderers=[p_2],tooltips=[('Aptamer/nM', '@a_0{0.000}'),( 'Aptamer Signal', '@av_per{0}' )],mode='vline')
        self.p.add_tools(hover_tool_2)
        self.p.circle('a_0','rv_per',source=self.raw_data,color='red',line_width=3,alpha=0.75,legend='Receptor')
        self.p.circle('a_0','av_per',source=self.raw_data,color='green',line_width=3,alpha=0.5,legend='Aptamer')
        self.p.legend.click_policy = 'hide'
        self.p.legend.location = 'center_left'
        self.p.legend.border_line_alpha = 0
        self.p.legend.background_fill_alpha = 0.1
        self.p.plot_height = 400
        self.p.plot_width = 600
        self.pn = figure(x_axis_label='Aptamer Concentration (nM)', y_axis_label='Normalized Receptor Aptamer Signal', x_axis_type='log',x_range=self.p.x_range)
        self.pn.title.text = 'Normalized Signal'
        self.pn.title_location='above'
        pn_1=self.pn.line('a_0','rv_norm', source=self.fit_data,line_color='green',line_width=2,alpha=0.75,legend='Receptor')
        hover_tool_1n = HoverTool(renderers=[pn_1],tooltips=[('Aptamer/nM', '@a_0{0.00}'),( 'RV signal', '@rv_norm{0.00}' )],mode='vline')
        self.pn.add_tools(hover_tool_1n)
        pn_2=self.pn.line('a_0','av_norm', source=self.fit_data,line_color='red',line_width=2,alpha=0.5,legend='Aptamer')
        hover_tool_2n = HoverTool(renderers=[pn_2],tooltips=[('Aptamer/nM', '@a_0{0.00}'),( 'AV signal', '@av_norm{0.00}' )],mode='vline')
        self.pn.add_tools(hover_tool_2n)
        self.pn.circle('a_0','rv_norm',source=self.raw_data,color='red',line_width=3,alpha=0.75,legend='Receptor')
        self.pn.circle('a_0','av_norm',source=self.raw_data,color='green',line_width=3,alpha=0.5,legend='Aptamer')
        self.pn.legend.click_policy = 'hide'
        self.pn.legend.location = 'center_left'
        self.pn.legend.border_line_alpha = 0
        self.pn.legend.background_fill_alpha = 0.1
        self.pn.plot_height = 400
        self.pn.plot_width = 600
        self.slider_kd_1.on_change('value',self.callback)
        self.slider_kd_2.on_change('value',self.callback)
        self.slider_r_0.on_change('value',self.callback)
        self.slider_v_0.on_change('value',self.callback)
        self.slider_Fminmax.on_change('value',self.callback)
        self.slider_conc_range.on_change('value',self.callback)
        self.slider_points.on_change('value',self.callback)
        self.slider_linear.on_change('value',self.callback)
        self.slider_proportional.on_change('value',self.callback)
        self.slider_sets.on_change('value',self.callback)
        refresh_button.on_click(self.refresh_button_cb)
        add_button.on_click(self.add_button_cb)
        self.layout =([self.p,self.pn],[para_input,para_slider,rand_opt])

    def rv_solver(self,a_0, r_0, v_0, kd_r, kd_a):
        a=kd_a-kd_r
        d=-v_0*kd_a*r_0**2
        b=kd_r*r_0-2*kd_a*r_0+kd_r**2-kd_r*kd_a-kd_r*a_0+kd_r*v_0-kd_a*v_0
        c=r_0**2*kd_a+kd_r*kd_a*r_0-v_0*r_0*kd_r+2*v_0*r_0*kd_a+a_0*kd_r*r_0
        if a!=0:
            b=b.astype('complex128')
            c=c.astype('complex128')
            p=(3*a*c-b**2)/(3*a**2)
            q=(27*a**2*d-9*a*b*c+2*b**3)/(27*a**3)
            middle=((q/2)**2+(p/3)**3)**(1/2)
            left=(-q/2+middle)**(1/3)
            right=(-q/2-middle)**(1/3)
            w = (3**(1/2)*1j-1)/2
            convert=b/(3*a)
            x_1=left+right-convert
            x_2=w*left+w**2*right-convert
            x_3=w**2*left+w*right-convert
        else:
            x_1=(-c+(c**2-4*b*d)**(1/2))/(2*b)
            x_2=(-c-(c**2-4*b*d)**(1/2))/(2*b)
            x_3=np.zeros(len(a_0))
        min_rv=min(r_0,v_0)
        result=np.array([])
        for root in zip(x_1,x_2,x_3):
            real_root = [i.real for i in root if 0<i.real<min_rv]
            if real_root:
                real_root = min(real_root)
            else:
                real_root=0
            result=np.append(result,real_root)
        return result


    def randomizer(self,signal,linear=0.001, proportional=0.001,seed=42):
        np.random.seed(seed)
        size = len(signal)
        max_sig = max(signal)
        return signal * np.random.normal(loc=1,scale=proportional,size=size) + max_sig*np.random.normal(loc=0,scale=linear,size=size)


    def generate_plotdata(self,r_0=1,v_0=1,kd_1=1,kd_2=1,Fmax=10000,Fmin=1000,start=None,end=None,point=100,randomize=False,proportional=0.001,linear=0.001,sets=1,**kwargs):
        if start == None and end == None:
            ic_50 = kd_2*(1+r_0/kd_1)
            a_0 = np.geomspace(ic_50 / 1000, ic_50 * 1000, 300)
        else:
            a_0= np.geomspace(start,end,point)
        rv = self.rv_solver(a_0,r_0,v_0,kd_1,kd_2)
        av =v_0-rv-kd_1*rv/(r_0-rv)
        rv_per = rv/v_0*(Fmax-Fmin)+Fmin
        av_per = av/v_0*(Fmax-Fmin)+Fmin
        if randomize:
            rv_per = self.randomizer(np.repeat(rv_per,sets),linear=linear,proportional=proportional,**kwargs)
            av_per = self.randomizer(np.repeat(av_per,sets),linear=linear,proportional=proportional,**kwargs)
            rv_norm = (rv_per-min(rv_per))/(max(rv_per)-min(rv_per))*100
            av_norm = (av_per-min(av_per))/(max(av_per)-min(av_per))*100
            result_={'a_0':np.repeat(a_0,sets), 'rv_per':rv_per,'av_per':av_per,'av_norm':av_norm,'rv_norm':rv_norm}
        else:
            rv_norm = (rv_per-min(rv_per))/(max(rv_per)-min(rv_per))*100
            av_norm = (av_per-min(av_per))/(max(av_per)-min(av_per))*100
            result_ = {'a_0':a_0, 'rv_per':rv_per, 'av_per':av_per,'av_norm':av_norm,'rv_norm':rv_norm}
        return result_


    def fetch_para(self,randomize):
        kd_1 = 10 ** self.slider_kd_1.value
        kd_2 = 10 ** self.slider_kd_2.value
        r_0 = 10 ** self.slider_r_0.value
        v_0 = 10 ** self.slider_v_0.value
        sets=self.slider_sets.value
        Fmax = self.slider_Fminmax.value[1]
        Fmin = self.slider_Fminmax.value[0]
        start = 10**self.slider_conc_range.value[0]
        end = 10**self.slider_conc_range.value[1]
        points = self.slider_points.value
        linear = self.slider_linear.value/100
        proportional = self.slider_proportional.value/100
        if randomize:
            result = dict(kd_1=kd_1,kd_2=kd_2,r_0=r_0,v_0=v_0,Fmax=Fmax,Fmin=Fmin,start=start,end=end,point=points,randomize=True,linear=linear,proportional=proportional,sets=sets)
        else:
            result = dict(kd_1=kd_1,kd_2=kd_2,r_0=r_0,v_0=v_0,Fmax=Fmax,Fmin=Fmin)
        return result


    def callback(self,attr,old,new):
        kd_1 = 10 ** self.slider_kd_1.value
        kd_2 = 10 ** self.slider_kd_2.value
        r_0 = 10 ** self.slider_r_0.value
        v_0 = 10 ** self.slider_v_0.value
        self.input_kd_1.value = str(kd_1)
        self.input_kd_2.value = str(kd_2)
        self.input_r_0.value = str(r_0)
        self.input_v_0.value = str(v_0)
        self.input_ic50.value = str(kd_2*(1+r_0/kd_1))
        self.fit_data.data = self.generate_plotdata(**self.fetch_para(False))
        self.raw_data.data = self.generate_plotdata(**self.fetch_para(True))
        self.p.title.text = self.title_generator(self.fetch_para(False))


    def refresh_button_cb(self):
        seed = np.random.randint(90000)
        self.raw_data.data = self.generate_plotdata(**self.fetch_para(True),seed=seed)
        self.p.title.text = self.title_generator(self.fetch_para(False))

    def title_generator(self,para_dict):
        kd_1 = para_dict['kd_1']
        kd_2 = para_dict['kd_2']
        a_0 = para_dict['r_0']
        v_0 = para_dict['v_0']
        Fmax = para_dict['Fmax']
        Fmin = para_dict['Fmin']
        title = 'kd_r:{:.3g}nM; kd_a:{:.3g}nM; r_0:{:.3g}; v_0:{:.3g}; Fmax:{:.3g}, Fmin:{:.3g}'.format(kd_1,kd_2, a_0,v_0,Fmax,Fmin)
        return title

    def add_button_cb(self):
        plojo_data=plojo_data_init()
        new_entry=plojo_data.next_exp(1)[0]
        new_index = False
        for key in plojo_data.index.keys():
            if 'Simuojo' in key:
                new_index = key
        if not new_index:
            new_index = plojo_data.new_index('Simuojo Data')
        concentration=list(self.raw_data.data['a_0'])
        signal=list(self.raw_data.data['rv_per'])
        name=self.name_input.value
        para = self.fetch_para(True)
        tags = 'IC50: {:.4g}nM'.format(float(self.input_ic50.value))
        notes = self.title_generator(para) + '; R_l:{:.3g}; R_p:{:.3g}'.format(self.fetch_para(True)['linear'],self.fetch_para(True)['proportional'])
        dict_tosave = dict(date='SimuDate',concentration=concentration,signal=signal,tag=tags,note=notes,fit_method='ic_50',author='simu',name=name)
        plojo_data.experiment[new_entry]=dict_tosave
        plojo_data.index[new_index].add(new_entry)
        plojo_data.experiment_to_save.update({new_entry:'sync'})
        plojo_data.save_experiment()
        self.p.title.text='Data Saved to plojo.'
