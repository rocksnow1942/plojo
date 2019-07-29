import os,glob
from os import path
# import shelve
# from bokeh.models import HoverTool,Slider,RangeSlider
from bokeh.models.widgets import Button, TextInput,Div,TextAreaInput,Select,Dropdown, CheckboxGroup,PreText,Tabs,Panel
from bokeh.layouts import widgetbox,column,row,layout
# from bokeh.models import CustomJS
# import numpy as np
# from simu_utils import file_save_location,file_name
from _structurepredict import Structure,SingleStructureDesign,StructurePerturbation,MultiStructureDesign,Multistrand
import datetime
from MSA import Alignment, buildMSA
cache_loc=path.join(path.dirname(__file__),'static','cache')


"""
known bug:
weirdly, NUPACK defect function won't work if server start from ~ folder on mac mini.
"""


class structure_prediction():
    def __init__(self):
        #parameters default
        width=125
        buttonwidth=100

        # define gadgets
        #inputs
        self.sequence =TextAreaInput(title="Enter Sequence:",value='ACCCTTGCTTGCGTAGCATTTTACCAGTGAGTCGGATCTCCGCATATCTGCG',
                            rows=7,cols=150,max_length=5000,width=1000,css_classes=['custom1'])
        self.name=TextInput(title='Sequence Name',value='NewSequence',width=width)
        self.inputa_backbone =Select(title='Backbone type:',value='rna',options=[('rna','RNA'),('dna','DNA')],width=width)
        self.inputb_SetTemperature = TextInput(title='Set Temperature (C):',value='37',width=width)
        self.inputc_percent= TextInput(title='Max % suboptimal',value='25',width=width)
        self.inputd_window = TextInput(title='Window',value='2',width=width)
        self.inpute_maxstr = TextInput(title='Max Structure NO',value='8',width=width)
        self.inputf_ForceSingleStranded =TextInput(title='Force Single Strand',value='None',width=width)
        self.inputg_ForceDoubleStranded = TextInput(title='Force Double Strand',value='None',width=width)
        self.inputh_ForcePair = TextInput(title='Force Pair',value='None',width=width)
        self.inputi_ForceProhibitPair = TextInput(title='Force No Pair',value='None',width=width)
        self.method = Select(title='Prediction Algorithm',value='RNAstructure',width=width,
                                    options=[('RNAstructure','RNAstructure'),('ViennaRNA','ViennaRNA'),('Nupack', 'Nupack'),('Compare_RVN','Compare Algorithms')])
        self.parainputs = [i for k,i in sorted(self.__dict__.items(),key=lambda x:x[0]) if k.startswith('input')]

        # align parameters
        self.align_gap = TextInput(title='Align Gap Penalty',value='4',width=width)
        self.align_gapext = TextInput(title='Align GapExt Penalty',value='2',width=width)
        self.align_options = CheckboxGroup(labels=["Align Offset", 'Sequence Order'], active=[1])
        self.align_widget = widgetbox(self.align_gap,self.align_gapext,self.align_options)


        # fold parameters
        self.fold_cluster = Select(title='Structure Filter Method',value='hamming',options=[('hamming','Hamming Distance'),('basepair','Base Pair Distance'),('tree','Tree Edit Distance')],width=150)
        self.fold_center = Select(title='Structure Pick Method',value='energy',options=[('energy','Lowest Energy'),('distance','Distance Centroid')],width=150)
        self.fold_widget = widgetbox(self.fold_cluster,self.fold_center)
        self.default_mode_para = row(self.align_widget,self.fold_widget)
        self.para_settings = column(PreText(text='Default Mode Parameter Settings'),self.default_mode_para)

        # NUPACK parameters
        self.dangles = Select(title='Set Dangle Energy',value='some',width=180,options=[('some','Some'),('all','All'),('none','None')])
        self.pseudo = CheckboxGroup(labels=["Allow Pseudo Knot"], active=[],width=180)
        self.energypara = Select(title='RNA Parameters',value='1995',width=180,options=[('1995','Serra and Turner, 1995'),('1999','Mathews et al., 1999')])
        self.sodium = TextInput(title='Sodium Concentration (M)',value='0.15',width=180)
        self.magnesium = TextInput(title='Magnesium Concentration (M)',value='0.002',width=180)
        self.nupack_para=widgetbox(self.energypara,self.dangles,self.pseudo,self.sodium,self.magnesium)

        # plot mode
        self.plot_mode_energy = TextInput(title='Enter Energy Value',value='-1',width=180)

        # structure exclusion parameters
        self.strexc_threshold = TextInput(title='Exclude structure threshold',value='10',width=180)
        self.strexc_method = Select(title='Structure Match Method',value='match',options=[('match','Exact Match'),
                        ('hamming','Hamming Distance'),('basepair','Base Pair Distance'),('tree','Tree Edit Distance'),('starmatch','Star Match')],width=180)
        self.strexc_strthreshold =  TextInput(title='Structure Match Threshold',value='4',width=180)
        self.strexc_mode_para = widgetbox(self.strexc_threshold,self.strexc_method,self.strexc_strthreshold)

        # perturbation mode parameters
        perturb_mode_explain="""<h2>Structure Perturbation</h2>
        <p>Mutate a number of sites concurrently and promote or inhibit a target structure.</p>
        <p>Use ViennaRNA Algorithm only.</p>
        """
        self.perturb_mode_explain = Div(text=perturb_mode_explain,width=700,height=100)
        self.perturb_range = TextInput(title='Mutation Region:',value='ALL',width=180)
        self.perturb_n = TextInput(title='Total Mutation Sites:',value='1',width=180)
        self.perturb_iteration = TextInput(title='Maximum Iterations:',value='1000',width=180)
        self.perturb_goal = Select(title="Perturbation Goal",value='inhibit',options=[('inhibit',"Inhibit Target Structure"),('promote',"Promote Target Structure")])
        self.perturb_mode_para=widgetbox(self.perturb_goal,self.perturb_range,self.perturb_n,self.perturb_iteration)


        # single strand parameters
        sgl_nexplaintext = """<h2>Single Strand Design</h2>
        <p>Enter IUPAC notation sequence and target structure. Then predict.</p>
        <p>Use ViennaRNA Algorithm is much faster than RNAstructure.</p>
        <p> <b>*</b> is allowed to match any structure. However, This will drastically slow down prediction.</p>
        <h3>IUPAC notation:</h3>
        <p> W=AT S=CG M=AC K=GT R=AG Y=CT B=CGT D=AGT H=ACT V=ACG N=AGCT</p>
        """
        self.sgld_n = TextInput(title='Maximum Iterations:',value='10',width=180)
        self.sgld_top = TextInput(title='Show Top:',value='50',width=180)
        self.sgl_nexplain = Div(text=sgl_nexplaintext,width=700,height=200)
        self.sgld_mdoe_para =widgetbox(self.sgld_n, self.sgld_top,self.strexc_method,self.strexc_strthreshold)


        #multistrand design parameters
        mul_nexplaintext = """<h2>Multiple Strand Design</h2>
        <p>Enter IUPAC notation sequence and target structure. Strands seprated by '+'.</p>
        <h3> Only support <b>NUPACK</b> algorithm. </h3>
        <h3>IUPAC notation:</h3>
        <p> W=AT S=CG M=AC K=GT R=AG Y=CT B=CGT D=AGT H=ACT V=ACG N=AGCT</p>
        """
        self.mul_nexplain=Div(text=mul_nexplaintext,width=700,height=200)

        # cofold parameters
        cofold_nexplaintext="""<h2>Co-Fold</h2>
        <p>Enter sequence and concentration in uM. Enter max allowed complex strands.</p>
        <h3> Only support <b>NUPACK</b> algorithm. </h3>
        """
        self.cofold_nexplain=Div(text=cofold_nexplaintext,width=700,height=100)
        self.maxcofoldstrand = TextInput(title='Max Complex Strands:',value='2',width=180)


        # action buttons
        self.predict = Button(label=chr(9193)+' Predict ',button_type='success',width=buttonwidth,css_classes=['button'])
        self.tools = Dropdown(label=chr(127922)+' Tools ',button_type='success',width=buttonwidth,css_classes=['button'],
        value=None,menu = [('Draw Structure','plot'),None,('Default mode','default'),('Structure Exclusion','exclusion'),('Single Strand Design','single'),
        ('Multi-Strand Design','multi'),('Structure Perturbation','perturb'),None,('Co-Fold','cofold'),('Reset Parameter','reset')]) #('Co-Fold','cofold'),
        self.align = Button(label = chr(128260)+ ' Align',button_type='success',width=buttonwidth,css_classes=['button'])
        self.settings = Dropdown(label=chr(128295)+' Settings', button_type='success',value='.png',css_classes=['button'],
                                    menu = [('As PNG','.png'),('As SVG','.svg'),None,('Help','help')],width=buttonwidth)



        # misc
        #
        self.header= Div(text=header.default,width=1000,height=40,css_classes=['custom1'])
        self.plot = Div(text="<h2>Welcome!<br><br>Today is {}.</h2>".format(datetime.datetime.now().strftime("%A, %B %-d")),width=800,height=900,)

        self.text = Div(text='',width=800,height=900)
        self.plottab = Tabs(active=0,width=810,height=930,tabs=[Panel(child=self.plot,title='Output1'),Panel(child=self.text,title='Output2'),
                Panel(child=self.para_settings,title="Parameters")])
        self.div=Div(width=50,css_classes=['divider'])
        self.div1=Div(width=50,css_classes=['divider'])
        self.div2=Div(width=50,css_classes=['divider'])


        # add callbacks
        self.predict.on_click(self.predict_cb)

        # self.predict.callback=CustomJS(args=dict(plot=self.plot),code="""document.getElementById("Test1234").innerHTML = "whatever" """)

        self.tools.on_change('value',self.tools_cb)
        self.method.on_change('value',self.method_cb)
        self.settings.on_change('value',self.settings_cb)
        self.align.on_click(self.align_cb)

        # layouts
        self.layout=layout([self.header],[self.sequence],[widgetbox(self.name,self.method,*self.parainputs,css_classes=['widgetbox']),
                    column(row(self.predict,self.div,self.align,self.div1,self.tools,self.div2,self.settings),self.plottab,css_classes=['plotarea'])],)


        # status attributes
        self.tool_mode='default'
        self.fold_status={}
        self.last_predict=''
        self.align_status={}
        self.align=Alignment()
        self.showingplot=True
        self.plotbackend = '.png'

    def align_cb(self):
        name=self.name.value
        backend=self.plotbackend
        mode = self.tool_mode
        if backend=='dot': backend='.png'
        try:
            self.clear_cache()
            gap=int(self.align_gap.value)
            gapext=int(self.align_gapext.value)
            offset=0 in self.align_options.active
            order=1 in self.align_options.active
            para = dict(gap=gap,gapext=gapext,offset=offset,order=order)
            sequence=self.parsesequence()
            if mode == 'default':
                if set(sequence)==set(self.align.seq) and para==self.align_status:
                    self.align.name = name
                else:

                    self.align = buildMSA(sequence,name=name,**para)

                    self.align_status=para
                self.sequence.value=self.align.format()
                savename= name +'_LOGO_'+ datetime.datetime.now().strftime("%m%d_%H%M%S") + backend
                save = path.join(cache_loc,savename)

                self.align.dna_logo(save=save,show=False)

                self.plot.text="""
                <h3>Alignment:</h3>
                <p style= "font-family:monospace">{}</p>
                <img src="foldojo/static/cache/{}"  hspace='20' height='200' width='700'>
                """.format(self.align.format(index=True,link=True).replace('\n','</br>'),savename)
            elif mode == 'exclusion':
                seq1,seq2=sequence
                ali1 = buildMSA(seq1,name=name,**para)
                ali2 = buildMSA(seq2,name=name,**para)
                ali1s,ali2s=ali1.format(),ali2.format()
                self.sequence.value=ali1s+'\n/\n'+ali2s
                ali1l,ali2l = list(zip(*ali1s.split('\n'))), list(zip(*ali2s.split('\n')))
                linker=''
                for i,j in zip(ali1l,ali2l):
                    if set(i)==set(j):
                        linker+='|'
                    else:
                        linker+='*'
                self.plot.text="""
                <h3>Alignment</h3>
                <p style= "font-family:monospace">{}<br>{}<br>{}</p>
                """.format(ali1.format(index=True).replace('\n','</br>'),linker,ali2.format().replace('\n','</br>'))

        except Exception as e:
            self.plot.text=str(e)

    def settings_cb(self,attr,old,new):
        if new in ['.png','.svg','dot']:
            self.plotbackend = new
            return 0
        elif new =='help':
            self.plot.text=helptext
        self.settings.value=old

    def viennainputtoggle(self,onoroff):
        self.inputa_backbone.disabled=onoroff
        self.inputi_ForceProhibitPair.disabled=onoroff

    def method_cb(self,attr,old,new):

        if new=='ViennaRNA':
            self.viennainputtoggle(True)
        else:
            self.viennainputtoggle(False)

        if new in ('Nupack','Compare_RVN'):
            self.para_settings.children.append(self.nupack_para)
        if old in ('Nupack','Compare_RVN'):
            self.para_settings.children.pop(-1)

    def tools_cb(self,attr,old,new):
        """
        reset layout
        """
        if new=='none':
            return 0
        elif new == 'reset':
            self.fold_status={}
            self.name.value='NewSequence'
            self.plot.text="<h2>Parameters reset.</h2>"
            self.align_gap.value='4'
            self.align_gapext.value='2'
            self.align_options.active=[1]
            self.fold_cluster.value='hamming'
            self.fold_center.value='energy'
        else:
            self.header.text=getattr(header,new)
            self.tool_mode=new
            self.change_parameter(new)
        self.tools.value='none'

    def change_parameter(self,mode):
        if mode == 'default':
            recommendmethod=self.method.value
            self.para_settings.children = [PreText(text='Default Mode Parameter Settings'),self.default_mode_para]
            if not self.fold_status:self.sequence.value = 'ACCCTTGCTTGCGTAGCATTTTACCAGTGAGTCGGATCTCCGCATATCTGCG'

        elif mode == 'plot':
            recommendmethod=self.method.value
            self.para_settings.children = [Div(text="<h2>Nothing To Set Here.</h2>"),self.plot_mode_energy]
            if not self.fold_status:self.sequence.value = 'WANTSUBWAYHAMBURGRSTHURSDAYWHACKYSHMUCHYUK\n((((((..((((....))))..((((....))))..))))))'


        elif mode == 'exclusion':
            recommendmethod=self.method.value
            self.para_settings.children = [PreText(text='Structure Exclusion Mode Parameter Settings'),self.strexc_mode_para]
            if not self.fold_status:self.sequence.value = 'ACCCTTGCTTGCGTAGCATTTTACCAGTGAGTCGGATCTCCGCATATCTGCG\n/\nACCCTTGCTTGCGTAGCATTTTGCCAGTGAGTCGGATCTCCGCATATCACGC'

        elif mode == 'single':
            recommendmethod='ViennaRNA'
            self.para_settings.children = [self.sgl_nexplain,self.sgld_mdoe_para]
            if not self.fold_status:self.sequence.value = 'SWRNATATNYWS\n(((......)))'

        elif mode == 'perturb':
            recommendmethod='ViennaRNA'
            self.para_settings.children = [self.perturb_mode_explain,self.perturb_mode_para]
            if not self.fold_status:self.sequence.value = 'ACCCTTGCTTGCGTAGCATTTTACCAGTGAGTCGGATCTCCGCATATCTGCG\n......(((.(((........))).)))(((......)))((((....))))'

        elif mode == 'multi':
            self.para_settings.children = [self.mul_nexplain,self.sgld_top]
            recommendmethod='Nupack'
            if not self.fold_status:self.sequence.value = 'ACGNNNNNN+NNNNNNNNNNNAGCACNNNNNN\n...((((((+))))))....((((...)))).'

        elif mode == 'cofold':
            self.para_settings.children = [self.cofold_nexplain,self.maxcofoldstrand]
            recommendmethod='Nupack'
            if not self.fold_status:self.sequence.value = 'ACCCTTGCTTGCGTAGCATTTTACCA-1uM\nTGAGTCGGATCTCCGCATATCTGCG-0.5uM'

        self.method_cb(1,'old',self.method.value)
        self.method.value=recommendmethod

    def _sequence_cleaner(self,seq):
        seq = seq.strip().split('\n')
        result=[]
        for i in seq:
            temp=''.join([j.upper() for j in i if j in "ATGCUWSMKRYBDHVN-(.){}*+" ])
            if temp:
                result.append(temp)
        return result

    def parsesequence(self,*args):
        mode =self.tool_mode
        seq = self.sequence.value
        if mode in ['default','plot']:
            return self._sequence_cleaner(seq)
        elif mode == 'exclusion':
            seq1,seq2 = seq.strip().split('/')
            seq1,seq2=self._sequence_cleaner(seq1),self._sequence_cleaner(seq2)
            return seq1,seq2
        elif mode in ['single','perturb','multi']:
            r = self._sequence_cleaner(seq)
            assert len(r)==2,('Must enter two row.')
            assert set(r[1])<=set("(.*){}+"),('Wrong structure format.')
            return r
        elif mode == 'cofold':
            r=seq.strip().split('\n')
            cl,sl=[],[]
            for i in r:
                conc=self.parse_conc(i.split('-')[1])
                s=self._sequence_cleaner(i.split('-')[0])
                cl.append(conc)
                sl.extend(s)
            return sl,cl

    def parse_conc(self,conc):
        amout=float(conc.lower().rstrip('nmupf'))
        unit=conc.lower().lstrip('0123456789.')
        unitdict={'um':1e-6,'m':1,'nm':1e-9,'mm':1e-3,'pm':1e-12,'fm':1e-15}
        return amout*unitdict[unit]


    def parsepara(self):
        mode =self.tool_mode
        para={'method':self.method.value}
        if self.method.value in ('Nupack','Compare_RVN'):
            pseudo=0 in self.pseudo.active
            sodium=float(self.sodium.value)
            magnesium=float(self.magnesium.value)
            para.update(dangles=self.dangles.value,pseudo=pseudo,energypara=self.energypara.value,sodium=sodium,magnesium=magnesium)

        para.update(cluster=self.fold_cluster.value,center=self.fold_center.value)
        for k,i in filter(lambda x:x[0].startswith('input'),self.__dict__.items()):
            key=k.split('_')[1]
            value=i.value.strip()
            if 'Force' in key:
                if value == 'None' or (not value):
                    continue
                elif '-' in value:
                    _v=[(int(i.split('-')[0]),int(i.split('-')[1])) for i in value.split(',')]
                else:
                    _v=[int(i) for i in value.split(',')]
                para.update({key:_v})
            elif key=='backbone':
                para.update(backbone=value)
            elif key=='SetTemperature':
                para.update({key:float(value)})
            else:
                para.update({key:float(value)})

        if mode == 'exclusion':
            para.update(exclusion={'threshold':float(self.strexc_threshold.value),
                        'method':self.strexc_method.value,'diffthreshold':int(self.strexc_strthreshold.value)})
        if mode in ['single','multi']:
            para.update(n=int(self.sgld_n.value),top=int(self.sgld_top.value),
                        strexc_method=self.strexc_method.value,strexc_threshold=int(self.strexc_strthreshold.value))

        if mode == 'cofold':
            para.update(maxcofoldstrand=int(self.maxcofoldstrand.value))

        if mode =='perturb':
            mutrange = self.perturb_range.value
            if '-' in mutrange:
                mutrange= [(int(i.split('-')[0]),int(i.split('-')[1])) for i in mutrange.split(',')]
            else:
                mutrange=None
            para.update(n=int(self.perturb_n.value),mutrange=mutrange,
                    iteration=int(self.perturb_iteration.value),goal=self.perturb_goal.value)

        return para


    def clear_cache(self):
        file_list = sorted(glob.glob(path.join(cache_loc,'*')),key=lambda x: path.getmtime(x))
        if len(file_list)>10:
            os.remove(file_list[0])
        return len(file_list)

    def predict_cb(self,):
        name=self.name.value
        backend=self.plotbackend
        mode=self.tool_mode
        try:
            self.clear_cache()
            sequence=self.parsesequence()

            para=self.parsepara()
            ct=para.copy()
            ct.update(sequence=sequence,name=name,backend=backend,mode=mode)

            if ct == self.fold_status:
                self.plot.text="""
                <img src="foldojo/static/cache/{}"  hspace='20' height='{:.0f}' width='{:.0f}'>
                """.format(*self.last_predict)
                return 0

            # only creaet new structure if sequence and key parameters changed.
            diff=[]
            for k,i in ct.items():
                if i != self.fold_status.get(k,None):
                    diff.append(k)
            plotnew = False
            if set(diff)<=set(['window','name','cluster','center','maxstr']):
                rna=self.structure
                rna.name=name
            else:
                plotnew=True

            #pop out unncessary parameters before pass to fold.
            window=para.pop('window',4)
            cluster = para.pop('cluster','hamming')
            center = para.pop('center','energy')
            maxstr = para.pop('maxstr',8)

            if plotnew:
                lengthlist=[len(i) for i in sequence]
                if mode == 'cofold':
                    return self.update_cofold(sequence,para)
                elif mode == 'default':
                    assert max(lengthlist) == min(lengthlist), ('input sequence not same length')
                    self.sequence.value='\n'.join(sequence)
                    align = Alignment(sequence,name=name)

                    rna=Structure(align,name,save_loc=cache_loc).fold(**para)

                elif mode == 'plot':
                    assert max(lengthlist) == min(lengthlist), ('input sequence not same length')
                    self.sequence.value='\n'.join(sequence)
                    ntseq=[]
                    for i in sequence:
                        if '.' in i or '(' in i:
                            dotbracket = i
                        else: ntseq.append(i)
                    rna=Structure(name=name,save_loc=cache_loc)
                    rna.plot_dot_bracket(dotbracket,seq=ntseq,energy=float(self.plot_mode_energy.value),**para)
                elif mode =='exclusion':
                    self.sequence.value='\n'.join(sequence[0])+'\n/\n'+'\n'.join(sequence[1])
                    align1 = Alignment(sequence[0])
                    align2 = Alignment(sequence[1])
                    exclusionpara=para.pop('exclusion',{})
                    rna=Structure(align1,name,save_loc=cache_loc).fold(**para)
                    b=Structure(align2,name,save_loc=cache_loc).fold(**para)
                    rna.subtract(b,**exclusionpara)
                elif mode in ['single','perturb','multi']:
                    title ={'single':'Single Strand Design',"perturb":"Structure Perturbation",'multi':"Multi-Strand Design"}[mode]
                    self.sequence.value='\n'.join(sequence)
                    seq,target=sequence
                    cls = {'single':SingleStructureDesign,'perturb':StructurePerturbation,'multi':MultiStructureDesign}[mode]
                    ssd=cls(seed=seq,target=target,**para)
                    _=ssd.generate(**para)
                    self.text.text="""
                    <table style="width:100%">
                    <caption style="text-align:center;font-family:cursive;font-size:150%">{} Results: ({}/{:.1e})</caption>
                    {}
                    """.format(title,len(ssd.result),ssd.totalmutation,ssd.to_html())
                    rna=Structure(ssd.result[0][0],name,save_loc=cache_loc).fold(**para)
                self.structure = rna


            save = name +'_'+ datetime.datetime.now().strftime("%m%d_%H%M%S")
            rna.restrict_fold(window,cluster,center)

            if mode not in ['single','perturb','multi']: # avoid update text display in other modes
                self.text.text="""
                <table style="width:100%">
                <caption style="text-align:center;font-family:cursive;font-size:150%">Structure Prediction Results: ({})</caption>
                {}
                """.format(len(rna.dot),rna.print())

            rna.init_dotgraph(maxstr)
            h,w=rna.plot_fold(save=save,plotbackend=backend,showpara=(mode!='plot'))
            figh,figw= 800*h/(max(h,w)),800*w/(max(h,w))
            self.plot.text="""
            <img src="foldojo/static/cache/{}"  hspace='20' height='{:.0f}' width='{:.0f}'>
            """.format(save+backend,figh,figw)
            self.last_predict=(save+backend,figh,figw)
            self.fold_status=ct
        except Exception as e:
            self.plot.text=str(e)

    def update_cofold(self,sequence,para):
        ms=Multistrand(*sequence)
        ms.fold(**para)
        title=' '.join(['S'+str(i+1)+':'+"{:.1e}".format(j) for i,j in enumerate(sequence[1])])
        self.text.text="""
        <table style="width:100%">
        <caption style="text-align:center;font-family:cursive;font-size:150%">Co-Fold Structure: {}</caption>
        {}
        """.format(title,ms.to_html())
        self.plot.text="""
        <table style="width:100%">
        <caption style="text-align:center;font-family:cursive;font-size:150%">Co-Fold Concentration: {}</caption>
        {}
        """.format(title,ms.to_html('complex_df'))


helptext="""
<h2>Secondary Structure Prediction (This is outdated help.)</h2>
<h3>- based on <a href='https://rna.urmc.rochester.edu/RNAstructure.html'>RNAstructure package</a> and
    <a href='https://www.tbi.univie.ac.at/RNA/'> ViennaRNA package</a> </h3>
<p>Predicts the lowest free energy structure and suboptimal structures.</p>
<h3>How To Use:</h3>
<p>1. <b>Enter Sequence</b> Only enter oligo sequence. Case doesn't matter. Letters other than "ATGCUatgcu" will be automatically filtered out.
That is, f[G]f[U]f{T}ome-t = GUTT</p>
<p>2. <b>Sequence Name</b> Sequence name is used for plot title and figure file name.</p>
<p>3. <b>Backbone Type</b> Select RNA or DNA. Thermodynamics paramters will be different for DNA or RNA.</p>
<p>4. <b>Set Temperature</b> This allows the user to specify folding temperatures other than 37deg C.  </p>
<p>5. <b>Max % suboptimal</b> is the maximum % difference in free energy in suboptimal structures from the lowest free energy structure. The higher the number, more suboptimal structures will be permitted. </p>
<p>6. <b>Window</b> is a parameter that specifies how different the suboptimal
structures should be from each other (0=no restriction and larger
integers require structures to be more different). </p>
<p>7. <b>Max Structure NO</b> is the maximum number of suboptimal structures to
generate.</p>
<p>8. <b>Force Single Strand</b> Force one or more nucleotide to be single stranded. format: n.t. position index separated by ",". e.g. 1,4,38 </p>
<p>9. <b>Force Double Strand</b> Force one or more nucleotide to be double stranded. format: same as #8.</p>
<p>10.<b>Force Pair</b> Force a pair between two nucleotides. format: use "-" to link paired n.t., separate multiple pairs by ",". e.g. 1-39,2-38,3-37 </p>
<p>11.<b>Force No Pair</b> Prohibit a pair between two nucleotides. format: same as #10.</p>
<p>12.<b>Modification Site</b> This paramater indicates a nucleotide that is accessible to chemical
modification. In subsequent structure prediction, this nucleotide will
be single stranded, at the end of a helix, or in or adjacent to a GU
pair. format: same as #8</p>
<p>13.<b>FMN Cleavage</b> Indicate a nucleotide that is accessible to FMN cleavage (a U in GU
pair) (<a href='https://pubs.acs.org/doi/10.1021/ja962918p'>Reference</a>).
In subsequent structure prediction, this nucleotide will be in a GU
pair.</p>
<p>Click <b>Predict</b> will use current paramters to predict secondary structures.</p>
<p>Click <b>Reset</b> will reset all parameters.</p>
"""


class Header():
    def __init__(self):
        self.template="""
        <style>
        h1 {{
          position: relative;
          animation: mymove 2s;
          animation-iteration-count: infinite;
        }}

        @keyframes mymove {{
        0%  {{ border:4px outset #81F7F3;  }}
        20% {{ border:4px outset #81F7F3; }}
        50% {{ border:4px outset #ff66ff; }}
        70% {{ border:4px outset #81F7F3; }}
        100%{{ border:4px outset #81F7F3;}}
        }}
        </style>
        <h1 style="width:1050px;height:50px;border: 4px outset #81F7F3;text-align:center;font-family:cursive;font-size:230%;color:#FF00BF;background-color:{color}"">
        &#128540
        <span style="color:#0000FF">F</span>
        <span style="color:red">O</span>
        <span style="color:#FFFF00">L</span>
        <span style="color:#31B404">D</span>
        <span style="color:#FF00BF">ojo </span>
        &#129322
        {subtitle}
        </h1>
        """

    @property
    def default(self):
        return self.template.format(subtitle='',color='#81F7F3')
    @property
    def plot(self):
        return self.template.format(subtitle='Plot Structure Mode',color='#088A85')
    @property
    def exclusion(self):
        return self.template.format(subtitle='Structure Exclusion Mode',color='#0404B4')
    @property
    def single(self):
        return self.template.format(subtitle='Single Strand Design Mode',color='#2EFE2E')
    @property
    def multi(self):
        # two strand design
        return self.template.format(subtitle='Multi-Strand Design',color='#100719')
    @property
    def cofold(self):
        return self.template.format(subtitle='Co-Fold of Multiple Sequences',color='#FE9A2E')
    @property
    def perturb(self):
        return self.template.format(subtitle='Structure Perturbation Mode',color='#F5A9D0')


header=Header()
