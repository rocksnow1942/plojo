import os,glob
from os import path
# import shelve
# from bokeh.models import HoverTool,Slider,RangeSlider
from bokeh.models.widgets import Button, TextInput,Div,TextAreaInput,Select,Dropdown
from bokeh.layouts import widgetbox,column,row
# import numpy as np
# from simu_utils import file_save_location,file_name
from _structurepredict import Structure
import datetime
from MSA import Alignment, buildMSA

cache_loc=path.join(path.dirname(__file__),'static','cache')

class structure_prediction():
    def __init__(self):
        #parameters default
        self.default = ['rna','37','50','1','8']+['None']*6
        width=125
        buttonwidth=100

        # define gadgets
        #inputs
        self.sequence =TextAreaInput(title="Enter Sequence:",value='ACCCTTGCTTGCGTAGCATTTTACCAGTGAGTCGGATCTCCGCATATCTGCG',rows=5,cols=145,max_length=5000,width=1000)
        self.name=TextInput(title='Sequence Name',value='NewSequence',width=width)
        self.inputa_backbone =Select(title='Backbone type:',value='rna',options=[('rna','RNA'),('dna','DNA')],width=width)
        self.inputb_SetTemperature = TextInput(title='Set Temperature (C):',value='37',width=width)
        self.inputc_percent= TextInput(title='Max % suboptimal',value='50',width=width)
        self.inputd_window = TextInput(title='Window',value='2',width=width)
        self.inpute_maxstr = TextInput(title='Max Structure NO',value='8',width=width)
        self.inputf_ForceSingleStranded =TextInput(title='Force Single Strand',value='None',width=width)
        self.inputg_ForceDoubleStranded = TextInput(title='Force Double Strand',value='None',width=width)
        self.inputh_ForcePair = TextInput(title='Force Pair',value='None',width=width)
        self.inputi_ForceProhibitPair = TextInput(title='Force No Pair',value='None',width=width)
        # self.inputj_ForceModification = TextInput(title='Modification Site',value='None',width=width)
        # self.inputk_ForceFMNCleavage = TextInput(title='FMN Cleavage',value='None',width=width)
        self.method = Select(title='Prediction Algorithm',value='RNAstructure',width=width,
                                    options=[('RNAstructure','RNAstructure'),('ViennaRNA','ViennaRNA'),('Compare_RV','Compare RS and VR')])
        self.parainputs = [i for k,i in sorted(self.__dict__.items(),key=lambda x:x[0]) if k.startswith('input')]


        # action buttons
        self.predict = Button(label=chr(9193)+' Predict ',button_type='success',width=buttonwidth)
        self.tools = Dropdown(label=chr(127922)+' Tools ',button_type='success',width=buttonwidth,
                            value=None,menu = [('Reset','reset'),None,('',''),('','')])
        self.align = Button(label = chr(128260)+ ' Align',button_type='success',width=buttonwidth)
        self.settings = Dropdown(label=chr(128295)+' Settings', button_type='success',value='.png',
                                    menu = [('As PNG','.png'),('As SVG','.svg'),None,('Parameters','para')('Help','help')],width=buttonwidth)


        # misc
        self.plot = Div(text="",width=800,height=900)
        bottom= Div(text="""
        <p align='center'><a href='http://www.aptitudemedical.com'>Aptitude Medical Systems, Inc.<a></p>
        """,width=800,height=50)
        self.div=Div(width=50)
        self.div1=Div(width=50)
        self.div2=Div(width=50)


        # add callbacks
        self.predict.on_click(self.predict_cb)
        self.tools.on_change('value',self.tools_cb)
        self.method.on_change('value',self.method_cb)
        self.settings.on_change('value',self.settings_cb)
        self.align.on_click(self.align_cb)

        # layouts
        self.layout=([self.sequence],[widgetbox(self.name,self.method,*self.parainputs),column(row(self.predict,self.div,self.align,self.div2,self.tools,self.div1,self.settings),self.plot)],[bottom])

        # status attributes
        self.status=0
        self.plotbackendstatus='.png'

    def align_cb(self):
        name=self.name.value
        backend=self.settings.value
        try:
            sequence=self.parsesequence(self.sequence.value)

            align = buildMSA(sequence,name=name)
            self.sequence.value=align.format()
            savename= name +'_LOGO_'+ datetime.datetime.now().strftime("%m%d_%H%M%S") + backend
            save = path.join(cache_loc,savename)

            align.dna_logo(save=save,show=False)

            self.plot.text="""
            <h3>Alignment:</h3>
            <p style= "font-family:monospace">{}</p>
            <img src="foldojo/static/cache/{}"  hspace='20' height='200' width='700'>
            """.format(align.format().replace('\n','</br>'),savename)

        except Exception as e:
            self.plot.text=str(e)


    def settings_cb(self,attr,old,new):
        if new in ['.png','.svg']:
            return 0
        elif new =='help':
            self.plot.text=helptext
        elif new == 'para':
            # switch to para adjustement interface.
            pass

        self.settings.value=old

    def viennainputtoggle(self,onoroff):
        self.inputa_backbone.disabled=onoroff
        self.inputi_ForceProhibitPair.disabled=onoroff

    def method_cb(self,attr,old,new):

        if new=='ViennaRNA':
            self.viennainputtoggle(True)
        else:
            self.viennainputtoggle(False)

    def tools_cb(self,attr,old,new):
        """
        reset layout
        """
        if new=='none':
            return 0
        elif new == 'reset':
            for i,j in zip(self.default,self.parainputs):
                j.value=i
            self.name.value='NewSequence'
            self.plot.text="<h2>Parameters reset.</h2>"
            self.predict.disabled=False
            self.predict.button_type='success'
        self.tools.value='none'

    def parsesequence(self,seq):
        seq = seq.strip().split('\n')
        result=[]
        for i in seq:
            temp=''.join([j.upper() for j in i if j in "ATGCUWSMKRYBDHVN-" ])
            if temp:
                result.append(temp)
        return result

    def parsepara(self):
        para={'method':self.method.value}
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
                para.update({key:float(value)+273.15})
            else:
                para.update({key:float(value)})

        return para

    def clear_cache(self):
        file_list = sorted(glob.glob(path.join(cache_loc,'*')),key=lambda x: path.getmtime(x))
        if len(file_list)>10:
            os.remove(file_list[0])
        return len(file_list)

    def predict_cb(self,):
        name=self.name.value
        backend=self.settings.value
        try:
            sequence=self.parsesequence(self.sequence.value)
            self.sequence.value='\n'.join(sequence)
            lengthlist=[len(i) for i in sequence]
            assert max(lengthlist) == min(lengthlist), ('input sequence not same length')
            para=self.parsepara()
            ct=para.copy()
            ct.update(sequence=sequence,name=name)
            if ct == self.status and backend==self.plotbackendstatus:
                return 0

            align = Alignment(sequence,name=name)

            self.status=ct
            self.plotbackendstatus=backend
            self.clear_cache()
            save = name +'_'+ datetime.datetime.now().strftime("%m%d_%H%M%S")
            print(para)
            rna=Structure(align,name,save_loc=cache_loc)

            window=para.pop('window',4)
            cluster = para.pop('cluster','hamming')
            center = para.pop('center','energy')
            maxstr = para.pop('maxstr',8)

            rna.fold(**para)
            rna.restrict_fold(window,cluster,center)
            rna.init_dotgraph(maxstr)
            self.structure = rna
            h,w=rna.plot_fold(save=save,plotbackend=self.plotbackendstatus)
            figh,figw= 800*h/(max(h,w)),800*w/(max(h,w))
            self.plot.text="""
            <img src="foldojo/static/cache/{}"  hspace='20' height='{:.0f}' width='{:.0f}'>
            """.format(save+self.plotbackendstatus,figh,figw)
        except Exception as e:
            self.plot.text=str(e)




helptext="""
<h2>Secondary Structure Prediction</h2>
<h3>- based on <a href='https://rna.urmc.rochester.edu/RNAstructure.html'>RNAstructure package<a> and
    <a href='https://www.tbi.univie.ac.at/RNA/'> ViennaRNA package<a> </h3>
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
pair) (<a href='https://pubs.acs.org/doi/10.1021/ja962918p'>Reference<a>).
In subsequent structure prediction, this nucleotide will be in a GU
pair.</p>
<p>Click <b>Predict</b> will use current paramters to predict secondary structures.</p>
<p>Click <b>Reset</b> will reset all parameters.</p>
"""
