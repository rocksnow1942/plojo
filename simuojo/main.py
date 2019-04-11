
from simuojo_module import ic50_simu,DR_5PL,Kd_simu,ric50_simu
from bokeh.io import curdoc
from bokeh.layouts import row, layout
from bokeh.models.widgets import Div, Select


simu_sele = Select(title='Select a simulation model',value='Select a model',options=['Select a model','kd simulation','ic_50 simulation','ric_50 simulation','Dose-Response 5 Para Logistic'])

model_display = Div(text="""<div style="text-align:center;"><h1>Welcome to Simuojo!</h1> <div>""",width=600,height=75)

def simu_sele_cb(attr,old,new):
    if new == 'kd simulation':
        kd_layout = Kd_simu().layout
        temp = layout([model_display,simu_sele],*kd_layout)
        select_layout.children = temp.children
        model_display.text = """ <div style="text-align:center;">
        <img src="simuojo/static/kd_formula.png" height='50' width="530" "><div>"""
    elif new=='ic_50 simulation':
        ic_50_layout = ic50_simu().layout
        temp = layout([model_display,simu_sele],*ic_50_layout)
        select_layout.children = temp.children
        model_display.text = """ <div style="text-align:center;">
        <img src="simuojo/static/ic_50.png" height='75' width="403" "><div>"""
    elif new == 'ric_50 simulation':
        ric_50_layout = ric50_simu().layout
        temp = layout([model_display,simu_sele],*ric_50_layout)
        select_layout.children = temp.children
        model_display.text = """ <div style="text-align:center;">
        <img src="simuojo/static/ric_50.png" height='75' width="485" "><div>"""
    elif new == 'Dose-Response 5 Para Logistic':
        dr5pl_layout = DR_5PL().layout
        temp = layout([model_display,simu_sele],*dr5pl_layout)
        select_layout.children = temp.children
        model_display.text = """ <div style="text-align:center;">
        <img src="simuojo/static/DR_5PL.png" height='60' width="503" "><div>"""


simu_sele.on_change('value',simu_sele_cb)
select_layout=layout([row(model_display,simu_sele)])

curdoc().add_root(select_layout)
