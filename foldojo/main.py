from foldojo_module import structure_prediction
from bokeh.io import curdoc
predictlayout=structure_prediction().layout


curdoc().add_root(predictlayout)
curdoc().title = chr(128540) + 'Foldojo'
