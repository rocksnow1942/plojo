from foldojo_module import structure_prediction
from bokeh.io import curdoc
from bokeh.layouts import row, layout
from bokeh.models.widgets import Div, Select

predictlayout=structure_prediction().layout

# predictlayout = layout([header],*structlayout)

curdoc().add_root(predictlayout)
curdoc().title = chr(128540) + 'Foldojo'
