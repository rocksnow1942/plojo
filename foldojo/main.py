from foldojo_module import structure_prediction
from bokeh.io import curdoc
from bokeh.layouts import row, layout
from bokeh.models.widgets import Div, Select

structlayout=structure_prediction().layout

temp = Div(text="""
<style>
h1 {
  position: relative;
  animation: mymove 2s;
  animation-iteration-count: infinite;
}

@keyframes mymove {
0%  { border:4px outset #81F7F3;             }
20%  { border:4px outset #81F7F3;             }
50%  { border:4px outset #ff66ff;        }
70%  { border:4px outset #81F7F3;             }
100%  { border:4px outset #81F7F3;      }
}
</style>
<h1 style="width:1050px;height:50px;border: 4px outset #81F7F3;text-align:center;font-family:cursive;font-size:230%;color:#FF00BF;background-color:#81F7F3"">
&#128540
<span style="color:#0000FF">F</span>
<span style="color:red">O</span>
<span style="color:#FFFF00">L</span>
<span style="color:#31B404">D</span>
<span style="color:#FF00BF">ojo </span>
&#129322
</h1>
""",width=1000,height=40)
predictlayout = layout([temp],*structlayout)

curdoc().add_root(predictlayout)
