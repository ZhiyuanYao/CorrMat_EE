#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
# itertools handles the color cycling
import itertools
import bokeh.io
# the high-level bokeh.plotting interface
from bokeh.plotting import ColumnDataSource, figure, output_file, show, save
# select a palette, the _10 or other number specify how number used, necessary
from bokeh.palettes import Category10_10 as palette
import smart
from bokeh.models import Label
from bokeh.util.compiler import TypeScript

L = 20
# ----------------------------------------------------------------------------#
#                         Bokeh basic setting                                 #
# ----------------------------------------------------------------------------#
# create a new plot
p = figure(
   tools="pan,wheel_zoom,box_zoom,reset,save",plot_width=1000,
   x_axis_label='E_n', y_axis_label='S_z', title="E-Sz plot for L=" + str(L) +
   " with (k, I) = 0+"
)
p.output_backend = "svg"

# ----------------------------------------------------------------------------#
#                         Read in the data                                    #
# ----------------------------------------------------------------------------#
filename = f"../data/L20_E_sz.dat"
xdat, ydat = smart.columndata(filename, [0,1])
source = ColumnDataSource(data=dict(
    x=xdat,
    y=ydat,
    ieigen=list(range(1, 1 + len(xdat)))
))
# Add the hover tool with annotation of value of data points
tooltips = [
    ("index", "@ieigen"),
    ("(x,y)", "(@x{0.2f}, @y{0.2f})"),
]
p.add_tools(bokeh.models.HoverTool(tooltips=tooltips))

# create a color iterator, use line_color=color if color cycle needed
# colors = itertools.cycle(palette)
# p.line(xdat, ydat, line_width=2, line_alpha=1.0, line_color='blue', legend="line")
p.circle('x', 'y', size=7, alpha=0.5, color='green', legend="E-s_z",
        source=source)

# k= 0 sector data
# filename = "../symmetry_method/E_szL" + str(L) + ".dat"
# xdat, ydat = smart.columndata(filename, [0,1])
# p.circle(xdat, ydat, size=5, alpha=0.9, color='blue', legend="k=0")

p.legend.location = "top_left"
p.legend.click_policy="hide"
# output to static HTML file
output_file(filename + ".html")
save(p)
show(p)
