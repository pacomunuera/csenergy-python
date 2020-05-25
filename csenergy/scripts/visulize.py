# -*- coding: utf-8 -*-
"""
Created on Sat May  2 11:05:52 2020

@author: paco
"""

import seaborn as sns
import pandas as pd
import sys


# df = pd.read_csv('./simulations_outputs/20200502 025558LS3_1.csv',
#                  decimal = ',',
#                  index_col=0)

# strfilename = sys.argv
strfilename = './simulations_outputs/20200502 025558LS3_1.csv'
df = pd.read_csv(strfilename,
                 decimal = ',',
                 sep =  ';',
                 index_col=0)

df.index = pd.to_datetime(df.index, infer_datetime_format=True)
df_filtrado = df[['DNI', 'REQ_MASSFLOW']]
g = sns.pairplot(df_filtrado)


"""
ts
I am trying to plot using Seaborn in Tkinter. My approaches so far were different variations of this and I could not get it to work.

I tried the matplotlib.use("Agg"), which works fine one the normal Matplotlib graphs on the page but doesn't seem to work on the Seaborn plots

matplotlib.use("TkAgg") # 'Agg' doesnt work either
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import seaborn as sns
import tkinter as tk

def graphspage():
    pf = tk.Tk()
    pf.geometry("1000x800")


### Works
    f = Figure(figsize=(5, 5), dpi=100)
    a = f.add_subplot(111)
    a.plot(df['date'],daily_drawdown, 'b-',df['date'], daily_drawdownbm, 'k-', linewidth=1)
    f.tight_layout()

    canvas = FigureCanvasTkAgg(f,pf)
    canvas.get_tk_widget().grid(row=1,column=1)

### Doesnt Work
    pct = diststats()[4]
    pctbm = diststats()[5]
    f = Figure(figsize=(5, 5), dpi=100)
    a = f.add_subplot(111)
    a.sns.distplot(pct,label = 'Portfolio')
    a.sns.distplot(pctbm,axlabel='Distribution of returns',label='Benchmark')

    canvas = FigureCanvasTkAgg(f,pf)
    canvas.get_tk_widget().grid(row=2,column=1)

graphspage()
python matplotlib tkinter seaborn sis
shareeditfollow
asked Apr 14 '19 at 23:06

Joan Arau
3588 bronze badges
did you see any example which use Seaborn in tkinter ? – furas Apr 14 '19 at 23:30
1
Instead of a.sns.distplot(pct,label = 'Portfolio') it needs to be sns.distplot(pct,label = 'Portfolio', ax=a) – ImportanceOfBeingErnest Apr 14 '19 at 23:36
1
create code which we could run. Now there is some df and diststats() so we can't run it. – furas Apr 14 '19 at 23:45
add a comment
2 Answers
Active
Oldest
Votes
¿No encuentras la respuesta? Pregunta en Stack Overflow en español.

✕

2

The OP's original code only lacked a canvas.draw(), if I'm not mistaken. This has also been indicated by furas. I recently found it difficult to find a full example of how to draw with Seaborn on a Tkinter GUI and especially, how to redraw on the canvas.

So, let me give you a fully working but mininmal snippet for a program that initially draws and redraws on every keypress.

import tkinter
from typing import Callable

import numpy as np
import seaborn as sns
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


def init_gui(root, update_function: Callable) -> FigureCanvasTkAgg:
    def event_key_press(event):
        print("you pressed {}".format(event.key))
        update_function()
        key_press_handler(event, canvas)

    # create empty figure and draw
    init_figure = create_figure()
    canvas = FigureCanvasTkAgg(init_figure, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
    # call key press event
    canvas.mpl_connect("key_press_event", event_key_press)
    return canvas


def create_figure() -> Figure:
    # generate some data
    matrix = np.random.randint(20, size=(10, 10))
    # plot the data
    figure = Figure(figsize=(6, 6))
    ax = figure.subplots()
    sns.heatmap(matrix, square=True, cbar=False, ax=ax)
    return figure


def redraw_figure():
    figure = create_figure()
    canvas.figure = figure
    canvas.draw()


sns.set()
root = tkinter.Tk()
canvas = init_gui(root, update_function=redraw_figure)

tkinter.mainloop()
"""