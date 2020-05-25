# -*- coding: utf-8 -*-
"""
Created on Sat May  2 18:28:14 2020

@author: paco
"""

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