import pandas as pd
import numpy as np
from astropy.time import Time
import plotly.graph_objects as go
import chart_studio.plotly as py
import random

beta = 1
gamma = 0.3
N = 1000
I0 = round(N*0.05)
T = 15
N_runs =100
delta_t = 0.5
resolution = 0.2
def stoach_mean(beta,gamma,N,I0,T,N_runs,resolution):
    plotter = 0
    stoach_df = stoach_model(beta,gamma,N,I0,T,plotter)
    for i in range (0,N_runs -1):
        next_df = stoach_model(beta,gamma,N,I0,T,plotter)
        stoach_df = pd.concat((stoach_df,next_df))

    stoach_sort_df = sort(stoach_df,resolution)
    return(stoach_sort_df)

def sort(data,n):
    minimum = np.floor((data.t.min()*100))/100
    maximum = np.ceil((data.t.max()*100))/100
    resolution = n
    bins = np.arange(minimum,maximum, resolution)
    bin_labels = bins[:-1] +(n/2)
    data.t = pd.cut(x=data.t, bins=bins, labels = bin_labels)
    mean = (data.groupby(['t']).mean())
    std = (data.groupby(['t']).std())
    data_1 = {'t': bin_labels,
              's' : mean.s.values,
              's_error' : std.s.values,
              'i' : mean.i.values,
              'i_error' : std.i.values,
              'r' : mean.r.values,
              'r_error' : std.r.values}
    bdf = pd.DataFrame(data=data_1)
    return(bdf)


def stoach_mean_plot(beta,gamma,N,I0,delta_t,T,N_runs):
    plotter = 0
    
    stoach_df = stoach_model(beta,gamma,N,I0,T,plotter)
    stoach_mean_df = stoach_mean(beta,gamma,N,I0,T,N_runs,delta_t)
    
    
    t = stoach_mean_df.t.to_list()
    
    stoach_mean_df.i_error = stoach_mean_df.i_error/N
    i_upper = stoach_mean_df.i/N + stoach_mean_df.i_error
    i_lower = stoach_mean_df.i/N - stoach_mean_df.i_error
    
    stoach_mean_df.s_error = stoach_mean_df.s_error/N
    s_upper = stoach_mean_df.s/N + stoach_mean_df.s_error
    s_lower = stoach_mean_df.s/N - stoach_mean_df.s_error
    
    stoach_mean_df.r_error = stoach_mean_df.r_error/N
    r_upper = stoach_mean_df.r/N + stoach_mean_df.r_error
    r_lower = stoach_mean_df.r/N - stoach_mean_df.r_error
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=stoach_df.t, y=stoach_df.i/N,
        hoverinfo='x+y',
        mode='lines',
        name = 'Infected',
        line=dict(width=1, color='hotpink'),
    ))

    fig.add_trace(go.Scatter(
        x=stoach_mean_df.t, y=stoach_mean_df.i/N,
        hoverinfo='x+y',
        mode='lines',
        name = 'Infected averaged 100 repeats',
        line=dict(width=1, color='red'),
    ))
    
    fig.add_trace(go.Scatter(
        x=t+t[::-1], 
        y=i_upper.to_list()+i_lower.to_list()[::-1],
        fill='toself',
        fillcolor='rgba(10,10,10,0.05)',
        line=dict(color='rgba(255,255,255,0)'),
        hoverinfo="skip",
        showlegend=False
    ))

    fig.update_layout(
    title="Applying Monte Carlo Method <br>to SIR stochastic simulation",
    xaxis_title="Time (days)",
    yaxis_title="Fraction of population",
    
    font=dict(
        family="Arial Black",
        size=20,  # Set the font size here
    ))
    fig.show()
    print(stoach_mean_df.i_error.mean())
stoach_mean_plot(beta,gamma,N,I0,delta_t,T,N_runs)
