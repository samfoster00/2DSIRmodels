import pandas as pd
import numpy as np
from astropy.time import Time
import plotly.graph_objects as go
import chart_studio.plotly as py
import random

def Eulers(beta,gamma,N,I_0,T,plotter):
    S_0 = N - I_0
    R_0 = 0


    delta_t= 0.01
    h = round(T/delta_t)


    S = np.empty(h)
    I = np.empty(h)
    R = np.empty(h)
    t = np.empty(h)

    S[0] = S_0
    I[0] = I_0 
    R[0] = R_0
    t[0] = 0

    for n in range (1,h):
        S[n] = S[n-1] - beta*S[n-1]*I[n-1]*(delta_t/N)
        I[n] = I[n-1] + ((beta*I[n-1]*S[n-1])/N -gamma*I[n-1])*delta_t
        R[n] = R[n-1] + gamma*I[n-1]*delta_t
        t[n] = n*delta_t

    df = pd.DataFrame(list(zip(t, S,I,R)),
                   columns =['t', 's','i','r'])

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df.t, y=df.i,
        hoverinfo='x+y',
        mode='lines',
        name = 'Infected',
        line=dict(width=1, color='red'),
        stackgroup='one'
    ))
    fig.add_trace(go.Scatter(
        x=df.t, y=df.r,
        hoverinfo='x+y',
        name = 'Recovered',
        mode='lines',
        line=dict(width=1, color='blue'),
        stackgroup='one'

    ))
    fig.add_trace(go.Scatter(
        x=df.t, y=df.s,
        hoverinfo='x+y',
        mode='lines',
        name = 'Susceptible',
        line=dict(width=1, color='green'),
        stackgroup='one'

    ))

    fig.update_layout(yaxis_range=(0, N))
    if plotter == 1:
        fig.show()
    else:
        return(df)
Eulers(1,0.3,1000,10,40,1)
