import pandas as pd
import numpy as np
from astropy.time import Time
import plotly.graph_objects as go
import chart_studio.plotly as py
import random

def stoach_model(beta,gamma,N,I0,total_T,plotter):

    S = [N-I0]
    I = [I0]
    R = [0]
    T = [0]
    while (T[-1] < total_T) and (I[-1]>0):
        v = S[-1],I[-1],R[-1]
        v_new = stoch(v,beta,gamma,N)
        S.append(v_new[0])
        I.append(v_new[1])
        R.append(v_new[2])
        T.append(v_new[3] + T[-1])
    df = pd.DataFrame(list(zip(S, I,R,T)),
               columns =['s','i','r','t'])
    if plotter == 1:
        plot_stoach(df,N)
    else:
        return (df)

def stoch(v,beta,gamma,N):
    u_1,u_2 = np.random.uniform(0,1,2)
    s,i,r = v
    lambd = s*i*beta/N +gamma*i
    
    dt = tau(u_1,lambd)
    
    p_1 = (s*i*beta/N)/lambd
    p_2 = (gamma*i)/lambd
    
    if u_2 < p_1:
        s -= 1
        i += 1
    else:
        i -= 1
    r = N - s - i
    return s,i,r,dt

def tau(u_1,lambd):
    return -np.log(u_1)/lambd

def plot_stoach(df,N):
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df.t, y=df.i/N,
        hoverinfo='x+y',
        mode='lines',
        name = 'Infected',
        line=dict(width=1, color='red'),

    ))
    fig.add_trace(go.Scatter(
        x=df.t, y=df.r/N,
        hoverinfo='x+y',
        name = 'Recovered',
        mode='lines',
        line=dict(width=1, color='blue'),


    ))

    fig.add_trace(go.Scatter(
        x=df.t, y=df.s/N,
        hoverinfo='x+y',
        mode='lines',
        name = 'Susceptible',
        line=dict(width=1, color='green'),
        


    ))

 
    fig.update_layout(yaxis_range=(0, 1))
    fig.update_layout(
        title="Stochastic simulation of the SIR model",
        xaxis_title="Time (days)",
        yaxis_title="Fraction of the population",
        legend_title="Compartment",
    )
    fig.show()

stoach_model(1,0.3,100,5,15,1)
