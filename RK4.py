import pandas as pd
import numpy as np
from astropy.time import Time
import plotly.graph_objects as go
import chart_studio.plotly as py
import random

def rk4(beta,gamma,N,I0,T,plotter):
    S0 = N - I0
    R0 = 0
    
    delta_t = 0.01
    h=round(T/delta_t)

    S,I,R,t = construct_lists(S0,I0,R0,h)
    S,I,R,t = iterate_rk4(S,I,R,t,h,beta,gamma,N,delta_t)
    df = pd.DataFrame(list(zip(t, S,I,R)),
               columns =['t', 's','i','r'])
    if plotter == 1:
        plot_rk4(df,N)
    else:
        return df

def plot_rk4(df,N):
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
    fig.update_layout(
    font=dict(
        family="Arial Black",
        size=20,  # Set the font size here
    ))
    fig.update_layout(yaxis_range=(0, N))
    fig.show()

def iterate_rk4(S,I,R,t,h,beta,gamma,N,delta_t):
    for i in range (0,h-1):
        Si = S[i]
        Ii = I[i]
        
        S_k1 = dSdt(Si,Ii,beta,N)
        I_k1 = dIdt(Si,Ii,beta,gamma,N)
        
        S_k2 = dSdt(Si + delta_t / 2 * S_k1, Ii + delta_t / 2 * I_k1,beta,N)
        I_k2 = dIdt(Si + delta_t / 2 * S_k1, Ii + delta_t / 2 * I_k1,beta,gamma,N)
        
        S_k3 = dSdt(Si + delta_t / 2 * S_k2, Ii + delta_t / 2 * I_k2,beta,N)
        I_k3 = dIdt(Si + delta_t / 2 * S_k2, Ii + delta_t / 2 * I_k2,beta,gamma,N)
        
        S_k4 = dSdt( Si + delta_t * S_k3, Ii + delta_t * I_k3,beta,N)
        I_k4 = dIdt( Si + delta_t * S_k3, Ii + delta_t * I_k3,beta,gamma,N)
        
        S[i + 1] = Si + delta_t / 6 * (S_k1 + 2 * S_k2 + 2 * S_k3 + S_k4)
        I[i + 1] = Ii + delta_t / 6 * (I_k1 + 2 * I_k2 + 2 * I_k3 + I_k4)
        t[i + 1] = t[i]+delta_t
        R[i + 1] = N - I[i+1] - S[i+1]
    return(S,I,R,t)
    
def construct_lists(S0,I0,R0,h):
    S = np.empty(h)
    I = np.empty(h)
    R = np.empty(h)
    t = np.empty(h)
    
    S[0] = S0
    I[0] = I0 
    R[0] = R0
    t[0] = 0
    
    return(S,I,R,t)

def dSdt(S,I,beta,N):
    return(-beta * S * I / N)
def dIdt(S,I,beta,gamma,N):
    return(beta * S * I / N - gamma * I)
def dRdt(I):
    return (gamma* I)

rk4(1,0.3,1000,10,40,1)
