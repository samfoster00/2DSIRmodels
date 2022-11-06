import pandas as pd
import numpy as np
from astropy.time import Time
import plotly.graph_objects as go
import chart_studio.plotly as py
import random

def twoD_trans(beta,gamma,N,I0,T,contacts,trans_prob,plotter):
    prob = beta/contacts
    ratio = (contacts+1)/25
    total_N = round(N/ratio)
    
    N_root = round(total_N**(1/2))
    I = I0
    S = N-I0
    df = create_map(N_root,S,I,total_N)
    data_h=[df_to_plotly(df)]
    model_df,frames = modeltrans(df,data_h,T,N_root,prob,gamma,contacts,trans_prob)
    if plotter == 1:
        animation(df,data_h,T)
        plot_model(model_df)
    else:
        return model_df
    

def df_to_plotly(df):
    return {'z': df.values.tolist(),
            'x': df.columns.tolist(),
            'y': df.index.tolist()}

def create_map(N_root,S,I,total_N):
    df = pd.DataFrame(np.random.randint(0,total_N,size=(N_root, N_root)))
    df = df.replace([range(0,I)],2)
    df = df.replace([range(I,S+I)],3)
    df = df.replace([range(S+I,total_N)],0)
    return(df)

def infecttrans(df,N_root,prob,gamma,contacts):
    i_list = [(df.index[x], df.columns[y]) for x, y in zip(*np.where(df.values == 2))]
    for n in range(0,len(i_list)):
        for i in [-2,-1,0,1,2]:
            for j in [-2,-1,0,1,2]:
                trail_x,trail_y = ([i_list[n][1]+i,i_list[n][0]+j])
                invalid = [-2,-1,N_root,N_root+1]
                if not(any(elem in [trail_x,trail_y] for elem in invalid)):
                    if (df.iat[trail_y,trail_x]) == 3:
                        if i == 1 or j == 1:
                            if random.random() < prob*3/2:
                                df.iat[trail_y, trail_x] = 2
                        else:
                            if random.random() < prob*3/4:
                                df.iat[trail_y, trail_x] = 2
        if random.random() < gamma:
            df.iat[i_list[n][0],i_list[n][1]] = 1
            
    
             
    return df

def transport(df,trans_prob):
    t_list = [(df.index[x], df.columns[y]) for x, y in zip(*np.where(df.values != 0))]
    for n in range(0,len(t_list)-1):
        if random.random() < trans_prob:
            point = t_list[n]
            point_value = df.iat[t_list[n][0],t_list[n][1]]
            
            select = random.randint(0, len(t_list)-1)
            trans_point = t_list[select]
            trans_value = df.iat[t_list[select][0],t_list[select][1]]
            df.iat[point[0],point[1]] = trans_value
            df.iat[trans_point[0],trans_point[1]] = point_value
    return df

def plot(df):
    fig= go.Figure(data=go.Heatmap(df_to_plotly(df), colorscale=["black", "blue","red", "green"]))
    fig.show()

def modeltrans(df,data_h,T,N_root,prob,gamma,contacts,trans_prob):
    s = np.empty(T)
    i = np.empty(T)
    r = np.empty(T)
    t = np.empty(T)
    
    count = df.stack().value_counts()
    s[0] = count[3]
    i[0] = count[2]
    r[0] = 0
    t[0] = 0
    for dt in range(1,T):
        df = infecttrans(df,N_root,prob,gamma,contacts)
        df = transport(df,trans_prob)
        data_h.append(df_to_plotly(df))
        count = df.stack().value_counts()
        try: 
            s[dt] = count[3]
        except:
            s[dt] = 0
        try: 
            i[dt] = count[2]
        except:
            i[dt] = 0
        try: 
            r[dt] = count[1]
        except:
            r[dt] = 0
        t[dt] = t[dt-1] +1
    model_df = pd.DataFrame(list(zip(t, s,i,r)),
                            columns =['t', 's','i','r'])
    return(model_df,data_h)
    
def plot_model(model_df):
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=model_df.t, y=model_df.i,
        hoverinfo='x+y',
        mode='lines',
        name = 'Infected',
        line=dict(width=1, color='red'),
        stackgroup='one'
    ))
    fig.add_trace(go.Scatter(
        x=model_df.t, y=model_df.r,
        hoverinfo='x+y',
        name = 'Recovered',
        mode='lines',
        line=dict(width=1, color='blue'),
        stackgroup='one'

    ))
    fig.add_trace(go.Scatter(
        x=model_df.t, y=model_df.s,
        hoverinfo='x+y',
        mode='lines',
        name = 'Susceptible',
        line=dict(width=1, color='green'),
        stackgroup='one'

    ))
    fig.show()

def animation(df,data_h,T):
    fig = go.Figure(data=go.Heatmap( df_to_plotly(df), colorscale=["black", "blue","red", "green"], zmin=0, zmax=3),
                    frames=[go.Frame(data = go.Heatmap(data_h[i], colorscale=["black", "blue","red", "green"], zmin=0, zmax=3))for i in range(T)])
    fig.update_layout(
    updatemenus=[
        dict(type="buttons", visible=True,
        buttons=[dict(label="Play", method="animate", args=[None])]
            )])
    fig.show()


twoD_trans(1,0.3,10000,10,50,24,0.01,1)
