import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import netCDF4
import os
from ..utils.constants import *

# function to plot the satellite orbit in 3D space with the earth at the origin 

def plot_groundtrack(Sat_ID, Coords, res=(800, 800), Mode='markers+lines', style='seaborn', 
                     marker_size=0.5, tracecolor='red', play_duration=300): 
    Lat = Coords[:, 0]
    Lon = Coords[:, 1]

    # converting latitude and longitude to equirectangular projection
    Lat_equi = 800 * np.cos(np.deg2rad(Lat)) * np.cos(np.deg2rad(Lon))
    Lon_equi = 400 * np.cos(np.deg2rad(Lat)) * np.sin(np.deg2rad(Lon))


    animation_steps = 40
    
    fig = make_subplots(rows=1, cols=2, subplot_titles=("Orthographic Projection", "Flat Projection"), specs=[[{"type": "geo"}, {"type": "geo"}]], 
                    print_grid=True)

    fig.add_trace(go.Scattergeo(lat=Lat, lon=Lon, mode=Mode, marker=dict(size=marker_size, color=tracecolor), name=Sat_ID), row=1, col=1)
    fig.add_trace(go.Scattergeo(lat=Lat_equi, lon=Lon_equi, mode=Mode, marker=dict(size=marker_size, color=tracecolor), showlegend=False), row=1, col=2)

    fig.update_geos(projection_type="orthographic", row=1, col=1, 
                    showland=True, landcolor="LightGreen", showocean=True, oceancolor="#179fff", showlakes=True, lakecolor="#179fff", 
                    showrivers=True, rivercolor="#179fff")
    fig.update_geos(projection_type="equirectangular", row=1, col=2, 
                    showland=True, landcolor="LightGreen", showocean=True, oceancolor="#179fff", showlakes=True, lakecolor="#179fff", 
                    showrivers=True, rivercolor="#179fff")
    fig.update_layout(template=style, width=res[0], height=res[1]//2, title_text="Groundtrack", showlegend=True)

    frames = []
    for i in range(1, animation_steps + 1):
        alpha = i / animation_steps
        frames.append(go.Frame(data=[go.Scattergeo(lat=Lat[:int(alpha*len(Lat))], lon=Lon[:int(alpha*len(Lon))],
                                                    mode=Mode, marker=dict(size=marker_size, color=tracecolor),
                                                    name=Sat_ID)],
                               layout=dict(geo=dict(projection_type='orthographic'))))

    fig.frames = frames
    fig.update_layout(updatemenus=[dict(type='buttons', showactive=False, buttons=[dict(label='Play',
                                                                                      method='animate',
                                                                                      args=[None, dict(frame=dict(duration=play_duration),
                                                                                                        fromcurrent=True,
                                                                                                        mode="immediate",
                                                                                                        transition=dict(duration=0))]),
                                                                                     dict(label='Pause',
                                                                                          method='animate',
                                                                                          args=[[None], dict(frame=dict(duration=0),
                                                                                                              mode='immediate',
                                                                                                              transition=dict(duration=0))])])])
    fig.show()