import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import netCDF4
import os
from ..utils.constants import *

# function to plot the satellite orbit in 3D space with the earth at the origin 

def plot_groundtrack(Coords, res=(800, 800), Mode='markers+lines', proj_type='equirectangular'): 

    Lat = Coords[:, 0]
    Lon = Coords[:, 1]
    types = ['orthographic', 'natural earth', 'kavrayskiy7', 'miller', 'robinson', 'equirectangular', 'mercator', 'orthographic', 'natural earth', 'kavrayskiy7', 'miller', 'robinson', 'equirectangular', 'mercator']
    if proj_type not in types:
        raise ValueError("Plot type must be picked from the plotly list of projections: \n" + str(types))
    
    fig = make_subplots(rows=1, cols=2, subplot_titles=("Orthographic Projection", "Flat Projection"))

    
    # Define the data for your plot
    data_ortho = [
        go.Scattergeo(
            lat=Lat, lon=Lon, mode=Mode, marker=dict(size=0, color='red'),
            hoverinfo="none"
        )
    ]

    data_flat = [
        go.Scattergeo(
            lat=Lat, lon=Lon, mode=Mode, marker=dict(size=0, color='red'),
            hoverinfo="none"
        )
    ]

    # Define the layout for your plot
    layout_ortho = go.Layout(
        geo=dict(
            showland=True, landcolor="Green",
            showocean=True, oceancolor="LightBlue",
            showcountries=False, showcoastlines=True,
            projection=dict(type=type, rotation=dict(lat=0, lon=0, roll=0)),
            bgcolor="black"
        ),
        width=res[0]//2, height=res[1], 
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectratio=dict(x=1, y=1, z=1), 
            domain=dict(x=[0, 0.6], y=[0, 1]))
    )

    layout_flat = go.Layout(
        xaxis=dict(visible=True, range=[-180, 180]),
        yaxis=dict(visible=True, range=[-90, 90]),
        width=res[0]//2, height=res[1],
        margin=dict(l=0, r=0, t=40, b=0),
        plot_bgcolor="black",
        paper_bgcolor="black",
    )
    
    # Add the traces and layouts to the figure
    fig.add_trace(data_ortho[0], row=1, col=1)
    fig.add_trace(data_flat[0], row=1, col=2)
    fig.update_layout(layout_ortho, row=1, col=1)
    fig.update_layout(layout_flat, row=1, col=2)
    
    # Update the flat projection subplot to use the chosen projection
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    fig.update_xaxes(scaleanchor="y", scaleratio=1)
    fig.update_layout(scene=dict(projection=dict(type=proj_type)))

    fig.show()


