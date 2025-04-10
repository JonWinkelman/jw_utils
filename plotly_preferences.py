from plotly import graph_objects as go
import plotly.express as px
import numpy as np


colors = ['rgba(132,137,145,1)',
'rgb(0,208,132)',
'rgb(0,122,255)',
'rgb(171,184,195)',
'rgb(255,105,0)',
'rgb(252,185,0)',
'rgb(123,220,181)',
'rgb(123,220,181)',
'rgb(142,209,252)',
'rgb(6,147,227)',
'rgb(155,81,224)',
'rgb(207,46,46)',
'rgb(247,141,167)',
]

def rgba_to_hex(r, g, b, a=1.0):
    a_int = round(a * 255)
    return f'#{r:02X}{g:02X}{b:02X}{a_int:02X}'

def rgb_to_hex(rgb_value = 'rgb(100,200,255)'):
    r,g,b = rgb_value.strip('rgb(').strip(')').split(',')
    r=int(r)
    g=int(g)
    b=int(b)
    return f'#{r:02X}{g:02X}{b:02X}'


def interpolate_colorscale(colorscale_name='earth', n=20, return_swatch_fig=False):
    """
    Interpolates a Plotly colorscale to have `n` equally spaced colors.
    
    Args:
        colorscale (list): Original Plotly colorscale.
        n (int): Number of points in the new colorscale.

        *use px.colors.named_colorscales() to list all colorscale name options    
    Returns:
        list: Interpolated colorscale with `n` entries.
    """
    colorscale = px.colors.get_colorscale(colorscale_name)
    colors = px.colors.sample_colorscale(colorscale, samplepoints=n, low=0, high=1)
    return colors 


def plot_colorscale(colorscale_name, n=20):
    """Generate a visualization of a Plotly colorscale.
    
    *use px.colors.named_colorscales() to list all colorscale name options    
    
    """
    colorscale = px.colors.sample_colorscale(px.colors.get_colorscale(colorscale_name), 
                                             [i / (n - 1) for i in range(n)])

    # Create a heatmap with one column and multiple rows for the colorscale
    fig = go.Figure(go.Heatmap(
        z=np.linspace(0, 1, n).reshape(n, 1),  # Gradient values from 0 to 1
        x=[0],  # Dummy x-axis
        y=np.linspace(1, 0, n),  # Reverse y-axis for top-to-bottom display
        colorscale=colorscale,
        showscale=False  # Hide colorbar
    ))

    fig.update_layout(
        title=f"Colorscale: {colorscale_name}",
        yaxis=dict(showticklabels=False),
        xaxis=dict(showticklabels=False)
    )

    return fig




def bargraph_layout(title, x_title, y_title, height=800, width=800, plot_bgcolor ='rgb(255,255,255)',
                    legend_fontsize=30, yaxis_fontsize=25, xaxis_fontsize = 25):
    return go.Layout(title ={'text':title,
                              'font':{'size':35}},
                      plot_bgcolor = plot_bgcolor,
                       yaxis = {'title':{'text':y_title,
                                        'font':{'size':yaxis_fontsize},
                                        },
                                'showline':True,
                                'linewidth':2,
                                'linecolor':'rgb(0,0,0)',
                                'ticks':'outside',
                                'tickwidth':2,
                                'tickfont':{'size':15}

                               },
                          xaxis = {'title':{'text':x_title, 
                                        'font':{'size':xaxis_fontsize},
                                        },
                                'showline':True,
                                'linewidth':2,
                                'linecolor':'rgb(0,0,0)',
                                'ticks':'outside',
                                'tickwidth':2,
                                'tickfont':{'size':15}

                               },
                       height = height,
                       width = width,
                       legend = {'font':{'size':legend_fontsize}}

                      )
                      
                      
                      
def scatter_layout(title=None, 
                   title_fsize=20,
                   x_title=None, 
                   x_title_fsize=15,
                   y_title=None,
                   y_title_fsize=15, 
                   height=600, width=1100, 
                   legend_fontsize=15):
    return go.Layout(title ={'text':title,
                              'font':{'size':title_fsize}},
                      plot_bgcolor = 'rgb(255,255,255)',
                       yaxis = {'title':{'text':y_title,
                                        'font':{'size':y_title_fsize},
                                        },
                                'showline':True,
                                'linewidth':2,
                                'linecolor':'rgb(0,0,0)',
                                'ticks':'outside',
                                'tickwidth':2,
                                'tickfont':{'size':15}

                               },
                          xaxis = {'title':{'text':x_title, 
                                        'font':{'size':x_title_fsize},
                                        },
                                'showline':True,
                                'linewidth':2,
                                'linecolor':'rgb(0,0,0)',
                                'ticks':'outside',
                                'tickwidth':2,
                                'tickfont':{'size':15}

                               },
                       height = height,
                       width = width,
                       legend = {'font':{'size':legend_fontsize}}

                      )