from plotly import graph_objects as go

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
                      
                      
                      
def scatter_layout(title=None, x_title=None, y_title=None, height=800, width=800):
    return go.Layout(title ={'text':title,
                              'font':{'size':35}},
                      plot_bgcolor = 'rgb(255,255,255)',
                       yaxis = {'title':{'text':y_title,
                                        'font':{'size':25},
                                        },
                                'showline':True,
                                'linewidth':2,
                                'linecolor':'rgb(0,0,0)',
                                'ticks':'outside',
                                'tickwidth':2,
                                'tickfont':{'size':15}

                               },
                          xaxis = {'title':{'text':x_title, 
                                        'font':{'size':25},
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
                       legend = {'font':{'size':30}}

                      )



