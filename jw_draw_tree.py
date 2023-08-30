from Bio import Phylo
from plotly import graph_objects as go
import pandas as pd

def get_x_coordinates(tree):
    """Associates to  each clade an x-coord.
       returns dict {clade: x-coord}
    """
    xcoords = tree.depths()
    # tree.depth() maps tree clades to depths (by branch length).
    # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth
    # is the distance from root to clade

    #  If there are no branch lengths, assign unit branch lengths
    if not max(xcoords.values()):
        xcoords = tree.depths(unit_branch_lengths=True)
    return xcoords


def get_y_coordinates(tree, dist=1.3):
    """
       returns  dict {clade: y-coord}
       The y-coordinates are  (float) multiple of integers (i*dist below)
       dist depends on the number of tree leafs
    """
    maxheight = tree.count_terminals()  # Counts the number of tree leafs.
    # Rows are defined by the tips/leafs
    ycoords = dict((leaf, maxheight - i * dist) for i, leaf in enumerate(reversed(tree.get_terminals())))

    def calc_row(clade):
        for subclade in clade:
            if subclade not in ycoords:
                calc_row(subclade)
        ycoords[clade] = (ycoords[clade.clades[0]] +
                          ycoords[clade.clades[-1]]) / 2

    if tree.root.clades:
        calc_row(tree.root)
    return ycoords


def get_clade_lines(orientation='horizontal', y_curr=0, x_start=0, x_curr=0, y_bot=0, y_top=0,
                    line_color='rgb(25,25,25)', line_width=0.5):
    """define a shape of type 'line', for branch
    """
    branch_line = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )
    if orientation == 'horizontal':
        branch_line.update(x0=x_start,
                           y0=y_curr,
                           x1=x_curr,
                           y1=y_curr)
    elif orientation == 'vertical':
        branch_line.update(x0=x_curr,
                           y0=y_bot,
                           x1=x_curr,
                           y1=y_top)
    else:
        raise ValueError("Line type can be 'horizontal' or 'vertical'")

    return branch_line



def draw_clade(clade, x_start, line_shapes, line_color='rgb(15,15,15)', line_width=1, x_coords=0, y_coords=0, 
               cl_to_highlight = None, highlight_line_width=2):
    """Recursively draw the tree branches, down from the given clade"""

    x_curr = x_coords[clade]
    y_curr = y_coords[clade]
    if cl_to_highlight:
        clade_names_to_highlight= [cl.name for cl in get_x_coordinates(cl_to_highlight)]
        if clade.name in clade_names_to_highlight:
            line_width = highlight_line_width
    # Draw a horizontal line from start to here
    branch_line = get_clade_lines(orientation='horizontal', y_curr=y_curr, x_start=x_start, x_curr=x_curr,
                                  line_color=line_color, line_width=line_width)

    line_shapes.append(branch_line)

    if clade.clades:
        # Draw a vertical line connecting all children
        y_top = y_coords[clade.clades[0]]
        y_bot = y_coords[clade.clades[-1]]

        line_shapes.append(get_clade_lines(orientation='vertical', x_curr=x_curr, y_bot=y_bot, y_top=y_top,
                                           line_color=line_color, line_width=line_width))

        # Draw descendants
        for child in clade:
            draw_clade(child, x_curr, line_shapes, x_coords=x_coords, y_coords=y_coords, cl_to_highlight=cl_to_highlight, highlight_line_width = highlight_line_width)



def read_treefile(filename):
    'create tree object from newick format using Bio.Phylo'
    tree = Phylo.read(filename, "newick")
    return tree



def create_tree_w_bargraphs(tree_obj, data_dict, colors=None, intern_node_label=None, node_color_dict=None,
                                in_node_size=5, t_node_size=10, labels=None, label_mode='markers', height=850,
                                node_size_dict=None, cl_to_highlight=None, highlight_line_width=1.5, hover_text=None):
    """Create a phylogenetic plotly tree figure with 1 or more bargraphs
    
    Parameters:
    tree_obj (Bio.Phylo newick tree object):

    colors (dict): Colors for the bars in each bargraph. Format = {'t_blue': 'rgba(0,102,153,255)', ...}
    intern_node_label (str): These are the clade attribute that will be assigned to the internal node hoverdata,
                            options are: 'name', 'confidence', 'branch_length'. If None, an empty string is assigned
                            to each internal node.
    data_dict (nested dict): {title:{leaf_name:int}} Can put multiple bargraphs onto tree
    label_mode (str): options 'markers', 'markers+text'. 'markers' - only show text on hover, 'markers+text' - show text always
    cl_to_highlight (Bio.Phylo clade object): Highlight all lines emanating from input clade
    highlight_line_width (int): size of the lines emanating from cl_to_highlight
    """
    if not colors:
        colors = {
        't_blue': 'rgba(0,102,153,255)',
        't_green': 'rgba(61,174,43,255)',
        't_red': 'rgb(255,20,20)',
        'seagreen':'#2c8d42',
        'orange':'#F9B257',
        'purple':'rgba(52, 30, 170, 0.25)',
        'lavender':'rgba(170, 30, 125, 0.25)',
        'greyish': 'rgba(80, 80, 69, 0.25)'}
    
    fig=create_tree(tree_obj, intern_node_label=intern_node_label, node_color_dict=node_color_dict,
                     in_node_size=in_node_size, t_node_size=t_node_size, label_mode=label_mode,height=height,
                     node_size_dict=node_size_dict, cl_to_highlight=cl_to_highlight, highlight_line_width=highlight_line_width,
                     hover_text=hover_text)
    fig_sp_tree = go.Figure(fig)
    shift = 0.01
    bar_thickness = 0.95
    scale = 0.02
    new_max_x = max(fig_sp_tree['data'][0]['x']) #initialize max_x
    for i, graph in enumerate(data_dict):
        l_colors = list(colors.values())
        color = l_colors[i]
        fig_sp_tree_added_traces, new_max_x = _add_bargraph_trace(fig_sp_tree,
                                                    bar_thickness, shift, scale, 
                                                    new_max_x, tree_obj, color, data_dict[graph],
                                                    bargraph_title = graph)
    return fig_sp_tree_added_traces



def _add_bargraph_trace(fig_sp_tree, bar_thickness, shift, 
                      scale, max_x, tree_obj, color, data_dict,
                      bargraph_title = None):
    """Add plotly bar graph trace to the input tree figure"""
    HOG = 'N1.HOG0001809'
    x_y_coords = get_bar_coords_dict(tree_obj, data_dict)
    new_max_x = max_x
    for acc, x_y in x_y_coords.items():
        copies = int(x_y[0])
        text = f'{acc}<br>{bargraph_title} copies: {copies}'
        bar_trace = get_bar_traces(max_x, copies, x_y[1], bar_thickness = bar_thickness, 
                                   shift=shift, scale=scale, text=text, 
                                   color=color)
        fig_sp_tree.add_trace(bar_trace)
        #find new max_x v_val
        if max(bar_trace['x']) > new_max_x:
            new_max_x = max(bar_trace['x'])
    fig_sp_tree.add_trace(make_baseline_trace(shift, x_y_coords, max_x, bar_thickness, title=bargraph_title))
    return fig_sp_tree, new_max_x


def get_bar_coords_dict(sp_tree, data_dict):
    """
    Parameters:
    data_dict (dict): {leaf_name:counts}
    """
    #sp_tree = read_treefile(path_to_species_tree)
    ycoords = get_y_coordinates(sp_tree)
    leaf_names  = [cl.name for cl in sp_tree.get_terminals()]
    ycoords_terminal = {clade.name:y_value for clade, y_value in ycoords.items() if clade.name in leaf_names}
    x_y_coords = {acc:[] for acc in data_dict.keys()}
    for acc, val in ycoords_terminal.items():
        x_y_coords[acc].append(data_dict[acc])
        x_y_coords[acc].append(val)
    return x_y_coords


def get_bar_traces(x_start, x_val, y_val, bar_thickness, shift=0,
                   scale=1, text='', color='rgba(100,100,100,0.5)'):
    """return a shape dict for a bar for a horizontal bargraph
    
    parameters:
    x_start (int): left_most x_coordinate of graph
    x (int or float): lenght of the bar
    y (int or float): y position on bargraph
    bar_thickness (float): thickness of a bar
    shift (float): shifts the left_hand side of graph 
    
    return (trace): trace containing one bar
    """
    x_left = x_start+shift
    x_right = x_left+(x_val*scale)
    y_low = y_val-(0.5*bar_thickness)
    y_high = y_val+(0.5*bar_thickness)
    x = [x_left, x_right, x_right, x_left, x_left]
    y = [y_low,  y_low,   y_high,  y_high, y_low]
    
    trace = go.Scatter(
        x = x,
        y = y,
        mode='lines',
        fill='toself',
        line={'color':color,
             'width':0.01},
        text = text,
        showlegend=False)
    return trace

def make_baseline_trace(shift, x_y_coords, max_x, bar_thickness, title=''):
    x_start = max_x + shift
    xy_coords_df = pd.DataFrame(x_y_coords).transpose()
    xy_coords_df.columns = [title+'_x',title+'_y' ]
    max_y = xy_coords_df[title+'_y'].max() + 0.5*bar_thickness
    min_y = xy_coords_df[title+'_y'].min() - 0.5*bar_thickness
    baseline_trace = go.Scatter(x=[x_start, x_start],
                                y=[min_y, max_y],
                                mode='lines',
                                line={'color':'rgba(100,100,100,0.5)',
                                     'width': 1})
    return baseline_trace



#############base tree functions###########################################################################################
def create_tree(tree, title=None, intern_node_label=None,  node_color_dict=None, in_node_size=6, t_node_size=10,
                 label_mode='markers', height=850, node_size_dict=None, cl_to_highlight=None, highlight_line_width=2,
                 hover_text = None):
    """Return a plotly tree  

    Parameters:  
    tree: Bio.Phylo newick tree object
    title (str): title to appear on graph 
    intern_node_label (str): clade attribute from which to name the internal nodes. 
                            options: None, 'confidence', 'name', 'branch_length'
    node_color_dict (dict): define color for nodes on tree: {node_name:rgb()}. If
                            node name is not in dict it will be given default color
    label_mode (str): options 'markers', 'markers+text'. 'markers' - only show text on hover, 'markers+text' - show text always
    """

    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    draw_clade(tree.root, 0, line_shapes, line_color='rgb(25,25,25)', line_width=1, x_coords=x_coords,
               y_coords=y_coords, cl_to_highlight=cl_to_highlight, highlight_line_width=highlight_line_width)
            
    my_tree_clades = x_coords.keys() #returns all clades in tree, {clade:dist from root to clade}

    X = []
    Y = []
    text = []
    node_sizes = []
    
    #hoverdata for ALL nodes
    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
        if cl in tree.get_terminals():
            if cl.name:
                text.append(cl.name)
            node_sizes.append(t_node_size)
        else:
            if intern_node_label == 'confidence':
                text.append(cl.confidence)
            elif intern_node_label == 'name':   
                text.append(cl.name)
            elif intern_node_label == 'branch_length':
                text.append(cl.branch_length)
            else:
                text.append('') 
            node_sizes.append(in_node_size)
    if hover_text:
        text=hover_text
    
    #make color dict:
    if node_color_dict:
        node_colors = make_tree_node_colorlist(my_tree_clades, node_color_dict)
    else:
        node_colors = []
    if node_size_dict:
        node_sizes = make_tree_node_sizelist(my_tree_clades, node_size_dict, default_size = 6)

    data = dict(type='scatter',
                x=X,
                y=Y,
                mode=label_mode,#'markers',
                textposition='middle right',
                textfont = dict(size=7.5),
                marker=dict(color=node_colors,
                            size=node_sizes),
                text=text,  # vignet information of each node
                hoverinfo='text',
                )
    layout = get_tree_layout(line_shapes=line_shapes,height=height)
    fig = dict(data=[data], layout=layout)
    return fig


def make_tree_node_colorlist(my_tree_clades, node_color_dict):
    node_colors = [] 
    for cl in my_tree_clades:
        color = node_color_dict.get(cl.name)
        if color:
            node_colors.append(color)
        else: 
            node_colors.append('rgb(25,25,25)')
    return node_colors


def make_tree_node_sizelist(my_tree_clades, node_size_dict, default_size = 2):
    node_sizes = [] 
    for cl in my_tree_clades:
        size = node_size_dict.get(cl.name)
        if size:
            node_sizes.append(size)
        else:
            node_sizes.append(default_size)
    return node_sizes



def get_tree_layout(title=None, paper_bgcolor ='rgb(248,248,248)', height=850, line_shapes=None,
                    plot_bgcolor='rgb(248,248,248)'):
    axis = dict(showline=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
            title=''  # y title
            )
    layout = dict(title=title,
                #uirevision='constant',
                paper_bgcolor=paper_bgcolor,
                dragmode="lasso",
                font=dict(family='Balto', size=14),
                #width=750,
                height=height,
                autosize=True,
                showlegend=False,
                xaxis=dict(showline=False,
                            zeroline=False,
                            showgrid=False,  # To visualize the vertical lines
                            ticklen=4,
                            showticklabels=False,
                            title=''),
                yaxis=axis,
                hovermode='closest',
                shapes=line_shapes,
                plot_bgcolor=plot_bgcolor,
                legend={'x': 0, 'y': 1},
                margin={'b': 0, 'l': 0, 'r': 0, 't': 0},
                )
    return layout



def prune_tree(path_to_tree, leaves_to_keep):
    tree = read_treefile(path_to_tree)
    pruned_tree = read_treefile(path_to_tree)
    for leaf in tree.get_terminals():
        if leaf.name not in leaves_to_keep:
            pruned_tree.prune(leaf.name)
    leaves_removed= tree.count_terminals() - pruned_tree.count_terminals()
    return pruned_tree



