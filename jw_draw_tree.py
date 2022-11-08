from Bio import Phylo

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


def draw_clade(clade, x_start, line_shapes, line_color='rgb(15,15,15)', line_width=1, x_coords=0, y_coords=0):
    """Recursively draw the tree branches, down from the given clade"""

    x_curr = x_coords[clade]
    y_curr = y_coords[clade]

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
            draw_clade(child, x_curr, line_shapes, x_coords=x_coords, y_coords=y_coords)


def read_treefile(filename):
    'create tree object from newick format using Bio.Phylo'
    tree = Phylo.read(filename, "newick")
    return tree



def create_plotly_tree(tree, title=None, t_nodes_color_dict=None, in_node_size=2, t_node_size=10,
                            in_node_color='rgb(100,100,100)'):
    """make plotly figure from newick tree

    parameters:
    tree_filepath (str or Phylo.Newick.Tree object): path to newick tree or tree object
    title (str): title of plotly graph
    t_nodes_color_dict (dict): differentiate groups of leaves on tree: 'rgb()':['leafname1', 'leafname2','leafname3'...]
    node_size_dict (dict): set the size of the different colors entered as the t_nodes_color_dict dict keys
                            this is in the form 'rgb()':['leafname1', 'leafname2','leafname3'...]

    return (Plotly.graph_objects.Figure)
    """
    if type(tree) != Phylo.Newick.Tree:
        tree = Phylo.read(tree, "newick")
    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    draw_clade(tree.root, 0, line_shapes, line_color='rgb(25,25,25)', line_width=1, x_coords=x_coords,
               y_coords=y_coords)
    
    X = []
    Y = []
    text = []
    node_sizes = []
    color_dict = {}
    t_node_names = [clade.name for clade in tree.get_terminals()]
    loops=0
    if type(t_nodes_color_dict)!= dict:
        print('this was triggered')
        t_nodes_color_dict={}
        t_nodes_color_dict['rgb(100,100,100)'] = t_node_names
    
    if tree.get_nonterminals()[0].name: 
        my_tree_clades = tree.depths().keys()
    else:
        my_tree_clades = tree.get_terminals()
    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
        #generate hover text and node colors
        i=0
        loops+=1
        for color in t_nodes_color_dict.keys():#check all groups to see if clade is included
            if cl.name in t_nodes_color_dict[color]:
                text.append(cl.name)
                color_dict[cl.name]=color
                i+=1
        if i==0:
            if cl.name:
                #print(cl.name, loops)
                text.append(cl.name)
                color_dict[loops]= in_node_color
                i+=1
    #node sizes
    for cl in my_tree_clades:
        if cl.name in t_node_names:
            node_sizes.append(t_node_size)
        else:
            node_sizes.append(in_node_size)


    axis = dict(showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''  # y title
                )

    data = dict(type='scatter',
                x=X,
                y=Y,
                mode='markers',
                marker=dict(color=list(color_dict.values()),
                            size=node_sizes
                ),
                text=text,  # vignet information of each node
                hoverinfo='text',
                )
    if title:
        title=title
    layout = dict(title=title,
                  paper_bgcolor='rgb(248,248,248)',
                  dragmode="lasso",
                  font=dict(family='Balto', size=14),
                  #width=750,
                  height=550,
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
                  plot_bgcolor='rgb(248,248,248)',
                  legend={'x': 0, 'y': 1},
                  margin={'b': 0, 'l': 0, 'r': 0, 't': 0}
                  )
    fig = dict(data=[data], layout=layout)
    return fig


def create_special_tree(tree_filepath, ingroup, outgroup, HOG, copies, title=None):
    missing_dict = a.genomes_missing_HOGs(ingroup, outgroup, HOG, copies)
    tree = read_treefile(tree_filepath)
    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    draw_clade(tree.root, 0, line_shapes, line_color='rgb(25,25,25)', line_width=1, x_coords=x_coords,
               y_coords=y_coords)
    my_tree_clades = x_coords.keys()
    X = []
    Y = []
    text = []
    node_sizes = []
    color_dict = {}
    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
        if cl.name in name_dict.keys():
            text.append(name_dict[cl.name])
        else:
            text.append(cl.name)
            
        if cl.name in ingroup:
            color_dict[cl.name] = colors['t_green']  #seagreen trestle
            node_sizes.append(10)
        elif cl.name in outgroup:
            color_dict[cl.name] = colors['t_blue']
            node_sizes.append(10)
        else:
            color_dict[cl.name] = '#92A0A9' #grey'
            node_sizes.append(5)
        if cl.name in missing_dict['ingroup']:
           color_dict[cl.name] = 'rgba(150,255,150,0.65)'
        if cl.name in missing_dict['outgroup']:
           color_dict[cl.name] = 'rgba(0,130,255,0.4)'

    axis = dict(showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''  # y title
                )

    data = dict(type='scatter',
                x=X,
                y=Y,
                mode='markers',
                marker=dict(color=list(color_dict.values()),
                            size=node_sizes),
                text=text,  # vignet information of each node
                hoverinfo='text',
                )
    if title:
        title=title
    layout = dict(title=title,
                  paper_bgcolor='rgb(248,248,248)',
                  dragmode="lasso",
                  font=dict(family='Balto', size=14),
                  #width=750,
                  height=550,
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
                  plot_bgcolor='rgb(248,248,248)',
                  legend={'x': 0, 'y': 1},
                  margin={'b': 0, 'l': 0, 'r': 0, 't': 0}
                  )
    fig = dict(data=[data], layout=layout)
    return fig 