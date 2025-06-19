from pathlib import Path
from typing import List, Dict, Optional, Union, Tuple


def write_internal_anc_node_names(leaf_names, tree):
    """returns tree, and node_dict, with all ancestral nodes of given branches named <leafname>__anc-<depth>
    
    *names all ancestral nodes in tree and adds to dict
    *replaces '_' with '-' in any leaf names because itol doesn't handle underscores well
    """
    node_d={name.replace('_', '-'):[] for name in leaf_names}
    new_tree=tree.copy()
    for leaf_name in leaf_names:
        for cl in new_tree:
            if cl.name==leaf_name:
                leaf_name=leaf_name.replace('_', '-')
                for i, anc in enumerate(cl.get_ancestors(), start=1):
                    anc.name = leaf_name+ "-anc" + str(-i) 
                    node_d[leaf_name].append(anc.name)
    return new_tree, node_d
    

def find_LCAs_of_monophyletic_clades(tree, target_leaves):
    """Return a list of LCA node names for each monophyletic clade withing the target leaves

    * if the group is monophyletic, it will return a single node
    * I have used this for highlighting clades in itol
    
    tree (Bio.Phylo tree):
    target_leaves (iterable): Leaves of interest where you want to find LCA nodes in tree, e.g. to annotate clades in itol...

    
    ------------------------
    return (list): List of node.names from Bio.Phylo tree
    """
    for cl in tree.get_nonterminals():
        if not cl.name:
            raise Exception(print('Internal nodes of tree must be named'))
    monophyletic_anc_nodes = []
    monophyletic_nodes = []
    for node in tree.find_clades(order='preorder'):
        # if node.is_terminal():
        #     continue
        tips = {leaf.name for leaf in node.get_terminals()}
        if tips.issubset(target_leaves):  #then this is monophyletic.
            t=node
            # we now need to see if this node is a child of a node already known to be monophyletic
            # since we are trraversing tree starting at root of tree (preorder), then we come accross most ancient nodes first
            if node.name not in monophyletic_nodes:
                monophyletic_anc_nodes.append(node.name)
            monophyletic_nodes = monophyletic_nodes + [child.name for child in node.clades]
    return monophyletic_anc_nodes

    

def itol_highlight_lineages(node_dict, out_path, colors=None, line_type='normal', line_width=2, 
                           label_type='bold', TYPE='branch', WHAT='node'):
    """
    *** use write_internal_anc_node_names() to generate tree with anc nodes labeled...
    node_dict (dict): dict with {leaf_name:[internal_node1, internal_node2,...]}
    
    TYPE: can be either 'branch' or 'label'. 'branch' will apply customizations to the tree branches, while 'labels' apply to the leaf text labels
    WHAT: can be either 'node' or 'clade', only relevant for internal tree nodes. 'Node' will apply the customization only to a single node, while 'clade' will apply to all child nodes as well.
    COLOR: can be in hexadecimal, RGB or RGBA notation. If RGB or RGBA are used, dataset SEPARATOR cannot be comma.
    WIDTH_OR_SIZE_FACTOR: for type 'branch', specifies the relative branch width, compared to the global branch width setting.
                           for type 'label', specifies the relative font size, compared to the global font size
    STYLE: for type 'branch', can be either 'normal' or 'dashed'
            for type 'label', can be one of 'normal', 'bold', 'italic' or 'bold-italic'
    BACKGROUND_COLOR (optional): only relevant for type 'label', specifies the color of the label background. The value is optional.
    label_type (str): 'normal', 'bold', 'italic' or 'bold-italic'
    """
    s = ''
    for l in node_dict:
        s  = s+l+','
    s=s.strip(',')
    if not colors:
        colors=['#58D68D','#F4D03F','#F5B041','#AAB7B8','#566573','#A93226','#EC7063', '#A569BD', '#5DADE2','#48C9B0'] 
    with open(out_path, 'w') as f:
        f.write('DATASET_STYLE\n')
        f.write('SEPARATOR COMMA\n')
        f.write('DATASET_LABEL,example style\n')
        f.write(f'LEGEND_LABELS,{s}\n')
        f.write('DATA\n')
        for i, (leaf, lineage) in enumerate(node_dict.items()):
            #label leaf name
            s = f"{leaf},label,{WHAT},{colors[i]},1,{label_type}"
            f.write(f'{s}\n')
            for int_node in lineage:
                int_node = int_node.replace(':', '_')
                int_node = int_node.replace('_', ' ')
                s=f"{int_node},{TYPE},{WHAT},{colors[i]},{line_width},{line_type}"
                f.write(f'{s}\n')
                

# def rgb_to_hex():
#     return {'rgba(132,137,145,1)':'#848991',
#                 'rgb(0,208,132)':'#00d084',
#                 'rgb(0,122,255)':'#007aff',
#                 'rgb(171,184,195)':'#abb8c3',
#                 'rgb(255,105,0)':'#ff6900',
#                 'rgb(252,185,0)':'#fcb900',
#                 'rgb(123,220,181)':'#7bdcb5',
#                 'rgb(142,209,252)':'#8ed1fc',
#                 'rgb(6,147,227)':'#0693e3',
#                 'rgb(155,81,224)':'#9b51e0',
#                 'rgb(207,46,46)':'#cf2e2e',
#                 'rgb(247,141,167)':'#f78da7'}


def rgba_to_hex(rgba, include_alpha=True):
    """
    Convert an (R, G, B, A) tuple to a HEX color string.

    Parameters:
    rgba (tuple): A tuple of (R, G, B, A), each ranging from 0-255.
    include_alpha (bool): Whether to include the alpha channel in the HEX code.

    Returns:
    str: The HEX color string.
    """
    r, g, b, a = rgba
    hex_color = "#{:02X}{:02X}{:02X}".format(r, g, b)
    if include_alpha:
        hex_color += "{:02X}".format(a)  # Append alpha as hex
    return hex_color



def _make_field_labels(name_list):
    """return itol line with field labels pulled from from name list"""
    
    field_labels='FIELD_LABELS,'
    for name in name_list:
        field_labels =f'{field_labels}{name},'
    field_labels=field_labels.strip(',')
    return field_labels


def _make_field_colors(name_list,hexcolors=None):
    """Return itol line with hexcolors for each field"""
    if not hexcolors:
        hexcolors = ['#58D68D','#F4D03F','#F5B041','#AAB7B8','#566573','#A93226','#EC7063', '#A569BD', '#5DADE2','#48C9B0']
    field_colors = 'FIELD_COLORS'
    for i, _ in enumerate(name_list):
        field_colors = f'{field_colors},{hexcolors[i]}'
    field_colors=field_colors.strip(',')
    return field_colors


def _make_legend_colors(name_list,hexcolors=None):
    """Return itol line with hexcolors for each field"""
    if not hexcolors:
        hexcolors = ['#58D68D','#F4D03F','#F5B041','#AAB7B8','#566573','#A93226','#EC7063', '#A569BD', '#5DADE2','#48C9B0']
    legend_colors = 'LEGEND_COLORS'
    for i, _ in enumerate(name_list):
        legend_colors = f'{legend_colors},{hexcolors[i]}'
    legend_colors=legend_colors.strip(',')
    return legend_colors


def _make_legend_labels(name_list):
    """return itol line with legend labels pulled from name list"""
    
    legend_labels='LEGEND_LABELS,'
    for name in name_list:
        legend_labels =f'{legend_labels}{name},'
    legend_labels=legend_labels.strip(',')
    return legend_labels


def make_itol_pie_dataset(outfile_path,
                            count_dict,
                            name_list,
                            dataset_label='dataset_label',
                            color='#848991',
                            legend_title='Dataset legend',
                            hexcolors=None,
                         ):
    
    """    

    Parameters
    ----------
    outfile_path
    count_dict :dict
        location can be 0 or -1, representing internal or external pie slices
        {term_node_name : [location, radius, val1, val2,...]}
    name_list (list): name for each value in the pie, must be same number of elements as the 
                      number of *values* in the pie, e.g. 2 names for val1 and val2 
                      in the count_dict {leaf_name1:[location, radius, *val1*, *val2*], ...}


    DATA:
    #the following fields are required for each node:
    #ID,position,radius,value1,value2,value3...
    #position defines the position of the pie chart on the tree:
    #  -1 = external pie chart
    #  a number between 0 and 1 = internal pie chart positioned at the specified value along the node branch (for example, position 0 is exactly at the start of node
    branch, position 0.5 is in the middle, and position 1 is at the end)
    DATA
    #Examples
    
    # node 9606 will have an external pie chart
    # 9606,-1,10,10000,15000,9000
    
    # node 9123 will have an internal pie chart directly over the node, and with radius 50 
    # (actual display radius will depend on other values in the dataset and the MAXIMUM_SIZE specified)
    #9132,0,50,11000,9000,120007
    """
    
    if not hexcolors:
        hexcolors = ['#F5B041','#AAB7B8','#566573','#A93226','#EC7063', '#A569BD', '#5DADE2','#48C9B0','#58D68D','#F4D03F']
    file_lst = []
    
    file_lst.append('DATASET_PIECHART\n')
    file_lst.append('SEPARATOR COMMA\n')
    file_lst.append(f'DATASET_LABEL,{dataset_label}\n')
    file_lst.append(f'COLOR,{color}\n')
    file_lst.append(f'{_make_field_colors(name_list, hexcolors=hexcolors)}\n')
    file_lst.append(f'{_make_field_labels(name_list)}\n')
    file_lst.append(f'LEGEND_TITLE,{legend_title}\n')
    file_lst.append(f"LEGEND_SHAPES,{','.join(['1' for _ in name_list])}\n")
    file_lst.append(f'{_make_legend_colors(name_list, hexcolors=hexcolors)}\n')
    file_lst.append(f'{_make_legend_labels(name_list)}\n')
    file_lst.append('ALIGN_FIELDS,1\n')
    file_lst.append('DATA\n')
    
    
    for name, counts in count_dict.items():
        line = f'{name},'
        for count in counts:
            line = line + str(count)+','
        line=line.strip(',') 
        file_lst.append((line+'\n'))
    file_lst
    with open(outfile_path, 'w') as f:
        for ele in file_lst:
            f.write(ele)


def make_simple_itol_bargraph_dataset(outfile_path, count_dict, name_list, dataset_label='dataset_label',
                                     color='#848991', legend_title='Dataset legend', hexcolors=None,):
    """
    make bargraph template for itol dataset adn write to file.

    Parameters
    ----------
    outfile_path
    name_list (list): same length (n) as [count1, count2,...countn]
    data_dict :dict
        {term_node_name : [count]}

    """
    if not hexcolors:
        hexcolors = ['#58D68D','#F4D03F','#F5B041','#AAB7B8','#566573','#A93226','#EC7063', '#A569BD', '#5DADE2','#48C9B0']

    with open(outfile_path, 'w') as f:
        f.write('DATASET_SIMPLEBAR\n')
        f.write('SEPARATOR COMMA\n')
        f.write(f'DATASET_LABEL,{dataset_label}\n')
        f.write(f'COLOR,{color}\n')
        f.write(f'{_make_field_colors(name_list, hexcolors=hexcolors)}\n')
        f.write(f'{_make_field_labels(name_list)}\n')
        f.write(f'LEGEND_TITLE,{legend_title}\n')
        f.write(f"LEGEND_SHAPES,{','.join(['1' for _ in name_list])}\n")
        f.write(f'{_make_legend_colors(name_list, hexcolors=hexcolors)}\n')
        f.write(f'{_make_legend_labels(name_list)}\n')
        f.write('ALIGN_FIELDS,1\n')
        f.write('DATA\n')

        for name, counts in count_dict.items():
            line = f'{name},'
            for count in counts:
                line = line + str(count)+','
            line=line.strip(',') 
            f.write(line+'\n')




def make_itol_multi_bargraph_dataset(outfile_path, count_dict, name_list, dataset_label='dataset_label',
                                     color='#848991', legend_title='Dataset legend', hexcolors=None,):
    """
    make bargraph template for itol dataset adn write to file.

    Parameters
    ----------
    outfile_path
    name_list (list): same length (n) as [count1, count2,...countn]
    data_dict :dict
        {term_node_name : [count1, count2,...countn]}

    """
    if not hexcolors:
        hexcolors = ['#58D68D','#F4D03F','#F5B041','#AAB7B8','#566573','#A93226','#EC7063', '#A569BD', '#5DADE2','#48C9B0']

    with open(outfile_path, 'w') as f:
        f.write('DATASET_MULTIBAR\n')
        f.write('SEPARATOR COMMA\n')
        f.write(f'DATASET_LABEL,{dataset_label}\n')
        f.write(f'COLOR,{color}\n')
        f.write(f'{_make_field_colors(name_list, hexcolors=hexcolors)}\n')
        f.write(f'{_make_field_labels(name_list)}\n')
        f.write(f'LEGEND_TITLE,{legend_title}\n')
        f.write(f"LEGEND_SHAPES,{','.join(['1' for _ in name_list])}\n")
        f.write(f'{_make_legend_colors(name_list, hexcolors=hexcolors)}\n')
        f.write(f'{_make_legend_labels(name_list)}\n')
        f.write('ALIGN_FIELDS,1\n')
        f.write('DATA\n')

        for name, counts in count_dict.items():
            line = f'{name},'
            for count in counts:
                line = line + str(count)+','
            line=line.strip(',') 
            f.write(line+'\n')

def relabel_itol_treeleafs(tree, relabel_dict, outfile_path):
    """Write itol annotation file to relabel terminal leaves in tree
    
    parameters:
    tree (Bio.Phylo.Newick.Tree): tree that is to be relabeled
    relabel_dict (dict): dict {old_leaf_name:new_leaf_name}
    outfile_path (str): path for new itol annotation file
    """
    import warnings
    tree_leafnames = [cl.name for cl in tree.get_terminals()]
    if len(tree_leafnames) != len(relabel_dict.keys()):
        warnings.warn('The number of tree leafs and the number of dict key names to be replaced are not equal')
    with open(outfile_path, 'w') as f:
        f.write('LABELS\n')
        f.write('SEPARATOR COMMA\n')
        f.write('DATA\n')
        for old_name, new_name in relabel_dict.items():
            f.write(f'{old_name},{new_name}\n')



def make_itol_binary_trait_dataset(outfile_path, count_dict, name_list, dataset_label='dataset_label',
                                     color='#848991',field_shapes=1, legend_title='Dataset legend', hexcolors=None,):
    """
    Add data to simple bargraph template for itol dataset adn write to file.

    Parameters
    ----------
    outfile_path : str

    name_list : str
        list of names for legend and hoverdata. Each name corresponds to a value in count_dict.values(), 
        should be same length (n) as count_dict[node_name] = [val1, val2,..valn]

    file_path : str
        path for new file that is to be created
    count_dict :dict
        (terminal_node_name: [val1, val2,..valn]} Possible values for each node are:   
                                            1 (filled shapes), 0 (empty shapes) and -1 (completely omitted)
    field_shapes: int       
            #1: square
            #2: circle
            #3: star
            #4: right pointing triangle
            #5: left pointing triangle
            #6: checkmark

    """
    if not hexcolors:
        hexcolors = ['#58D68D','#F4D03F','#F5B041','#AAB7B8','#566573','#A93226','#EC7063', '#A569BD', '#5DADE2','#48C9B0']

    with open(outfile_path, 'w') as f:
        f.write('DATASET_BINARY\n')
        f.write('SEPARATOR COMMA\n')
        f.write(f'DATASET_LABEL,{dataset_label}\n')
        f.write(f'COLOR,{color}\n')
        f.write(f"FIELD_SHAPES,{','.join([str(field_shapes) for _ in name_list])}\n")
        f.write(f'{_make_field_labels(name_list)}\n')
        f.write(f'{_make_field_colors(name_list, hexcolors=hexcolors)}\n')

        f.write(f'LEGEND_TITLE,{legend_title}\n')
        f.write(f"LEGEND_SHAPES,{','.join([str(field_shapes) for _ in name_list])}\n")
        f.write(f'{_make_legend_colors(name_list, hexcolors=hexcolors)}\n')
        f.write(f'{_make_legend_labels(name_list)}\n')
        f.write('DATA\n')

        for name, counts in count_dict.items():
            line = f'{name},'
            for count in counts:
                line = line + str(count)+','
            line=line.strip(',') 
            f.write(line+'\n')



def make_itol_colorstrip_dataset(
    outfile_path: Union[str, Path],
    data: List[Dict[str, str]],
    dataset_label: str = 'Phyla',
    dataset_label_color: str = '#848991',
    color_branches: int = 0,
    strip_width: int = 25,
    legend_title: str = 'Legend',
    margin: int = 15, 
    border_width: int = 0, 
    border_color: str = "#0000ff", 
    show_strip_labels: int = 0,
    strip_label_position: str = 'center', 
  
):
    """
    Write an iTOL color strip annotation file from structured data.

    Parameters
    ----------
    outfile_path : Path or str
        Path to output annotation file.
    data : list of dicts
        Each dict should contain:
            - 'node': node name (str)
            - 'color': hex color (e.g., '#FF0000')
            - 'label': label for strip and legend
            - 'shape': legend shape code (e.g., '1' for square)
            e.g.
            data = [
                    {'node': 'NodeA', 'label': 'Proteobacteria', 'color': '#FF0000', 'shape': '1'},
                    {'node': 'NodeB', 'label': 'Firmicutes', 'color': '#00FF00', 'shape': '2'},
                    {'node': 'NodeC', 'label': 'Actinobacteria', 'color': '#0000FF', 'shape': '3'},
                    {'node': 'NodeD', 'label': 'Firmicutes', 'color': '#00FF00', 'shape': '2'},
                    ]
    """
    outfile_path = Path(outfile_path)

    def dedup_in_order(items: List[Tuple]) -> List[Tuple]:
        seen = set()
        out = []
        for item in items:
            if item not in seen:
                out.append(item)
                seen.add(item)
        return out

    # Deduplicate legend entries by (label, color, shape)
    legend_entries = dedup_in_order([
        (entry['label'], entry['color'], entry.get('shape', '1')) for entry in data
    ])

    with open(outfile_path, 'w') as f:
        f.write('DATASET_COLORSTRIP\n')
        f.write('SEPARATOR SPACE\n')
        f.write(f'DATASET_LABEL {dataset_label}\n')
        f.write(f'COLOR {dataset_label_color}\n')
        f.write(f'COLOR_BRANCHES {color_branches}\n')
        f.write(f'LEGEND_TITLE {legend_title}\n')
        f.write(f'LEGEND_SHAPES {" ".join(entry[2] for entry in legend_entries)}\n')
        f.write(f'LEGEND_LABELS {" ".join(entry[0] for entry in legend_entries)}\n')
        f.write(f'LEGEND_COLORS {" ".join(entry[1] for entry in legend_entries)}\n')
        f.write(f'STRIP_WIDTH {strip_width}\n')
        f.write(f'MARGIN {margin}\n')
        f.write(f'BORDER_WIDTH {border_width}\n')
        f.write(f'BORDER_COLOR {border_color}\n')
        f.write(f'SHOW_STRIP_LABELS {show_strip_labels}\n')
        f.write(f'STRIP_LABEL_POSITION {strip_label_position}\n')
        f.write('DATA\n')
        for entry in data:
            f.write(f"{entry['node']} {entry['color']} {entry['label']}\n")
