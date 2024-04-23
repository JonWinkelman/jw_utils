def rgb_to_hex():
    return {'rgba(132,137,145,1)':'#848991',
                'rgb(0,208,132)':'#00d084',
                'rgb(0,122,255)':'#007aff',
                'rgb(171,184,195)':'#abb8c3',
                'rgb(255,105,0)':'#ff6900',
                'rgb(252,185,0)':'#fcb900',
                'rgb(123,220,181)':'#7bdcb5',
                'rgb(142,209,252)':'#8ed1fc',
                'rgb(6,147,227)':'#0693e3',
                'rgb(155,81,224)':'#9b51e0',
                'rgb(207,46,46)':'#cf2e2e',
                'rgb(247,141,167)':'#f78da7'}



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


def make_itol_colorstrip_dataset(outfile_path, data_lists, dataset_label='Phyla',
                                     color='#848991', color_branches=1, strip_width=25,
                                     legend_title='Dataset legend', hexcolors=None,):
    """
    make bargraph template for itol dataset adn write to file.

    Parameters
    ----------
    outfile_path   : (str)
    data_lists     : list of f strings, with each word separated by a space. First
                     word corresponds to the name of the internal node, second is the color
                     of colorstrip and branches, third is the displayed name.
                    [['p__Firmicutes_A rgb(100,100,100,0.5) Firmicutes_A'], [...], [...], ...]
    color_branches : (int) 0 or 1. Determines whether branches descending from named nodes will be colored.
    strip_width    : (int)  Width of the annotation strip displaying node name
    

    """
    if not hexcolors:
        hexcolors = ['#58D68D','#F4D03F','#F5B041','#AAB7B8','#566573','#A93226','#EC7063', '#A569BD', '#5DADE2','#48C9B0']

    with open(outfile_path, 'w') as f:
        f.write('DATASET_COLORSTRIP\n')
        f.write('SEPARATOR COMMA\n')
        f.write(f'DATASET_LABEL,{dataset_label}\n')
        f.write(f'COLOR,{color}\n')
        f.write(f'COLOR_BRANCHES,{1}\n')
        # f.write(f'{ita._make_field_colors(name_list, hexcolors=hexcolors)}\n')
        # f.write(f'{ita._make_field_labels(name_list)}\n')
        f.write(f'LEGEND_TITLE,{legend_title}\n')
        f.write(f"LEGEND_SHAPES,{','.join(['1' for _ in name_list])}\n")
        # f.write(f'{_make_legend_colors(name_list, hexcolors=hexcolors)}\n')
        # f.write(f'{_make_legend_labels(name_list)}\n')
        f.write(f'STRIP_WIDTH,{25}\n')
        f.write('DATA\n')

        for annot_line in data_lists:
            f.write(annot_line+'\n')
