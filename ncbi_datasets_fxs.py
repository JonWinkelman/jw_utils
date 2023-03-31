import os

def make_summary_dict(ncbi_summary):
    summ_dict = {}
    if not ncbi_summary.get('assemblies'):
        raise Exception('There are no assemblies in the summary Dict')
    for ass in ncbi_summary['assemblies']:
        acc = ass['assembly']['assembly_accession']
        summ_dict[acc] = ass['assembly']
    return summ_dict
    
    
    
    
def sci_strain_names(summ_dict):
    name_d = {}
    for accession, d in summ_dict.items():
        sci_name_lst = d['org']['sci_name'].split(' ')
        strain = d['org'].get('strain')
        if strain:
            sci_name = ' '.join(sci_name_lst[:2])
            full_name = f'{sci_name} {strain}'

            name_d[accession] = full_name
        else:
            full_name = ' '.join(sci_name_lst)
            name_d[accession] = full_name
    return name_d
    
    
    
def check_for_file(data_folder_path, prefix = 'GC', suffix = '.faa'):
    """return folder names that do not contain a file with the given extension
    
    parameters:
    data_folder_path (str): Path to the folder named data, direct parent of the assembly folders
    prefix (str): A string that all of the assembly folders begin with
    suffix (str): A string on the end of the file you are querying for the presence of 
    """
    folders = [f for f in os.listdir(data_folder_path) if f.startswith(prefix)]
    d = {}
    for f in folders:
        fp = os.path.join(data_folder_path, f)
        files = os.listdir(fp)
        foi =  f'file ending in {suffix} Not present'
        for file in files:
            if file.endswith(suffix):
                foi = 1
        if foi == f'file ending in {suffix} Not present':
            d[f] = foi
    return d
    
 
    
def copy_proteomes_to_own_folder(path_to_ncbi_data):
    """Copy proteomes from their ncbi assembly folders into ./Proteomes and rename them
    
    parameters:
    path_to_ncbi_data (str): This is the folder named 'data' inside the 'ncbi folder'
    return (list): accessions that were moved
    """
    import os
    import shutil
    path=path_to_ncbi_data
    os.makedirs('./Proteomes')
    folders =  os.listdir(path)
    folders = [fold for fold in folders if fold.startswith('GC')]
    print(f'{len(folders)} assemblies found in {path}')
    failed = []
    for folder in folders:
        proteome_file = [file for file in os.listdir(os.path.join(path,folder)) if file.endswith('.faa')]
        if len(proteome_file) == 1:
            path_to_proteome = os.path.join(path,folder,proteome_file[0])
            shutil.copy(path_to_proteome, f'./Proteomes/{folder}.faa')
        else:
            print(f'{folder} did not have 1 and only 1 .faa file')
            failed.append(folder)
    proteomes = os.listdir('./Proteomes')
    accessions = [acc.replace('.faa', '') for acc in proteomes]
    print(f'{len(folders)} folders were processed and there are {len(proteomes)} proteomes  in "./Proteomes" ')
    return accessions