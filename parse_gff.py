#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 13:34:22 2022

@author: jonwinkelman
"""
import pandas as pd

class _seq_attributes:
    'make a sequence object from a partially processed line in a GFF file'
    def __init__(self, line, anot_lst, attributes_dict, *args, **kwargs):
        """ make a sequence object for a each feature of a given type
        
        parameters
        line (str): the entire feature line in gff file
        annot_lst (list): beginning of feature line, tab-delimited
        attribute_dict (dict): second part of feature line. Field are variable depending
                            how it was annotated

        """
        self.line = line
        self.ID = attributes_dict['ID']
        self.Parent = attributes_dict.get('Parent')
        self.Dbxref = attributes_dict.get('Dbxref')
        self.Name = attributes_dict.get('Name')
        self.gbkey = attributes_dict.get('gbkey')
        self.inference = attributes_dict.get('inference')
        self.locus_tag = attributes_dict.get('locus_tag')
        self.product = attributes_dict.get('product')
        self.protein_id = attributes_dict.get('protein_id')
        self.transl_table = attributes_dict.get('transl_table')
        self.gene_biotype = attributes_dict.get('gene_biotype')
        self.Note = attributes_dict.get('Note')
        self.gene = attributes_dict.get('gene')
        self.transcript_id = attributes_dict.get('transcript_id')
        self.coverage = attributes_dict.get('cov')
        self.FPKM = attributes_dict.get('FPKM')
        self.TPM = attributes_dict.get('TPM')


        self.chromosome = anot_lst[0]
        self.source = anot_lst[1]
        self.feature_type = anot_lst[2]
        self.start = int(anot_lst[3])
        self.end = int(anot_lst[4])
        self.score = anot_lst[5]
        self.strand = anot_lst[6]   
        self.phase = anot_lst[7]
        
        
def make_full_seq_obj_dict(path_to_gff):
    """Return a dict {feature_ID:seq_obj} for all lines in the gff file"""
    
    obj_dict  = {}
    with open(path_to_gff, 'r') as f:
        for line in f:
            if line[0] != '#' and line[0] != '':
                anot_lst = line.split('\t')
                att_dict = _make_attribut_dict(line)
                seq_obj = _seq_attributes(line,anot_lst, att_dict)
                obj_dict[seq_obj.ID] = seq_obj
    return obj_dict
        
        

def _make_attribut_dict(gff_line):
    """return a dict of attributes parsed from a NCBI-formatted gff3 line"""
    
    gff_line = gff_line.strip(';\n')
    attr_list = gff_line.split('\t')[-1].split(';') #list of ['key=val','key=val','key2=val2',...]
    temp = [key_val.split('=') for key_val in attr_list]
    return {key_val[0]:key_val[1].strip() for key_val in temp}





def make_seq_object_dict(path_to_gff, feature_type = 'gene'):               
    """
    Parse GFF and make dict containing feature objects.
    
    feature_types = 'gene', 'CDS', 'tRNA', 'exon', 'pseudogene', 'SRP_RNA', 'sequence_feature',...
    
    Parameters:
    path_to_gff3 (str): path to gff3 file
    feature_type (str): 
    
    Returns:
    dict: dictionary key is the feature name, values are feature objects that contain
    various attributes such as start and end position in the genome, product description,
    what strand of DNA feature is encoded on, parent gene name, 
    """
    seq_dict = {}
    with open(path_to_gff, 'r') as f:
        for gff_line in f:
            anot_lst = _get_line_list(gff_line, feature_type)
            if anot_lst:
                attributes_dict = _make_attribut_dict(gff_line)
                seq_dict[attributes_dict['ID']] = _seq_attributes(gff_line, anot_lst, attributes_dict)
    return seq_dict


def get_geneObjects_on_strand(path_to_gff, strand):
    """return dict {ID:object} of all gene objects that exist on given strand
    
    strand (str): '+' or '-'
    """

    seq_obj_dict = make_seq_object_dict(path_to_gff)
    return {obj.ID:obj for obj in seq_obj_dict.values() if obj.strand==strand}



def get_contig_names(gff_path):
    "return set of contigs present in the gff file"
    with open(gff_path, 'r') as f:
        contigs =  [line.split('\t')[0] for line in f if line[0] != '#']
        return set(contigs)



def _get_line_list(gff_line, feature_type):
    """parse gff line and return a list of annotations if line is annot for feature type"""

    if gff_line[0] == '#' or gff_line[0] == ' ':
        return None

    anot_lst = gff_line.split('\t')
    if len(anot_lst) <6:
        raise Exception(f'line {gff_line} in gff does not contain standard annotation format')
    if feature_type=='all':
        return anot_lst
    if anot_lst[2] != feature_type:
        return None
    if anot_lst[-1].find('ID=') == -1:
        return None
    return anot_lst
        
        
        
def make_gene2prot_dict(path_to_gff3):               
    """
    Return a dictionary {geneID:proteinID}
    
    
    Parameters:
    path_to_gff3 (str): path to gff3 file
    
    Returns:
    dict
    """
    seq_dict = make_seq_object_dict(path_to_gff3, feature_type = 'CDS')
    return {seq_dict[protein].Parent:protein for protein in seq_dict.keys()}
   


def make_prot2gene_dict(path_to_gff3):
    """
    Return a dictionary {proteinID:geneID}
    
    
    Parameters:
    path_to_gff3 (str): path to gff3 file
    
    Returns:
    dict
    """
    seq_dict = make_seq_object_dict(path_to_gff3, feature_type = 'CDS')
    return {protein:seq_dict[protein].Parent for protein in seq_dict.keys()}



def write_simple_annot_file(path_to_gff, filepath):
    '''
    write tsv file with protein accession, parent gene, common name and function in each row
    
    parameters:
    
        path_to_gff (str): path to gff3 annotation file.
        
        filename/path (str): name you want to give file with ".txt" extension added. If
        path is given, then file will be written to that path, else it will be written in 
        current working directory.
    '''
    jb_bug = make_seq_object_dict(path_to_gff, feature_type='CDS')
    jb_bug_gene = make_seq_object_dict(path_to_gff, feature_type='gene')
    with open(f'{filepath}.txt', 'w') as f:
        for prot in jb_bug.keys():
            gene = jb_bug[prot].Parent.strip()
            if jb_bug_gene.get(gene):
                common = jb_bug_gene[gene].Name
            else: common = None
            product = jb_bug[prot].product
            f.write(f'{prot}\t{gene}\t{common}\t{product}\n')
            
            
                     
def make_simple_annot_df(path_to_gff, start_end=False, contig=False):
    """generate simple dataframe from gff

    index: gene_ID
    columns: 'protein_ID', 'common_name', 'product'

    parameters:
    path_to_gff (str): path to the gff file
    start_end (bool): include start and end of gene as individual columns
    contig (bool): if True, add contig names to df
    
    """
    jb_bug = make_seq_object_dict(path_to_gff, feature_type='CDS')
    jb_bug_gene = make_seq_object_dict(path_to_gff, feature_type='gene')
    contig_name = []
    genes = []
    proteins = []
    common_name = []
    product = []
    starts = []
    end = []
    strand = []
    for prot in jb_bug.keys():
        proteins.append(prot)
        gene = jb_bug[prot].Parent.strip()
        genes.append(gene)
        if start_end:
            starts.append(jb_bug[prot].start)
            end.append(jb_bug[prot].end)
            strand.append(jb_bug[prot].strand)
        if jb_bug_gene.get(gene):
            common_name.append(jb_bug_gene[gene].Name)
        else: 
            common_name.append('')
        if contig:
            contig_name.append(jb_bug[prot].chromosome)
        if jb_bug[prot].product:
            product.append(jb_bug[prot].product.strip())
        else:
            product.append(None)
    df = pd.DataFrame()
    df['gene_ID'] = genes
    df['protein_ID'] = proteins
    df['common_name'] = common_name
    df['product'] = product
    if start_end:
        df['start'] = starts
        df['end'] = end
        df['strand'] = strand
    if contig:
        df['contig'] = contig_name
        
    return df.set_index('gene_ID')


def make_saf_annot_df(local_gff_fp):
    df = make_simple_annot_df(local_gff_fp, start_end=True, contig=True)[['contig_name', 'start', 'end', 'strand']].reset_index()
    #GeneID Chr	Start	End	Strand
    df.columns = ['GeneID', 'Chr', 'Start', 'End', 'Strand']
    #remove 'gene-'
    geneids = [id.replace('gene-','') for id in list(df['GeneID'])]
    df['GeneID'] = geneids
    df = df.set_index('GeneID')
    return df


def get_gff_summary(path_to_gff):
    """NEEDS WORK!!!! return summary of genome assembly derived from gff
    
    parameters
    path_to_gff (str):

    return (dict): nested dict of dicts parsed from the 'region' line of the gff file 
    """ 

    summary = {}
    with open(path_to_gff, 'r') as f:
        for line in f:
            if line[0] == '#':
                if line.startswith('#!genome-build '):
                    summary['genome-build']=line.split(' ')[-1].strip()
                elif line.startswith('#!genome-build-accession'):
                    summary['accession'] = line.split(' ')[-1].strip()
                elif line.startswith('##species '):
                    summary['NCBI_taxon_url'] = line.split(' ')[-1].strip()
                elif line.startswith('#!processor '):
                    summary['processor'] = line.split(' ')[-1].strip()

            else:
                line_lst = line.split('\t')
                if line_lst[2] == 'region':
                    key_val_lst =line_lst[-1].split(';') #[key=val1,key=val2,key=val3,..]
                    key_vals_dict = {key_val.split('=')[0]:key_val.split('=')[1] for key_val in key_val_lst}
                    contig_name=line_lst[0]
                    length = line_lst[4]
                    summary['contigs']={contig_name:{'length':length}}
                    for key, value in key_vals_dict.items():
                        summary['contigs'][contig_name].update({key:value.strip()})
    return summary



def change_gff_attribute(path_to_gff, feature_val_dict, attribute, new_file_path=None, feature_type='gene'):
    """Change value of attribute in GFF file either inplace or return new file.
    
    **make sure key (feature ID) in feature_val_dict matches  is in the correct format. 
    e.g. if the feature_type == 'gene', then the  key should be a gene ID. if the 
    feature type =='CDS', then the key should be a protein ID. 
    
    Parameters:
    path_to_gff (str):
    feature_val_dict (dict): {feature_ID,attribute_value}
    attribute (str): attribute dictionary key for which the value will be changed.
            Available keys vary in gffs. 
            Example of an attribute dict:
            E.g. {'ID': 'gene-PA14_73420','Name': 'rnpA','gbkey': 'Gene','gene': 'rnpA',
                    'gene_biotype': 'protein_coding','locus_tag': 'PA14_73420'}

            

    in_place (bool): if True, edits and saves changes to the input gff. If False, saves a copy
    new_file_path (str): if in_place=False, saves copy of edited GFF here
    feature_type (str): feature type fo which to change attribute value

    """ 
    for key in feature_val_dict:
        feature_val_dict = {feature_type.lower()+'-'+feature if not feature.startswith(feature_type.lower()+'-') else feature:val for feature,val in feature_val_dict.items()}
    with open(path_to_gff, 'r') as f:
        lines = []
        for line in f:
            if line[0]=='#':
                lines.append(line)
            else:
                line_lst = line.split('\t')
                if line_lst[2]==feature_type:
                    attribute_list = line_lst[-1].split(';')
                    att_dict = {ele.split('=')[0]:ele.split('=')[1].strip() for ele in attribute_list}
                    att_dict = {ele.split('=')[0]:ele.split('=')[1].strip() for ele in attribute_list}
                    #change attribute value
                    if att_dict['ID'] in feature_val_dict.keys():
                        att_dict[attribute]=feature_val_dict[att_dict['ID']]
                        new_line_beg = '\t'.join(line_lst[:-1])+'\t'
                        new_line_end = ';'.join([key+'='+val for key, val in att_dict.items()])
                        print(new_line_end)
                        lines.append(new_line_beg + new_line_end + '\n')
                    else:
                        new_line_beg = '\t'.join(line_lst[:-1])+'\t'
                        new_line_end = ';'.join([key+'='+val for key, val in att_dict.items()])
                        lines.append(new_line_beg + new_line_end + '\n')
                else:
                    lines.append(line)
        if not new_file_path:
            new_file_path = path_to_gff.replace('.gff', '_edit.gff') 
        with open(new_file_path, 'w') as f:
            for line in lines:
                f.write(line)




### newer ###############
from dataclasses import dataclass

@dataclass
class Feature:
    ID: str
    parent: str
    feature_type: str
    start: int
    end: int
    strand: str
    attributes: dict

def parse_gff(path, feature_type='gene'):
    feature_dict = {}
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # skip bad lines
            if feature_type != 'all' and fields[2] != feature_type:
                continue

            attr_dict = {kv.split('=')[0]: kv.split('=')[1] for kv in fields[8].split(';') if '=' in kv}
            feature_id = attr_dict.get('ID')
            if not feature_id:
                continue  # skip features without ID

            feature = Feature(
                ID=feature_id,
                parent=attr_dict.get('Parent'),
                feature_type=fields[2],
                start=int(fields[3]),
                end=int(fields[4]),
                strand=fields[6],
                attributes=attr_dict
            )
            feature_dict[feature_id] = feature
    return feature_dict



if __name__ == '__main__':
    pass
    # parser=argparse.ArgumentParser()
    # parser.add_argument('--path_to_gff', type=str, required=True)
    # args=parser.parse_args()
    # gff_obj = make_seq_object_dict(args.path_to_gff)
    # get_gff_summary(args.path_to_gff)