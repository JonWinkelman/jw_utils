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
        self.ID = attributes_dict['gene_id']
        self.Parent = attributes_dict.get('Parent')
        self.Dbxref = attributes_dict.get('Dbxref')
        self.Name = attributes_dict.get('gene_name')
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



        self.chromosome = anot_lst[0]
        self.source = anot_lst[1]
        self.feature_type = anot_lst[2]
        self.start = int(anot_lst[3])
        self.end = int(anot_lst[4])
        self.score = anot_lst[5]
        self.strand = anot_lst[6]   
        self.phase = anot_lst[7]


def _make_attribut_dict(gtf_line):
    """return a dict of attributes parsed from a NCBI-formatted gff3 line"""
    attr_list = gtf_line.strip(';\n').split('\t')[-1].split('; ')
    temp = [key_val.split(' ') for key_val in attr_list]
    return  {key_val[0].strip():key_val[1].strip('"') for key_val in temp}


def _get_line_list(gtf_line, feature_type):
    """parse gff line and return a list of annotations if line is annot for feature type"""

    if gtf_line[0] == '#' or gtf_line[0] == ' ':
        return None

    anot_lst = gtf_line.split('\t')
    if len(anot_lst) <6:
        raise Exception(f'line {gtf_line} in gff does not contain standard annotation format')
    if feature_type=='all':
        return anot_lst
    if anot_lst[2] != feature_type:
        return None
    return anot_lst

def make_seq_object_dict(path_to_gtf, feature_type = 'gene'):               
    """
    Parse GTF and make dict containing feature objects.
    
    
    Parameters:
    path_to_gtf (str): path to gtf file
    feature_type (str): 
    
    Returns:
    dict: dictionary key is the feature name, values are feature objects that contain
    various attributes such as start and end position in the genome, product description,
    what strand of DNA feature is encoded on, parent gene name, 
    """
    seq_dict = {}
    with open(path_to_gtf, 'r') as f:
        for line in f:
            anot_lst = _get_line_list(line, feature_type)
            if anot_lst:
                attributes_dict = _make_attribut_dict(line)
                seq_dict[attributes_dict['gene_id']] = _seq_attributes(line, anot_lst, attributes_dict)
    return seq_dict


def make_SAF_annotation(annotation_fp, feature_type):
    """"""

    GeneID=[]
    Chr=[]	
    Start=[]
    End=[]
    Strand=[]
    
    seq_obj_d = make_seq_object_dict(annotation_fp, feature_type=feature_type)
    for gene_id, seq_obj in seq_obj_d.items():
        GeneID.append(gene_id)
        Chr.append(seq_obj.chromosome)
        Start.append(seq_obj.start)
        End.append(seq_obj.end)
        Strand.append(seq_obj.strand)
        
    df = pd.DataFrame()
    df['GeneID']=GeneID
    df['Chr']=Chr
    df['Start']=Start
    df['End']=End
    df['Strand']=Strand
    return df.set_index('GeneID')