import random
from jw_utils import dna_utils as du
import logomaker

def get_numlines_in_file(fp):
    "return number of lines in a file (1-based)"
    
    lines= 0
    with open(fp, 'r') as f:
        for line in f:
            lines +=1
    return lines
    
    
def _get_4th_nums(max_num):
    """return list of indeces corresponding to fastq ID lines, 0-based
    [0,4,8,12,...]
    """
    nums = []
    for i in range(max_num):
        if i%4 == 0:
            nums.append(i)
    return nums


def _sample_fastq(sample_size, fp):
    """Return a random sample of the fastq id line indeces
        
    """
    id_line_indeces = _get_4th_nums(get_numlines_in_file(fp))
    if sample_size >= len(id_line_indeces):
        sample=id_line_indeces
    else:
        sample = random.sample(population = id_line_indeces, k=sample_size)
    sample.sort()
    return sample
    
    
def _get_lines_from_fastq(id_line_nums, fp, single_or_paired_end='single'):
    """Return a list containing the sampled lines from the full fastq
    
    """
    line_list = []
    with open(fp, 'r') as f:
        list_index = 0
        other_lines = []
        t = 0
        for i, line in enumerate(f):
            if list_index>=len(id_line_nums):
                #get next three lines below ID and add to line_list
                if i in other_lines:
                    line_list.append(line)
                    t+=1
                    if t<3:
                        continue
                    else:
                        break
            if i == id_line_nums[list_index]:
                id_line_num = i
                other_lines = [id_line_num+1, id_line_num+2, id_line_num+3]
                list_index+=1
                line_list.append(line)
            if i in other_lines:
                line_list.append(line)

    return line_list


def sample_fastq(fp,fp2=None, sample_size=10000):
    """Return fastq file(s) containing random sample of sequence reads
    
    *if two filepaths are entered then fastqs are assumed to be from paired end sequencing,
    a second file with mates is written.
    """
    if fp2:
        fps = [fp, fp2]
    else:
        fps = [fp]
    sample = _sample_fastq(sample_size, fp)
    for fp in fps:
        new_fp = fp.replace('.fastq', '_sample.fastq')
        with open(new_fp, 'w') as f:
            line_lst = _get_lines_from_fastq(sample, fp)
            for ele in line_lst:
                f.write(ele)
        print(f'file written to "{new_fp}"')


def read_fastq(file_path):
    """
    A generator function that reads a FASTQ file and yields a sequence object
    for each record consisting of four lines.

    Parameters:
    - file_path (str): The path to the FASTQ file.

    Yields:
    - FastQ_Seq_Obj :  Object containing the header, sequence, and quality score of each read
                       and associated methods
    .
    """
    with open(file_path, 'r') as file:
        while True:
            # Read the next four lines as one block
            header = file.readline().strip()
            if not header:
                break  # Stop iteration if we reach the end of the file
            sequence = file.readline().strip()
            file.readline()  # This reads the '+' separator line and discards it
            quality = file.readline().strip()

            # Yield a dictionary containing the components of the FASTQ record
            yield FastQ_Seq_Obj(header,len(sequence), sequence, quality)



class FastQ_Seq_Obj:
    def __init__(self, id, length, seq, qual):
        self.seq = seq
        self.length = length
        self.id = id
        self.quality = qual

    def rev_comp(self):
        return du.rev_comp(self.seq)


class FastQ_File_Sample_Obj:
    """Create a object from a sample of a fastq file.

    fp : str
        path to fastq file.
    samplesize : int
        number of sequences to randomly sample from the fastq file.
    
    """ 
    
    def __init__(self, fp, samplesize):
        self.filepath = fp
        self.samplesize = samplesize
        self.seq_obj_dict = self.make_dict()
        self.seqs = list(self.seq_obj_dict.values())
        
    
    
    def get_freq_matrix(self, indeces = [0,15], type='information'):
        """
        indeces : [start:end]
            
        """
        seqs = self.get_seqs()
        seqs = [s[indeces[0]:indeces[1]] for s in seqs]
        return logomaker.alignment_to_matrix(seqs,to_type=type )

    def make_logo(self,indeces=[0,15], type='information',):
        """Return a logomaker logo from the indeces entered"""
        
        return logomaker.Logo(self.get_freq_matrix(indeces, type))
    
    
    def get_seqs(self):
        return [obj.seq for obj in list(self.seq_obj_dict.values())]
        

    def make_dict(self, ):
        fastq_line_list = _get_lines_from_fastq(_sample_fastq(self.samplesize, self.filepath), self.filepath, single_or_paired_end='single')
        
        fastq_obj_d = {}
        for i, ele in enumerate(fastq_line_list):
            if i%4 == 0:
                line_lst = ele.split(' ')
                id = line_lst[0]
                continue
            if i%4 == 1:
                seq=ele.strip()
                seq_length = int(len(seq))
                continue
            if i%4 == 2:
                continue
            if i%4 == 3:
                qual=ele
            fastq_obj_d[id]= FastQ_Seq_Obj(id,seq_length,seq,qual)
        return fastq_obj_d