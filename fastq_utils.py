import random


def get_numlines_in_file(fp):
    "return number of lines in a file (1-based)"
    
    lines= 0
    with open(fp, 'r') as f:
        for line in f:
            lines +=1
    return lines
    
    
def _get_4th_nums(max_num):
    """return list of indeces corresponding to fastq ID lines"""
    nums = []
    for i in range(max_num):
        if i%4 == 0:
            nums.append(i)
    return nums


def _sample_fastq(sample_size, fp):
    """Return a random sample of the fastq id line indeces
        
    """
    id_line_indeces = _get_4th_nums(get_numlines_in_file(fp))
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