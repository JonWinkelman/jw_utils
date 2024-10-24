import os
import subprocess

# Define common suffix patterns for paired-end reads
SUFFIX_PATTERNS = [
    ('_1', '_2'),   # Standard Illumina suffix
    ('.R1', '.R2'),  # R1/R2 suffix with dots
    ('_R1', '_R2'),  # R1/R2 suffix with underscores
    ('-1', '-2'),    # Hyphen-separated suffix
    ('.1', '.2'),    # Dot-separated numeric suffix
    ('_read1', '_read2')  # "_read1" and "_read2" suffixes
]

def find_paired_files(base_name, input_dir, suffix):
    """
    Finds paired-end FASTQ files matching a set of suffix patterns.

    Args:
        base_name (str): The base name of the input FASTQ files.
        input_dir (str): The directory containing the FASTQ files.
        suffix (str): The suffix for the FASTQ files (e.g., '.fastq' or '.fastq.gz').

    Returns:
        tuple: A tuple containing the paths to the R1 and R2 files, or (None, None) if not found.
    """
    for pattern1, pattern2 in SUFFIX_PATTERNS:
        f1_fp = os.path.join(input_dir, base_name + f'{pattern1}{suffix}')
        f2_fp = os.path.join(input_dir, base_name + f'{pattern2}{suffix}')
        if os.path.exists(f1_fp) and os.path.exists(f2_fp):
            return f1_fp, f2_fp
    return None, None

def bowtie2_align(base_name, input_dir, genome_index, suffix='.fastq', 
                  convert_to_bam=False, no_unal=False, threads=None):
    """
    Aligns sequencing reads using Bowtie2. Optionally converts the SAM output to BAM format
    and removes the SAM file after conversion.

    Args:
        base_name (str): The base name of the input FASTQ files.
        input_dir (str): The directory containing the FASTQ files.
        genome_index (str): Path to the Bowtie2 genome index.
        suffix (str): The suffix of the FASTQ files (default: '.fastq' or '.fastq.gz').
        convert_to_bam (bool): If True, converts the SAM output to BAM format and deletes
                               the original SAM file (default: False).
        no_unal (bool): If True, excludes unaligned reads from the SAM/BAM output using the
                        `--no-unal` option in Bowtie2 (default: False).
        threads (str or int): number of threads for bowtie to use, default=1

    Raises:
        FileNotFoundError: If the required FASTQ files are not found.
        subprocess.CalledProcessError: If the Bowtie2 or samtools commands fail.
    """
    # Ensure the suffix starts with a dot
    if not suffix.startswith('.'):
        suffix = '.' + suffix

    # Attempt to find paired-end files
    f1_fp, f2_fp = find_paired_files(base_name, input_dir, suffix)
    if f1_fp and f2_fp:
        paired = True
        dest_sam = os.path.join(input_dir, f"{base_name}_paired.sam")
        dest_bam = os.path.join(input_dir, f"{base_name}_paired.bam")
        log_out = os.path.join(input_dir, f"{base_name}_paired.log")
    else:
        paired = False
        f1_fp = os.path.join(input_dir, base_name + suffix)
        if not os.path.exists(f1_fp):
            raise FileNotFoundError(f'File not found: {f1_fp}')
        dest_sam = os.path.join(input_dir, f"{base_name}.sam")
        dest_bam = os.path.join(input_dir, f"{base_name}.bam")
        log_out = os.path.join(input_dir, f"{base_name}.log")

    # Construct Bowtie2 command based on alignment type
    if paired:
        print(f'Processing paired-end data:\nR1: {f1_fp}\nR2: {f2_fp}\nSAM: {dest_sam}')
        bowtie2_cmd = ["bowtie2", "-x", genome_index, "-1", f1_fp, "-2", f2_fp, "-S", dest_sam]
    else:
        print(f'Processing single-end data:\nR1: {f1_fp}\nSAM: {dest_sam}')
        bowtie2_cmd = ["bowtie2", "-x", genome_index, "-U", f1_fp, "-S", dest_sam]

    # Add --no-unal option if specified
    if no_unal:
        bowtie2_cmd.append("--no-unal")
    if threads:
        bowtie2_cmd.append("-p")
        bowtie2_cmd.append(str(threads))
    # Run Bowtie2 and log output
    with open(log_out, 'w') as log_file:
        try:
            subprocess.run(bowtie2_cmd, stderr=log_file, check=True)
        except subprocess.CalledProcessError as e:
            print(f'Error running Bowtie2: {e}')
            return

    # Optional: Convert SAM to BAM and remove SAM file
    if convert_to_bam:
        try:
            subprocess.run(["samtools", "view", "-o", dest_bam, "-b", dest_sam], check=True)
            print(f'BAM file created: {dest_bam}')
            os.remove(dest_sam)
            print(f'SAM file removed: {dest_sam}')
        except subprocess.CalledProcessError as e:
            print(f'Error converting SAM to BAM: {e}')
        except OSError as e:
            print(f'Error removing SAM file: {e}')

# Usage example
if __name__ == "__main__":
    bowtie2_align(
        base_name="sample",
        input_dir="/path/to/fastq/files",
        genome_index="/path/to/genome/index",
        convert_to_bam=True,
        no_unal=True
    )
