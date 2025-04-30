import subprocess
from pathlib import Path

def run_prodigal(input_file, output_dir=".", **kwargs):
    """
    Run Prodigal with flexible arguments, using a required input FASTA file
    and optional output directory.

    Parameters:
    - input_file (str or Path): Input FASTA file (required)
    - output_dir (str or Path): Directory to store output files (default: current directory)
    - **kwargs: Prodigal flags as keyword arguments (e.g., p='meta', f='gff', q=True)

    If -o, -a, or -d are not provided, they are auto-generated from the input filename.
             -a:  Write protein translations to the selected file.
             
     ### available arguments in prodigal:
     -c:  Closed ends.  Do not allow genes to run off edges.
     -d:  Write nucleotide sequences of genes to the selected file.
     -f:  Select output format (gbk, gff, or sco).  Default is gbk.
     -g:  Specify a translation table to use (default 11).
     -h:  Print help menu and exit.
     -i:  Specify FASTA/Genbank input file (default reads from stdin).
     -m:  Treat runs of N as masked sequence; don't build genes across them.
     -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
     -o:  Specify output file (default writes to stdout).
     -p:  Select procedure (single or meta).  Default is single.
     -q:  Run quietly (suppress normal stderr output).
     -s:  Write all potential genes (with scores) to the selected file.
     -t:  Write a training file (if none exists); otherwise, read and use
          the specified training file.
     -v:  Print version number and exit.
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_name = input_path.stem

    # Add required input flag
    kwargs["i"] = str(input_path)

    # Set defaults for output filenames if not explicitly passed
    default_outputs = {
        "o": output_dir / f"{base_name}.gff",
        "a": output_dir / f"{base_name}.faa",
        "d": output_dir / f"{base_name}.fna"
    }
    for flag, default_path in default_outputs.items():
        if flag not in kwargs:
            kwargs[flag] = default_path

    # Define valid Prodigal flags and whether they take a value
    valid_flags = {
        'a': True, 'c': False, 'd': True, 'f': True, 'g': True, 'h': False,
        'i': True, 'm': False, 'n': False, 'o': True, 'p': True, 'q': False,
        's': True, 't': True, 'v': False
    }

    # Build the command
    cmd = ["prodigal"]
    for k, v in kwargs.items():
        if k not in valid_flags:
            raise ValueError(f"Invalid Prodigal option: -{k}")
        if valid_flags[k]:
            cmd += [f"-{k}", str(v)]
        elif v:  # Boolean flags
            cmd += [f"-{k}"]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)