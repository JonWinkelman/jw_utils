from pathlib import Path
import subprocess

def run_fasttree(input_fasta, output_tree, **kwargs):
    """
    Runs FastTree on the given alignment file with user-specified flags.

    Parameters:
        input_fasta (str or Path): Path to the input concatenated alignment file.
        output_tree (str or Path): Path to save the output Newick tree file.
        kwargs (dict): Optional FastTree flags. Pass them as keyword arguments, e.g.:
                       run_fasttree("input.fasta", "output.nwk", quiet=True, log="logfile.txt")

    Returns:
        str: The standard output from FastTree execution.
    """
    input_fasta = Path(input_fasta).resolve()
    output_tree = Path(output_tree).resolve()

    # Base command
    cmd = ["fasttree"]

    # Detect “help” request
    if kwargs.pop("help", False) or kwargs.pop("h", False):
        cmd = ["fasttree", "-help"]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # FastTree may print help to stderr instead of stdout:
        help_text = result.stdout or result.stderr
        print(help_text)
        return help_text

    # Add user-defined flags
    for flag, value in kwargs.items():
        if isinstance(value, bool):  # Flags that are just toggles
            if value:  
                cmd.append(f"-{flag}")
        else:  # Flags that require a value (e.g., -log logfile.txt)
            cmd.append(f"-{flag}")
            cmd.append(str(value))

    # Add input and output
    cmd.append(str(input_fasta))

    try:
        # Run FastTree and capture output
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)

        # Write output to the specified tree file
        with open(output_tree, "w") as tree_file:
            tree_file.write(result.stdout)

        print(f"FastTree successfully generated: {output_tree}")
        return result.stdout

    except subprocess.CalledProcessError as e:
        print(f"Error running FastTree: {e.stderr}")
        return None



