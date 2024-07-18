import subprocess
import os

def run_fastqc(for_fp, rev_fp, outdir):
    
    os.makedirs(outdir, exist_ok=True)
    cmd_f = f"""fastqc -o {outdir} {for_fp}"""
    cmd_r = f"""fastqc -o {outdir} {rev_fp}"""
    result = subprocess.run(cmd_f, shell=True,
                        capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error executing command {cmd_f}: {result.stderr}")
    result = subprocess.run(cmd_r, shell=True,
                        capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error executing command {cmd_r}: {result.stderr}")


def run_multiqc(input_dir, title, outdir):
    multiqc_cmd = (
        f"multiqc --title {title} --outdir {outdir} {input_dir} "
    )
    result = subprocess.run(multiqc_cmd, shell=True,
                            capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error executing command: {result.stderr}")