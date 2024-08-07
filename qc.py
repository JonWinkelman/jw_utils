import subprocess
import os

def run_fastqc(fp, outdir):
    
    os.makedirs(outdir, exist_ok=True)
    cmd = f"""fastqc -o {outdir} {fp}"""

    result = subprocess.run(cmd, shell=True,
                        capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error executing command {cmd_f}: {result.stderr}")


def run_multiqc(input_dir, title, outdir):
    multiqc_cmd = (
        f"multiqc --title {title} --outdir {outdir} {input_dir} "
    )
    result = subprocess.run(multiqc_cmd, shell=True,
                            capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error executing command: {result.stderr}")



