"""
fastqc
"""

def fastqc(a, b):
    cmd = "Fastqc "
    cmd+= input('fichier(s).fastq')

    return cmd
p = subprocess.popen(cmd )#PIpe