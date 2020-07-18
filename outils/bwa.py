"""
fonction qui imprime la oommande cmd dans fichier name.sh qui est capable de rouler sur slurm
"""
import re
import sys
import glob
import os
import argparse

def init():
    parser = argparse.ArgumentParser(
        description = 'this scrpit runs bwa on gzipped FASTQC files'
    )
#finir de reecrire la fonction pour les arguments


def bwa(sample, output, read1, read2=None, index =None,):
    if read2 is None: read2 =''

    cmd = 'bwa mem -M -t4'
    cmd+= "{} {} {}>{}.sam\n".format(index, read1, read2, output)
    cmd+= "samtools sort -@4 -0 ban -o {0}.bam -T {0} {0}.sam\n".format(output)
    cmd+= "samtools index {}.bam\n".format(output)
    cmd+= "rm {}.sam\n.format(ouput)"

    name = 'bwa_'+ sample
    return cmd, name
def scrpitBWA(cmdScript, fh):

   fh = open("name.sh", "w")
   cmdScript = bwa()
   fh = print(cmdScript, file = fh)
   fh.close()
