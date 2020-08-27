import re
import sys
import glob
import os
import argparse
import subprocess

def init():
    parser = argparse.ArgumentParser(
        description='Runs the GATK module HaplotypeCaller ')
    parser.add_argument('-s', '--source', required=True, help='folder containing bam files')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-r', '--reference',required=True, help=' FASTA reference sequence to which BAM was aligned')
    parser.add_argument('-d', '--dbsnp',required=True,
                        help='current dbSNP reference in VCF format (if sample on multiple lanes and need to merge)')
    parser.add_argument('-g', '--bed',required=True,
                        help='BED file containing intervals of interest (if sample on multiple lanes and need to merge)'
                        )
    options = parser.parse_args()
    return options


def print2file(cmd, filename):
    f = open(filename, "w")
    f.write("#!/bin/bash\n")
    f.write("#SBATCH -J {} \n".format(filename))
    f.write("#SBATCH -p all\n")
    f.write("#SBATCH -n 1 -c 8\n")
    f.write("#SBATCH --mem=8gb\n")
    f.write("#SBATCH -t 72:00:00\n")
    f.write("#SBATCH -o %x-%j.out\n")
    f.write(cmd)
    f.close()

def haplotype(bam, outputdir, ref=None, bed=None, dbsnp=None,  ext='.HaplotypeC.raw.vcf'):
    tmpdir = "haplotype_tmp"
    cmd = '\nmkdir -p {}\nmkdir -p {}\n'.format(tmpdir,outputdir)

    match = re.match('(.+?)(\.processed)*\.bam', os.path.basename(bam))
    if not match:
        print('ERROR: bam file name should end with ".bam" or ".processed.bam"', file=sys.stderr)
        exit(1)
    name = match.group(1)
    args=""
    args += " --intervals " + bed + " --interval_padding 100"


    raw = os.path.join(outputdir, name + ext)

    cmd += "java -Xmx4g -Djava.io.tmpdir={} -jar $gatk_dir/GenomeAnalysisTK.jar -T HaplotypeCaller" \
          " {} -nct 4 -R {} -I {} -o {} -D {}  -stand_emit_conf 10 -rf BadCigar ".\
        format(tmpdir, args, ref, bam, raw, dbsnp)
    cmd += '\nrm -rf {}'.format(tmpdir)

    return cmd,name


if __name__ == '__main__':
    arguments = init()
    source = arguments.source
    outputdir = arguments.output
    dbsnp = arguments.dbsnp
    ref = arguments.reference
    bed = arguments.bed

    bams = glob.glob(os.path.join(source, '*.bam'))
    for bam in bams:
        cmd,name=haplotype(bam=bam, outputdir=outputdir, ref=ref, bed=bed, dbsnp=dbsnp)
        print2file(cmd,"haplotype_"+name+".sh")
        job_id = subprocess.check_output('sbatch {}'.format(script).split())
        print(job_id.decode("utf-8").strip())
