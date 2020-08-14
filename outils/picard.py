#!/bin/bash
import re
import sys
import glob
import os
import argparse
import subprocess 

def init():
    parser = argparse.ArgumentParser(
        description = 'MultiQC submodule to parse output from Picard MarkDuplicates')
    parser.add.argument('-s','--source',required = True ,help ='Source directory')
    parser.add.argument('-o', '--output', required=True, help='Source directory')
    parser.add.argument('-r', '--markduplicates',dest ='picard-markduplicates ', default = None, help='Indexed')
    args = parser.parse_args()
    return args


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



# Mark duplicates and sort
def markDuplicates(sample, output, read1, read2=None, index =None,):
   cmd ='java $JAVA_OPT'
   cmd +=  '-jar $MRKDUP'
   cmd += "{} {} {}>{}.bam\n".format(index, read1, read2, output)

   cmd +=  'INPUT=$OUTDIR/${read1}_${read2}_${output}-s.bam '.format(index, read1, read2, output)
   cmd +=  'OUTPUT=$OUTDIR/${read1}_${read2}_${output}-smd.bam '
   cmd +=  'METRICS_FILE=$OUTDIR/${read1}_${read2}_${output}-smd.metrics'
   cmd +=  'AS=TRUE '
   cmd +=  'VALIDATION_STRINGENCY=LENIENT'
   cmd +=  'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 '
   cmd +=   'REMOVE_DUPLICATES=TRUE'
   cmd +=  'samtools sort $OUTDIR/${read1}_${read2}_${output}-smd.bam $OUTDIR/${read1}_${read2}_${output}-smds'
   cmd += 'samtools index $OUTDIR/${read1}_${read2}_${output}-smds.bam'
   cmd +=  'samtools flagstat $OUTDIR/${read1}_${read2}_${output}-smds.bam > $OUTDIR/${read1}_${read2}_${output}-smds.flagstat'.format(index, read1, read2, output)
   name = 'markduplicates'+ sample+".bam"
   return cmd , name



if __name__ == '__main__':
    args = init()
    pattern = re.compile('.*/(.+)_([A-Z0-9]+)_([LS]\d+)_R1(_001)*.bam*')
    sample_list = set()
    for bam in glob.glob(os.path.join(args.source, "*markduplicates.bam*"))) +
glob.glob(os.path.join(args.source,"*markduplicates.bam*")):
        matched = pattern(bam)
        if matched :
            sample = matched.group(1)

        else :
            splitor = 'markduplicates.bam' if 'markduplicates.bam' in fastq else 'markduplicates.bam'
            sample = os.path.basename(bam).split(splitor)[0]
        sample_list.add(sample)
    for sample in sample_list:
        cmd,script = bam_sample(source=args.source,
outputdir=args.output,samplename=sample, index=args.bamindex)

        print2file(cmd,script)
        job_id = subprocess.check_output('sbatch {}'.format(script).split())
        print(job_id.decode("utf-8").strip())
