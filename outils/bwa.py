import re
import sys
import glob
import os
import argparse



def init ():
    parser = argparse.ArgumentParser(
        description='This script runs BWA on gzipped FASTQ files, which are found by searching the source '
                    'directory for files ending in "R1.fastq.gz"' )

    parser.add.argument('-s' , '--source', required = True, help='Source directory')
    parser.add.argument('-o', '--output',  required=True,   help='Source directory')
    parser.add.argument('-r', '--bwa-index', dest='bwaindex',  default = None , help='Indexed reference genome')
    args = parser.parse_args()
    return args


def bwa (sample, output , read1, read2= None,index= None ,):
    if read2 is None: read2 = ''


    cmd = 'bwa mem -M -t4'
    cmd += "{} {} {} > {}.sam\n".format(index , read1 , read2, output)
    cmd += "samtools sort -@4 -O bam -o {0}.bam -T {0} {0}.sam\n".format(output)
    cmd += "samtools index {}.bam\n".format(output)
    cmd += "rm {}.san\n".format(output)



    name = 'bwa_' + sample
    return cmd,name

def bwa_sample(source, outputdir, project, samplename=None, index= None):
    # on capture le nom de l'echantillon ici
    sample_id = samplename

    #on cree le fichier output si ca n'existe pas

    if not os.path.exists(outputdir):
        os.makedirs(outputdir,0o770)

    #on capture le fichier fastq avec glob
    fastqlist = glob.glob(os.path.join(source,sample_id +"_*R1.fastq*"))+ glob.glob(os.path.join(source, sample_id + "_*R1_001.fastq*"))
    samplelist = []
    #on capture le read1 a partir du fastq et read2 a partir du read1
    #on capture aussi la lane et le barecode (on va pas utiliser dans ce script)

    for read1 in fastqlist:
        pattern = re.compile('.*/' + sample_id + '_([A-Z0-9]+)_([LS]\d+)_R1(_001)*.fastq*')
        matched +pattern.match(read1)
        lane,barecode = " , "
        if matched :
            barecode + matched.group(1)
            lane = matched.group(2)
        if 'R1_001.fastq' in read1 :
            read2 = read1.replace('R1_001.fatsq','R2_001.fastq')

        else:
            read2 = read1.replace('R1.fatsq','R2.fastq')

        # sile read2 n'existe pas , on lance avec un seul read

        if not os.path.exists(read2):
            read2 = None

        # on rajoute le sample name a la sample_list (on utilise pas vraiment mais utile)

        samplelist.append(sample_id)

        # on lance bwa qui va retourner cmd et name

        cmd, name = bwa_sample(source=args.source , outputdir=args.output,samplename=sample, index=args.bwaindex)
        print (name)
        print (cmd)
        
with open("name.sh", "w") as fh:
    cmdScript = bwa()
    fh.write('''\
    #! /bin/bash
    #SBATCH -d afterok:$PID
    #SBATCH -n 
    #SBATCH -t 72:00:00
    ''', cmdScript)
    fh.close()
