import argparse
import os
import glob
import subprocess



def init():
    parser = argparse.ArgumentParser(
        description='Finds all the BAM in the source directory and run the recommended GATK pre-processing steps')
    parser.add_argument('-s', '--source', required=True, help='Source directory')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-d', '--dbsnp', required=True,help='current dbSNP reference in VCF format')
    parser.add_argument('-b', '--bed', required=True,help='BED file containing intervals of interest')
    parser.add_argument('-r', '--ref',required=True, help='reference fasta file')
    parser.add_argument('-p', '--picard_dir',required=True, help='picard directory')
    parser.add_argument('-g', '--gatk_dir',required=True, help='gatk directory')
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


def mark_duplicate(source, outputdir, sample,picard_dir):
    tmpdir="picard_tmp"
    dedup = os.path.join(outputdir, sample + '.dedup.bam')
    metrics = os.path.join(outputdir, sample + '.dedup')

    cmd = '\nmkdir -p {}\nmkdir -p {}\n'.format(tmpdir,outputdir)
    cmd += "java -Xmx12g -Djava.io.tmpdir={} -jar {}/picard.jar MarkDuplicates".format(tmpdir, picard_dir)
    cmd += " INPUT={} OUTPUT={} METRICS_FILE={} ASSUME_SORTED=true MAX_RECORDS_IN_RAM=100000 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true".format(os.path.join(source, sample + '.bam'), dedup,
                                                                           metrics)
    cmd += "\nrm -rf {}\n\n".format(tmpdir)

    return cmd


def indel_realignment(cmd,source, output, sample, bed, dbsnp,ref, extension='.dedup.bam'):
    source = os.path.abspath(source)

    coverage = "-dt None"
    if isinstance(sample, str):
        sample = [sample]
    input=""
    input += ' -I '.join([os.path.join(source, i + extension) for i in sample])
    name = '_'.join(sample)[:200]
    intervals = name + '.intervals'
    region = "--intervals " + bed + " --interval_padding 100"    

    tmpdir="indel_realignment_tmp"
    cmd += 'cd {}\nmkdir {}\n '.format(output,tmpdir)

    cmd += "java -Xmx8g -Djava.io.tmpdir={} -jar {}/GenomeAnalysisTK.jar -T RealignerTargetCreator --disable_auto_index_creation_and_locking_when_reading_rods -nt 4 -R {}" \
           " {} -o {} {} {} ".format(tmpdir,gatk_dir, ref, input,  intervals, coverage, region)
    cmd += "\njava -Xmx4g -Djava.io.tmpdir={} -jar {}/GenomeAnalysisTK.jar --disable_auto_index_creation_and_locking_when_reading_rods -T IndelRealigner {} -nWayOut " \
           ".realigned.bam -targetIntervals {} -R {} {}  -compress 0".format(tmpdir,gatk_dir, input, intervals, ref, coverage)

    # remove intermediate files
    cmd += '\nrm -rf {}'.format(tmpdir)
    cmd += '\nrm {}\n'.format(intervals)
    cmd += 'cd ..\n\n'

    return cmd


def bqsr(cmd, source, output, sample, bed, dbsnp, ref,extension='.dedup.realigned.bam'):
   
    tmpdir="bqsr_tmp" 
    if bed:
        region = "--intervals " + bed + " --interval_padding 100"
    coverage = "-dt None"
    recaldata = os.path.join(source, sample + ".recal_data.grp")
    input = os.path.join(source, sample + extension)
    cmd += 'mkdir {}\n'.format(tmpdir)
    cmd += "java -Xmx4g -Djava.io.tmpdir={}  -jar {}/GenomeAnalysisTK.jar -T BaseRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods -nct 8 -I {} -o {}" \
          " -R {} -knownSites {} -rf BadCigar -cov ReadGroupCovariate -cov ContextCovariate -cov CycleCovariate -cov" \
          " QualityScoreCovariate {} {}".format(tmpdir,gatk_dir, input, recaldata, ref, dbsnp, coverage, region)
    recal = os.path.join(source, sample + '.processed.bam')
    cmd += "\n\njava -Xmx4g -Djava.io.tmpdir={} -jar {}/GenomeAnalysisTK.jar -T PrintReads --disable_auto_index_creation_and_locking_when_reading_rods -nct 8 -I {} -R {}" \
           " -BQSR {} -o {} -rf BadCigar {}".format(tmpdir,gatk_dir, input, ref, recaldata, recal, coverage)
    cmd += '\nrm {} {} {}\n'.format(recaldata, input, input.replace('.bam', '.bai'))
    cmd += '\nrm -rf {}'.format(tmpdir)

    return cmd



if __name__ == '__main__':
    args = init()
    source, outputdir,dbsnp, bed, ref, picard_dir, gatk_dir = args.source, args.output, args.dbsnp, args.bed, args.ref, args.picard_dir, args.gatk_dir

    bams = glob.glob(os.path.join(source, '*.bam'))
    for bamfile in bams:
        sample = os.path.basename(bamfile)[:-4]  # remove .bam
        cmd=mark_duplicate(source=source, outputdir=outputdir, sample=sample,picard_dir=picard_dir)
        cmd=indel_realignment(cmd=cmd,source=outputdir, output=outputdir, sample=sample, bed=bed, dbsnp=dbsnp, ref=ref,extension='.dedup.bam')        
        cmd=bqsr(cmd=cmd,source=outputdir, output=outputdir, sample=sample, bed=bed, dbsnp=dbsnp,ref=ref,extension='.dedup.realigned.bam')
        print2file(cmd,"process_"+sample+".sh")
        job_id = subprocess.check_output('sbatch {}'.format(script).split())
        print(job_id.decode("utf-8").strip())
