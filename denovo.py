#!/usr/bin/python3

"""Pipeline for ATAC seq and ChIP seq

Usage:
 denovo.py align-pe <fqfile1> <fqfile2> REF [-t <x>] -o OUTPUT
 denovo.py align-se <fqfile> REF [-t <x>] -o OUTPUT
 denovo.py assembly-pe <home_dir> <fqfile1> <fqfile2> <out_dir> [REF]
 denovo.py assembly-se <home_dir> <fqfile1> <out_dir> [REF]
 denovo.py (-h | --help)


Options:
  -h --help                 Show this screen.
   REF                      Full path to reference genome file.
  -t --threads=<x>          Number of threads to use for alignment [default: 2].
  -o --output=OUTPUT        Output file in BAM format for arguments align-pe, extension '.txt' when argument annotate (IMPORTANT: only use '_' in file names).
  <home_dir>                Full path to home directory of stored tools.
  <out_dir>                 Name of directory to store assembler output.
"""


from __future__ import print_function, division
from  subprocess import *
#from statistics import *
#from numpy import *
import sys, os, time, docopt, errno, itertools, gzip, re


dir_list = ['spades.py']

########################################################################################################################

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as ex:
        if ex.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def dir_names(homdir, *args):
    tdir = []
    ndir_list = dir_list
    try:
        for i in ndir_list:
            dc = Popen(["find", homdir, "-name", i], stdout=PIPE)
            work_dir = dc.stdout.read().rstrip()
            work_dir = work_dir.decode('utf-8')
            tdir.append(work_dir)
            dc.kill()
        tdir = tdir + list(args)
        #print(tdir)
        for i in tdir:
            if len(i) == 0:
                print("Cannot find %s, Please enter the correct parent directory for the mentioned tools: %s" % (
                str(i), str(ndir_list[0])))
                sys.exit(-1)
        return tdir

    except CalledProcessError as cpe:
        print('Error %d.' % cpe.returncode)
        sys.exit(-1)



########################
# ALIGN READS WITH BWA #
########################


def bwape(fq1, fq2, ref, th, bamfile):
    bam_sort = bamfile.rstrip('.bam') + "_sorted.bam"
    stat_file = bamfile.rstrip('.bam') + "_stat_mem_pe.txt"
    align_mem = Popen(
        "bwa mem -M -t %d %s %s %s | samtools view -@ %d -bh - >%s" % (int(th), ref, fq1, fq2, int(th), bamfile),
        bufsize=-1, shell=True, executable='/bin/bash')
    align_mem.wait()
    samstat = Popen("samtools sort -@ 5 %s -o %s" % (bamfile, bam_sort), bufsize=-1, shell=True, executable='/bin/bash')
    samstat.wait()
    print("Indexing ...")
    Popen("samtools index %s" % (bam_sort), bufsize=-1, shell=True, executable='/bin/bash').wait()
    samstat = Popen("samtools flagstat %s > %s" % (bam_sort, stat_file), bufsize=-1, shell=True, executable='/bin/bash')
    samstat.wait()

    ##Extract reads
    Popen("samtools view -bh -f 4 %s > %s " % (bam_sort, ''.join([bam_sort.rstrip('.bam'), '_unmapped.bam'])), bufsize=-1, shell=True, executable='/bin/bash').wait()
    Popen("bamToFastq  -i %s -fq %s -fq2 %s " % (''.join([bam_sort.rstrip('.bam'), '_unmapped.bam']),
                                         ''.join([bam_sort.rstrip('.bam'), '_unmapped_R1.fastq']),''.join([bam_sort.rstrip('.bam'), '_unmapped_R2.fastq'])),
          bufsize=-1, shell=True, executable='/bin/bash').wait()





########################################################################################################################



def bwase(fq, ref, th, bamfile):
    bam_sort = bamfile.rstrip('.bam') + "_sorted.bam"
    stat_file = bamfile.rstrip('.bam') + "_stat_mem_se.txt"

    align_mem = Popen(
        "bwa mem -t %d %s %s | samtools view -@ %d -bh - > %s" % (int(th), ref, fq, int(th), bamfile),
        bufsize=-1, shell=True, executable='/bin/bash')
    align_mem.wait()
    samstat = Popen("samtools sort -@ 5 %s -o %s" % (bamfile, bam_sort), bufsize=-1, shell=True, executable='/bin/bash')
    samstat.wait()
    print("Indexing ...")
    Popen("samtools index %s" % (bam_sort), bufsize=-1, shell=True, executable='/bin/bash').wait()
    samstat = Popen("samtools flagstat %s > %s" % (bamfile, stat_file), bufsize=-1, shell=True, executable='/bin/bash')
    samstat.wait()

    #Extract reads
    Popen("samtools view -bh -f 4 %s > %s " % (bam_sort, ''.join([bam_sort.rstrip('.bam'), '_unmapped.bam'])),
          bufsize=-1, shell=True, executable='/bin/bash').wait()
    Popen("bamToFastq  -i %s -fq %s " % (''.join([bam_sort.rstrip('.bam'), '_unmapped.bam']),
                                                 ''.join([bam_sort.rstrip('.bam'), '_unmapped_R1.fastq'])),
          bufsize=-1, shell=True, executable='/bin/bash').wait()


def spadepe(homdir,fq1, fq2,outdir, ref):
    ndir = dir_names(homdir)

    if args["REF"]:
        asm = Popen("%s --careful --only-assembler --trusted-contigs %s -1 %s -2 %s -o %s" % (
               ndir[0],ref, fq1, fq2,outdir), bufsize=-1, shell=True,
            executable='/bin/bash')
        asm.wait()

    else:
        asm = Popen("%s --careful --only-assembler -1 %s -2 %s -o %s" % (
            ndir[0], fq1, fq2, outdir), bufsize=-1, shell=True,
                    executable='/bin/bash')
        asm.wait()

def spadese(homdir,fq1, outdir, ref):
    ndir = dir_names(homdir)

    if args["REF"]:
        asm = Popen("%s --careful --only-assembler --trusted-contigs %s -s %s -o %s" % (
               ndir[0],ref, fq1, outdir), bufsize=-1, shell=True,
            executable='/bin/bash')
        asm.wait()

    else:
        asm = Popen("%s --careful --only-assembler -s %s -o %s" % (
            ndir[0], fq1, outdir), bufsize=-1, shell=True,
                    executable='/bin/bash')
        asm.wait()







#####################################
#############    MAIN   #############
#####################################


if __name__ == '__main__':
    try:
        args = docopt.docopt(__doc__, version='denovo-v1.0')
        if args['align-pe']:
            bwape(args['<fqfile1>'], args['<fqfile2>'], args['REF'], args['--threads'], args['--output'])

        elif args['align-se']:
            bwase(args['<fqfile1>'], args['REF'], args['--threads'], args['--output'])

        elif args['assembly-pe']:
            spadepe(args['<home_dir>'], args['<fqfile1>'], args['<fqfile2>'], args['--out_dir'],args['REF'])

        elif args['assembly-se']:
            spadese(args['<home_dir>'], args['<fqfile1>'], args['--out_dir'],  args['REF'])

    except docopt.DocoptExit:
        print("Operation Failed (Check Arguments) !! for help use: python ATAC-seq.py -h or --help")
