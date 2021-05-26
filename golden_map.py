#!/usr/bin/env python3
# coding: utf-8
# 2020-05-24, by Huihui Yu


import sys
import os
import glob
import time
import functools
print = functools.partial(print, flush=True)


def map(base, r1, r2, pattern, length, fasta, index, Nbases, junction, ref, assemble, expresion, threads):
    if pattern[-3:] == '.gz':
        readFilesCommand = "--readFilesCommand zcat"
    elif pattern[-4:] == '.bz2':
        readFilesCommand = "--readFilesCommand bzip2 -c"
    else: 
        readFilesCommand = ""

    ### Step 1: first round of mappping
    print ("### Step 1: first round of mappping")
    if not os.path.isdir("map1"):
        os.mkdir("map1") and exit()
    for i in range(len(base)):
        if os.path.isfile("map1/%s_Log.final.out" % base[i]): continue
        print ("-" * 50)
        print ("***Start to analyze sample: " + base[i])
        cmd = "STAR --runThreadN %d --genomeDir %s --readFilesIn %s %s %s --outFileNamePrefix map1/%s_ --alignIntronMin 20 --alignIntronMax 20000 \
--outSAMtype None --outSJfilterReads Unique --outSJfilterCountUniqueMin 10 3 3 3 --outSJfilterCountTotalMin 10 3 3 3" \
% (threads, index, r1[i], r2[i], readFilesCommand, base[i])
        print("$ " + cmd)
        res=os.system(cmd)
        if res!=0:
            if os.path.isfile("map1/%s_Log.final.out" % base[i]): os.remove("map1/%s_Log.final.out" % base[i])
            print ("Error found! Please fix the error and rerun to continue.")
            exit() 
    
    ### Step 2: generate new genome index
    print("#" * 50 + "\n### Step 2: generate new genome index")
    if not os.path.isfile("star_merge_idx/SAindex"):
        cmd = "cat map1/*_SJ.out.tab > star_merge_SJ"
        print("$ " + cmd)
        os.system(cmd) and exit()
        if not os.path.isdir("star_merge_idx"):
            os.mkdir("star_merge_idx") and exit()
        cmd = "STAR --runThreadN %d --runMode genomeGenerate --genomeDir star_merge_idx --genomeFastaFiles %s --genomeSAindexNbases %d \
--sjdbFileChrStartEnd star_merge_SJ %s --sjdbOverhang %d --limitSjdbInsertNsj 2000000" % (threads, fasta, Nbases, junction, length-1)
        print("$ " + cmd)
        res=os.system(cmd)
        if res!=0:
            print ("Error found! Please fix the error and rerun to continue.")
            exit() 
        os.rename("Log.out","star_merge_idx_log.out") and exit()
 
    ### Step 3: second round of mappping
    print("#" * 50 + "\n### Step 3: second round of mappping")
    if not os.path.isdir("map2"):
        os.mkdir("map2") and exit()
    for i in range(len(base)):
        if os.path.isfile("map2/%s_Log.final.out" % base[i]): continue
        print ("-" * 50 )
        print ("***Start to analyze sample: " + base[i])
        cmd = "STAR --runThreadN %d --genomeDir star_merge_idx --readFilesIn %s %s %s --outFileNamePrefix map2/%s_ --alignIntronMin 20 --alignIntronMax 20000  \
--limitBAMsortRAM 5000000000 --outSAMstrandField intronMotif --alignSJoverhangMin 20 --outSAMtype BAM SortedByCoordinate"  % (threads, r1[i], r2[i], readFilesCommand, base[i])
        print("$ " + cmd)
        res=os.system(cmd)
        if res!=0:
            if os.path.isfile("map2/%s_Log.final.out" % base[i]): os.remove("map2/%s_Log.final.out" % base[i])
            print ("Error found! Please fix the error and rerun to continue.")
            exit() 
        os.rename("map2/%s_Aligned.sortedByCoord.out.bam" % base[i],"map2/%s.bam" % base[i]) and exit()

    ### Step 4: Assemble novel transcripts
    if assemble == True:
        print("#" * 50 + "\n### Step 4: assemble novel transcripts")
        for i in range(len(base)):
            if os.path.isfile("map2/%s.gtf" % base[i]): continue
            print("$ stringtie -j 5 -c 5 -g 10 -p %d -G %s map2/%s.bam -o map2/%s.gtf" % (threads, ref, base[i], base[i]))
            res=os.system("stringtie -j 5 -c 5 -g 10 -p %d -G %s map2/%s.bam -o map2/%s.gtf" % (threads, ref, base[i], base[i]) )
            if res!=0:
                if os.path.isfile("map2/%s.gtf" % base[i]): os.remove("map2/%s.gtf" % base[i])
                print("Error found! Please fix the error and rerun to continue.")
                exit()
        if  not os.path.isfile("merged.gtf"):
            gtf_files=glob.glob("map2/*.gtf")
            for x in gtf_files:
                if os.popen("grep -v '^#' %s" % x).readline() == "":
                    print("WARNING: no transcripts were found in GTF file " + x)
                    gtf_files.remove(x)
            if len(gtf_files)==0:
                sys.stderr.write("Error: No non-zero size GTF files!\n"); exit()
            cmd="stringtie --merge -G %s -o merged.gtf %s" % (ref," ".join(gtf_files)) 
            print("-" * 50+"\n"+cmd)
            res=os.system(cmd)
            if res!=0:
                if os.path.isfile("merged.gtf"): os.remove("merged.gtf")
                print("Error found! Please fix the error and rerun to continue.")
                exit()
 
    ### Step 5: estimate gene expression
    if expresion == True:
        print ("#" * 50 + "\n### Step %d: estimate gene expression" % (4 + assemble))
        for i in range(len(base)):
            if os.path.isfile("ballgown/%s/%s.gtf" % (base[i],base[i])): continue
            if assemble == True:
                cmd = "stringtie -e -B -p %d -G merged.gtf -o ballgown/%s/%s.gtf map2/%s.bam" % (threads, base[i], base[i], base[i])
            else:
                cmd = "stringtie -e -B -p %d -G %s -o ballgown/%s/%s.gtf map2/%s.bam" % (threads, ref, base[i], base[i], base[i])            
            print("$ " + cmd)
            res=os.system(cmd)
            if res!=0:
                if os.path.isfile("ballgown/%s/%s.gtf" % (base[i],base[i])): os.remove("ballgown/%s/%s.gtf" % (base[i],base[i]))
                print("Error found! Please fix the error and rerun to continue.")
                exit()


def main(options):
    # chech input files
    r1 = []; r2 = []; base = []
    if os.path.isdir(options.input):
        fs = [os.path.join(options.input, "*" + options.pattern),os.path.join(options.input, "*" + options.pattern.replace("1","2"))] 
        r1 = glob.glob(fs[0])
        r2 = glob.glob(fs[1])
        for x in r1+r2:
            if os.path.getsize(x)==0:
                sys.stderr.write("Error: The file is empty:\n\t" + x + "\n"); exit()
    elif os.path.isfile(options.input):
        fs = open (options.input, 'r')
        for line in fs.read().splitlines():
            line = line.strip()
            if line == "": continue
            if not os.path.isfile(line):
                sys.stderr.write("Error: The file does not exist:\n\t" + line + "\n"); exit()
            elif os.path.getsize(line) == 0:
                sys.stderr.write("Error: The file is empty:\n\t" + line + "\n"); exit()
            if line[-len(options.pattern):] ==  options.pattern:  r1.append(line)
            if line[-len(options.pattern):] ==  options.pattern.replace("1","2"): r2.append(line)
        fs.close()
    else:
        sys.stderr.write("Error: No such file or directory:\n\t" + options.input + "\n"); exit()
   
    if len(r1) == 0 or len(r2) ==0 or len(r1) != len(r2):
        if ("1" in options.pattern):
            p="*" + options.pattern + "; *" + options.pattern.replace("1","2") 
        else:
            p="*"+ options.pattern
        sys.stderr.write("Error: No such pattern files or some of the files unpaired:\n\t\t" + p + '\n');
        exit()

    r1.sort(); r2.sort()
    for i in range(len(r1)):
        p1 = r1[i][:-len(options.pattern)]
        p2 = r2[i][:-len(options.pattern)]
        if (p1 != p2):
            sys.stderr.write("Error: Some of the files unpaired!\n" ); exit()
        else:
            base.append(os.path.basename(p1))
    if r1 == r2:
        for i in range(len(r2)): r2[i] = ""

    map (base, r1, r2, options.pattern, options.length, options.fasta, options.index, options.Nbases, \
options.junction, options.gtf, options.assemble, options.expresion, options.threads)



if __name__ == "__main__":
    from optparse import OptionParser    
    parser = OptionParser()
    parser.add_option("-i", dest="input", default=".", help="input file directory or file list (.)")
    parser.add_option("-p", dest="pattern", default="_1P.fq.gz", help="input file pattern (_1P.fq.gz)")
    parser.add_option("-l", dest="length", type=int, default=101, help="read length (101)")
    parser.add_option("-f", dest="fasta", help="reference genome fasta file")
    parser.add_option("-x", dest="index", help="STAR genome index directory")
    parser.add_option("-n", dest="Nbases", type=int, default=14, help="length of the STAR SA indexing string (14)")
    parser.add_option("-j", dest="junction", default="", help="file for the splice junction introns (\"\")")
    parser.add_option("-g", dest="gtf", help="reference annotation GTF file")
    parser.add_option("-@", dest="threads", type=int, default=30, help="number of threads (30)")
    parser.add_option("-e", "--expresion", action="store_true", default=False, help="estimate gene expression")
    parser.add_option("-a", "--assemble", action="store_true",default=False, help="assemble novel transcripts")
    
    (options, args) = parser.parse_args()

    print('options:')
    for (key, value) in vars(options).items():
        print("\t--%s"%key, "=", value )

    if options.index == None:
        sys.stderr.write("Error: No STAR genome index directory specified!\n"); exit()
    elif not os.path.isdir(options.index):
        sys.stderr.write("Error: No such directory:\n\t" + options.index + "\n"); exit()
    if not (options.junction == "" or os.path.isfile(options.junction)):
        sys.stderr.write("Error: No such file:\n\t" + options.junction + "\n"); exit()
    if options.fasta == None:
        sys.stderr.write("Error: No genome fasta file specified!\n"); exit()
    elif not os.path.isfile(options.fasta):
        sys.stderr.write("Error: No such file:\n\t" + options.fasta + "!\n"); exit()
    if (options.expresion == True) or (options.assemble == True): 
        if options.gtf == None:
            sys.stderr.write("Error: No reference GTF file specified!\n"); exit()
        elif not os.path.isfile(options.gtf):
            sys.stderr.write("Error: No such file:\n\t" + options.gtf + "!\n"); exit()
    
    print("\nStart time is: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime()))
    print ("#"*50)
    main(options)
    print ("#"*50)
    print ("Stop time is: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime()))

