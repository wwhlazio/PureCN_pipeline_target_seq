#!/usr/bin/python
import os
import subprocess
import re

def run_fastqc_raw():                                                             #run fastqc on the fastq files
    datadir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/data/"
    workdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/workdir/fastqc/"
    outdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/fastqc/"
    func = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/bin/imple_fastqc.py"
    acc = "premium"
    mem1 = "4000"
    time1 = "1:00"
 
    inputfiles = workdir + "file_list"
    print(inputfiles)
    f = open(inputfiles,'r')
    samp = []
    samf = []
    ns = 0
    for line in f:
            tmp = line.split('\n')
            samf.append(tmp[0])
            tmp1 = line.split('.fastq.gz')
            samp.append(tmp1[0])
            ns = ns+1
    f.close()
    print("number of samples:%d\n" % ns)
    print(os.path.exists(workdir))
    if(not os.path.exists(workdir)):
        cmd="mkdir " + workdir
        subprocess.Popen(cmd.split())
    sh_file = workdir + "run_fastqc_raw.sh"
    fx = open(sh_file,'w')
    for i in range(ns):
        inputfile = datadir + samf[i]
        samp_lsf = workdir + "fastqc_" + samp[i] + ".lsf"
        f=open(samp_lsf,'w')
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -J fastqc_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %sfastqc_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %sfastqc_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")
        f.write("module load fastqc\n")
        f.write("python %s \"%s\" \"%s\"\n" % (func,outdir,inputfile))
        f.close()
        fx.write("bsub < %s\n" % samp_lsf)
    fx.close()
def run_cutadpt_filter():                                                         #filter reads contaminated by adapter, low quality reads, reads totally the same
    datadir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/data/"
    workdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/workdir/filter/"
    outdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/cut_adapter/"
    bindir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/bin/"
    bin2 = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/bin/reads_clear_hash_v5x.pl"
    inputfiles = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/workdir/filter/file_list_redo"
    thr = "39"
    acc = "premium"
    mem1 = "50000"
    time1 = "6:00"
    
    adapter_seq1 = "AGATCGGAAGAG"  #experiment specific adapter sequence
    ad_len = 7                                           #cutadapt paramters to claim a adapter contaminated sequence
    n_freq = 0.1                                    #freqency of N to filter a read  
 
    #inputfiles = datadir + "file_list"
    print(inputfiles)
    f = open(inputfiles,'r')
    samp = []
    samf = []
    ns = 0
    for line in f:
            tmp = line.split('\n')
            samf.append(tmp[0])
            tmp1 = line.split('.fastq.gz')
            tmp1 = re.split('.R1.fastq.gz|.R2.fastq.gz',line)
            if tmp1[0] not in samp:
                samp.append(tmp1[0])
                ns = ns+1
    f.close()
    print("number of samples:%d\n" % ns)
    print(os.path.exists(workdir))
    if(not os.path.exists(workdir)):
        cmd="mkdir " + workdir
        subprocess.Popen(cmd.split())
    sh_file = workdir + "run_cutadapt.sh"
    fx = open(sh_file,'w')
    for i in range(ns):
        inputfile1 = datadir + samp[i] + ".R1.fastq.gz"
        inputfile2 = datadir + samp[i] + ".R1.fastq.gz"
        outputfile11 = outdir + samp[i] + ".R1_cutadapt.fastq.gz"
        outputfile12 = outdir + samp[i] + ".R2_cutadapt.fastq.gz"
        midfile1 = outdir + samp[i] + ".R1_cutadapt.fastq"
        midfile2 = outdir + samp[i] +  ".R2_cutadapt.fastq"
        outputfile21 = outdir + samp[i] + ".R1_cutadapt.not_pass.fastq"
        outputfile22 = outdir + samp[i] + ".R2_cutadapt.not_pass.fastq"
        outputfile31 = outdir + samp[i] + ".R1_cutadapt.pass.fastq"
        outputfile32 = outdir + samp[i] + ".R2_cutadapt.pass.fastq"
        logfile = workdir + samp[i] + "_cutadapt.log"
        samp_lsf = workdir + "cutadapt_" + samp[i] + ".lsf"
        f=open(samp_lsf,'w')
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -J cutadapt_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %scutadapt_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %scutadapt_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")
        f.write("cd %s\n" % outdir)
        f.write("cutadapt -b \"A{100}\" -b \"T{100}\" -b \"%s\" -B \"A{100}\" -B \"T{100}\" -B \"%s\" --discard-trimmed -O %s --max-n %f -o %s -p %s %s %s > %s\n" %(adapter_seq1,adapter_seq1,ad_len,n_freq,outputfile11,outputfile12,inputfile1,inputfile2,logfile))
        f.write("gunzip %s\n" % outputfile11)
        f.write("gunzip %s\n" % outputfile12)
        f.write("perl %s %s %s %s %s %s %s %s\n" %(bin2,midfile1,midfile2,thr,outputfile21,outputfile22,outputfile31,outputfile32))
        f.write("rm %s\n" % midfile1)
        f.write("rm %s\n" % midfile2)
        f.close()
        fx.write("bsub < %s\n" % samp_lsf)
    fx.close()


def run_fastqc(datadir,workdir,outdir):
    #datadir = "/sc/orga/projects/zhuj05a/Wenhui/single_cell_RNAseq/data/"
    #workdir = "/sc/orga/projects/zhuj05a/Wenhui/single_cell_RNAseq/workdir/fastqc_raw/"
    #outdir = "/hpc/users/wangm08/Wenhui/single_cell_RNAseq/fastqc_raw/"
    func = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/bin/imple_fastqc.py"
    acc = "premium"
    mem1 = "4000"
    time1 = "1:00"
 
    inputfiles = datadir + "file_list"
    print(inputfiles)
    f = open(inputfiles,'r')
    samp = []
    samf = []
    ns = 0
    for line in f:
            tmp = line.split('\n')
            samf.append(tmp[0])
            tmp1 = re.split('_cutadapt.pass.fastq',line)
            samp.append(tmp1[0])
            ns = ns+1
    f.close()
    print("number of samples:%d\n" % ns)
    print(os.path.exists(workdir))
    if(not os.path.exists(workdir)):
        cmd="mkdir " + workdir
        subprocess.Popen(cmd.split())
    sh_file = workdir + "run_fastqc.sh"
    fx = open(sh_file,'w')
    for i in range(ns):
        inputfile = datadir + samf[i]
        samp_lsf = workdir + "fastqc_" + samp[i] + ".lsf"
        f=open(samp_lsf,'w')
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J fastqc_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -o fastqc_%s.stdout\n" % samp[i])
        f.write("#BSUB -eo fastqc_%s.stderr\n" % samp[i])
        f.write("#BSUB -L /bin/bash\n")
        f.write("module load fastqc\n")
        f.write("python %s \"%s\" \"%s\"\n" % (func,outdir,inputfile))
        f.close()
        fx.write("bsub < %s\n" % samp_lsf)
    fx.close()    

    
def collect_fastqc(datadir,workdir):                   #collect fastqc results into a matrix of metric
    
    output=datadir + "collected_summary";
    
    inputfiles = datadir + "file_list"
    print(inputfiles)
    f = open(inputfiles,'r')
    samp = []
    samf = []
    ns = 0
    for line in f:
            tmp = line.split('\n')
            samf.append(tmp[0])
            tmp1 = line.split('_fastqc.zip')
            samp.append(tmp1[0])
            ns = ns+1
    f.close()
    print("number of samples:%d\n" % ns)
    
    
    fx = open(output,'w')
    for i in range(ns):
        tfile = datadir + samf[i]
        cmd = "unzip "+ inputfile + " -d " + workdir
        #print(cmd)
        os.system(cmd)
        inputfile2 = workdir + samp[i] + "_fastqc/summary.txt"
        cmd = "cat " + inputfile2
        a=subprocess.check_output(cmd,shell=True).decode("utf-8")
        b=a.split('\n')
        for i in range(len(b)-1):
             c=b[i].split('\t')[0]
             fx.write("%s\t" % c)
        fx.write("\n")
    fx.close()  

def run_bwa(datadir,workdir,outdir):     #run reads alignment with bwa, filter duplicate reads, reads quality recalibration
    
    inputfiles=datadir + "file_list"
    samp=[]
 
    ns=0
    f=open(inputfiles,'r')
    
    for line in f:
        tmp1 = line.split(".")
        if tmp1[0] not in samp:
            samp.append(tmp1[0])
            ns=ns+1
    f.close()
    shfile=workdir + 'run_bwa.sh'
  
    acc = "premium"
    mem1 = "25000"
    time1 = "4:00"
    ref = "/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/genome_dir/b37/human_g1k_v37_decoy.fasta"
    dbgap = "/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/genome_dir/b37/dbsnp_138.b37.vcf"
    mills = "/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/genome_dir/b37/Mills_and_1000G_gold_standard.indels.b37.vcf"

    fx=open(shfile,'w')
    
    for i in range(ns):
        input1=datadir + samp[i] + ".R1_cutadapt.pass.fastq"
        input2=datadir + samp[i] + ".R2_cutadapt.pass.fastq"
        samp_lsf = workdir + "bwa_" + samp[i] + ".lsf"
        f=open(samp_lsf,'w')
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J bwa_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -o bwa_%s.stdout\n" % samp[i])
        f.write("#BSUB -eo bwa_%s.stderr\n" % samp[i])
        f.write("#BSUB -L /bin/bash\n")
    
        f.write("cd %s\n" % outdir)
        f.write("module load samtools\n")
        f.write("module load bwa\n")
        f.write("module load java\n")
        
        f.write("bwa mem -M -t 8 -p %s %s %s > %s.sam\n" % (ref,input1,input2,samp[i]))
        f.write("samtools view -bS %s.sam > %s.bam\n\n" % (samp[i],samp[i]))
        

        f.write("java -jar /hpc/packages/minerva-common/picard/1.93/bin/AddOrReplaceReadGroups.jar INPUT=%s.bam OUTPUT=%s.addgroup.bam RGID=%s RGSM=%s RGLB=DNA RGPL=illumina RGPU=none SO=coordinate\n" % (samp[i],samp[i],samp[i],samp[i]))
        f.write("java -Xmx8G -jar /hpc/packages/minerva-common/picard/1.93/bin/MarkDuplicates.jar INPUT=%s.addgroup.bam OUTPUT=%s.addgroup.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=%s.output.metrics\n" % (samp[i],samp[i],samp[i]))
        f.write("java -jar /hpc/packages/minerva-common/gatk/3.6-0/src/GenomeAnalysisTK.jar -T BaseRecalibrator -R %s -I %s.addgroup.dedupped.bam -knownSites %s -knownSites %s -o %s.recal.table\n" % (ref,samp[i],dbgap,mills,samp[i]))
        f.write("java -jar /hpc/packages/minerva-common/gatk/3.6-0/src/GenomeAnalysisTK.jar -T PrintReads -R %s -I %s.addgroup.dedupped.bam -BQSR %s.recal.table -o %s.addgroup.dedupped.recal.bam\n" % (ref,samp[i],samp[i],samp[i]))
        f.write("rm %s.sam\n" % samp[i])
        f.write("rm %s.bam\n" % samp[i])
        f.write("rm %s.addgroup.bam\n" % samp[i])
        f.write("rm %s.addgroup.dedupped.bam\n" % samp[i])
        f.close()
        fx.write("bsub < %s\n" % samp_lsf)
    fx.close()

def run_mutect(datadir,workdir,outdir):   #run mutect to detect snps

    ref="/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/genome_dir/b37/human_g1k_v37_decoy.fasta"
    dbgap="/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/genome_dir/b37/dbsnp_138.b37.vcf"
    cos="/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/genome_dir/b37/cosmic_coding_and_noncoding_sorted2.vcf"
    
    inputfiles=datadir + "file_list"
    
    f=open(inputfiles,'r')
    samp=[]
    ns=0
    
    for line in f:
        tmp1 = line.split(".")
        if tmp1[0] not in samp:
            samp.append(tmp1[0])
            ns=ns+1
    f.close()
    

    shfile=workdir + 'run_mutect.sh'

    acc = "premium"
    mem1 = "25000"
    time1 = "4:00"
    

    fx=open(shfile,'w')

    for i in range(ns):
        input1=datadir + samp[i] + ".addgroup.dedupped.recal.bam"
        samp_lsf = workdir + "mutect_" + samp[i] + ".lsf"
        f=open(samp_lsf,'w')
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J mutect_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -o mutect_%s.stdout\n" % samp[i])
        f.write("#BSUB -eo mutect_%s.stderr\n" % samp[i])
        f.write("#BSUB -L /bin/bash\n")
        
        f.write("module load java/1.7.0_51\n")
        f.write("cd %s\n" % outdir )
        f.write("java -jar /hpc/packages/minerva-centos7/mutect/1.1.6/muTect-1.1.6-10b1ba92.jar \\\n")
        f.write("    --analysis_type MuTect \\\n")
        f.write("    -R %s \\\n" % ref)
        f.write("    --dbsnp %s \\\n" % dbgap)
        f.write("    --cosmic %s \\\n" % cos)
        f.write("    -I:tumor %s \\\n" % input1)
        f.write("    -o %s%s_mutect_stats.txt \\\n" % (outdir,samp[i]))
        f.write("    -vcf %s%s_mutect.vcf \n" % (outdir,samp[i]))
        f.close()
        fx.write("bsub < %s\n" % samp_lsf)
    fx.close()
        
def run_samtools(datadir,workdir,outdir):                 #run samtools mpileup on the bam files
    
    ref="/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/genome_dir/b37/human_g1k_v37_decoy.fasta"
    
    inputfiles=datadir + "bam_file_list"

    f=open(inputfiles,'r')
    samp=[]
    ns=0

    for line in f:
        tmp1 = line.split(".")
        if tmp1[0] not in samp:
            samp.append(tmp1[0])
            ns=ns+1
    f.close() 


    shfile=workdir + 'run_samtools.sh'

    acc = "premium"
    mem1 = "15000"
    time1 = "1:00"
    
    fx=open(shfile,'w')

    for i in range(ns):
        input1=datadir + samp[i] + ".addgroup.dedupped.recal.bam"
        samp_lsf = workdir + "samtools_" + samp[i] + ".lsf"
        f=open(samp_lsf,'w')
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J satmools_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -o samtools_%s.stdout\n" % samp[i])
        f.write("#BSUB -eo samtools_%s.stderr\n" % samp[i])
        f.write("#BSUB -L /bin/bash\n")

        f.write("module load samtools\n")
        f.write("cd %s\n" % outdir )
        f.write("samtools mpileup -f %s %s > %s/%s.mpileup\n" % (ref,input1,outdir,samp[i]))
        f.write("samtools view -c -F 4 %s > %s/%s.readnumber\n" % (input1,outdir,samp[i]))
        f.close()
        fx.write("bsub < %s\n" % samp_lsf)
    fx.close()

def extract_coverage(datadir,workdir,outdir):    #collect covrage information for cnv detection
   
    inputfiles=datadir + "bam_file_list"

    f=open(inputfiles,'r')
    samp=[]
    ns=0

    for line in f:
        tmp1 = line.split(".")
        if tmp1[0] not in samp:
            samp.append(tmp1[0])
            ns=ns+1
    f.close()


    shfile=workdir + 'run_coverage.sh'
    
    acc = "premium"
    mem1 = "15000"
    time1 = "1:00"

    fx=open(shfile,'w')

    for i in range(ns):
        input1=datadir + samp[i] + ".addgroup.dedupped.recal.bam"
        samp_lsf = workdir + "coverage_" + samp[i] + ".lsf"
        f=open(samp_lsf,'w')
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J coverage_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -o coverage_%s.stdout\n" % samp[i])
        f.write("#BSUB -eo coverage_%s.stderr\n" % samp[i])
        f.write("#BSUB -L /bin/bash\n")

        f.write(" export PURECN=\"/hpc/packages/minerva-centos7/rpackages/bioconductor/3.8/PureCN/extdata\"\n")
        f.write("module load R\n")
        f.write("Rscript $PURECN/Coverage.R \\\n")
        f.write("--outdir %s \\\n" % outdir)
        f.write("--bam %s \\\n" % input1)
        f.write("--intervals /sc/arion/projects/zhuj05a/Wenhui/bladder/purity/purecn/baits_g1k_intervals.txt\n")
        f.close()
        fx.write("bsub < %s\n" % samp_lsf)
    fx.close()

 
def run_purecn(datadir,datadir2,workdir,outdir):  #run purecn 
    
    inputfiles=datadir + "file_list"

    f=open(inputfiles,'r')
    samp=[]
    ns=0

    for line in f:
        tmp1 = line.split(".")
        if tmp1[0] not in samp:
            samp.append(tmp1[0])
            ns=ns+1
    f.close()


    shfile=workdir + 'run_purecn.sh'

    acc = "premium"
    mem1 = "15000"
    time1 = "4:00"

    fx=open(shfile,'w')

    for i in range(ns):
        input1=datadir + samp[i] + ".addgroup.dedupped.recal_coverage_loess.txt"
        samp_lsf = workdir + "purecn_" + samp[i] + ".lsf"
        f=open(samp_lsf,'w')
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J purecn_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -o purecn_%s.stdout\n" % samp[i])
        f.write("#BSUB -eo purecn_%s.stderr\n" % samp[i])
        f.write("#BSUB -L /bin/bash\n")

        f.write("export PURECN=\"/hpc/packages/minerva-centos7/rpackages/bioconductor/3.8/PureCN/extdata\"\n")
        f.write("module load R\n")
        f.write("Rscript $PURECN/PureCN.R \\\n")
        f.write("--out %s%s \\\n" % (outdir,samp[i]))
        f.write("--tumor %s \\\n" % input1)
        f.write("--sampleid %s \\\n" % samp[i])
        f.write("--vcf %s%s_mutect.vcf \\\n" % (datadir2,samp[i]))
        f.write("--statsfile %s%s_mutect_stats.txt \\\n" % (datadir2,samp[i])) 
        f.write("--intervals /sc/arion/projects/zhuj05a/Wenhui/bladder/purity/purecn/baits_g1k_intervals.txt \\\n")
        f.write("--normaldb /sc/arion/projects/zhuj05a/Wenhui/bladder/purity/purecn/normal/normalDB_hg19.rds \\\n")
        f.write("--normal_panel /sc/arion/projects/zhuj05a/Wenhui/bladder/purity/purecn/normal/mapping_bias_hg19.rds \\\n")
        f.write("--genome hg19 \\\n")
        f.write("--force --postoptimize --seed 123\n")
        f.close()
        fx.write("bsub < %s\n" % samp_lsf)
    fx.close()


    

if __name__=="__main__":
    run_fastqc_raw()
    datadir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/fastqc/";
    workdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/workdir/collect_fastqc/";
    collect_fastqc(datadir,workdir)
    run_cutadpt_filter()
    fastqc_filter_datadir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/cut_adapter/"
    fastqc_filter_workdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/workdir/fastqc_filter/"
    fastqc_filter_outdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/fastqc_filter/"
    run_fastqc(fastqc_filter_datadir,fastqc_filter_workdir,fastqc_filter_outdir)
    datadir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/fastqc_filter/";
    workdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/workdir/collect_fastqc/";
    collect_fastqc(datadir,workdir)
    datadir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/cut_adapter/"
    workdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/workdir/align/"
    outdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/alignment/"
    run_bwa(datadir,workdir,outdir) 
    datadir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/alignment/" 
    workdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/workdir/mutect/"
    outdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/mutect/"
    run_mutect(datadir,workdir,outdir)
    datadir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/alignment/"
    workdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/workdir/mpileup/"
    outdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/mpileup/"
    run_samtools(datadir,workdir,outdir)

    datadir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/alignment/"
    workdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/workdir/purecn/"
    outdir = "/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/purecn/coverage/"    

    extract_coverage(datadir,workdir,outdir)
    
    datadir="/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/purecn/coverage/"
    datadir2="/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/mutect/"
    workdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/workdir/purecn/purecn/"
    outdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/purity/purecn/tumor_only/"
    run_purecn(datadir,datadir2,workdir,outdir)
    
