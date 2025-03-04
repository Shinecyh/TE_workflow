#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 09:24:45 2020

@author: chenyunhao
"""


import argparse
import os
import subprocess
import sys
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
import logging
#from alive_progress import alive_bar    
def ensure_directory_exists(directory):
    """确保目录存在"""
    os.makedirs(directory, exist_ok=True)
    return directory

def check_and_exit_if_not_exists(directory):
    """检查目录是否存在，不存在则退出程序"""
    if not os.path.exists(directory):
        print(f"{directory} does not exist")
        sys.exit(2)
def run_command(command):
    """运行命令并处理异常"""
    try:
        res = subprocess.run(command, shell=True, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
        logging.info(f"Command executed successfully: {command}")
        logging.info(f"Command executed successfully: {res.stderr}")
        return res.returncode == 0 
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running command: {command}")
        logging.error(e)
        print(f"Error running command: {command}")
        print(e)
        sys.exit(2)
def time_out():
    return time.strftime( '%Y-%m-%d %H:%M:%S',time.localtime())

def file_name(file_dir):
    """获取目录下的文件"""
    for root, dirs, files in os.walk(file_dir):
        return files, dirs, root 


def get_fastqgz_filename_in_wd(wd):
    """获取当前目录下的fastq.gz文件"""
    files, dirs, root = file_name(wd)
    fastqgz_filename=[]
    for i in files:
        if(i.endswith('fastq.gz')):
            fastqgz_filename.append(i)
            fq_pattern = 'fastq.gz'
        if(i.endswith('fq.gz')):
            fastqgz_filename.append(i)
            fq_pattern = 'fq.gz'
    return fastqgz_filename, fq_pattern

def get_samp_names(fqfiles, fq_pattern):
    """获取样本名称"""
    names=[]
    for i in fqfiles:
        if(i.endswith('_R1.'+fq_pattern)):
            names.append(i[0:i.rfind("_R1."+fq_pattern)])
            pattern = 0
        if(i.endswith('_1.'+fq_pattern)):
            names.append(i[0:i.rfind("_1."+fq_pattern)]) 
            pattern = 1
        if(i.endswith('.R1.'+fq_pattern)):
            names.append(i[0:i.rfind(".R1."+fq_pattern)])
            pattern = 2
    return names  , pattern

def samp_pairs(i, pattern, fq_pattern):
    """获取样本对"""
    if pattern == 0:
        R1 = i + "_R1."+fq_pattern
        R2 = i + "_R2."+fq_pattern
    elif pattern == 1:
        R1 = i + "_1."+fq_pattern
        R2 = i + "_2."+fq_pattern
    elif pattern == 2:
        R1 = i + ".R1."+fq_pattern
        R2 = i + ".R2."+fq_pattern
    return R1, R2

def fastp_commands(sampnames, threads, pattern, fq_pattern, wd):
    """生成fastp命令"""
    commands = []
    for i in sampnames:
        R1, R2 = samp_pairs(i, pattern, fq_pattern)
        com = ['fastp --detect_adapter_for_pe -f 15 -F 15 -t 15 -T 15',
               '-R', i,
               '--html', i + '-fastp.html',
               '--json', i + '-fastp.json',
               '-i', wd + R1,
               '-I', wd + R2,
               '-o', os.path.join(trim_dir, i + '_1_trimmed.fastq.gz'),
               '-O', os.path.join(trim_dir, i + '_2_trimmed.fastq.gz'),
               '--thread=' + str(threads),
               '--compression=6']
        fastp_command = " ".join(com)
        commands.append(fastp_command)
    return commands

def creat_folder(path):
    """创建目录"""
    if not os.path.exists(path):
        os.mkdir(path)

def samp_pairs_trimmed(i):
    """获取修剪后的样本对"""
    R1 = i + "_1_trimmed.fastq.gz"
    R2 = i + "_2_trimmed.fastq.gz"
    return R1, R2

def map_commands(sampnames, threads):
    """生成映射命令"""
    commands = []
    for i in sampnames:
        R1, R2 = samp_pairs_trimmed(i)
        com = ['hisat2 --rna-strandness RF -p', str(threads),
               '-x', index_loc,
               '-1', os.path.join(trim_dir, R1),
               '-2', os.path.join(trim_dir, R2),
               '| samtools sort -@', str(threads),
               '-o', i + ".bam"]
        map_command = " ".join(com)
        commands.append(map_command)
    return commands
def index_bam_commands(sampnames):
    commands = []
    for i in sampnames:
        com = ['samtools index', os.path.join(map_dir, i + ".bam")]
        index_command = " ".join(com)
        commands.append(index_command)
    return commands
def count_commands(sampnames, gtf_loc, feature_type, threads):
    """生成计数命令"""
    commands = []
    for i in sampnames:
        com = ['featureCounts -s 2 -p -T', str(threads),
               '-t', feature_type,
               '-M -a', gtf_loc,
               '-o', i + ".out",
               os.path.join(map_dir, i + ".bam")]
        count_command = " ".join(com)
        commands.append(count_command)
    return commands

def index_commands(reference, proj_name):
    """生成索引命令"""
    commands = []
    com = ['hisat2-build', reference, proj_name]
    index_command = " ".join(com)
    commands.append(index_command)
    return commands

def rpkm_calculate(counts, lengths):
    """Calculate reads per kilobase transcript per million reads.

    RPKM = (10^9 * C) / (N * L)
    Where:
    C = Number of reads mapped to a gene
    N = Total mapped reads in the experiment
    L = Exon length in base pairs for a gene

    Parameters
    ----------
    counts: array, shape (N_genes,)
        RNAseq (or similar) count data where columns are individual samples
        and rows are genes.
    lengths: array, shape (N_genes,)
        Gene lengths in base pairs in the same order
        as the rows in counts.

    Returns
    -------
    normed : array, shape (N_genes,)
        The RPKM normalized counts matrix.
    """
    N = np.sum(counts, axis=0) # 将每一列相加以得到每个样本的总read数
    L = lengths
    C = counts


    normed = (1e9 * C) / (N* L)

    return(normed)
def count_percent_calculate(counts):
    """Calculate reads per kilobase transcript per million reads.

    RPKM = (10^9 * C) / (N * L)
    Where:
    C = Number of reads mapped to a gene
    N = Total mapped reads in the experiment
    L = Exon length in base pairs for a gene

    Parameters
    ----------
    counts: array, shape (N_genes,)
        RNAseq (or similar) count data where columns are individual samples
        and rows are genes.
    lengths: array, shape (N_genes,)
        Gene lengths in base pairs in the same order
        as the rows in counts.

    Returns
    -------
    normed : array, shape (N_genes,)
        The RPKM normalized counts matrix.
    """
    N = np.sum(counts, axis=0) # 将每一列相加以得到每个样本的总read数
    C = counts
    normed = (100 * C)/N

    return(normed)

def rpkm_df(name,count_dir):
    
    #file_list = [x+'.out' for x in sampnames]
    filename = os.path.join(count_dir,name+'.out')

    with open(filename, 'rt') as f:
        data_table = pd.read_table(f, header =1 ,index_col=0)
    index = data_table.keys() #获取列名
    lengths = np.array(data_table[index[4]])
    
    counts = np.array(data_table[index[5]])
    rpkm = rpkm_calculate(counts,lengths)
    rpkm_df = pd.DataFrame(rpkm,index=data_table.index, columns=[name])
    return rpkm_df
def count_percent_df(name,count_dir):
    
    #file_list = [x+'.out' for x in sampnames]
    filename = os.path.join(count_dir,name+'.out')

    with open(filename, 'rt') as f:
        data_table = pd.read_table(f, header =1 ,index_col=0)
    index = data_table.keys() #获取列名
    lengths = np.array(data_table[index[4]])
    
    counts = np.array(data_table[index[5]])
    count_percent = count_percent_calculate(counts)
    count_percent_df = pd.DataFrame(count_percent,index=data_table.index, columns=[name])
    return count_percent_df

def count_df(name,count_dir):
    
    #file_list = [x+'.out' for x in sampnames]
    filename = os.path.join(count_dir,name+'.out')

    with open(filename, 'rt') as f:
        data_table = pd.read_table(f, header =1 ,index_col=0)
    index = data_table.keys() #获取列名
    lengths = np.array(data_table[index[4]])
    
    counts = np.array(data_table[index[5]])
    count_df = pd.DataFrame(counts,index=data_table.index, columns=[name])
    return count_df
def sort_bam(map_dir, name, fna, threads):
    """
    Sorts BAM files and performs various samtools operations.
    
    Args:
        map (str): Path to the map directory.
        name (str): Base name for output files.
        fna (str): Path to the reference genome in FASTA format.
        threads (int): Number of threads to use.
        wd (str): Working directory.
    """
    # Index the reference genome
    

    # Sort BAM file by name
    #com5 = ['samtools', 'sort', '-n', '-@', str(threads), os.path.join(wd, 'map', f'{name}.bam'), '-o', os.path.join(wd, 'map', f'{name}.1sorted.bam')]
    #com5 = ['samtools', 'index', os.path.join(map_dir, f'{name}.bam')]

    # Fix mates
    com55 = ['samtools', 'fixmate', '-@', str(threads), '-m', os.path.join(map_dir, f'{name}.bam'), os.path.join(map_dir, f'{name}.fix.1sorted.bam')]
    run_command(com55)

    # Sort fixed BAM file
    com56 = ['samtools', 'sort', '--verbosity','1','-@', str(threads), os.path.join(map_dir, f'{name}.fix.1sorted.bam'), '-o', os.path.join(map_dir, f'{name}.sorted.bam')]
    run_command(com56)

    # Generate depth file
    com7 = ['samtools', 'depth', os.path.join(map_dir, f'{name}.sorted.bam'), '>>', os.path.join(map_dir, f'{name}_depth.txt')]
    run_command(com7)

    # Generate statistics file
    com8 = ['samtools', 'flagstat', os.path.join(map_dir, f'{name}.sorted.bam'), '>>', os.path.join(map_dir, f'{name}_stat.txt')]
    run_command(com8)

    # Mark duplicates
    com9 = ['samtools', 'markdup', '-r', '-@', str(threads), os.path.join(map_dir, f'{name}.sorted.bam'), os.path.join(map_dir, f'{name}_nopcr.bam')]
    run_command(com9)

    # Index final BAM file
    com6 = ['samtools', 'index', os.path.join(map_dir, f'{name}_nopcr.bam')]
    run_command(com6)
    #com_rm0 = ['rm -rf',os.path.join(map_dir, f'{name}.sam')]
    #com_rm1 = ['rm -rf',os.path.join(map_dir, f'{name}.1sorted.bam')]
    com_rm2 = ['rm -rf',os.path.join(map_dir, f'{name}.fix.1sorted.bam')]
    com_rm3 = ['rm -rf',os.path.join(map_dir, f'{name}.sorted.bam')]
    run_command(com_rm2)
    run_command(com_rm3)
    #res = subprocess.Popen(" ".join(com4), shell=True)
    #res.wait()
    #res = subprocess.Popen(" ".join(com5), shell=True)
    #res.wait()
    # res = subprocess.Popen(" ".join(com55), shell=True)

    # res.wait()
    # res = subprocess.Popen(" ".join(com56), shell=True)

    # res.wait()
    # res = subprocess.Popen(" ".join(com7), shell=True)

    # res.wait()
    # res = subprocess.Popen(" ".join(com8), shell=True)

    # res.wait()
    # res = subprocess.Popen(" ".join(com9), shell=True)

    # res.wait()
    # res = subprocess.Popen(" ".join(com6), shell=True)

    # res.wait()
    #res = subprocess.Popen(" ".join(com_rm0), shell=True)
    #res.wait()
    #res = subprocess.Popen(" ".join(com_rm1), shell=True)
    #res.wait()
    # res = subprocess.Popen(" ".join(com_rm2), shell=True)
    # res.wait()
    # res = subprocess.Popen(" ".join(com_rm3), shell=True)
    # res.wait()

    # res.wait()

def split_bam(name, map_dir, split_dir,):
    #print(f'Split ' + name + '.bam')
    
    # 定义命令
    comf1 = ['samtools', 'view', '--threads', '8', '-b', '-f', '128', '-F', '16', f'{map_dir}/{name}.bam', '-o', f'{split_dir}/{name}.fwd1.bam']
    comf2 = ['samtools', 'view', '--threads', '8', '-b', '-f', '80', f'{map_dir}/{name}.bam', '-o', f'{split_dir}/{name}.fwd2.bam']
    comf3 = ['samtools', 'index', f'{split_dir}/{name}.fwd1.bam']
    comf4 = ['samtools', 'index', f'{split_dir}/{name}.fwd2.bam']
    comf5 = ['samtools', 'merge', '--threads', '8', '-f','-o', f'{split_dir}/{name}.fwd.bam', f'{split_dir}/{name}.fwd1.bam', f'{split_dir}/{name}.fwd2.bam']
    comf6 = ['samtools', 'index', f'{split_dir}/{name}.fwd.bam']
    com_frm1 = ['rm', '-rf', f'{split_dir}/{name}.fwd1.bam']
    com_frm2 = ['rm', '-rf', f'{split_dir}/{name}.fwd2.bam']
    com_frm3 = ['rm', '-rf', f'{split_dir}/{name}.fwd1.bam.bai']
    com_frm4 = ['rm', '-rf', f'{split_dir}/{name}.fwd2.bam.bai']
    comr1 = ['samtools', 'view', '--threads', '8', '-b', '-f', '144', f'{map_dir}/{name}.bam', '-o', f'{split_dir}/{name}.rev1.bam']
    comr2 = ['samtools', 'view', '--threads', '8', '-b', '-f', '64', '-F', '16', f'{map_dir}/{name}.bam', '-o', f'{split_dir}/{name}.rev2.bam']
    comr3 = ['samtools', 'index', f'{split_dir}/{name}.rev1.bam']
    comr4 = ['samtools', 'index', f'{split_dir}/{name}.rev2.bam']
    comr5 = ['samtools', 'merge', '--threads', '8', '-f','-o', f'{split_dir}/{name}.rev.bam', f'{split_dir}/{name}.rev1.bam', f'{split_dir}/{name}.rev2.bam']
    comr6 = ['samtools', 'index', f'{split_dir}/{name}.rev.bam']
    com_rrm1 = ['rm', '-rf', f'{split_dir}/{name}.rev1.bam']
    com_rrm2 = ['rm', '-rf', f'{split_dir}/{name}.rev2.bam']
    com_rrm3 = ['rm', '-rf', f'{split_dir}/{name}.rev1.bam.bai']
    com_rrm4 = ['rm', '-rf', f'{split_dir}/{name}.rev2.bam.bai']

    # 执行命令
    try:    
        resf1 = subprocess.Popen(comf1)
        resf1.wait()
        resf2 = subprocess.Popen(comf2)
        
        resf2.wait()
        resf3 = subprocess.Popen(comf3)
        resf3.wait()
        resf4 = subprocess.Popen(comf4)
        
        resf4.wait()
        resf5 = subprocess.Popen(comf5)
        resf5.wait()
        resf6 = subprocess.Popen(comf6)
        resf6.wait()

        resrm1 = subprocess.Popen(com_frm1)
        resrm1.wait()
        resrm2 = subprocess.Popen(com_frm2)
        resrm2.wait()
        resrm3 = subprocess.Popen(com_frm3)
        resrm3.wait()
        resrm4 = subprocess.Popen(com_frm4)
        resrm4.wait()

        resr1 = subprocess.Popen(comr1)
        resr2 = subprocess.Popen(comr2)
        resr1.wait()
        resr2.wait()
        resr3 = subprocess.Popen(comr3)
        resr4 = subprocess.Popen(comr4)
        resr3.wait()
        resr4.wait()
        resr5 = subprocess.Popen(comr5)
        resr5.wait()
        resr6 = subprocess.Popen(comr6)
        resr6.wait()

        resrm1 = subprocess.Popen(com_rrm1)
        resrm1.wait()
        resrm2 = subprocess.Popen(com_rrm2)
        resrm2.wait()
        resrm3 = subprocess.Popen(com_rrm3)
        resrm3.wait()
        resrm4 = subprocess.Popen(com_rrm4)
        resrm4.wait()
    except Exception as e:
        print(f"Error executing command: {e}")
def bam2bw(name,split_dir,bigwig_dir,threads):
    com1 = f"bamCoverage --normalizeUsing RPKM -p {threads} --bam {split_dir}/{name}.fwd.bam --outFileName {bigwig_dir}/{name}.fwd.bw"
    com2 = f"bamCoverage --normalizeUsing RPKM -p {threads} --bam {split_dir}/{name}.rev.bam --outFileName {bigwig_dir}/{name}.rev.bw"
    # com2 = [
    #     "bamCoverage",
    #     "--normalizeUsing", "RPKM",
    #     "-p", "16",
    #     "--bam", os.path.join(split_dir, f"{name}.rev.bam"),
    #     "--outFileName", os.path.join(bigwig_dir, f"{name}.rev.bw")
    # ]
    run_command(command=com1)
    run_command(command=com2)
    # res = subprocess.Popen(com1, shell=True)
    # res.wait()
    # res = subprocess.Popen(com2, shell=True)
    # res.wait()
    
def bw2cov(name,bigwig_dir,bed_path,threads):
    com3 = f"computeMatrix reference-point --referencePoint TES --numberOfProcessors {threads} -R {bed_path} --binSize 50 --beforeRegionStartLength 100 --afterRegionStartLength 100 --missingDataAsZero --outFileNameMatrix {bigwig_dir}/{name}.fwd.TES.transterm.matrix.txt -o {bigwig_dir}/{name}.fwd.TES.matrix.mat.gz -S {bigwig_dir}/{name}.fwd.bw"
    #     "",
    #     "", "TES",
    #     "--numberOfProcessors", "16",
    #     "-R", bed_path,
    #     "--binSize", "50",
    #     "--beforeRegionStartLength", "100",
    #     "--afterRegionStartLength", "100",
    #     "--missingDataAsZero",
    #     "--outFileNameMatrix", os.path.join(bigwig_dir, f"{name}.fwd.TES.transterm.matrix.txt"),
    #     "-o", os.path.join(bigwig_dir, f"{name}.fwd.TES.matrix.mat.gz"),
    #     "-S", os.path.join(bigwig_dir, f"{name}.fwd.bw")
    # ]
    com4 = f"computeMatrix reference-point --referencePoint TES --numberOfProcessors {threads} -R {bed_path} --binSize 50 --beforeRegionStartLength 100 --afterRegionStartLength 100 --missingDataAsZero --outFileNameMatrix {bigwig_dir}/{name}.rev.TES.transterm.matrix.txt -o {bigwig_dir}/{name}.rev.TES.matrix.mat.gz -S {bigwig_dir}/{name}.rev.bw"
    # com4 = [
    #     "computeMatrix",
    #     "reference-point",
    #     "--referencePoint", "TES",
    #     "--numberOfProcessors", "16",
    #     "-R", bed_path,
    #     "--binSize", "50",
    #     "--beforeRegionStartLength", "100",
    #     "--afterRegionStartLength", "100",
    #     "--missingDataAsZero",
    #     "--outFileNameMatrix", os.path.join(bigwig_dir, f"{name}.rev.TES.transterm.matrix.txt"),
    #     "-o", os.path.join(bigwig_dir, f"{name}.rev.TES.matrix.mat.gz"),
    #     "-S", os.path.join(bigwig_dir, f"{name}.rev.bw")
    # ]
    com5 = f"pigz -d -p {threads} {bigwig_dir}/{name}.fwd.TES.matrix.mat.gz"
    # com5 = [
    #     "pigz",
    #     "-d",
    #     "-p", {threads},
    #     os.path.join(bigwig_dir, f"{name}.fwd.TES.matrix.mat.gz")
    # ]
    com6 = f"pigz -d -p {threads} {bigwig_dir}/{name}.rev.TES.matrix.mat.gz"
    # com6 = [
    #     "pigz",
    #     "-d",
    #     "-p", {threads},
    #     os.path.join(bigwig_dir, f"{name}.rev.TES.matrix.mat.gz")
    # ]
    # res = subprocess.Popen(com3, shell=True)
    # res.wait()
    # res = subprocess.Popen(com4, shell=True)
    # res.wait()
    # res = subprocess.Popen(com5, shell=True)
    # res.wait()
    # res = subprocess.Popen(com6, shell=True)
    # res.wait()
    run_command(com3)
    run_command(com4)
    run_command(com5)
    run_command(com6)

def calcu_TE(name, bigwig_dir, TE_dir):
    fwd_TES_matrix_file = os.path.join(bigwig_dir, f"{name}.fwd.TES.matrix.mat")
    #fwd_TES_matrix_file = os.path.abspath(f"{bigwig_dir}/{name}.fwd.TES.matrix.mat")
    fwd_coverage_information=pd.read_table(fwd_TES_matrix_file,header=None,index_col=[3,0,1,2,4,5],comment="@",names=["chr","start","end","name","score","strand","plus_U100","plus_U50","plus_D50","plus_D100"])
    fwd_coverage_information[f"{name}_plus_strand_TE%"]=(1-fwd_coverage_information["plus_D50"].div(fwd_coverage_information["plus_U100"]))*100
    rev_TES_matrix_file = os.path.join(bigwig_dir, f"{name}.rev.TES.matrix.mat")
    #rev_TES_matrix_file = os.path.abspath(f"{bigwig_dir}/{name}.rev.TES.matrix.mat")
    rev_coverage_information=pd.read_table(rev_TES_matrix_file,header=None,index_col=[3,0,1,2,4,5],comment="@",names=["chr","start","end","name","score","strand","minus_U100","minus_U50","minus_D50","minus_D100"])
    rev_coverage_information[f"{name}_minus_strand_TE%"]=(1-rev_coverage_information["minus_D50"].div(rev_coverage_information["minus_U100"]))*100
    
    cov_FR_list=[fwd_coverage_information,rev_coverage_information]
    cov_FR = pd.concat(cov_FR_list,axis=1, join='outer')
    cov_FR_simp = cov_FR[[f"{name}_plus_strand_TE%",f"{name}_minus_strand_TE%"]]
    cov_F = cov_FR[f"{name}_plus_strand_TE%"]
    cov_R = cov_FR[f"{name}_minus_strand_TE%"]
    
    cov_FR.to_csv(os.path.join(TE_dir,f"{name}_UD_TE.csv.csv"))
    #cov_FR.to_csv(f"{TE_dir}/{name}_UD_TE.csv")
    cov_FR_simp.to_csv(os.path.join(TE_dir,f"{name}_TE.csv"))
    #cov_FR_simp.to_csv(f"{TE_dir}{name}_TE.csv")
    cov_F.to_csv(os.path.join(TE_dir,f"{name}_F_TE.csv"))
    #cov_F.to_csv(f"{TE_dir}/{name}_F_TE.csv")
    cov_R.to_csv(os.path.join(TE_dir,f"{name}_R_TE.csv"))
    #cov_R.to_csv(f"{TE_dir}/{name}_R_TE.csv")
    
    
def combind_TE_FR(TE_dir,sampnames):
    fwd_TE_list=[]
    rev_TE_list = []
    for name in sampnames:
        fwd_TE = pd.read_csv(os.path.join(TE_dir,f"{name}_F_TE.csv"),header=0,index_col=0)
        #fwd_TE = pd.read_csv(f"{TE_dir}/{name}_F_TE.csv",header=0,index_col=0)
        fwd_TE = fwd_TE[f"{name}_plus_strand_TE%"]
        fwd_TE_list.append(fwd_TE)
        rev_TE = pd.read_csv(os.path.join(TE_dir,f"{name}_R_TE.csv"),header=0,index_col=0)
        #rev_TE = pd.read_csv(f"{TE_dir}/{name}_R_TE.csv",header=0,index_col=0)
        rev_TE = rev_TE[f"{name}_minus_strand_TE%"]
        rev_TE_list.append(rev_TE)        
    FWD = pd.concat(fwd_TE_list,axis=1, join='outer')
    REV = pd.concat(rev_TE_list,axis=1, join='outer')
    FWD.to_csv(os.path.join(TE_dir,f"TE_of_terminaotrs_in_plus_strand_FR.csv"))
    #FWD.to_csv(f"{TE_dir}/TE_of_terminaotrs_in_plus_strand_FR.csv")
    REV.to_csv(os.path.join(TE_dir,f"TE_of_terminaotrs_in_minus_strand_FR.csv"))
    #REV.to_csv(f"{TE_dir}/TE_of_terminaotrs_in_minus_strand_FR.csv")

def combind_TE(TE_dir, sampnames):
    fwd_TE_list=[]
    rev_TE_list = []
    for name in sampnames:    
        fwd_TE = pd.read_csv(os.path.join(TE_dir,f"{name}_F_TE.csv"),header=0,index_col=0)
        #fwd_TE = pd.read_csv(f"{TE_dir}/{name}_F_TE.csv",header=0,index_col=0)
        fwd_TE = fwd_TE[fwd_TE["strand"]=="+"]
        fwd_TE_list.append(fwd_TE[[f"{name}_plus_strand_TE%"]])
        rev_TE = pd.read_csv(os.path.join(TE_dir,f"{name}_R_TE.csv"),header=0,index_col=0)
        #rev_TE = pd.read_csv(f"{TE_dir}/{name}_R_TE.csv",header=0,index_col=0)
        rev_TE = rev_TE[rev_TE["strand"]=="-"]
        rev_TE_list.append(rev_TE[[f"{name}_minus_strand_TE%"]])
    FWD = pd.concat(fwd_TE_list,axis=1, join='outer')
    REV = pd.concat(rev_TE_list,axis=1, join='outer')
    FWD.to_csv(os.path.join(TE_dir,f"TE_of_terminaotrs_in_plus_strand.csv"))
    #FWD.to_csv(f"{TE_dir}/TE_of_terminaotrs_in_plus_strand.csv")
    REV.to_csv(os.path.join(TE_dir,f"TE_of_terminaotrs_in_minus_strand.csv"))
    #REV.to_csv(f"{TE_dir}/TE_of_terminaotrs_in_minus_strand.csv")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This workflow is designed to calculate RPKM and calculate terminator efficiency by RNA-seq data.')
    parser.add_argument('-d', '--work_dir', type=str,help='Folder where store your sequencing data.')
    parser.add_argument('-t', '--threads', type=int, default=8,help='How much threads you need.')
    parser.add_argument('-p', '--project_name', type=str, help='Name for this project.')
    parser.add_argument('-r','--reference', type=str, help='Reference sequence in FASTA format.')    
    parser.add_argument('--gtf', type=str, help='Annotation file in GTF format.')
    parser.add_argument('--bed', type=str, help='Terminaotr location file in BED format')
    parser.add_argument('--type', type=str, default = 'exon', help='Feature type(gene,exon,CDS)', choices=['gene','exon','CDS'])
    
    group = parser.add_argument_group('Debug option',description='When you need to skip certain step, you can use these options. But be attention, every step of this program requires the result of the previous step as input.' )
    group.add_argument('--trim', type=str,default = 'True',choices=['True','False'])
    group.add_argument('--map', type=str,default = 'True',choices=['True','False'])
    #group.add_argument('--sort', type=str,default = 'True',choices=['True','False'])
    group.add_argument('--split', type=str,default = 'True',choices=['True','False'])
    group.add_argument('--count', type=str,default = 'True',choices=['True','False'])
    group.add_argument('--bigwig', type=str,default = 'True',choices=['True','False'])
    group.add_argument('--TE', type=str,default = 'True',choices=['True','False'])
    
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(file = sys.stderr)
        sys.exit(1)
    
    wd = os.path.abspath(args.work_dir)
    
    if not os.path.isdir(wd):
        print("The path specified does not exist")
        sys.exit(2)
    if not(wd.endswith("/")) :
        wd=wd+'/'
    
    logging.basicConfig(
    filename=os.path.join(wd,'app.log'),  # 日志文件名
    level=logging.DEBUG,  # 设置日志级别
    format='%(asctime)s - %(levelname)s - %(message)s'  # 日志格式
    )
    
    group =parser.parse_args()
    trim_bool = group.trim
    map_bool = group.map
    #sort_bool = group.sort
    split_bool = group.split
    count_bool = group.count
    bigwig_bool = group.bigwig
    TE_bool = group.TE
    
    proj_name = args.project_name
    threads = args.threads
    feature_type = args.type
    gtf_loc =  os.path.abspath(args.gtf)
    reference = os.path.abspath(args.reference)
    bed_path = os.path.abspath(args.bed)


    # 使用 os.path.join 来安全地拼接路径，并确保所有路径都以斜杠结尾
    trim_dir = ensure_directory_exists(os.path.join(os.path.abspath(wd), "trimmed"))
    map_dir = ensure_directory_exists(os.path.join(os.path.abspath(wd), "hisat2_map_RF"))
    count_dir = ensure_directory_exists(os.path.join(os.path.abspath(wd), "counts_hisat2_RF"))
    anno_dir = ensure_directory_exists(os.path.join(os.path.abspath(wd), "anno"))
    split_dir = ensure_directory_exists(os.path.join(os.path.abspath(wd), "split"))
    bigwig_dir = ensure_directory_exists(os.path.join(os.path.abspath(wd), "bigwig"))
    TE_dir = ensure_directory_exists(os.path.join(os.path.abspath(wd), "TE"))

    # 确保所有目录存在
    check_and_exit_if_not_exists(trim_dir)
    check_and_exit_if_not_exists(map_dir)
    check_and_exit_if_not_exists(count_dir)
    check_and_exit_if_not_exists(anno_dir)
    check_and_exit_if_not_exists(split_dir)
    check_and_exit_if_not_exists(bigwig_dir)
    check_and_exit_if_not_exists(TE_dir)

    '''
    Building Index  
    '''
    
    if map_bool == 'True'  :
        res = run_command(index_commands(reference, proj_name))
        #print("\nBuilding complete.")
        index_loc = wd+proj_name
        
    os.chdir(wd)
    

    
    
    fqfiles, fq_pattern = get_fastqgz_filename_in_wd(wd)
    #print(fqfiles)
    sampnames, pattern = get_samp_names(fqfiles,fq_pattern)
    #print(sampnames)
    
    """
    Enter Trimprocess
    """
       
    if trim_bool == 'True' :
        
        os.chdir(trim_dir)
        fastp_comds = fastp_commands(sampnames,threads,pattern,fq_pattern,wd)
        check_and_exit_if_not_exists(trim_dir)
        print("Trimming reads. Start at "+time_out())
        with tqdm(total=len(fastp_comds), desc="Trimming reads") as pbar:
            
            for com in fastp_comds:
                res = run_command(com)
                if res:
                    pbar.update(1)
                    #print("\ndone")
                else:
                    print('fastp Error!\n')
                    sys.exit(2)
                
    #trimmed_files,dirs,root = file_name(trim_dir)
    #print(trimmed_files)

    if map_bool == 'True' :
        
        map_coms = map_commands(sampnames, threads)
        check_and_exit_if_not_exists(map_dir)
        os.chdir(map_dir)
        """
        Enter Submap process
        """
        print("Mapping reads to reference. Start at "+time_out())
        with tqdm(total=len(map_coms), desc="Mapping reads to reference") as pbar:
            
            for com in map_coms:
                run_command(com)
                pbar.update(1)
                #print("\ndone")

        index_bam_coms = index_bam_commands(sampnames)
        check_and_exit_if_not_exists(map_dir)
        os.chdir(map_dir)
        """
        Enter index process
        """
        # print("Index BAM filess. Start at "+time_out())
        # with tqdm(total=len(index_bam_coms), desc="Index BAM filess") as pbar:
            
        #     for com in index_bam_coms:
        #         run_command(com)
        #         pbar.update(1)
        #         #print("\ndone")
        
        
    # if sort_bool == 'True' :
    #     print("Sort mapped reads and remove PCR duplicates. Start at "+time_out())
    #     com4 = ['samtools', 'faidx', reference]
    #     run_command(com4)
    #     with tqdm(total=len(sampnames), desc="Sort reads") as pbar:
    #         for name in sampnames:
    #             sort_bam(map_dir=map_dir,name=name,fna=reference,threads=threads)
    #             pbar.update(1)
    #     os.chdir(wd)
        

    if split_bool == 'True' :
        print("Splitting BAM files. Start at "+time_out())
        with tqdm(total=len(sampnames), desc="Splitting BAM files") as pbar:
            
            for name in sampnames:
                split_bam(name, map_dir, split_dir)
                pbar.update(1)
                #print("\ndone")
    """
    Enter Featurecounts process
    """    
    if count_bool == 'True'  :
        check_and_exit_if_not_exists(count_dir)
        os.chdir(count_dir)
        count_coms = count_commands(sampnames,gtf_loc,feature_type,threads)
        with tqdm(total=len(count_coms), desc="Featurecounts process") as pbar:  # 使用`tqdm`
            for com in count_coms:
                run_command(com)
                pbar.update(1)
            
        rpkm_list=[]
        count_list=[]
        count_percent_list=[]
        for name in sampnames:
            rpkm_list.append(rpkm_df(name, count_dir))
            count_list.append(count_df(name, count_dir))
            count_percent_list.append(count_percent_df(name, count_dir))
        rf = pd.concat(rpkm_list,axis=1, join='outer')
        cf = pd.concat(count_list,axis=1, join='outer')
        cpf = pd.concat(count_percent_list,axis=1, join='outer')

        rf.to_csv(os.path.join(count_dir,'RPKMs.csv'))
        cf.to_csv(os.path.join(count_dir,'counts.csv'))
        cpf.to_csv(os.path.join(count_dir,'count_percent.csv'))
    os.chdir(wd)

    """
    Enter TE calculation process
    """   
    if bigwig_bool == 'True' :
        print("Convert splited BAM files into BigWig file and calculate Terminator Efficiency. Start at "+time_out())
        with tqdm(total=len(sampnames), desc="Convert BAM to BigWig") as pbar:  # 使用`tqdm`
            
            for name in sampnames:
                bam2bw(name,split_dir,bigwig_dir,threads=threads)
                #bw2cov(name,bigwig_dir,bed_path,threads=threads)
                pbar.update(1)  # 更新进度条

    if TE_bool == 'True' :
        print("Calculate TE. Start at "+time_out())
        with tqdm(total=len(sampnames), desc="Calculate TE") as pbar:  # 使用`tqdm`
            
            for name in sampnames:
                bw2cov(name,bigwig_dir,bed_path,threads=threads)
                calcu_TE(name,bigwig_dir,TE_dir)
                pbar.update(1)  # 更新进度条
        combind_TE(TE_dir,sampnames)
        combind_TE_FR(TE_dir,sampnames)
        