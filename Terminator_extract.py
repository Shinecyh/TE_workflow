#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 20:21:36 2023

@author: Chen Yunhao
"""
import argparse
import sys
import os
import csv
from textwrap import dedent
from Bio.Seq import Seq
from Bio import SeqIO
import primer3
from alive_progress import alive_bar
def rev_comp(seq):
    """
    返回给定序列的反向互补序列。

    参数:
    seq (str): 输入的核酸序列。

    返回:
    str: 输入序列的反向互补序列。

    示例:
    rev_comp("GATCGATGGGCCTATATAGG")
    'CCCGGTACCGCGATCGCGATTAGCGCGA'
    """
    try:
        # 将输入序列转换为Bio.Seq.Seq对象
        my_seq = Seq(seq)
        # 生成并返回反向互补序列的字符串表示
        return str(my_seq.reverse_complement())
    except Exception as e:
        # 异常处理
        print(f"Error processing sequence: {e}")
        return None

# 测试代码
# if __name__ == "__main__":
#     test_seq = "GATCGATGGGCCTATATAGG"
#     result = rev_comp(test_seq)
#     print(result)

def find_next_gene(i,data):
    """
    寻找下一个基因的ID和方向。

    该函数通过分析存储在'data'中的基因信息，从给定的索引'i'开始寻找。
    它主要通过分割字符串并检查特定关键字来确定基因的ID和方向。

    参数:
    i -- 开始搜索的索引位置。

    返回:
    gene_id -- 找到的基因ID。
    now_gene_sense -- 基因的方向（'sense'或'antisense'）。
    """
    while True:   
        # 检查索引是否越界
        if i >= len(data):
            return " ", " "
        
        # 获取当前行
        current_line = data[i]
        
        # 分割当前行
        current_split = current_line.split()
        if not current_split:
            continue
        
        
        # 当前行以TERM开头，向下索引i+2
        if current_split[0].startswith("TERM"):
            if i + 2 < len(data):
                i += 2
            else:
                break
        
        # 行以locus开头，找到基因ID和方向，返回
        elif current_split[0].startswith(locus):
            gene_id = current_split[0]
            now_gene_sense = current_split[4]
            return gene_id, now_gene_sense
    return " ", " "
def pick_primers(sequence: str, start: int, end: int, sense: str):
    """
    使用primer3工具设计引物对给定DNA序列。

    参数:
    sequence (str): DNA序列。
    start (int): 引物设计的起始位置。
    end (int): 引物设计的结束位置。
    sense (str): DNA链的方向，'+'表示正向，'-'表示反向。

    返回:
    tuple: 包含左引物和右引物的元组。

    引物设计遵循指定的起始和结束位置，以及链的方向。该函数适用于任何符合要求的DNA序列。
    """

    # 类型检查
    if not isinstance(sequence, str) or not isinstance(start, int) or not isinstance(end, int) or not isinstance(sense, str):
        raise TypeError("Invalid argument types")
   
    # 处理边界条件
    if start == end:
        raise ValueError("Start and end positions must be different")

    # 根据起始和结束位置，构建primer3所需的序列参数
    seq_args = {
        'SEQUENCE_ID': 'Terminaotor',
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_FORCE_LEFT_START': min(int(start), int(end))-16,
        'SEQUENCE_FORCE_RIGHT_START': max(int(start), int(end))+14,
    }

    # 设置全局参数，用于控制引物设计的条件
    global_args = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_OPT_SIZE': 30,
        'PRIMER_MIN_SIZE': 20,
        'PRIMER_MAX_SIZE': 36,
        'PRIMER_OPT_TM': 62.0,
        'PRIMER_MIN_TM': 58.0,
        'PRIMER_MAX_TM': 68.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 3,
        'PRIMER_MAX_SELF_END_TH': 35.00,
        'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45,
        'PRIMER_PAIR_MAX_COMPL_END_TH': 35,
        'PRIMER_NUM_RETURN': 1,
        'PRIMER_PRODUCT_SIZE_RANGE': [40, 100],
        'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
        'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 1,
        'PRIMER_PICK_ANYWAY': 1
    }

    # 尝试设计引物
    try:
        result = primer3.design_primers(seq_args, global_args)
        # 根据链的方向返回相应的引物序列
        if sense == '-':
            return result["PRIMER_RIGHT_0_SEQUENCE"], result["PRIMER_LEFT_0_SEQUENCE"], result["PRIMER_RIGHT_0_TM"], result["PRIMER_LEFT_0_TM"]
        return result["PRIMER_LEFT_0_SEQUENCE"], result["PRIMER_RIGHT_0_SEQUENCE"], result["PRIMER_LEFT_0_TM"], result["PRIMER_RIGHT_0_TM"]
    except primer3.PrimerDesignError as e:
        raise RuntimeError("Failed to design primers") from e
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise RuntimeError("An unexpected error occurred during primer design") from e
    finally:
        # 可以在这里做一些清理工作，例如关闭文件等
        pass
    
    
def pick_primers_for_upstream_gene(sequence: str, start: int, end: int, sense: str, up_length:int,ter_right_primer: str):
    """
    使用primer3工具设计引物对给定DNA序列。

    参数:
    sequence (str): DNA序列。
    start (int): 引物设计的起始位置。
    end (int): 引物设计的结束位置。
    sense (str): DNA链的方向，'+'表示正向，'-'表示反向。

    返回:
    tuple: 包含左引物和右引物的元组。

    引物设计遵循指定的起始和结束位置，以及链的方向。该函数适用于任何符合要求的DNA序列。
    """

    # 类型检查
    if not isinstance(sequence, str) or not isinstance(start, int) or not isinstance(end, int) or not isinstance(sense, str):
        raise TypeError("Invalid argument types")
   
    # 处理边界条件
    if start == end:
        raise ValueError("Start and end positions must be different")

    # 根据起始和结束位置，构建primer3所需的序列参数
    seq_args_plus = {
        'SEQUENCE_ID': 'Terminaotor',
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_FORCE_LEFT_START': min(int(start), int(end))-16-up_length,
        
    }
    seq_args_minus = {
        'SEQUENCE_ID': 'Terminaotor',
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_FORCE_RIGHT_START': max(int(start), int(end))+14+up_length,
        
    }

    # 设置全局参数，用于控制引物设计的条件
    global_args_plus = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 0,
        'SEQUENCE_PRIMER_REVCOMP':(ter_right_primer),
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_OPT_SIZE': 30,
        'PRIMER_MIN_SIZE': 20,
        'PRIMER_MAX_SIZE': 36,
        'PRIMER_OPT_TM': 62.0,
        'PRIMER_MIN_TM': 58.0,
        'PRIMER_MAX_TM': 68.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 3,
        'PRIMER_MAX_SELF_END_TH': 35.00,
        'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45,
        'PRIMER_PAIR_MAX_COMPL_END_TH': 35,
        'PRIMER_NUM_RETURN': 1,
        'PRIMER_PRODUCT_SIZE_RANGE': [40, up_length+100],
        'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
        'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 1,
        'PRIMER_PICK_ANYWAY': 1
    }
    
    global_args_minus = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 0,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'SEQUENCE_PRIMER':(ter_right_primer),
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_OPT_SIZE': 30,
        'PRIMER_MIN_SIZE': 20,
        'PRIMER_MAX_SIZE': 36,
        'PRIMER_OPT_TM': 62.0,
        'PRIMER_MIN_TM': 58.0,
        'PRIMER_MAX_TM': 68.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 3,
        'PRIMER_MAX_SELF_END_TH': 35.00,
        'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45,
        'PRIMER_PAIR_MAX_COMPL_END_TH': 35,
        'PRIMER_NUM_RETURN': 1,
        'PRIMER_PRODUCT_SIZE_RANGE': [40, up_length+100],
        'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
        'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 1,
        'PRIMER_PICK_ANYWAY': 1
    }
    try:
        # 根据链的方向返回相应的引物序列
        if sense == '-':
            result = primer3.design_primers(seq_args_minus, global_args_minus)
            return result["PRIMER_RIGHT_0_SEQUENCE"], result["PRIMER_RIGHT_0_TM"]
        result = primer3.design_primers(seq_args_plus, global_args_plus)
        return result["PRIMER_LEFT_0_SEQUENCE"], result["PRIMER_LEFT_0_TM"]
    except primer3.PrimerDesignError as e:
        raise RuntimeError("Failed to design primers") from e
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise RuntimeError("An unexpected error occurred during primer design") from e
    finally:
        # 可以在这里做一些清理工作，例如关闭文件等
        pass
def get_up_gene_seq(sequence: str, start: int, end: int, sense: str, up_length: int):
    """
    Get the upstream gene sequence from the terminator region.
    
    This function extracts the upstream gene sequence based on the provided gene sequence, start and end positions, 
    sense strand marker, and the required upstream sequence length.
    
    Parameters:
    - sequence (str): The full gene sequence.
    - start (int): The start position of the gene.
    - end (int): The end position of the gene.
    - sense (str): The sense strand marker of the gene, '+' for the positive strand, '-' for the negative strand.
    - up_length (int): The length of the upstream sequence to be extracted.
    
    Returns:
    - str: The extracted upstream gene sequence.
    """
    # 检查输入参数的合法性
    if not sequence or start < 0 or end < 0 or start > len(sequence) or end > len(sequence):
        raise ValueError("Invalid input parameters")
    
    if sense not in ['+', '-']:
        raise ValueError("Invalid sense strand marker")
    
    # 计算上游序列的起始和结束位置
    if sense == '+':
        # 正链上的上游序列
        start_pos = max(0, start - 16 - up_length)
        end_pos = start - 16
        upstream_gene = sequence[start_pos:end_pos]
    else:
        # 负链上的上游序列
        start_pos = start + 15
        end_pos = min(len(sequence), start + 15 + up_length)
        upstream_gene = rev_comp(sequence[start_pos:end_pos])
    
    
    return upstream_gene
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Terminator_extract.py',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=dedent('''\
                                    ----------------
                                    Simple Usage:
                                    %(prog)s -i ./Transterm_out -o ./Terminator_info.csv
                                        Extract terminator info predictied by TransTermHP.
                                    ----------------
                                    FORMAT OF THE TRANSTERM OUTPUT

                                    The organism's genes are listed sorted by their end coordinate and terminators
                                    are output between them. A terminator entry looks like this:
                                    
                                        TERM 19  15310 - 15327  -      F     99      -12.7 -4.0 |bidir
                                        (name)   (start - end)  (sense)(loc) (conf) (hp) (tail) (notes)
                                    
                                    where 'conf' is the overall confidence score, 'hp' is the hairpin score, and
                                    'tail' is the tail score. 'Conf' (which ranges from 0 to 100) is what you
                                    probably want to use to assess the quality of a terminator. Higher is better.
                                    The confidence, hp score, and tail scores are described in the paper cited
                                    above.  'Loc' gives type of region the terminator is in:
                                    
                                        'G' = in the interior of a gene (at least 50bp from an end),
                                        'F' = between two +strand genes,
                                        'R' = between two -strand genes,
                                        'T' = between the ends of a +strand gene and a -strand gene,
                                        'H' = between the starts of a +strand gene and a -strand gene,
                                        'N' = none of the above (for the start and end of the DNA)
                                    
                                    Because of how overlapping genes are handled, these designations are not
                                    exclusive. 'G', 'F', or 'R' can also be given in lowercase, indicating that
                                    the terminator is on the opposite strand as the region.  Unless the
                                    --all-context option is given, only candidate terminators that appear to be in
                                    an appropriate genome context (e.g. T, F, R) are output. 
                                    
                                    Following the TERM line is the sequence of the hairpin and the 5' and 3'
                                    tails, always written 5' to 3'.
                                        TTCCTAAAGGTTCCA  GCG CAAAA TGC  CATAAGCACCACATT
                                        (left context)     (hairpin)    (right contenxt)
                                     '''))
    parser.add_argument('-i', '--input_file', type=str,required=True,help='Path to TransTermHP output file')
    parser.add_argument('-o', '--output_csv', type=str,default='Terminator_info.csv',help='output file name. Default = Terminator_info.csv')
    parser.add_argument('-v','--version',action = 'version',version = '%(prog)s V2.2 20240916')
    parser.add_argument('-l','--locus',type=str,default='ZMO',help='Header of locus_tag, default = ZMO')
    parser.add_argument('-c','--chr',type=str,required=True,help='Chromesome accession')
    parser.add_argument('-t','--site_list',type=str,help='A file which contains Term-seq detacted TTS sites in format as Name\tStrand\tSite')
    parser.add_argument('-s','--sequence',type=str,help='Template sequence for terminator cloning, usually the genome sequence of the organism')
    parser.add_argument('--upstream_length',type=int,default=100,help='The length of upstream sequence to be amplified, default = 100')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(file = sys.stderr)
        sys.exit(1)
    root = os.getcwd()
    file = args.input_file
    output_path = args.output_csv
    locus = args.locus
    chromesome = args.chr
    
    sequence_file = args.sequence
    upstream_length = args.upstream_length
    genome=SeqIO.read(os.path.abspath(sequence_file), "fasta")

    
    # root = "G:/实验数据/终止子新/文章数据/终止子预测"
    # file = "G:/实验数据/终止子新/文章数据/终止子预测/out_genome_rr.tt"
    # output_path = "G:/实验数据/终止子新/文章数据/终止子预测/sss.csv"
    # locus = "ZMO"
    # site_list = "G:实验数据/终止子新/文章数据/终止子预测/Term-seq得到的终止子位置/CP023715_Term-seq.txt"
    # genome = SeqIO.read(os.path.abspath("G:/实验数据/终止子新/文章数据/终止子预测/Test/ZMNP.fasta"), "fasta")
    # upstream_length = 100
    # chromesome = "CP023715"
    
    class Terminator:
        def __init__(self, name, location, start,end, sense, loc, conf, hp, tail, notes, A_tract, Left_stem, loop, Right_stem, U_tract,Left_primer,Right_primer, Left_primer_TM,Right_primer_TM,up_gene_primer=" ",up_gene_TM=" ",up_gene=" ", down_gene=" ",up_gene_seq_50=" ", up_gene_seq_length=""):
            
            self.name = name
            self.location = location
            self.sense = sense
            self.loc = loc
            self.conf = conf
            self.hp = hp
            self.tail = tail
            self.notes = notes
            self.A_tract = A_tract
            self.Left_stem = Left_stem
            self.loop = loop
            self.Right_stem = Right_stem
            self.U_tract = U_tract
            self.up_gene = up_gene
            self.down_gene = down_gene
            self.Left_primer = Left_primer
            self.Right_primer = Right_primer
            self.Left_primer_TM = Left_primer_TM
            self.Right_primer_TM = Right_primer_TM
            self.up_gene_primer = up_gene_primer
            self.up_gene_TM = up_gene_TM
            self.up_gene_seq_50 = up_gene_seq_50
            self.up_gene_seq_length = up_gene_seq_length
             # 尝试将 start 和 end 转换为整数，并处理异常
            try:
                self.start = int(start.strip())
                self.end = int(end.strip())
            except ValueError:
                raise ValueError("Start and end must be convertible to integers.")
                
            self.Term_list = []
        #     drop_seq = "".join([i.A_tract,i.Left_stem, i.loop, i.Right_stem])
        #     full_seq_length = "".join([i.up_gene_seq_length,i.terminator_seq()])
        #     full_seq_50= "".join([i.up_gene_seq_50,i.terminator_seq()])
        #     drop_seq_length = "".join([i.up_gene_seq_length,drop_seq])
        #     drop_seq_50 = "".join([i.up_gene_seq_50,drop_seq])
        #     hair_pin = "".join([i.Left_stem, i.loop, i.Right_stem])
        # def get_drop_seq(self):
        #     if self.sense == "+":
        #         drop_seq = self.A_tract + self.Left_stem + self.loop + self.Right_stem
        #     else:
        #         drop_seq = self.A_tract + self.Right_stem + self.loop + self.Left_stem
        #     return drop_seq
        # def get_full_seq_length(self):
        #     if self.sense == "+":
        #         full_seq_length = self.up_gene_seq_length + self.terminator_seq()
        #     else:
        #         full_seq_length = self.terminator_seq() + self.up_gene_seq_length
        #     return full_seq_length
        # def get_full_seq_50(self):
        #     if self.sense == "+":
        #         full_seq_50 = self.up_gene_seq_50 + self.terminator_seq()
        #     else:
        #         full_seq_50 = self.terminator_seq() + self.up_gene_seq_50
        #     return full_seq_50
        # def get_drop_seq_length(self):
        #     if self.sense == "+":
        #         drop_seq_length = self.up_gene_seq_length + self.get_drop_seq()
        #     else:
        #         drop_seq_length = self.get_drop_seq() + self.up_gene_seq_length
        #     return drop_seq_length
        # def get_drop_seq_50(self):
        #     if self.sense == "+":
        #         drop_seq_50 = self.up_gene_seq_50 + self.get_drop_seq()
        #     else:
        #         drop_seq_50 = self.get_drop_seq() + self.up_gene_seq_50
        #     return drop_seq_50
        # def get_hairpin(self):
        #     hair_pin = self.Left_stem + self.loop + self.Right_stem
        #     return hair_pin
        
        

        def terminator_seq(self):
            A_tract = self.A_tract
            Left_stem = self.Left_stem
            loop = self.loop
            Right_stem = self.Right_stem
            U_tract = self.U_tract
            seq = A_tract+Left_stem+loop+Right_stem+U_tract
            return seq
        def region_terms(self):
            return "/".join(self.Term_list)
    # 读取文件
    try:
        with open(file, encoding="utf-8") as f:
            data = f.readlines()
    except FileNotFoundError:
        print(f"文件 {file} 未找到")
        data = []
    # 计算终止子个数
    term_num = 0
    for line in data:
        data_split = line.split()  # 默认按空白字符分割
        if data_split and data_split[0].startswith("TERM"):
            term_num += 1

    print(f"Total terminaotr number: {term_num}")
    
    gene_id = ""
    term_list={}
    with alive_bar(len(range(term_num))) as bar:
        i=0
        while i < (len(data)-2):
            
            now = data[i]
            down = data[i+1]
            now_split = now.split()
            down_split = down.split()
            if len(now_split) > 0 and len(down_split) > 0:
                # 如果当前行和下一行都以locus开头，则跳过当前行
                if now_split[0].startswith(locus):
                    
                    gene_id=now_split[0]
                    now_gene_sense = now_split[4]
                    
                    i+=1
                    continue
                elif now_split[0].startswith("TERM") :
                    
                    term = data[i+1]
                    term_split = term.split()
                    name = f"TransTerm_{now_split[0]}_{now_split[1]}"
                    location = now_split[2]+now_split[3]+now_split[4]
                    start = now_split[2]
                    end = now_split[4]
                    sense = now_split[5]
                    loc = now_split[6]
                    conf = now_split[7]
                    hp = now_split[8]
                    tail = now_split[9]
                    notes = " ".join(now_split[10:])
                    A_tract = term_split[0]
                    Left_stem = term_split[1]
                    loop = term_split[2]
                    Right_stem = term_split[3]
                    U_tract = term_split[4]
                    # 调用 pick_primers 函数，并处理可能出现的异常
                    try: 
                        #print(f"Pick_primers:{name}")
                        left_primer, right_primer, left_primer_TM, right_primer_TM= pick_primers(str(genome.seq._data,encoding = 'UTF-8'), int(start), int(end), sense)
                        up_gene_primer,  up_gene_TM = pick_primers_for_upstream_gene(str(genome.seq._data,encoding = 'UTF-8'), int(start), int(end), sense,upstream_length,ter_right_primer=right_primer)
                        #print(f"Pick_primers:{name} done")
                    except Exception as e:
                        raise RuntimeError(f"Failed to pick primers: {name}")
                    next_name, next_sense = find_next_gene(i,data)
                    up_gene_seq_50 = get_up_gene_seq(str(genome.seq._data,encoding = 'UTF-8'), int(start), int(end), sense,50)
                    up_gene_seq_length = get_up_gene_seq(str(genome.seq._data,encoding = 'UTF-8'), int(start), int(end), sense,upstream_length)
                    if sense=='+' and now_gene_sense=="+":
                        term_list[name]=Terminator( name, location, start,end, sense, loc, conf, hp, tail, notes, A_tract, Left_stem, loop, Right_stem, U_tract,left_primer,right_primer,left_primer_TM,right_primer_TM,up_gene_primer=up_gene_primer, up_gene_TM=up_gene_TM, up_gene=gene_id,down_gene=next_name,up_gene_seq_50 = up_gene_seq_50,up_gene_seq_length=up_gene_seq_length)
                        #term_list.append(terminator( name, location, start,end, sense, loc, conf, hp, tail, notes, A_tract, Left_stem, loop, Right_stem, U_tract))
                    elif sense=='-'and next_sense =="-" :
                        term_list[name]=Terminator( name, location, start,end, sense, loc, conf, hp, tail, notes, rev_comp(U_tract), rev_comp(Right_stem), rev_comp(loop),rev_comp(Left_stem), rev_comp(A_tract),left_primer,right_primer,left_primer_TM,right_primer_TM,up_gene_primer=up_gene_primer, up_gene_TM=up_gene_TM,up_gene=next_name,down_gene=gene_id,up_gene_seq_50 = up_gene_seq_50,up_gene_seq_length=up_gene_seq_length)
                        #term_list.append(terminator( name, location, start,end, sense, loc, conf, hp, tail, notes, rev_comp(U_tract), rev_comp(Right_stem), rev_comp(loop),rev_comp(Left_stem), rev_comp(A_tract)))
                    elif sense=='+':
                        term_list[name]=Terminator( name, location, start,end, sense, loc, conf, hp, tail, notes, A_tract, Left_stem, loop, Right_stem, U_tract,left_primer,right_primer,left_primer_TM,right_primer_TM,up_gene_primer=up_gene_primer, up_gene_TM=up_gene_TM,up_gene_seq_50 = up_gene_seq_50,up_gene_seq_length=up_gene_seq_length)
                    elif sense=='-':
                        term_list[name]=Terminator( name, location, start,end, sense, loc, conf, hp, tail, notes, rev_comp(U_tract), rev_comp(Right_stem), rev_comp(loop),rev_comp(Left_stem), rev_comp(A_tract),left_primer,right_primer,left_primer_TM,right_primer_TM,up_gene_primer=up_gene_primer, up_gene_TM=up_gene_TM,up_gene_seq_50 = up_gene_seq_50,up_gene_seq_length=up_gene_seq_length)
                    i+=2
                    bar()
                    continue
                
            i+=1
    if args.site_list != None:
        site_list = args.site_list
        if len(site_list)>0:
            # 读取文件内容
            try:
                with open(site_list, 'r', encoding="utf-8") as f:
                    terms = f.readlines()
            except FileNotFoundError:
                print(f"文件 {site_list} 不存在")
                
            except IOError as e:
                print(f"读取文件时发生错误: {e}")
                

            # 处理每一行数据
            for line in terms:
                term_split = line.split("\t")
                try:
                    site = int(term_split[2].replace("\n", ""))
                except ValueError:
                    print(f"无效的 site 值: {term_split[2]}")
                    continue
                
                name = term_split[0]
                strand = term_split[1]
                
                # 遍历 term_list 并添加匹配项
                for key, terminator in term_list.items():
                    if int(min(terminator.start, terminator.end))-15 <= site <= int(max(terminator.start, terminator.end))+15:
                        terminator.Term_list.append("_".join([name, strand, str(site)]))
                
    
    
    with open(output_path,'w',newline='') as file:
        writer = csv.writer(file)
        header_list = ["name", "sequence","location", "start", "end","sense", "loc", "conf", "hp", "tail", "notes", "A_tract", "Left_stem", "loop", "Right_stem","U_tract", "Left_primer", "Right_primer", "Left_primer_TM", "Right_primer_TM","up_gene_primer", "up_gene_TM", "up_gene", "down_gene","Term-seq"]
        writer.writerow(header_list)
        for i in term_list.values():
            writer.writerow([i.name, i.terminator_seq(),i.location, i.start, i.end, i.sense, i.loc, i.conf, i.hp, i.tail, i.notes, i.A_tract, i.Left_stem, i.loop, i.Right_stem, i.U_tract,i.Left_primer,i.Right_primer,i.Left_primer_TM,i.Right_primer_TM, i.up_gene_primer, i.up_gene_TM, i.up_gene, i.down_gene,"/".join(i.Term_list)])
    
    with open(os.path.join(os.path.dirname(output_path),chromesome+"_Ter.bed"),'w',newline='') as file:
        writer = csv.writer(file,delimiter="\t")
        for i in term_list.values():
            writer.writerow([chromesome, int(min(i.start,i.end))-15-1,int(max(i.start,i.end))+15-1,i.name,"0",i.sense])
    
    with open(os.path.join(os.path.dirname(output_path),chromesome+"_Ter_seqs.csv"),'w',newline="") as file:
        writer = csv.writer(file)
        header_list = ["name", "Full_seq","Drop_seq", f"Full_seq_{upstream_length}", "Full_seq_50",f"Drop_seq_{upstream_length}", "Drop_seq_50", "Hair_pin"]
        writer.writerow(header_list)
        for i in term_list.values():
            drop_seq = "".join([i.A_tract,i.Left_stem, i.loop, i.Right_stem])
            full_seq_length = "".join([i.up_gene_seq_length,i.terminator_seq().strip('-')])
            full_seq_50= "".join([i.up_gene_seq_50,i.terminator_seq().strip('-')])
            drop_seq_length = "".join([i.up_gene_seq_length,drop_seq])
            drop_seq_50 = "".join([i.up_gene_seq_50,drop_seq])
            hair_pin = "".join([i.Left_stem, i.loop, i.Right_stem])
            writer.writerow([i.name, i.terminator_seq().strip('-'),drop_seq,full_seq_length,full_seq_50,drop_seq_length,drop_seq_50,hair_pin])
            
            

    sys.exit(0)
            
            