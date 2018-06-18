#!/usr/bin/env python
import argparse
import os
import sys
import subprocess
from Bio import SeqIO
from collections import defaultdict
import pandas as pd


def run_cmd(cmd):
    cmd = "t1=`date +%s`;"+cmd + \
        """;t2=`date +%s`;tdiff=`echo 'scale=3;('$t2'-'$t1')/60' | bc`;
        echo '##### Total time:  '$tdiff' mins'"""
    p = subprocess.Popen(cmd, bufsize=-1, shell=True, universal_newlines=True,
                         stdout=subprocess.PIPE, executable='/bin/bash')
    output = p.communicate()[0]
    return output


def is_tool(name):
    """Check whether `name` is on PATH."""
    from distutils.spawn import find_executable
    return find_executable(name) is not None


def pipeline_create(args):
    # step 1: read the genome with IR and the IR position info to create sequences with putative invertible region inverted

    flanksize = args.flanksize
    fastafile = args.fastafile
    einvertedtab = args.einvertedtab
    invertedfile = args.invertedfile

    outseq = list()
    seq_dict = SeqIO.to_dict(SeqIO.parse(fastafile, "fasta"))
    lines = [x.rstrip().split("\t") for x in open(einvertedtab)]
    f = open(invertedfile+".info.tab", "w")
    for each_line in lines:
        each_seq = seq_dict[each_line[0]]
        left_pos2 = int(each_line[2])
        right_pos1 = int(each_line[3])
        left_pos1 = max(int(each_line[1]) - flanksize, 0)
        right_pos2 = min(int(each_line[4]) + flanksize, len(each_seq.seq))

        left_seq = each_seq[left_pos1:left_pos2]
        right_seq = each_seq[right_pos1:right_pos2]
        midfwd_seq = each_seq[left_pos2:right_pos1]
        # inverted the sequence between pos2 and pos3
        midrev_seq = midfwd_seq.reverse_complement()
        Fversion = left_seq+midfwd_seq+right_seq
        Rversion = left_seq+midrev_seq+right_seq
        name = each_line[0]+":"+"-".join(each_line[1:5])

        outputpos = list(map(int, each_line[1:5]))
        outputpos.append(right_pos2)
        outputpos = list(map(lambda x: x-left_pos1, outputpos))

        print >>f, name+"\t"+"\t".join(map(str, outputpos))

        Fversion.id = name+"_F"
        Rversion.id = name+"_R"
        Fversion.description = ""
        Rversion.description = ""
        outseq.append(Fversion)
        outseq.append(Rversion)

    # write the inverted sequence and index with bowtie-build
    SeqIO.write(outseq, invertedfile, "fasta")
    f.close()
    cmd = ''' bowtie-build {genome} {genome} '''.format(genome=invertedfile)
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    print run_cmd(cmd)


def process(infile):
    df = pd.read_table(infile, sep="\t", names=("ID", "dir", "count"))
    df = df.pivot(index='ID', columns='dir', values='count').reset_index().rename_axis(None, axis=1)
    df.columns.ID = None
    df.reset_index().fillna(0)
    return df


def pipeline_align(args):
    # step 2: align reads to the inverted sequence and identify reads supporting either R or F orientations

    invertedfile = args.invertedfile
    fq1 = args.fq1
    fq2 = args.fq2
    oversize = 10 # require 10 base pairs spanning the invertible region and surrounding genome
    core = args.core
    output = args.output

    cmd = '''
    bowtie -p {core}  -a --best --strata {genome} -1 {fq1} -2 {fq2} -S|\
    samtools view -@ {core} -F 4 -h  |sam2bed -d|sortBed |cut -f 1-4,7 > {output}.bed '''.format(genome=invertedfile, fq1=fq1, fq2=fq2, output=output, core=core)
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    print run_cmd(cmd)

    cmd = '''
    awk -v out={out} 'BEGIN{{OFS="\\t"}}{{print $1"_F",$2,$3,$4,$5"\\n"$1"_R",$2,$3,$4,$5 > out".bed" ;print $1"_F", $6 "\\n" $1"_R",$6 > out".info" }}'  {out}
    awk '{{print $1"\\t"$2"\\t"$3"\\t1\\n"$1"\\t"$4"\\t"$5"\\t1\\n"$1"\\t"$2"\\t"$5"\\t-1"}}' {out}.bed |slopBed -b {oversize} -g {out}.info |\
    sortBed|intersectBed -c -f 1 -a - -b {output}.bed |awk '{{a[$1]+=$4*$5}}END{{for(i in a){{print i"\\t"a[i]}}}}'|sort -k1,1|sed 's/_F\\t/\\tF\\t/;s/_R\\t/\\tR\\t/' >{output}.span.count'''.format(out=invertedfile+".info.tab", output=output, oversize=oversize)
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    print run_cmd(cmd)

    cmd = '''
    cat {output}.bed|awk 'BEGIN{{OFS="\\t"}}{{print $4,$5,$1}}'|sed 's/_\(.\)$/\\t\\1/g'|awk '{{if(and(64,$2)){{P=1}}else{{P=2}};print $1"\\t"P"\\t"$3"\\t"$4}}' > {output}.tab
    cut -f 1-3 {output}.tab|sort|uniq -u|fgrep -f - {output}.tab|cut -f 3-4|sort|uniq -c|awk '{{print $2"\\t"$3"\\t"$1}}' >{output}.pe.count '''.format(output=output)
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    print run_cmd(cmd)

    dfspan=process("{output}.span.count".format(output=output))
    dfpe=process("{output}.pe.count".format(output=output))
    if 'F' not in list(dfpe): dfpe['F'] = 0
    if 'R' not in list(dfpe): dfpe['R'] = 0
    dfmerge = pd.merge(dfpe,dfspan, on = "ID", how='outer').fillna(0)
    df = pd.DataFrame({"ID":dfmerge["ID"],
        "Pe_F": dfmerge["F_x"].astype(int),
        "Pe_R": dfmerge["R_x"].astype(int),
        "Span_F": dfmerge["F_y"].astype(int),
        "Span_R": dfmerge["R_y"].astype(int)})
    df.to_csv(output+".ratio.txt",sep="\t",index=False)

    cmd = '''rm {output}.bed {output}.tab {output}.span.count {output}.pe.count {out}.bed {out}.info'''.format(
        output=output, out=invertedfile+".info.tab")
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    print run_cmd(cmd)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='Identify putative insertions in reference genomes from bam alignments')
    subparsers = parser.add_subparsers(title="actions")

    parser_create = subparsers.add_parser(
        "create",  help="create inverted fasta file")
    parser_create.add_argument('-f', '--fasta', help='input genome sequence file in fasta format',
                               required=True, dest='fastafile', metavar='')
    parser_create.add_argument('-t', '--tab', help='table with inverted repeats coordinates',
                               required=True, dest='einvertedtab', metavar='')
    parser_create.add_argument('-s', '--flanksize', help='flanked base pair',
                               type=int, required=True, dest='flanksize', metavar='')
    parser_create.add_argument('-i', '--inv', help='output path of the inverted fasta file',
                               required=True, dest='invertedfile', metavar='')
    parser_create.set_defaults(func=pipeline_create)

    parser_align = subparsers.add_parser(
        "ratio",  help="align reads to inverted fasta file")
    parser_align.add_argument('-i', '--inv', help='input path of the inverted fasta file',
                              required=True, dest='invertedfile', metavar='')
    parser_align.add_argument('-1', '--fastq1', help='first pair in fastq',
                              required=True, dest='fq1', metavar='')
    parser_align.add_argument('-2', '--fastq2', help='second pair in fastq',
                              required=True, dest='fq2', metavar='')
    parser_align.add_argument('-p', '--threads', help='number of threads used',
                              type=int, default=1, required=False, dest='core', metavar='')
    parser_align.add_argument('-o', '--output', help='output prefix',
                              required=True, dest='output', metavar='')
    parser_align.set_defaults(func=pipeline_align)

    for i in ["bowtie", "samtools", "sam2bed", "bc"]:
        if not is_tool(i):
            print "tool {i} is not installed".format(i=i)
            sys.exit(0)

    args = parser.parse_args()
    args.func(args)
