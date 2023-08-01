#!/usr/bin/env python3
import click
import os
import sys
import subprocess
import re
import tempfile
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def run_cmd(cmd):
    p = subprocess.Popen(
        cmd,
        bufsize=-1,
        shell=True,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        executable="/bin/bash",
    )
    output = p.communicate()[0]
    return output


def is_tool(name):
    """Check whether `name` is on PATH."""
    from distutils.spawn import find_executable
    return find_executable(name) is not None


def process(infile):
    df = pd.read_table(infile, sep="\t", names=("ID", "dir", "count"))
    df = (
        df.pivot(index="ID", columns="dir", values="count")
        .reset_index()
        .rename_axis(None, axis=1)
    )
    df.columns.ID = None
    df.reset_index().fillna(0)
    return df


@click.group()
def main():
    pass


@main.command(help="Locate putative inverted regions")
@click.option(
    "-f",
    "--fasta",
    help="Input genome sequence file in fasta format",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-t",
    "--tab",
    help="Output table with inverted repeats coordinates",
    required=True,
    type=click.Path(),
)
@click.option(
    "-e",
    "--einv",
    help="Einverted parameters, if unspecified run with PhaseFinder default pipeline",
    required=False,
    type=str,
)
@click.option(
    "-m",
    "--mismatch",
    help="Max number of mismatches allowed between IR pairs, used with -einv (default:3)",
    type=int,
    required=False,
    default=3,
)
@click.option(
    "-r",
    "--IRsize",
    help="Max size of the inverted repeats, used with -einv (default:50)",
    type=int,
    required=False,
    default=50,
)
@click.option(
    "-g",
    "--gcRatio",
    help="The minimum and maximum value of GC ratio",
    nargs=2,
    type=float,
    required=False,
    metavar="MIN MAX",
)
@click.option(
    "-p", "--polymer", help="Remove homopolymer inverted repeats", is_flag=True
)
def locate(fasta, tab, einv, mismatch, irsize, gcratio, polymer):
    # step 1: identify the IR position in the genome sequence
    reffile = fasta
    invtab = tab
    if gcratio is not None:
        minGC = gcratio[0]
        maxGC = gcratio[1]
    einvertedParam = einv
    homopolymer = polymer
    maxmis = mismatch
    maxIR = irsize

    f = tempfile.NamedTemporaryFile(mode="w+b", delete=False)
    tmpout = f.name
    f.close()

    if einvertedParam is None:
        # if the einverted parameter is unspecified
        cmd = """
        einverted -maxrepeat 750  -gap 100 -threshold 51 -match 5 -mismatch -9 -outfile {out}.51.outfile -outseq {out}.51.outseq -sequence {ref}
        einverted -maxrepeat 750  -gap 100 -threshold 75 -match 5 -mismatch -15 -outfile {out}.75.outfile -outseq {out}.75.outseq -sequence {ref}

        awk 'BEGIN{{OFS="\\t";ORS="";pass=0}}{{
            if(NR%5==2){{
                split($4,a,"/");
                if(a[2] <= 45 && (a[1]==a[2]) || (a[1]+1==a[2] && a[2] >=13) || (a[1]+2==a[2] && a[2] >=19)){{pass=1}}else{{pass=0}}
                sub(":","",$1);
                if(pass){{print $1"\\t"}}
            }}else if(NR%5==3 && pass ){{print $1-1,$3"\\t"}} else if(NR%5==0 && pass ){{print $3-1,$1"\\n"}}
        }}' {out}.51.outfile | awk '$4-$3>30 ' >{out}.pos.51.tab

        awk 'BEGIN{{OFS="\\t";ORS="";pass=0}}{{
            if(NR%5==2){{
                split($4,a,"/");
                if(a[2] <= 45 && (a[1]==a[2]) || (a[1]+1==a[2] && a[2] >=13) || (a[1]+2==a[2] && a[2] >=19)){{pass=1}}else{{pass=0}}
                sub(":","",$1);
                if(pass){{print $1"\\t"}}
            }}else if(NR%5==3 && pass ){{print $1-1,$3"\\t"}} else if(NR%5==0 && pass ){{print $3-1,$1"\\n"}}
        }}' {out}.75.outfile | awk '$4-$3>30 ' >{out}.pos.75.tab


        awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$5,$0}}' {out}.pos.51.tab |sortBed  > {out}.a.bed
        awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$5,$0}}' {out}.pos.75.tab |sortBed  > {out}.b.bed
        intersectBed   -a {out}.a.bed  -b {out}.b.bed  -v|cat - {out}.b.bed|cut -f 4- > {out}.pos.tab

        rm -rf {out} {out}.a.bed {out}.b.bed {out}.pos.51.tab {out}.pos.75.tab {out}.51.outfile  {out}.51.outseq  {out}.75.outfile  {out}.75.outseq 
        """.format(
            out=tmpout, ref=reffile
        )
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        run_cmd(cmd)
    else:
        # if the einverted parameter is specified
        cmd = """
        einverted  {einvertedParam} -outfile {out}.outfile -outseq {out}.outseq -sequence {ref}
        awk 'BEGIN{{OFS="\\t";ORS="";pass=0}}{{
            if(NR%5==2){{
                split($4,a,"/");
                if(a[2] <= {maxIR} && (a[2] - a[1] <= {maxmis} )){{pass=1}}else{{pass=0}}
                sub(":","",$1);
                if(pass){{print $1"\\t"}}
            }}else if(NR%5==3 && pass ){{print $1-1,$3"\\t"}} else if(NR%5==0 && pass ){{print $3-1,$1"\\n"}}
        }}' {out}.outfile  >{out}.pos.tab

        rm -rf {out} {out}.outfile  {out}.outseq
        """.format(
            out=tmpout,
            ref=reffile,
            maxIR=maxIR,
            maxmis=maxmis,
            einvertedParam=einvertedParam,
        )
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        run_cmd(cmd)

    seq_dict = SeqIO.to_dict(SeqIO.parse(reffile, "fasta"))
    lines = [x.rstrip().split("\t") for x in open(tmpout + ".pos.tab")]
    with open(invtab, "w+") as outfile:
        for each_line in lines:
            accept = 1
            each_seq = seq_dict[each_line[0]]
            posA = int(each_line[1])
            posB = int(each_line[2])
            posC = int(each_line[3])
            posD = int(each_line[4])

            left_seq = each_seq[posA:posB]
            right_seq = each_seq[posC:posD]
            mid_seq = each_seq[posB:posC]

            Lgc = gc_fraction(left_seq.seq) * 100
            Rgc = gc_fraction(right_seq.seq) * 100

            if (
                homopolymer
                and len(re.findall(r"([ACGT])\1{4,}", str(left_seq.seq))) > 0
                and len(re.findall(r"([ACGT])\1{4,}", str(right_seq.seq))) > 0
            ):
                # if homopolymer filter is specified
                accept = 0

            if gcratio is not None and (
                Lgc < minGC or Rgc < minGC or Lgc > maxGC or Rgc > maxGC
            ):
                # if GC ratio filter is specified
                accept = 0

            if accept:
                print(
                    "\t".join(each_line),
                    left_seq.seq,
                    mid_seq.seq,
                    right_seq.seq,
                    sep="\t",
                    file=outfile,
                )
    os.remove(tmpout + ".pos.tab")


@main.command(help="Create inverted fasta file")
@click.option(
    "-f",
    "--fasta",
    help="Input genome sequence file in fasta format",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-t",
    "--tab",
    help="Table with inverted repeat coordinates",
    required=True,
    type=click.Path(),
)
@click.option(
    "-s",
    "--flanksize",
    help="Base pairs of flanking DNA on both sides of the identified inverted repeats",
    type=int,
    required=True,
    default=500,
)
@click.option(
    "-i",
    "--inv",
    help="Output path of the inverted fasta file",
    required=True,
    type=click.Path(),
)
def create(fasta, tab, flanksize, inv):
    # step 2: read the genome with IR and the IR position info to create sequences with putative invertible region inverted

    outseq = list()
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    lines = [x.rstrip().split("\t") for x in open(tab)]
    with open(inv + ".info.tab", "w") as f:
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
            Fversion = left_seq + midfwd_seq + right_seq
            Rversion = left_seq + midrev_seq + right_seq
            name = each_line[0] + ":" + "-".join(each_line[1:5])

            outputpos = list(map(int, each_line[1:5]))
            outputpos.append(right_pos2)
            outputpos = list(map(lambda x: x - left_pos1, outputpos))

            f.write(name + "\t" + "\t".join(map(str, outputpos)) + "\n")

            Fversion.id = name + "_F"
            Rversion.id = name + "_R"
            Fversion.description = ""
            Rversion.description = ""
            outseq.append(Fversion)
            outseq.append(Rversion)

    # write the inverted sequence and index with bowtie-build
    with open(inv, "w") as fasta_out:
        SeqIO.write(outseq, fasta_out, "fasta")
    cmd = """ bowtie-build {genome} {genome} """.format(genome=inv)
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    run_cmd(cmd)


@main.command(help="Align reads to inverted fasta file")
@click.option(
    "-i",
    "--inv",
    help="Input path of the inverted fasta file",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-1",
    "--fastq1",
    help="First pair in fastq",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-2",
    "--fastq2",
    help="Second pair in fastq",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-p", "--threads", help="Number of threads", type=int, default=1, required=False
)
@click.option("-o", "--output", help="Output prefix", required=True, type=str)
def ratio(inv, fastq1, fastq2, threads, output):
    # step 3: align reads to the inverted sequence and identify reads supporting either R or F orientations

    oversize = 10  # require 10 base pairs spanning the invertible region and surrounding genome

    cmd = """ bowtie -p {core}  -a --best --strata {genome} -1 {fq1} -2 {fq2} -S|\
    samtools view -@ {core} -F 4 -h  |sam2bed -d|sortBed |cut -f 1-4,7 > {output}.bed """.format(
        genome=inv, fq1=fastq1, fq2=fastq2, output=output, core=threads
    )
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    run_cmd(cmd)

    cmd = """
    awk -v out={out} 'BEGIN{{OFS="\\t"}}{{print $1"_F",$2,$3,$4,$5"\\n"$1"_R",$2,$3,$4,$5 > out".bed";print $1"_F", $6 "\\n" $1"_R",$6 > out".info" }}'  {out}
    awk '{{print $1"\\t"$2"\\t"$3"\\t1\\n"$1"\\t"$4"\\t"$5"\\t1\\n"$1"\\t"$2"\\t"$5"\\t-1"}}' {out}.bed |slopBed -b {oversize} -g {out}.info |\
    sortBed|intersectBed -c -f 1 -a - -b {output}.bed |awk '{{a[$1]+=$4*$5}}END{{for(i in a){{print i"\\t"a[i]}}}}'|sort -k1,1|sed 's/_F\\t/\\tF\\t/;s/_R\\t/\\tR\\t/' >{output}.span.count""".format(
        out=inv + ".info.tab", output=output, oversize=oversize
    )
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    run_cmd(cmd)

    cmd = """
    cat {output}.bed|awk 'BEGIN{{OFS="\\t"}}{{print $4,$5,$1}}'|sed 's/_\(.\)$/\\t\\1/g'|awk '{{if(and(64,$2)){{P=1}}else{{P=2}};print $1"\\t"P"\\t"$3"\\t"$4}}' > {output}.tab
    cut -f 1-3 {output}.tab|sort|uniq -u|fgrep -f - {output}.tab|cut -f 3-4|sort|uniq -c|awk '{{print $2"\\t"$3"\\t"$1}}' >{output}.pe.count """.format(
        output=output
    )
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    run_cmd(cmd)

    dfspan = process("{output}.span.count".format(output=output))
    dfpe = process("{output}.pe.count".format(output=output))
    if "F" not in list(dfpe):
        dfpe["F"] = 0
    if "R" not in list(dfpe):
        dfpe["R"] = 0
    dfmerge = pd.merge(dfpe, dfspan, on="ID", how="outer").fillna(0)
    df = pd.DataFrame(
        {
            "ID": dfmerge["ID"],
            "Pe_F": dfmerge["F_x"].astype(int),
            "Pe_R": dfmerge["R_x"].astype(int),
            "Pe_ratio": (dfmerge["R_x"] / (dfmerge["R_x"] + dfmerge["F_x"]))
            .astype(float)
            .round(2),
            "Span_F": dfmerge["F_y"].astype(int),
            "Span_R": dfmerge["R_y"].astype(int),
            "Span_ratio": (dfmerge["R_y"] / (dfmerge["R_y"] + dfmerge["F_y"]))
            .astype(float)
            .round(2),
        }
    ).fillna("NA")
    df.to_csv(output + ".ratio.txt", sep="\t", index=False)

    cmd = """rm {output}.bed {output}.tab {output}.span.count {output}.pe.count {out}.bed {out}.info""".format(
        output=output, out=inv + ".info.tab"
    )
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    run_cmd(cmd)


if __name__ == "__main__":
    for i in ["bowtie", "samtools", "sam2bed", "einverted"]:
        if not is_tool(i):
            print("Tool {i} is not installed".format(i=i))
            sys.exit(0)

    main()
