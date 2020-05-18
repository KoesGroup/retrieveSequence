from optparse import OptionParser
import glob

print("\n#=========================#\n|       tijs bliek        |\n| University of Amsterdam |\n#=========================#\n")
# get the command line arguments
parser = OptionParser(description="Script that retrieves a fragment (containing a genomic feature like a gene) of a scaffold or contig from a multifasta file",
                    usage="usage: python %prog -g genomename -s scaffold name [options] -o output filename",
                    version="%prog 1.0\nTijs Bliek\nUvA Amsterdam\nMail: Bliek@uva.nl")
parser.add_option('-g', '--genome', 
                    type=str,
                    default="None",
                    metavar="",
                    help = "Name of of the multi fasta file containing the genome scaffolds and/or contigs")
parser.add_option('-s', '--scaffold_name', 
                    type=str,
                    default="None",
                    metavar="",
                    help = "Name of of the scaffold containing your fragment of interest.")
parser.add_option('-f', '--start', 
                    type=int,
                    default=-1,
                    metavar="",
                    help = "starting position of the requested fragment. If skipped beginning of the scaffold will be used")
parser.add_option('-t', '--stop', 
                    type=int,
                    default=-1,
                    metavar="",
                    help = "end position of the requested fragment. If skipped end of the scaffold will be used")
parser.add_option('-d', '--orientation', 
                    type=str,
                    default="fw",
                    metavar="",
                    help = "orientation of the fragment fw (=defauld) or rv.")
parser.add_option('-o', '--output_filename', 
                    type=str,
                    default="None",
                    metavar="",
                    help = "name of the output file.")

(options, args) = parser.parse_args()

g = options.genome
s = options.scaffold_name
f = options.start
t = options.stop
d = options.orientation
o = options.output_filename


def revComp(seq):
    seq = seq.upper()
    revSeq = ""
    base = ""
    for step in range(len(seq)-1,-1,-1):
        if seq[step] == 'A': base = 'T'
        if seq[step] == 'C': base = 'G'
        if seq[step] == 'G': base = 'C'
        if seq[step] == 'T': base = 'A'
        if seq[step] == 'N': base = 'N'
        revSeq = revSeq + base
    return(revSeq)


if s == "None":
    parser.error("name of target contig or scaffold is missing.\nRun \"python Getfasta.py -h\" for more info.")
elif g == "None":
    parser.error("name of genome or other multi fasta file is missing.\nRun \"python Getfasta.py -h\" for more info.")

if len(glob.glob(g)) == 0:
    parser.error("genome file (" + g + ") is not found\nRun \"python Getfasta.py -h\" for more info.")
    g = "None"


if g != "None" and s != "None":
    fas = open(g)
    go = False
    seq = "None"
    afbreken = False
    for line in fas:
        if line.startswith(">"):
            if s in line:
                name = line
                seq=""
                go = True
            else:
                if go:break
                go = False
        elif go:
            seq += line.rstrip()
    if seq == "None":
        print("scaffold name " + s + " is not found in the genome file (" + g + ")." )
    else:
        if f == -1: f = 0
        if t == -1: t = len(seq)
        if t < f: t,f = f,t
        if f > len(seq):
            parser.error("starting position (" + str(f) + ") is bigger then the length (" + str(len(seq)) + ") of the scaffold/contig.\nRun \"python Getfasta.py -h\" for more info.")
            afbreken = True
        if afbreken == False:
            if t > len(seq):
                t = len(seq)
            seq = seq[f:t]
            d = str(d.upper())
            if d == "RV" or d  == "REV" or d == "REVERSE":
                seq = revComp(seq)
                #print("sequence requested in reverse orientation.")
            fas.close()
            if o == "None":
                o = s + "(" + str(f) + "-" + str(t) + ").fa"
            outFile = open(o, "w")
            outFile.write(name)
            outFile.write(seq + "\n")
            outFile.close()
