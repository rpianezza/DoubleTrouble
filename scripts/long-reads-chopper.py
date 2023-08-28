import argparse
import gzip

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Chops long reads into short reads.")
parser.add_argument("fastq", help="Path to fastq.gz file")
parser.add_argument("length", help="Length of the chooped reads", type=int)
parser.add_argument("output", help="Path to the output file")
args = parser.parse_args()

'''

'''
def chop(long_read, rl):
    bases = len(long_read)
    reads = []
    for i in range(0, bases, rl):
        read = long_read[i:(i+rl)]
        reads.append(read)
    return reads

with gzip.open(args.fastq, 'rt') as fastq:
        with gzip.open(args.output, 'wt') as output:
            to_write = "no"
            quality = ""
            for line in fastq:
                if line.startswith("@"):
                    name = line.strip().split(' ')[0]
                elif line.startswith("+"):
                    quality = "yes"
                elif quality == "yes":
                    to_write = "yes"
                else:
                    chopped = chop(line[0:-2], args.length)
                if to_write == "yes":
                    for i in range(0, len(chopped)):
                        print(len(chopped))
                        read_number = i+1
                        cigar = "I"*len(chopped[i])
                        output.write(name + "_" + str(read_number) + "\n")
                        output.write(chopped[i] + "\n")
                        output.write("+\n")
                        output.write(cigar + "\n")
                    to_write = "no"
