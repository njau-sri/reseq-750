import os
import sys

def add_index(infile,outfile):
    with open(infile)as f1:
        with open(outfile,"w")as f2:
            n = f1.readlines()
            for i in n:
                f2.write(i.strip("\r\n")+"\t"+str(n.index(i)+1))
                f2.write("\n")

def find_blocks_haploview_result(infile,outfile):
    with open(infile)as f1:
        with open(outfile,"w")as f2:
            for line in f1:
                if line.startswith("BLOCK"):
                    n = line.strip("\n").split()
                    f2.write("\t".join(x for x in n[3:]))
                    f2.write("\n")

if __name__ == "__main__": 
    add_index(sys.argv[1],sys.argv[2])
    find_blocks_haploview_result(sys.argv[3],sys.argv[4])
    
            
