import os
import sys

with open(sys.argv[1])as f1:
    with open(sys.argv[2],"w")as f2:
        for line in f1:
            n = line.strip("\r\n").split("_")
            f2.write(line.strip("\r\n")+"\t"+n[1])
            f2.write("\n")
