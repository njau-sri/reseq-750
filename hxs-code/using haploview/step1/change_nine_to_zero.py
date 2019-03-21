import os
import sys

with open(sys.argv[1])as f1:
    with open(sys.argv[2],"w")as f2:
        for line in f1:
            n = line.strip("\r\n").split()
            n1 = n[0:5]
            n2 = n[6:]
            f2.write("\t".join(x for x in n1)+"\t"+"0"+"\t"+"\t".join(y for y in n2))
            f2.write("\n")
        
