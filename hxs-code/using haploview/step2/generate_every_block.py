import os
import sys
with open(sys.argv[1])as f1:
    s = f1.readlines()
    n = s[int(sys.argv[2])]
    n1 = n.strip("\n").split()
    with open(sys.argv[3])as f2:
        with open(sys.argv[4],"w")as f3:
            for line in f2:
                n2 = line.strip("\n").split()
                n21 = n2[1]
                for i in n1:
                    if int(n21) == int(i):
                        f3.write(n2[0]+"\t")

with open(sys.argv[4])as f1:
    with open(sys.argv[5],"w")as f2:
        for line in f1:
            n = line.strip("\n").split()
            f2.write(str(int(n[0].split("_")[0][3:]))+"\t"+n[0].split("_")[1]+"\t"+n[-1].split("_")[1])
            f2.write("\n")
            
            
            
