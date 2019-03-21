import os
import sys
with open(sys.argv[1])as f1:
    s = f1.readlines()
    n = s[int(sys.argv[2])]
    n1 = n.strip("\n").split()
    n10 = n1[0]
    n11 = n1[1]
    n12 = n1[2]
    with open(sys.argv[3])as f2:
        with open(sys.argv[4],"w")as f3:
            for line in f2:
                n2 = line.strip("\n").split("_")
                n20 = n2[0]
                n21 = n2[1]
                if int(n20[3:]) == int(n10):
                    if int(n21) >= int(n11) and int(n21) <= int(n12):
                        f3.write(line)
                
            
            
            
