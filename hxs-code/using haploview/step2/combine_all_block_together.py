import os
import sys

with open(sys.argv[1])as f1:
    with open(sys.argv[2],"a")as f2:
        for line in f1:
            f2.write(line)
