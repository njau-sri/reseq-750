import sys

MAX_HET = float(sys.argv[2])

with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith("##fileformat="):
            print line.strip()
            continue
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            print line.strip()
            break
    for line in f:
        v = line.split()[9:]
        m, n = 0, 0
        for e in v:
            a, b = e[:3].split("/")
            if a != "." and b != ".":
                n += 1
                if a != b:
                    m += 1
        if float(m)/n <= MAX_HET:
            print line.strip()
