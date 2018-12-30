import sys

f1 = open(sys.argv[1])
f2 = open(sys.argv[2])

for line in f1:
    if line.startswith("##"):
        continue
    if line.startswith("#CHROM"):
        ind1 = line.split()[9:]
        break

for line in f2:
    if line.startswith("##"):
        continue
    if line.startswith("#CHROM"):
        ind2 = line.split()[9:]
        break

ind = sorted(ind1)
if sorted(ind2) != ind:
    raise IOError

f = open(sys.argv[3],"w")
f.write("##fileformat=VCFv4.2\n")
f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
f.write("\t".join(ind))
f.write("\n")

idx1 = [ind1.index(e) for e in ind]
idx2 = [ind2.index(e) for e in ind]

for line1 in f1:
    line2 = f2.next()

    v1, v2 = line1.split(), line2.split()
    if v1[:5] != v2[:5]:
        raise IOError

    chr = v1[0]
    if chr[:3].lower() == "chr":
        chr = chr[3:]
    while chr.startswith("0"):
        chr = chr[1:]

    v = [chr, v1[1], "_".join(v1[:2]) if v1[2] == "." else v1[2], v1[3], v1[4], ".", ".", ".", "GT"]

    for i in range(len(ind)):
        i1, i2 = idx1[i] + 9, idx2[i] + 9
        if v1[i1][:3] == v2[i2][:3]:
            v.append(v1[i1][:3])
        else:
            v.append("./.")

    f.write("\t".join(v))
    f.write("\n")
