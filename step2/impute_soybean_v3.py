import sys
import gzip
import subprocess
from multiprocessing import Pool


def split_vcf(infile):
    if infile.endswith(".vcf"):
        ifs = open(infile)
    elif infile.endswith(".vcf.gz"):
        ifs = gzip.open(infile)
    else:
        raise IOError
    header = []
    for line in ifs:
        if line.startswith("##fileformat="):
            header.append(line.strip())
            continue
        if line.startswith("#CHROM"):
            header.append(line.strip())
            break
    ofs = None
    chrid = ""
    chrlst = []
    for line in ifs:
        v = line.split(None,1)
        if v[0] != chrid:
            if ofs:
                ofs.close()
                ofs = None
            if v[0] in chrlst:
                raise IOError
            chrid = v[0]
            ofs = open(chrid+".vcf","w")
            chrlst.append(chrid)
            ofs.write("\n".join(header))
            ofs.write("\n")
        ofs.write(line.strip())
        ofs.write("\n")
    ofs.close()
    return chrlst, header


def convert_npute(chrid):
    ifs = open(chrid+".vcf")
    ofs = open(chrid+".npute.in","w")
    for line in ifs:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            break
    for line in ifs:
        v = line.split()
        g = []
        for e in v[9:]:
            if e[0] == "." or e[2] == "." or e[0] != e[2]:
                g.append("?")
            else:
                g.append(v[3] if e[0] == "0" else v[4])
        ofs.write(",".join(g))
        ofs.write("\n")
    ifs.close()
    ofs.close()


def read_best_window(infile):
    x, y = 0, 0.0
    with open(infile) as f:
        for line in f:
            v = line.strip().split(",")
            yi = float(v[1])
            if yi > y:
                y = yi
                x = int(v[0])
    return x


def worker(chrid):
    convert_npute(chrid)
    subprocess.call(["./npute-test-lnx64", chrid+".npute.in", chrid+".windows.csv"] + [str(e) for e in range(10,101,10)])
    w = read_best_window(chrid+".windows.csv")
    subprocess.call(["./npute-lnx64", chrid+".npute.in", chrid+".npute.out", str(w)])


def main():
    chrlst, header = split_vcf(sys.argv[1])
    pool = Pool(8)
    pool.map(worker, chrlst)
    ofs = open("imputed.vcf","w")
    ofs.write("\n".join(header))
    ofs.write("\n")
    for chrid in chrlst:
        f1 = open(chrid+".vcf")
        for line in f1:
            if line.startswith("#CHROM"):
                break
        f2 = open(chrid+".npute.out")
        for l1 in f1:
            l2 = f2.next()
            v1 = l1.split(None,9)
            v2 = l2.strip().upper().split(",")
            g = []
            for e in v2:
                if e == v1[3]:
                    g.append("0/0")
                elif e == v1[4]:
                    g.append("1/1")
                else:
                    print v1[3], v1[4], e
                    raise RuntimeError
            ofs.write("\t".join(v1[:9]))
            ofs.write("\t")
            ofs.write("\t".join(g))
            ofs.write("\n")


if __name__ == "__main__":
    main()
