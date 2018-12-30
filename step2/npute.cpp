#include <string>
#include <vector>
#include <limits>
#include <numeric>
#include <fstream>
#include <iostream>
#include <algorithm>


// C++ implementation of NPUTE (http://compgen.unc.edu/software)

// Roberts A, McMillan L, Wang W, Parker J, Rusyn I, Threadgill D (2007) Inferring missing
//     genotypes in large SNP panels using fast nearest-neighbor searches over sliding windows.
//     Bioinformatics 23:i401-407


using std::size_t;


template<typename T>
void order(const std::vector<T> &v, std::vector<int> &z)
{
    z.resize(v.size());
    std::iota(z.begin(), z.end(), 0);
    std::sort(z.begin(), z.end(), [&v](int i, int j) { return v[i] < v[j]; });
}


template<typename T>
void add(const std::vector<T> &x, const std::vector<T> &y, std::vector<T> &z)
{
    int n = x.size();
    for (int i = 0; i < n; ++i)
        z[i] = x[i] + y[i];
}


template<typename T>
void subtract(const std::vector<T> &x, const std::vector<T> &y, std::vector<T> &z)
{
    int n = x.size();
    for (int i = 0; i < n; ++i)
        z[i] = x[i] - y[i];
}


// CircularQueue.py

class CircularQueue
{
public:

    void init(int L, int width)
    {
        mid_ = -1;
        L_ = half_ = L;
        height_ = 2*L+2;
        queue.assign(height_, std::vector<int>(width, 0));
    }

    void getEnds(std::vector<int> &top, std::vector<int> &bottom) const
    {
        auto i = (mid_ + L_) % height_;
        top = queue[i];

        auto j = (mid_ - L_ - 1) % height_;
        if (j < 0)
            j += queue.size();
        bottom = queue[j];
    }

    void enqueue(const std::vector<int> &e)
    {
        mid_ = (mid_ + 1) % height_;
        auto nextIn = (mid_ + half_) % height_;
        queue[nextIn] = e;
    }

public:
    std::vector< std::vector<int> > queue;

private:
    int L_ = 0;
    int height_ = 0;
    int half_ = 0;
    int mid_ = 0;
};


// NPUTE.py, SNPData.py

class NPUTE
{
public:

    void readInData(const std::string &inFile)
    {
        std::ifstream ifs(inFile);
        if ( ! ifs ) {
            std::cerr << "ERROR: can't open file: " << inFile << "\n";
            return;
        }
        size_t n = 0;
        for (std::string line; std::getline(ifs,line); ) {
            removeExtraChars(line);
            if (n == 0)
                n = line.size();
            if (line.size() != n) {
                std::cerr << "ERROR: inconsistent number of samples\n";
                return;
            }
            auto p = getAlleles(line);
            if (p.first != 0)
                std::replace(line.begin(), line.end(), p.first, '0');
            if (p.second != 0)
                std::replace(line.begin(), line.end(), p.second, '1');
            snps_.push_back(line);
            nucs_.push_back(p);
        }
        numSamps_ = n;
    }

    void genMismatchVectors(const std::string &g, std::vector<int> &d) const
    {
        int n = numSamps_;

        d.clear();
        d.reserve(n*(n-1)/2);

        for (int i = 0; i < n; ++i) {
            auto a = g[i];
            for (int j = i + 1; j < n; ++j) {
                auto b = g[j];
                if (a != '?' && b != '?')
                    d.push_back(a == b ? 0 : 2);
                else
                    d.push_back(1);
            }
        }
    }

    void genExtractIndices()
    {
        int n = numSamps_;
        extractIndices_.assign(n, std::vector<int>(n, -1));
        for (int i = 0, k = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j, ++k)
                extractIndices_[i][j] = extractIndices_[j][i] = k;
    }

    void extractRow(const std::vector<int> &mmv, int rowNum, std::vector<int> &row)
    {
        row.clear();
        for (auto i : extractIndices_[rowNum])
            row.push_back(i < 0 ? std::numeric_limits<int>::max() : mmv[i]);
    }

    char getMinIMp(const std::string &snp, const std::vector<int> &sA, const std::vector<int> &score) const
    {
        int lastM = 0;
        int points = 0;
        char winner = '0';
        for (auto i : sA) {
            if (snp[i] == '?')
                continue;
            auto m = score[i];
            if (m != lastM && points > 0)
                return winner;
            if (snp[i] == winner)
                points += 1;
            else {
                if (points == 0) {
                    winner = snp[i];
                    points = 1;
                }
                else
                    points -= 1;
            }
            lastM = m;
        }
        return winner;
    }

    void imputeSNP(const std::vector<int> &mmv, const std::string &in, std::string &out)
    {
        out = in;
        int n = numSamps_;
        std::vector<int> score, sA;
        for (int i = 0; i < n; ++i) {
            if (in[i] != '?')
                continue;
            extractRow(mmv, i, score);
            order(score, sA);
            sA.pop_back();
            out[i] = getMinIMp(in, sA, score);
        }
    }

    void impute(int L, const std::string &infile, const std::string &outfile)
    {
        readInData(infile);
        genExtractIndices();

        int n = numSamps_;
        int m = snps_.size();

        std::cerr << "INFO: " << n << " individuals and " << m << " SNPs were observed\n";

        if (n < 2 && m < L) {
            std::cerr << "ERROR: not enough data for NPUTE\n";
            return;
        }

        std::cerr << "INFO: imputing with window size " << L << "\n";

        int nn = n*(n-1)/2;
        CircularQueue vectorQueue;
        vectorQueue.init(L, nn);

        std::vector<int> mmv, acc(nn, 0);

        for (int i = 0; i < L; ++i) {
            genMismatchVectors(snps_[i], mmv);
            vectorQueue.queue[i] = mmv;
            add(acc, mmv, acc);
        }

        std::ofstream ofs(outfile);
        if ( ! ofs ) {
            std::cerr << "ERROR: can't open file: " << outfile << "\n";
            return;
        }

        std::string snp, out;
        std::vector<int> mmv0(nn,0), top, bottom;

        int step = m / 100;

        for (int i = 0; i < m; ++i) {
            if (i + L < m) {
                genMismatchVectors(snps_[i+L], mmv);
                vectorQueue.enqueue(mmv);
            }
            else
                vectorQueue.enqueue(mmv0);

            vectorQueue.getEnds(top, bottom);
            add(acc, top, acc);
            subtract(acc, bottom, acc);

            snp = snps_[i];
            if (snp.find('?') != std::string::npos)
                imputeSNP(acc, snps_[i], snp);

            auto p = nucs_[i];
            out.clear();
            for (auto a : snp) {
                out.push_back(a == '0' ? p.first : p.second);
                out.push_back(',');
            }
            out.pop_back();
            ofs << out << "\n";

            if (step > 100 && ((i+1) % step == 0 || i == m - 1))
                std::cerr << "INFO: imputed/total: " << i+1 << "/" << m << "\n";
        }
    }

private:

    void removeExtraChars(std::string &s) const
    {
        s.erase(s.find_last_not_of(" \t\r\n") + 1);
        auto end = std::remove(s.begin(), s.end(), ',');
        end = std::remove(s.begin(), end, ' ');
        end = std::remove(s.begin(), end, '\t');
        s.erase(end, s.end());
    }

    std::pair<char,char> getAlleles(const std::string &s)
    {
        int q = (s.size() - std::count(s.begin(), s.end(), '?') + 1) / 2;
        char major = 0, minor = 0;
        for (auto c : s) {
            if (c == '?' || c == major || c == minor)
                continue;
            if (major == 0 && std::count(s.begin(),s.end(),c) >= q) {
                major = c;
                continue;
            }
            if (minor == 0 && std::count(s.begin(),s.end(),c) <= q) {
                minor = c;
                continue;
            }
            std::cerr << "WARNING: SNP is not ternary\n";
        }
        return std::make_pair(major, minor);
    }

private:
    std::vector<std::string> snps_;
    std::vector< std::pair<char,char> > nucs_;
    std::vector< std::vector<int> > extractIndices_;
    std::size_t numSamps_ = 0;
};


int main(int argc, char *argv[])
{
    if (argc < 4) {
        std::cerr << "usage: " << argv[0] << " infile outfile window\n";
        return 1;
    }

    std::string infile = argv[1];
    std::string outfile = argv[2];
    int L = std::stoi(argv[3]);

    NPUTE npute;
    npute.impute(L, infile, outfile);

    return 0;
}
