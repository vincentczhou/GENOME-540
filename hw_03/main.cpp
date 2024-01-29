#include <bits/stdc++.h>

using namespace std;

template <typename T>
struct CDSite {
    int pos;
    double score;
    T seq;
    int strand;
};

void printData(map<char, int>& nuFreq, map<char, double>& nuPer,
               map<int, map<char, int>>& countMatrix,
               map<int, map<char, double>>& frequencyMatrix,
               map<int, map<char, double>>& weightMatrix, double maxScore,
               map<int, vector<CDSite<string>>>& cscoreHist,
               map<int, vector<CDSite<string_view>>>& gscoreHist,
               vector<CDSite<string_view>>& uniqueSites) {
    // std::map sorts keys!!!
    cout << "Nucleotide Histogram:" << '\n';
    for (const auto& [k, v] : nuFreq) {
        if (k == 'N') {
            continue;
        }
        cout << k << '=' << v << '\n';
    }
    cout << "N=" << nuFreq['N'] << '\n';
    cout << '\n';
    cout << fixed << setprecision(4);
    cout << "Background Frequency:" << '\n';
    for (const auto& [k, v] : nuPer) {
        if (k == 'N') {
            continue;
        }
        cout << k << '=' << v << '\n';
    }
    cout << '\n';
    cout << "Count Matrix:" << '\n';
    for (auto& [k, v] : countMatrix) {
        cout << k - 10;
        for (const auto& [n, c] : v) {
            if (n == 'N') {
                continue;
            }
            cout << ' ' << c;
        }
        cout << '\n';
    }
    cout << '\n';
    cout << "Frequency Matrix:" << '\n';
    for (auto& [k, v] : frequencyMatrix) {
        cout << k - 10;
        for (const auto& [n, c] : v) {
            if (n == 'N') {
                continue;
            }
            cout << ' ' << c;
        }
        cout << '\n';
    }
    cout << '\n';
    cout << "Weight Matrix:" << '\n';
    for (auto& [k, v] : weightMatrix) {
        cout << k - 10;
        for (const auto& [n, c] : v) {
            if (n == 'N') {
                continue;
            }
            cout << ' ' << c;
        }
        cout << '\n';
    }
    cout << '\n';
    cout << setprecision(10);
    cout << "Maximum Score: " << maxScore << '\n';
    cout << '\n';
    cout << "Score Histogram CDS:" << '\n';
    bool printLT = false;
    for (const auto& [k, v] : cscoreHist) {
        if (k == numeric_limits<int>::min()) {
            printLT = true;
        } else {
            cout << k << ' ' << v.size() << '\n';
        }
    }
    if (printLT) {
        cout << "lt-50 " << cscoreHist[numeric_limits<int>::min()].size()
             << '\n';
        printLT = false;
    }
    cout << '\n';
    cout << "Score Histogram All:" << '\n';
    for (const auto& [k, v] : gscoreHist) {
        if (k == numeric_limits<int>::min()) {
            printLT = true;
        } else {
            cout << k << ' ' << v.size() << '\n';
        }
    }
    if (printLT) {
        cout << "lt-50 " << gscoreHist[numeric_limits<int>::min()].size()
             << '\n';
        printLT = false;
    }
    cout << '\n';
    cout << setprecision(4);
    cout << "Position List:" << '\n';
    for (const auto& g : uniqueSites) {
        cout << g.pos << ' ' << g.strand << ' ' << g.score << '\n';
    }
    cout << '\n';
}

vector<tuple<deque<array<int, 2>>, int>> locRanges(vector<string>& locs) {
    vector<tuple<deque<array<int, 2>>, int>> lrs;
    for (string& loc : locs) {
        if (loc.find('<') != -1 || loc.find('>') != -1) {
            continue;
        }

        int rStart = loc.rfind('(') + 1;
        int rCount = loc.find(')') - rStart;
        string r = loc.substr(rStart, rCount);

        if (!isdigit(r[0])) {
            continue;
        }

        deque<array<int, 2>> ranges;
        while (r.find(',') != -1) {
            array<int, 2> startend;
            startend[0] = stoi(r);
            r.erase(0, r.find('.') + 2);
            startend[1] = stoi(r);
            r.erase(0, r.find(',') + 1);
            ranges.push_back(startend);
        }
        array<int, 2> startend;
        startend[0] = stoi(r);
        r.erase(0, r.find('.') + 2);
        startend[1] = stoi(r);
        ranges.push_back(startend);

        int strand = 0;
        if (loc.find("complement") != -1) {
            strand = 1;
        }

        tuple<deque<array<int, 2>>, int> curr{ranges, strand};
        lrs.push_back(curr);
    }
    return lrs;
}

vector<CDSite<string>> locCds(
    vector<tuple<deque<array<int, 2>>, int>>& rangestrands, string& sequence) {
    vector<CDSite<string>> cds;
    for (auto& rangestrand : rangestrands) {
        deque<array<int, 2>> ranges = get<0>(rangestrand);
        int strand = get<1>(rangestrand);

        string cd;
        int pos = 0;
        array<int, 2> currRange;
        int offset = 0;

        bool valid = true;

        if (strand == 1) {
            currRange = ranges.back();
            pos = currRange[1];
            for (int i = -10; i < 11; i++) {
                int currIdx = -1 * currRange[1] + i + offset;
                if (currIdx > -1 * currRange[0] && ranges.size() > 1) {
                    ranges.pop_back();
                    currRange = ranges.back();
                    offset = -1 * i;
                    currIdx = -1 * currRange[1] + i + offset;
                }

                if (-1 * currIdx - 1 < 0 ||
                    -1 * currIdx - 1 >= sequence.length()) {
                    valid = false;
                    break;
                }
                cd += sequence.at(-1 * currIdx - 1);
            }
            transform(cd.begin(), cd.end(), cd.begin(), [](char ch) {
                switch (ch) {
                    case 'A':
                        return 'T';
                    case 'T':
                        return 'A';
                    case 'G':
                        return 'C';
                    case 'C':
                        return 'G';
                }
                return ch;
            });
        } else {
            currRange = ranges.front();
            pos = currRange[0];
            for (int i = -10; i < 11; i++) {
                int currIdx = currRange[0] + i + offset;
                if (currIdx > currRange[1] && ranges.size() > 1) {
                    ranges.pop_front();
                    currRange = ranges.front();
                    offset = -1 * i;
                    currIdx = currRange[0] + i + offset;
                }

                if (currIdx - 1 < 0 || currIdx - 1 >= sequence.length()) {
                    valid = false;
                    break;
                }
                cd += sequence.at(currIdx - 1);
            }
        }

        if (valid) {
            cds.push_back(CDSite<string>{pos, 0, cd, strand});
        }
    }
    return cds;
}

map<char, int> nFreq(string& seq) {
    map<char, int> nuFreq;
    int len = seq.length();
    for (int i = 0; i < len; i++) {
        char nucleotide = seq[i];
        switch (nucleotide) {
            case 'A':
                nuFreq['A']++;
                break;
            case 'T':
                nuFreq['T']++;
                break;
            case 'G':
                nuFreq['G']++;
                break;
            case 'C':
                nuFreq['C']++;
                break;
            case 'N':
                nuFreq['N']++;
                break;
        }
    }
    return nuFreq;
}

map<char, double> nPer(map<char, int>& nuFreq, int total) {
    map<char, double> nuPer;
    for (const auto& [k, v] : nuFreq) {
        nuPer[k] = (double)v / total;
    }
    return nuPer;
}

map<int, map<char, int>> cMatrix(vector<CDSite<string>>& cds) {
    map<int, map<char, int>> countMatrix;
    for (CDSite<string>& c : cds) {
        string& s = c.seq;
        for (int i = 0; i < s.length(); i++) {
            countMatrix[i][s.at(i)]++;
        }
    }
    for (auto& [k, v] : countMatrix) {
        v['A'];
        v['C'];
        v['G'];
        v['T'];
        v['N'];
    }
    return countMatrix;
}

map<int, map<char, double>> fMatrix(map<int, map<char, int>>& countMatrix) {
    map<int, map<char, double>> frequencyMatrix;
    for (auto& [k, v] : countMatrix) {
        int total = accumulate(v.begin(), v.end(), 0,
                               [](const int previous, const auto& element) {
                                   return previous + element.second;
                               });
        frequencyMatrix[k] = nPer(v, total);
    }
    return frequencyMatrix;
}

map<int, map<char, double>> wMatrix(
    map<int, map<char, double>>& frequencyMatrix, map<char, double>& nuPer) {
    map<int, map<char, double>> weightMatrix;
    for (auto& [k, v] : frequencyMatrix) {
        map<char, double> pweightMatrix;
        for (auto& [n, c] : v) {
            double per = -99;
            if (c > 0) {
                per = log2(c / nuPer[n]);
            }
            pweightMatrix[n] = per;
        }
        weightMatrix[k] = pweightMatrix;
    }
    return weightMatrix;
}

double mScore(map<int, map<char, double>>& weightMatrix) {
    double maxScore = 0;
    for (auto& [k, v] : weightMatrix) {
        double currMax = 0;
        for (const auto& [n, c] : v) {
            if (n == 'N') {
                continue;
            }
            if (c > currMax) {
                currMax = c;
            }
        }
        maxScore += currMax;
    }
    return maxScore;
}

template <typename T>
map<int, vector<CDSite<T>>> sHist(vector<T>& stringlikes,
                                  map<int, map<char, double>>& weightMatrix) {
    map<int, vector<CDSite<T>>> cscoreHist;
    for (int i = 0; i < stringlikes.size(); i++) {
        CDSite<T> currGS = {i + 1 + 10};
        T& s = stringlikes[i];
        double currScore = 0;
        for (int j = 0; j < s.length(); j++) {
            if (s.at(j) != 'N') {
                currScore += weightMatrix[j][s.at(j)];
            }
        }
        currGS.seq = s;
        currGS.score = currScore;
        int score = (int)floorf(currScore);
        if (score < -50) {
            score = numeric_limits<int>::min();
        }
        cscoreHist[score].push_back(currGS);
    }
    return cscoreHist;
}

vector<string_view> kmer(string& sequence, int n) {
    vector<string_view> mer;
    const char* seq = sequence.data();
    for (int i = 0; i < sequence.length() - n + 1; i++) {
        mer.push_back(string_view(seq + i, n));
    }
    return mer;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    freopen("output.txt", "w", stdout);
    // string input = "input/test1.gbff";
    // string input = "input/test2.gbff";
    string input = "input/s_pyogenes.gbff";
    string inputName = input.substr(input.find('/') + 1);

    freopen(input.c_str(), "r", stdin);

    string sequence;
    vector<string> locs;
    string currLoc;
    string line;
    bool secRec = false;
    bool locRec = false;
    while (cin >> line) {
        if (secRec && line.find_first_of("atcgnATGCN") != -1) {
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            sequence += line;
        } else if (line == "ORIGIN") {
            secRec = true;
        } else if (locRec && line.find('/') == 0) {
            locRec = false;
            locs.push_back(currLoc);
            currLoc.clear();
        } else if (locRec) {
            currLoc += line;
        } else if (line == "CDS") {
            locRec = true;
        }
    }

    string revSeq = sequence;
    reverse(revSeq.begin(), revSeq.end());
    transform(revSeq.begin(), revSeq.end(), revSeq.begin(), [](char ch) {
        switch (ch) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'G':
                return 'C';
            case 'C':
                return 'G';
        }
        return ch;
    });

    // Nucleotide Histogram
    map<char, int> nuFreq = nFreq(sequence);

    // Background Frequency
    map<char, int> nunuFreq = nuFreq;
    map<char, char> rc = {
        {'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}, {'N', 'N'}};
    for (const auto& [k, v] : nuFreq) {
        nunuFreq[k] = v + nuFreq[rc[k]];
    }
    map<char, double> nuPer =
        nPer(nunuFreq, 2 * (sequence.length() - nuFreq['N']));

    // Count Matrix
    vector<tuple<deque<array<int, 2>>, int>> rangestrands = locRanges(locs);
    vector<CDSite<string>> cds = locCds(rangestrands, sequence);
    map<int, map<char, int>> countMatrix = cMatrix(cds);

    // Frequency Matrix
    map<int, map<char, double>> frequencyMatrix = fMatrix(countMatrix);

    // Weight Matrix
    map<int, map<char, double>> weightMatrix = wMatrix(frequencyMatrix, nuPer);

    // Maximum Score
    double maxScore = mScore(weightMatrix);

    // Score Histogram CDS
    vector<string> cdsseqs;
    for (CDSite<string>& g : cds) {
        cdsseqs.push_back(g.seq);
    }
    map<int, vector<CDSite<string>>> cscoreHist =
        sHist<string>(cdsseqs, weightMatrix);

    // Score Histogram All
    vector<string_view> mers = kmer(sequence, 21);
    map<int, vector<CDSite<string_view>>> gscoreHist =
        sHist<string_view>(mers, weightMatrix);
    for (auto& [k, v] : gscoreHist) {
        for (auto& g : v) {
            g.strand = 0;
        }
    }
    vector<string_view> revmers = kmer(revSeq, 21);
    map<int, vector<CDSite<string_view>>> revgscoreHist =
        sHist<string_view>(revmers, weightMatrix);
    for (auto& [k, v] : revgscoreHist) {
        for (auto& g : v) {
            g.strand = 1;
            g.pos = sequence.length() + 1 - g.pos;
        }
        gscoreHist[k].insert(gscoreHist[k].end(), v.begin(), v.end());
    }

    // Position List
    vector<CDSite<string_view>> uniqueSites;
    for (auto& [k, v] : gscoreHist) {
        if (k >= 10) {
            for (auto& g : v) {
                bool unique = true;
                for (CDSite<string>& c : cds) {
                    if (c.pos == g.pos && c.strand == g.strand) {
                        unique = false;
                        break;
                    }
                }
                if (unique) {
                    uniqueSites.push_back(g);
                }
            }
        }
    }
    sort(uniqueSites.begin(), uniqueSites.end(),
         [](auto& ga, auto& gb) { return ga.pos < gb.pos; });

    printData(nuFreq, nuPer, countMatrix, frequencyMatrix, weightMatrix,
              maxScore, cscoreHist, gscoreHist, uniqueSites);

    return 0;
}