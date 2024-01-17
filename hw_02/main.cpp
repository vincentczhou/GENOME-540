#include <bits/stdc++.h>

using namespace std;

void printData(int num, string inputName, int nonAlpha, string header,
               int total, map<char, int>& nuFreq, map<char, float>& nuPer,
               map<string, int>& diFreq, map<string, float>& diPer,
               map<string, float>& coPer) {
    // std::map sorts keys!!!
    cout << "Fasta " << num << ": " << inputName << '\n';
    cout << "Non-alphabetic characters: " << nonAlpha << '\n';
    cout << header << '\n';
    cout << "*=" << total << '\n';
    for (const auto& [k, v] : nuFreq) {
        cout << k << '=' << v << '\n';
    }
    cout << "N=" << nuFreq['N'] << '\n';
    cout << '\n';
    cout << "Nucleotide Frequencies:" << '\n';
    for (const auto& [k, v] : nuPer) {
        cout << k << "=" << fixed << setprecision(4) << v << '\n';
    }
    cout << '\n';
    cout << "Dinucleotide Count Matrix:";
    char prevChar;
    char currChar;
    for (const auto& [k, v] : diFreq) {
        currChar = k[0];
        if (currChar != prevChar) {
            prevChar = currChar;
            cout << '\n';
            cout << currChar << '=' << v;
        } else {
            cout << ' ' << v;
        }
    }
    cout << '\n';
    cout << '\n';
    cout << "Dinucleotide Frequency Matrix:";
    for (const auto& [k, v] : diPer) {
        currChar = k[0];
        if (currChar != prevChar) {
            prevChar = currChar;
            cout << '\n';
            cout << currChar << '=' << fixed << setprecision(4) << v;
        } else {
            cout << ' ' << fixed << setprecision(4) << v;
        }
    }
    cout << '\n';
    cout << '\n';
    cout << "Conditional Frequency Matrix:";
    for (const auto& [k, v] : coPer) {
        currChar = k[0];
        if (currChar != prevChar) {
            prevChar = currChar;
            cout << '\n';
            cout << currChar << '=' << fixed << setprecision(4) << v;
        } else {
            cout << ' ' << fixed << setprecision(4) << v;
        }
    }
    cout << '\n';
    cout << '\n';
}

string conGen(int length, map<string, float>& coPer) {
    map<char, discrete_distribution<>> prob;
    vector<float> currV;
    char prevChar;
    char currChar;
    for (const auto& [k, v] : coPer) {
        currChar = k[0];
        if (currChar != prevChar) {
            if (prevChar) {
                prob[prevChar] =
                    discrete_distribution(currV.begin(), currV.end());
            }
            currV.clear();
            currV.push_back(v);
            prevChar = currChar;
        } else {
            currV.push_back(v);
        }
    }
    prob[prevChar] = discrete_distribution(currV.begin(), currV.end());

    map<int, char> cMap = {{0, 'A'}, {1, 'C'}, {2, 'G'}, {3, 'T'}};

    string gen;
    random_device rd;
    mt19937 ranGen(1337);
    discrete_distribution<> initial({0.25, 0.25, 0.25, 0.25});

    gen += cMap[initial(ranGen)];

    for (int i = 0; i < length - 1; i++) {
        gen += cMap[prob[gen[i]](ranGen)];
        // vector<double> currProb = prob[gen[i]].probabilities();
        // for (auto n : currProb) {
        //     cout << n << ' ';
        // }
        // cout << '\n';
    }

    return gen;
}

map<string, float> cPer(map<string, float>& diPer, int dtotal) {
    map<string, float> coPer;
    map<char, float> cTotal;
    float currTotal;
    char prevChar;
    char currChar;
    for (const auto& [k, v] : diPer) {
        currChar = k[0];
        if (currChar != prevChar) {
            cTotal[currChar] = v;
            prevChar = currChar;
        } else {
            cTotal[currChar] += v;
        }
    }
    for (const auto& [k, v] : diPer) {
        coPer[k] = (float)v / cTotal[k[0]];
    }
    return coPer;
}

map<string, float> dPer(map<string, int>& diFreq, int dtotal) {
    map<string, float> diPer;
    for (const auto& [k, v] : diFreq) {
        diPer[k] = (float)v / dtotal;
    }
    return diPer;
}

map<string, int> dFreq(string& seq, int len) {
    map<string, int> diFreq;
    for (int i = 0; i < len - 1; i++) {
        diFreq[string() + seq[i] + seq[i + 1]]++;
    }
    return diFreq;
}

map<char, float> nPer(map<char, int>& nuFreq, int total) {
    map<char, float> nuPer;
    for (const auto& [k, v] : nuFreq) {
        nuPer[k] = (float)v / total;
    }
    return nuPer;
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
        }
    }
    return nuFreq;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    freopen("output.txt", "w", stdout);
    // const char* input = "input/CP003913.fna";
    const char* input = "input/mm39_chr19_35Mb_45Mb.fa";
    string inputName = string(input);
    inputName = inputName.substr(inputName.find('/') + 1);

    freopen(input, "r", stdin);

    string header;
    getline(cin, header, '\n');
    string sequence;
    int seqLen = 0;
    int n = 0, f = 0;
    map<char, int> nuFreq;
    string line;
    while (cin >> line) {
        int linLen = line.length();
        transform(line.begin(), line.end(), line.begin(), ::toupper);
        if (isdigit(line[0])) {
            f += linLen;
        } else {
            for (int i = 0; i < linLen; i++) {
                char nucleotide = line[i];
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
                        n++;
                        break;
                }
            }
            sequence += line;
        }
    }
    seqLen = sequence.length();
    int total = accumulate(nuFreq.begin(), nuFreq.end(), 0,
                           [](const int previous, const auto& element) {
                               return previous + element.second;
                           }) +
                n;
    map<char, float> nuPer = nPer(nuFreq, total);
    map<string, int> diFreq = dFreq(sequence, seqLen);
    int dtotal = accumulate(diFreq.begin(), diFreq.end(), 0,
                            [](const int previous, const auto& element) {
                                return previous + element.second;
                            });
    map<string, float> diPer = dPer(diFreq, dtotal);
    map<string, float> coPer = cPer(diPer, dtotal);
    printData(1, inputName, f, header, total, nuFreq, nuPer, diFreq, diPer,
              coPer);

    map<string, float> ePer;
    map<string, float> zeroPer;
    for (const auto& [k, v] : coPer) {
        ePer[k] = (float)0.25;
        zeroPer[k] = nuPer[k[1]];
    }

    string eName = "simulated_equal_freq.fa";
    string eHeader = ">simulated_equal_freq";
    string eSeq = conGen(seqLen, ePer);
    ofstream eFile("input/" + eName);
    eFile << eHeader << '\n';
    eFile << eSeq;
    eFile.close();
    map<char, int> enuFreq = nFreq(eSeq);
    map<char, float> enuPer = nPer(enuFreq, total);
    map<string, int> ediFreq = dFreq(eSeq, seqLen);
    map<string, float> ediPer = dPer(ediFreq, dtotal);
    map<string, float> ecoPer = cPer(ediPer, dtotal);
    printData(2, eName, 0, eHeader, total, enuFreq, enuPer, ediFreq, ediPer,
              ecoPer);

    string zeroName = "simulated_markov_0.fa";
    string zeroHeader = ">simulated_markov_0";
    string zeroSeq = conGen(seqLen, zeroPer);
    ofstream zeroFile("input/" + zeroName);
    zeroFile << zeroHeader << '\n';
    zeroFile << zeroSeq;
    zeroFile.close();
    map<char, int> zeronuFreq = nFreq(zeroSeq);
    map<char, float> zeronuPer = nPer(zeronuFreq, total);
    map<string, int> zerodiFreq = dFreq(zeroSeq, seqLen);
    map<string, float> zerodiPer = dPer(zerodiFreq, dtotal);
    map<string, float> zerocoPer = cPer(zerodiPer, dtotal);
    printData(3, zeroName, 0, zeroHeader, total, zeronuFreq, zeronuPer,
              zerodiFreq, zerodiPer, zerocoPer);

    string cName = "simulated_markov_1.fa";
    string cHeader = ">simulated_markov_1";
    string cSeq = conGen(seqLen, coPer);
    ofstream cFile("input/" + cName);
    cFile << cHeader << '\n';
    cFile << cSeq;
    cFile.close();
    map<char, int> cnuFreq = nFreq(cSeq);
    map<char, float> cnuPer = nPer(cnuFreq, total);
    map<string, int> cdiFreq = dFreq(cSeq, seqLen);
    map<string, float> cdiPer = dPer(cdiFreq, dtotal);
    map<string, float> ccoPer = cPer(cdiPer, dtotal);
    printData(4, cName, 0, cHeader, total, cnuFreq, cnuPer, cdiFreq, cdiPer,
              ccoPer);

    return 0;
}
