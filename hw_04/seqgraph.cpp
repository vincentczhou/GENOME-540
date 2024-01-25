#include <bits/stdc++.h>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    freopen("output.txt", "a", stdout);
    // string input = "input/CP003508.fna";
    string input = "input/s_pyogenes.fa";
    string inputName = input.substr(input.find('/') + 1);
    string sinput = "input/scoringscheme.txt";
    string sinputName = sinput.substr(sinput.find('/') + 1);

    freopen(input.c_str(), "r", stdin);
    string header;
    getline(cin, header, '\n');

    string sequence;
    int seqLen = 0;
    map<char, int> nuFreq;
    int f = 0;
    string line;
    while (cin >> line) {
        int linLen = line.length();
        transform(line.begin(), line.end(), line.begin(), ::toupper);
        if (isdigit(line[0])) {
            f += linLen;
        } else {
            for (char& nucleotide : line) {
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
            sequence += line;
            seqLen += line.length();
        }
    }

    cin.clear();
    freopen(sinput.c_str(), "r", stdin);

    map<char, double> ss;
    while (getline(cin, line)) {
        ss[line.front()] = stod(line.substr(line.find(' ') + 1));
    }

    string sqName = "seqgraph.txt";
    ofstream sqFile("input/" + sqName);
    for (int i = 0; i < seqLen; i++) {
        sqFile << "V " << i << '\n';
    }
    for (int i = 0; i < seqLen - 1; i++) {
        sqFile << "E " << sequence[i] << ' ' << i << ' ' << i + 1 << ' '
               << ss[sequence[i]] << '\n';
    }
    sqFile << "E " << sequence[seqLen - 1] << ' ' << seqLen - 1 << ' ' << seqLen
           << ' ' << ss[sequence[seqLen - 1]];
    sqFile.close();

    cout << "Fasta: " << inputName << '\n';
    cout << "Non-alphabetic characters: " << f << '\n';
    cout << header << '\n';
    cout << "*=" << seqLen << '\n';
    for (auto& [k, v] : nuFreq) {
        if (k == 'N') {
            continue;
        }
        cout << k << '=' << v << '\n';
    }
    cout << "N=" << nuFreq['N'] << '\n';
    cout << '\n';

    return 0;
}