#include <bits/stdc++.h>

using namespace std;

struct Suffix {
    string_view s;
    int i;
    int group;
};

bool lexico(const Suffix &a, const Suffix &b) { return a.s < b.s; }

template <size_t N>
vector<Suffix> sufArr(array<int, N> &seqLens, string &fSeq) {
    vector<Suffix> sa;
    int fSeqLen = fSeq.length();
    sa.reserve(fSeqLen);
    int gNum = seqLens[0];
    int group = 0;
    const char *pfSeq = fSeq.data();
    for (int i = 0; i < fSeqLen; i++) {
        string_view suffix = string_view(pfSeq + i, fSeqLen - i);
        Suffix curr = {suffix, i, group};
        sa.push_back(curr);
        if (i == gNum) {
            gNum += seqLens[++group] + 1;
        }
    }
    sort(sa.begin(), sa.end(), lexico);
    return sa;
}

int cMatch(string_view a, string_view b) {
    int match = 0;
    for (int i = 0; i < min(a.length(), b.length()); i++) {
        if (a.at(i) == b.at(i) && a.at(i) != '$') {
            match++;
        } else {
            break;
        }
    }
    return match;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    freopen("output.txt", "w", stdout);
    // const char *inputs[2] = {"input/CP001872.fna", "input/CP003913.fna"};
    const char *inputs[2] = {"input/hg38_chr10_90Mb_100Mb.fa",
                             "input/mm39_chr19_35Mb_45Mb.fa"};
    array<string, size(inputs)> inputNames;
    for (int i = 0; i < inputNames.size(); i++) {
        string inputString = string(inputs[i]);
        inputNames[i] = inputString.substr(inputString.find('/') + 1);
    }

    // Concatenate Genomic Sequences
    int a = 0, t = 0, g = 0, c = 0, n = 0, f = 0;
    array<string, size(inputs) + 1> sequences;
    array<int, size(inputs) + 1> seqLens;
    int counter = 0;
    for (const char *input : inputs) {
        cin.clear();
        freopen(input, "r", stdin);
        // cin.ignore(numeric_limits<streamsize>::max(), '\n');
        string header;
        getline(cin, header, '\n');
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
                            a++;
                            break;
                        case 'T':
                            t++;
                            break;
                        case 'G':
                            g++;
                            break;
                        case 'C':
                            c++;
                            break;
                        case 'N':
                            n++;
                            break;
                    }
                }
                sequences[counter] += line;
            }
        }
        seqLens[counter] = sequences[counter].length();

        cout << "Fasta " + to_string(counter + 1) + ": " + inputNames[counter]
             << '\n';
        cout << "Non-alphabetic characters: " << f << '\n';
        cout << header << '\n';
        cout << "*=" << a + c + g + t + n << '\n';
        cout << "A=" << a << '\n';
        cout << "C=" << c << '\n';
        cout << "G=" << g << '\n';
        cout << "T=" << t << '\n';
        cout << "N=" << n << '\n';
        cout << '\n';
        a = 0, t = 0, g = 0, c = 0, n = 0, f = 0;

        counter++;
    }

    string rRev(sequences[counter - 1]);
    reverse(rRev.begin(), rRev.end());
    transform(rRev.begin(), rRev.end(), rRev.begin(), [](char ch) {
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
    sequences[counter] = rRev;
    seqLens[counter] = rRev.length();

    // Generate Suffix Array
    string fSeq;
    for (string seq : sequences) {
        fSeq += seq + '$';
    }
    vector<Suffix> sa = sufArr(seqLens, fSeq);

    // Find Common Substrings
    map<string_view, int> maxSuf;
    map<int, int> frequency;
    int maxMatch = 0;
    int topIdx = -1;
    int botIdx = size(sequences);
    for (int i = size(sequences); i < sa.size(); i++) {
        int currMatch = 0;
        int topMatch = 0;
        int botMatch = 0;
        if (sa[i].group != 0) {
            topIdx = i;
        }
        while (botIdx < sa.size() && (sa[botIdx].group == 0 || botIdx <= i)) {
            botIdx++;
        }
        if (sa[i].group == 0 && topIdx != -1) {
            topMatch = cMatch(sa[i].s, sa[topIdx].s);
        }
        if (sa[i].group == 0 && botIdx != sa.size()) {
            botMatch = cMatch(sa[i].s, sa[botIdx].s);
        }
        if (topMatch != 0 || botMatch != 0) {
            currMatch = max(topMatch, botMatch);
            frequency[currMatch]++;
        }
        if (currMatch > maxMatch) {
            maxMatch = currMatch;
            maxSuf.clear();
            maxSuf[string_view(sa[i].s.data(), maxMatch)] = i;
        } else if (currMatch == maxMatch) {
            maxSuf[string_view(sa[i].s.data(), maxMatch)] = i;
        }
    }

    cout << "Match Length Histogram:" << '\n';
    for (const auto &[k, v] : frequency) {
        cout << k << ' ' << v << '\n';
    }
    cout << '\n';

    cout << "The longest match length: " << maxMatch << '\n';
    cout << "Number of match strings: " << maxSuf.size() << '\n';
    cout << '\n';

    for (const auto &[maxSeq, idx] : maxSuf) {
        cout << "Match string: " << maxSeq << '\n';
        cout << "Description: This sequence comes from [look up entry in .gbff "
                "annotation file using the position information below]"
             << '\n';
        cout << '\n';

        int topMaxIdx = idx - 1;
        int botMaxIdx = idx + 1;
        while (topMaxIdx >= 0 &&
               sa[topMaxIdx].s.substr(0, maxMatch) == maxSeq) {
            topMaxIdx--;
        }
        while (botMaxIdx < sa.size() &&
               sa[botMaxIdx].s.substr(0, maxMatch) == maxSeq) {
            botMaxIdx++;
        }
        topMaxIdx++;
        botMaxIdx--;

        for (int i = topMaxIdx; i <= botMaxIdx; i++) {
            int position = sa[i].i + 1;
            for (int j = 0; j < sa[i].group; j++) {
                int tosub = seqLens[j];
                position -= seqLens[j] + 1;
            }

            string strand = "forward";
            int nIdx = sa[i].group;
            if (nIdx == size(inputs)) {
                nIdx--;
                strand = "reverse";
            }

            cout << "Fasta: " << inputNames[nIdx] << '\n';
            cout << "Position: " << position << '\n';
            cout << "Strand: " << strand << '\n';
            cout << '\n';
        }
    }

    return 0;
}