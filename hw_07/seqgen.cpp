#include <bits/stdc++.h>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    freopen("output.txt", "a", stdout);
    string input = "input/hw_06.txt";

    freopen(input.c_str(), "r", stdin);

    int n_count = 8422401;

    string line;
    map<int, int> bgHist;
    map<int, int> targHist;
    int record = 0;
    while (getline(cin, line, '\n')) {
        if (line.empty()) {
            record = 0;
        } else if (line ==
                   "Read start histogram for non-elevated copy-number "
                   "segments:") {
            record = 1;
        } else if (record == 1) {
            bgHist[stoi(line.substr(line.find_first_of("0123456789"), 1))] =
                stoi(line.substr(line.rfind('=') + 1));
        } else if (line ==
                   "Read start histogram for elevated copy-number segments:") {
            record = 2;
        } else if (record == 2) {
            targHist[stoi(line.substr(line.find_first_of("0123456789"), 1))] =
                stoi(line.substr(line.rfind('=') + 1));
        }
    }
    bgHist[0] -= n_count;
    assert(bgHist.size() == targHist.size());

    int bgSum = accumulate(
        bgHist.begin(), bgHist.end(), 0,
        [](const int prev, const auto& v) { return prev + v.second; });
    int targSum = accumulate(
        targHist.begin(), targHist.end(), 0,
        [](const int prev, const auto& v) { return prev + v.second; });
    int bgTargSum = bgSum + targSum;

    map<int, double> bgFreq;
    map<int, double> targFreq;
    for (auto& [k, v] : bgHist) {
        bgFreq[k] = static_cast<double>(v + targHist[k]) / (bgTargSum);
    }
    for (auto& [k, v] : targHist) {
        targFreq[k] = static_cast<double>(v) / targSum;
    }

    map<int, double> score;
    for (auto& [k, v] : targFreq) {
        score[k] = log2(v / bgFreq[k]);
    }

    string scoreName = "scoring.scheme.txt";
    ofstream scoreFile("input/" + scoreName);
    for (auto& [k, v] : score) {
        if (k == next(score.rend())->first) {
            scoreFile << ">=";
        }
        scoreFile << k << ' ' << v << '\n';
    }
    scoreFile << "D -5" << '\n';
    scoreFile << "S 5" << '\n';
    scoreFile.close();

    string seqName = "chm13.chr16.txt";
    ofstream seqFile("input/" + seqName);
    string gen;
    random_device rd;
    mt19937 ranGen(1337);
    vector<double> bgFreqVal;
    for (auto& [k, v] : bgFreq) {
        bgFreqVal.push_back(v);
    }
    discrete_distribution<> dd(bgFreqVal.begin(), bgFreqVal.end());
    for (int i = 0; i < bgTargSum; i++) {
        seqFile << "16" << '\t' << i + 1 << '\t' << dd(ranGen) << '\n';
    }
    seqFile.close();

    cout << fixed << setprecision(4);
    cout << "Background frequencies:" << '\n';
    for (auto& [k, v] : bgFreq) {
        if (k == next(bgFreq.rend())->first) {
            cout << ">=";
        }
        cout << k << '=' << v << '\n';
    }
    cout << '\n';
    cout << "Target frequencies:" << '\n';
    for (auto& [k, v] : targFreq) {
        if (k == next(targFreq.rend())->first) {
            cout << ">=";
        }
        cout << k << '=' << v << '\n';
    }
    cout << '\n';
    cout << "Scoring scheme:" << '\n';
    for (auto& [k, v] : score) {
        if (k == next(score.rend())->first) {
            cout << ">=";
        }
        cout << k << '=' << v << '\n';
    }
    cout << '\n';

    return 0;
}