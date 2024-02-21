#include <bits/stdc++.h>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    freopen("output.txt", "a", stdout);
    string input = "output_main_real.txt";
    string sinput = "output_main_simulated.txt";

    freopen(input.c_str(), "r", stdin);

    string line;
    vector<double> rSegs;
    bool record = false;
    while (getline(cin, line, '\n')) {
        if (record && line.empty()) {
            break;
        } else if (line == "Segment List:") {
            record = true;
        } else if (record) {
            rSegs.push_back(stod(line.substr(line.rfind(' ') + 1)));
        }
    }

    cin.clear();
    freopen(sinput.c_str(), "r", stdin);

    vector<double> sSegs;
    record = false;
    while (getline(cin, line, '\n')) {
        if (record && line.empty()) {
            break;
        } else if (line == "Segment List:") {
            record = true;
        } else if (record) {
            sSegs.push_back(stod(line.substr(line.rfind(' ') + 1)));
        }
    }

    map<int, int> rNum;
    map<int, int> sNum;
    for (int i = 5; i < 31; i++) {
        rNum[i] = 0;
        sNum[i] = 0;
    }
    for (double seg : rSegs) {
        for (int i = 5; i < 31; i++) {
            if (seg >= i) {
                rNum[i]++;
            } else {
                break;
            }
        }
    }
    for (double seg : sSegs) {
        for (int i = 5; i < 31; i++) {
            if (seg >= i) {
                sNum[i]++;
            } else {
                break;
            }
        }
    }

    cout << fixed << setprecision(2);
    cout << "Real data:" << '\n';
    for (auto& [k, v] : rNum) {
        cout << k << ' ' << v << '\n';
    }
    cout << '\n';
    cout << '\n';
    cout << "Simulated data:" << '\n';
    for (auto& [k, v] : sNum) {
        cout << k << ' ' << v << '\n';
    }
    cout << '\n';
    cout << '\n';
    cout << "Ratios of simulated data:" << '\n';
    for (int i = 5; i < 30; i++) {
        cout << "N_seg(" << i << ")/N_seg(" << i + 1 << ") ";
        if (sNum[i + 1] == 0) {
            cout << -1.00;
        } else {
            cout << static_cast<double>(sNum[i]) / sNum[i + 1];
        }
        cout << '\n';
    }
    cout << '\n';
    cout << '\n';

    return 0;
}