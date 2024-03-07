#include <bits/stdc++.h>

using namespace std;

struct Aln {
    int start;
    int end;
    string human;
    string dog;
    string mouse;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    freopen("output.txt", "w", stdout);
    // string input = "input/ENm006_short.aln";
    string input = "input/ENm010.aln";
    string input_neu = "input/STATE1_anc_rep_counts.txt";
    string input_con = "input/STATE2_codon1_2_counts.txt";
    freopen(input.c_str(), "r", stdin);

    vector<Aln> alns;
    string line;
    while (getline(cin, line)) {
        if (line.front() == '#') {
            int d_a = line.find(':');
            int d_b = line.find('-');
            alns.push_back({stoi(line.substr(d_a + 1, d_b - d_a)),
                            stoi(line.substr(d_b + 1))});
        } else if (line.find("hg18") != -1) {
            alns.back().human = line.substr(line.find('\t') + 1);
        } else if (line.find("canFam2") != -1) {
            alns.back().dog = line.substr(line.find('\t') + 1);
        } else if (line.find("mm9") != -1) {
            alns.back().mouse = line.substr(line.find('\t') + 1);
        }
    }

    cin.clear();
    freopen(input_neu.c_str(), "r", stdin);

    map<string, int> counts_1;
    while (getline(cin, line)) {
        stringstream linestream(line);
        string seq;
        string val;
        linestream >> seq;
        linestream >> val;
        counts_1[seq] = stoi(val);
    }

    cin.clear();
    freopen(input_con.c_str(), "r", stdin);

    map<string, int> counts_2;
    while (getline(cin, line)) {
        stringstream linestream(line);
        string seq;
        string val;
        linestream >> seq;
        linestream >> val;
        counts_2[seq] = stoi(val);
    }

    int tot_1 = accumulate(
        counts_1.begin(), counts_1.end(), 0,
        [](const int prev, const auto& v) { return prev + v.second; });
    int tot_2 = accumulate(
        counts_2.begin(), counts_2.end(), 0,
        [](const int prev, const auto& v) { return prev + v.second; });
    map<string, double> e_1;
    for (auto& [k, v] : counts_1) {
        e_1[k] = static_cast<double>(v) / tot_1;
    }
    map<string, double> e_2;
    for (auto& [k, v] : counts_2) {
        e_2[k] = static_cast<double>(v) / tot_2;
    }

    // State 1
    double init_1 = 0.95;
    double a_11 = 0.95;
    double a_12 = 0.05;

    // State 2
    double init_2 = 0.05;
    double a_21 = 0.10;
    double a_22 = 0.90;

    string human;
    string dog;
    string mouse;
    for (auto& a : alns) {
        human += a.human;
        dog += a.dog;
        mouse += a.mouse;
    }

    vector<string> seqs;
    for (int i = 0; i < human.length(); i++) {
        string seq;
        seq += human.at(i);
        seq += dog.at(i);
        seq += mouse.at(i);
        seqs.push_back(seq);
    }

    vector<double> dp_1(seqs.size());
    vector<int> back_1(seqs.size());
    vector<double> dp_2(seqs.size());
    vector<int> back_2(seqs.size());

    // Initialization Case
    // "Transition" and Emission
    double curr_1 = log(init_1) + log(e_1.at(seqs[0]));
    double curr_2 = log(init_2) + log(e_2.at(seqs[0]));
    // DP
    int prev_1 = -1;
    int prev_2 = -1;
    dp_1[0] = curr_1;
    dp_2[0] = curr_2;

    for (int i = 1; i < seqs.size(); i++) {
        // Transition and Emission
        double h_11 = dp_1[i - 1] + log(a_11) + log(e_1.at(seqs[i]));
        double h_21 = dp_2[i - 1] + log(a_21) + log(e_1.at(seqs[i]));
        double h_12 = dp_1[i - 1] + log(a_12) + log(e_2.at(seqs[i]));
        double h_22 = dp_2[i - 1] + log(a_22) + log(e_2.at(seqs[i]));
        if (h_11 > h_21) {
            prev_1 = 1;
            curr_1 = h_11;
        } else {
            prev_1 = 2;
            curr_1 = h_21;
        }
        if (h_12 > h_22) {
            prev_2 = 1;
            curr_2 = h_12;
        } else {
            prev_2 = 2;
            curr_2 = h_22;
        }
        // DP
        dp_1[i] = curr_1;
        back_1[i - 1] = prev_1;
        dp_2[i] = curr_2;
        back_2[i - 1] = prev_2;
    }
    back_1[back_1.size() - 1] = 1;
    back_2[back_2.size() - 1] = 2;

    vector<int> viterbi;
    vector<int>* back;
    if (dp_1.back() > dp_2.back()) {
        back = &back_1;
    } else {
        back = &back_2;
    }
    for (int i = back->size() - 1; i >= 0; i--) {
        viterbi.push_back(back->at(i));
        if (back->at(i) == 1) {
            back = &back_1;
        } else {
            back = &back_2;
        }
    }
    reverse(viterbi.begin(), viterbi.end());

    vector<pair<int, int>> segList_1;
    vector<pair<int, int>> segList_2;

    int v_start = alns[0].start;
    int offset = 0;
    int prev = -1;

    map<int, int> s;
    for (int i = 0; i < viterbi.size(); i++) {
        s[viterbi[i]]++;
        if (viterbi[i] != prev && i != 0) {
            if (prev == 1) {
                segList_1.push_back({v_start + offset, v_start + i - 1});
            } else if (prev == 2) {
                segList_2.push_back({v_start + offset, v_start + i - 1});
            }
            offset = i;
        }
        prev = viterbi[i];
    }
    if (prev == 1) {
        segList_1.push_back({v_start + offset, v_start + viterbi.size() - 1});
    } else if (prev == 2) {
        segList_2.push_back({v_start + offset, v_start + viterbi.size() - 1});
    }

    sort(segList_2.begin(), segList_2.end(), [](const auto& a, const auto& b) {
        return (a.second - a.first) > (b.second - b.first);
    });

    cout << "State Histogram:" << '\n';
    for (auto& [k, v] : s) {
        cout << k << '=' << v << '\n';
    }
    cout << '\n';
    cout << "Segment Histogram:" << '\n';
    cout << "1=" << segList_1.size() << '\n';
    cout << "2=" << segList_2.size() << '\n';
    cout << '\n';
    cout << fixed << setprecision(5);
    cout << "Initial State Probabilities:" << '\n';
    cout << "1=" << init_1 << '\n';
    cout << "2=" << init_2 << '\n';
    cout << '\n';
    cout << "Transition Probabilities:" << '\n';
    cout << "1,1=" << a_11 << '\n';
    cout << "1,2=" << a_12 << '\n';
    cout << "2,1=" << a_21 << '\n';
    cout << "2,2=" << a_22 << '\n';
    cout << '\n';
    cout << "Emission Probabilities:" << '\n';
    for (auto& [k, v] : e_1) {
        cout << "1," << k << '=' << v << '\n';
    }
    for (auto& [k, v] : e_2) {
        cout << "2," << k << '=' << v << '\n';
    }
    cout << '\n';
    cout << "Longest Segment List:" << '\n';
    cout << '\n';
    for (int i = 0; i < 10; i++) {
        cout << segList_2[i].first << ' ' << segList_2[i].second << '\n';
    }
    cout << '\n';
    cout << "Annotations:" << '\n';
    cout << '\n';
    for (int i = 0; i < 5; i++) {
        cout << "Start: " << segList_2[i].first << '\n';
        cout << "End: " << segList_2[i].second << '\n';
        cout << '\n';
        cout << '\n';
    }
    cout << '\n';

    return 0;
}