#include <bits/stdc++.h>

using namespace std;

// struct Edge {
//     array<int, 3> src;
//     array<int, 3> dest;
//     string label;
//     int weight;
// };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    freopen("output.txt", "a", stdout);
    string input = "input/input.fa";
    string inputName = input.substr(input.find('/') + 1);
    string sinput = "input/scoringscheme.txt";
    string sinputName = sinput.substr(sinput.find('/') + 1);

    freopen(input.c_str(), "r", stdin);

    vector<string> sequences;

    string line;
    while (getline(cin, line, '\n')) {
        if (line.front() != '>' && !line.empty()) {
            sequences.push_back(line);
        }
    }

    if (sequences.size() != 3) {
        throw runtime_error("Only Three Sequences Are Allowed.");
    }

    cin.clear();
    freopen(sinput.c_str(), "r", stdin);

    vector<char> headerMap;
    string header;
    getline(cin, header, '\n');
    for (char c : header) {
        if (c != ' ') {
            headerMap.push_back(c);
        }
    }

    map<char, map<char, int>> sm;
    while (getline(cin, line, '\n')) {
        if (line.empty()) {
            break;
        }
        stringstream linestream(line);
        string row;
        linestream >> row;
        int colIdx = 0;
        string val;
        while (linestream >> val) {
            sm[row.front()][headerMap[colIdx++]] = stoi(val);
        }
    }
    getline(cin, line, '\n');
    string gp = "Gap penalty: ";
    int gpVal = stoi(line.substr(line.find(gp) + gp.length()));
    for (auto& [k, v] : sm) {
        v['-'] = gpVal;
    }
    for (char c : headerMap) {
        sm['-'][c] = gpVal;
    }
    sm['-']['-'] = 0;

    // Three Sequence Limit Logic Below
    string sqName = "seqgraph.txt";
    ofstream sqFile("input/" + sqName);

    int zeroL = sequences[0].length();
    int oneL = sequences[1].length();
    int twoL = sequences[2].length();

    vector<array<int, 3>> vertices;
    vertices.reserve((zeroL + 1) * (oneL + 1) * (twoL + 1));
    for (int i = 0; i < zeroL + 1; i++) {
        for (int j = 0; j < oneL + 1; j++) {
            for (int k = 0; k < twoL + 1; k++) {
                vertices.push_back({i, j, k});
                sqFile << "V (" << i << ',' << j << ',' << k << ')' << '\n';
            }
        }
    }

    // vector<Edge> edges;
    // edges.reserve((7 * zeroL * oneL * twoL) +
    //               ((3 * zeroL * oneL) + zeroL + oneL) + (3 * zeroL * twoL) +
    //               ((3 * oneL * twoL) + twoL));
    map<string, int> labelCounts;
    map<string, int> labelWeights;
    for (auto& src : vertices) {
        auto& [x, y, z] = src;
        array<array<int, 3>, 7> d;
        d[0] = {x + 1, y, z};
        d[1] = {x, y + 1, z};
        d[2] = {x, y, z + 1};
        d[3] = {x + 1, y + 1, z};
        d[4] = {x, y + 1, z + 1};
        d[5] = {x + 1, y, z + 1};
        d[6] = {x + 1, y + 1, z + 1};
        for (auto& dest : d) {
            auto& [dx, dy, dz] = dest;
            if (dx > zeroL || dy > oneL || dz > twoL) {
                continue;
            }
            string label;
            for (int i = 0; i < dest.size(); i++) {
                if (dest[i] - src[i] == 0) {
                    label += '-';
                } else if (dest[i] - src[i] == 1) {
                    label += sequences[i][src[i]];
                } else {
                    assert(dest[i] - src[i] == 1 || dest[i] - src[i] == 0);
                }
            }
            int weight = 0;
            if (labelCounts[label] > 0) {
                weight = labelWeights[label];
            } else {
                weight += sm[label[0]][label[1]];
                weight += sm[label[0]][label[2]];
                weight += sm[label[1]][label[2]];
                labelWeights[label] = weight;
            }
            labelCounts[label]++;
            sqFile << "E " << label << " (" << x << ',' << y << ',' << z
                   << ") (" << dx << ',' << dy << ',' << dz << ") " << weight
                   << '\n';
            // edges.push_back(Edge{src, dest, label, weight});
        }
    }
    sqFile.close();

    cout << "Edge weights: " << '\n';
    for (auto& [k, v] : labelWeights) {
        cout << k << " = " << v << '\n';
    }
    cout << '\n';
    cout << "Edge counts: " << '\n';
    for (auto& [k, v] : labelCounts) {
        cout << k << " = " << v << '\n';
    }
    cout << '\n';

    return 0;
}