#include <bits/stdc++.h>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    freopen("output.txt", "w", stdout);
    // string input = "input/example.input.txt";
    // string sinput = "input/example.scoring.scheme.txt";
    string input = "input/chm13.chr16.txt";
    string sinput = "input/scoring.scheme.txt";

    freopen(input.c_str(), "r", stdin);

    vector<int> reads;
    string line;
    while (getline(cin, line, '\n')) {
        reads.push_back(stoi(line.substr(line.rfind('\t') + 1)));
    }

    cin.clear();
    freopen(sinput.c_str(), "r", stdin);

    map<int, double> scoreMap;
    int thresh = 0;
    double threshMap = 0;
    double s = 0;
    double d = 0;
    while (getline(cin, line, '\n')) {
        stringstream linestream(line);
        string col;
        linestream >> col;
        if (col.find(">=") == 0) {
            thresh = stoi(col.substr(2));
            linestream >> col;
            threshMap = stod(col);
        } else if (col.front() == 'S') {
            linestream >> col;
            s = stod(col);
        } else if (col.front() == 'D') {
            linestream >> col;
            d = stod(col);
        } else {
            linestream >> scoreMap[stoi(col)];
        }
    }

    vector<tuple<int, int, double>> segList;
    double total = 0;
    double max = 0;
    int start = 0;
    int end = 0;
    int rs = reads.size();
    for (int i = 0; i < rs; i++) {
        double val = 0;
        int read = reads[i];
        if (read >= thresh) {
            val = threshMap;
        } else {
            val = scoreMap.at(read);
        }
        total += val;

        if (total >= max) {
            max = total;
            end = i;
        }
        if ((total <= 0) || (total <= max + d) || (i == rs - 1)) {
            if (max >= s) {
                segList.push_back({start, end, max});
            }
            total = 0;
            max = 0;
            start = i + 1;
            end = i + 1;
        }
    }

    vector<tuple<int, int>> bgsegList;
    int bgstart = 0;
    int maxEnd = rs - 1;
    for (auto& [start, end, max] : segList) {
        if (start == 0) {
            bgstart = end + 1;
            continue;
        }

        bgsegList.push_back({bgstart, start - 1});
        bgstart = end + 1;
    }
    if (bgstart <= maxEnd) {
        bgsegList.push_back({bgstart, maxEnd});
    }

    map<int, int> hist;
    for (auto& [start, end, max] : segList) {
        for (int i = start; i <= end; i++) {
            int read = reads[i];
            if (read >= thresh) {
                read = thresh;
            }
            hist[read]++;
        }
    }

    map<int, int> bgHist;
    for (auto& [start, end] : bgsegList) {
        for (int i = start; i <= end; i++) {
            int read = reads[i];
            if (read >= thresh) {
                read = thresh;
            }
            bgHist[read]++;
        }
    }

    cout << fixed << setprecision(2);
    cout << "Segment Histogram:" << '\n';
    cout << "Non-Elevated CN Segments=" << bgsegList.size() << '\n';
    cout << "Elevated CN Segments=" << segList.size() << '\n';
    cout << '\n';
    cout << "Segment List:" << '\n';
    for (auto& [start, end, max] : segList) {
        cout << start + 1 << ' ' << end + 1 << ' ' << max << '\n';
    }
    cout << '\n';
    cout << "Annotations:" << '\n';
    cout << '\n';
    int as = min(3, static_cast<int>(segList.size()));
    sort(segList.begin(), segList.end(),
         [](const auto& sega, const auto& segb) {
             return get<2>(sega) > get<2>(segb);
         });
    for (int i = 0; i < as; i++) {
        auto& [start, end, max] = segList[i];
        cout << "Start: " << start + 1 << '\n';
        cout << "End: " << end + 1 << '\n';
        cout << "Description: Something interesting (e.g., \"Overlaps with "
                "exon5 of the protein coding gene cMyc\")"
             << '\n';
        cout << '\n';
    }
    cout << "Read start histogram for non-elevated copy-number segments:"
         << '\n';
    for (auto& [read, count] : bgHist) {
        if (read == thresh) {
            cout << ">=";
        }
        cout << read << '=' << count << '\n';
    }
    cout << '\n';
    cout << "Read start histogram for elevated copy-number segments:" << '\n';
    for (auto& [read, count] : hist) {
        if (read == thresh) {
            cout << ">=";
        }
        cout << read << '=' << count << '\n';
    }
    cout << '\n';

    return 0;
}