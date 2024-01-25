#include <bits/stdc++.h>

using namespace std;

struct Edge {
    int u;
    int v;
    double w;
    string label;
};

struct Graph {
    vector<vector<Edge>> adj;

    void addVertex() { adj.push_back(vector<Edge>()); }

    void addEdge(int u, int v, double w, string l) {
        adj[u].push_back(Edge{u, v, w, l});
    }
};

vector<int> ideg(const Graph &g) {
    vector<int> indg(g.adj.size(), 0);

    for (vector<Edge> ve : g.adj) {
        for (Edge e : ve) {
            indg[e.v]++;
        }
    }

    return indg;
}

vector<int> topo(const Graph &g) {
    vector<int> sorted;
    vector<bool> visited(g.adj.size(), false);
    deque<int> dq;

    for (int i = 0; i < visited.size(); i++) {
        if (!visited[i]) {
            dq.push_back(i);
            visited[i] = true;
        }

        // DFS
        while (!dq.empty()) {
            int curr = dq.back();
            int numAdds = 0;
            for (Edge e : g.adj[curr]) {
                if (!visited[e.v]) {
                    numAdds++;
                    dq.push_back(e.v);
                    visited[e.v] = true;
                }
            }
            if (numAdds == 0) {
                sorted.push_back(curr);
                dq.pop_back();
            }
        }
    }

    reverse(sorted.begin(), sorted.end());
    return sorted;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    freopen("output.txt", "a", stdout);
    // string input = "input/example_graph.txt";
    // string input = "input/constraint_example_graph.txt";
    // string input = "input/graph.txt";
    // string input = "input/constraint_graph.txt";
    string input = "input/seqgraph.txt";
    string inputName = input.substr(input.find('/') + 1);
    freopen(input.c_str(), "r", stdin);

    // Create Graph
    Graph g = Graph();

    map<string, int> nti;
    vector<string> itn;
    int nameIdx = 0;

    // Constraint Case
    int startIdx = -1;
    int endIdx = -1;

    string line;
    while (getline(cin, line)) {
        if (line.front() == 'V') {
            if (line.find("START") != -1) {
                startIdx = nameIdx;
                line = line.substr(0, line.find("START") - 1);
            } else if (line.find("END") != -1) {
                endIdx = nameIdx;
                line = line.substr(0, line.find("END") - 1);
            }
            string vName = line.substr(line.find(' ') + 1);
            nti[vName] = nameIdx++;
            itn.push_back(vName);
            g.addVertex();
        } else if (line.front() == 'E') {
            line = line.substr(line.find(' ') + 1);
            string l = line.substr(0, line.find(' '));
            line = line.substr(line.find(' ') + 1);
            string u = line.substr(0, line.find(' '));
            line = line.substr(line.find(' ') + 1);
            string v = line.substr(0, line.find(' '));
            line = line.substr(line.find(' ') + 1);
            string w = line;
            g.addEdge(nti[u], nti[v], stod(w), l);
        }
    }

    // Topo Sort
    vector<int> t = topo(g);

    // In-Degree Vector
    // vector<int> indg = ideg(g);

    // Modified DP for SSSP on DAG

    // Constraint Case
    double dpinit = (startIdx != -1 && endIdx != -1)
                        ? -1 * numeric_limits<double>::infinity()
                        : 0;

    vector<double> dp(g.adj.size(), dpinit);
    vector<int> dpprev(g.adj.size());
    iota(dpprev.begin(), dpprev.end(), 0);

    double maxScore = 0;
    int maxEnd = 0;

    // Constraint Case
    if (startIdx != -1 && endIdx != -1) {
        dp[startIdx] = 0;
    }

    for (int i : t) {
        // if (indg[i] == 0) {
        //     dp[i] = 0;
        // }

        // Constraint Case
        if (dp[i] == -1 * numeric_limits<double>::infinity()) {
            continue;
        }

        for (Edge e : g.adj[i]) {
            double vw = dp[i] + e.w;
            // Only Accounts For First Instance Of Path Weight Score -> More
            // Than One Can Exist - Does Not Account For Multiple Instances
            if (vw > dp[e.v]) {
                dp[e.v] = vw;
                dpprev[e.v] = i;
            }
            if (vw > maxScore) {
                maxScore = vw;
                maxEnd = e.v;
            }
        }
    }

    // Constraint Case
    if (startIdx != -1 && endIdx != -1) {
        maxEnd = endIdx;
        maxScore = dp[endIdx];
    }

    // Generate Path
    deque<int> pathgen{maxEnd};
    while (dpprev[pathgen.front()] != pathgen.front()) {
        pathgen.push_front(dpprev[pathgen.front()]);
    }
    int maxStart = pathgen.front();

    // Generate Path Label
    string pathlabel;
    while (pathgen.size() != 1) {
        int u = pathgen.front();
        pathgen.pop_front();
        for (Edge e : g.adj[u]) {
            if (e.v == pathgen.front()) {
                pathlabel += e.label;
                break;
            }
        }
    }

    // Print
    cout << "Input: " << inputName << '\n';
    cout << "Score: " << maxScore << '\n';
    cout << "Begin: " << itn[maxStart] << '\n';
    cout << "End: " << itn[maxEnd] << '\n';
    cout << "Path: " << pathlabel << '\n';
    cout << '\n';

    return 0;
}