// Bounded TabuCol candidate generator. Failure never claims uncolorability.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  try {
    if (argc != 9) throw std::runtime_error("usage: graph init best stats seed seconds max_moves sample");
    std::ifstream input(argv[1]), initial(argv[2]);
    int n, m, k;
    if (!(input >> n >> m >> k) || n <= 0 || m < 0 || k < 2)
      throw std::runtime_error("invalid graph header");
    const uint64_t seed = std::stoull(argv[5]);
    const double seconds = std::stod(argv[6]);
    const uint64_t max_moves = std::stoull(argv[7]);
    const int sample_size = std::stoi(argv[8]);
    if (seconds <= 0 || max_moves == 0 || sample_size <= 0) throw std::runtime_error("invalid limit");
    std::mt19937_64 rng(seed);
    std::vector<std::array<int,2>> edges(m);
    std::vector<int> offsets(n+1, 0);
    for (auto &e : edges) {
      if (!(input >> e[0] >> e[1]) || e[0]<0 || e[0]>=n || e[1]<0 || e[1]>=n || e[0]==e[1])
        throw std::runtime_error("invalid edge");
      ++offsets[e[0]+1]; ++offsets[e[1]+1];
    }
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
    std::vector<int> adjacency(2*m), cursor(offsets);
    for (const auto &e : edges) { adjacency[cursor[e[0]]++]=e[1]; adjacency[cursor[e[1]]++]=e[0]; }
    std::vector<int> color(n), count(size_t(n)*k,0), uncolored;
    for (int v=0;v<n;++v) {
      if (!(initial >> color[v]) || color[v]<-1 || color[v]>=k) throw std::runtime_error("invalid initial color");
      if (color[v]<0) uncolored.push_back(v);
    }
    for (int v=0;v<n;++v) if (color[v]>=0)
      for (int j=offsets[v];j<offsets[v+1];++j) ++count[size_t(adjacency[j])*k+color[v]];
    std::stable_sort(uncolored.begin(),uncolored.end(),[&](int u,int v) {
      return offsets[u+1]-offsets[u] > offsets[v+1]-offsets[v];
    });
    for (int v : uncolored) {
      int chosen=0, value=std::numeric_limits<int>::max(), ties=0;
      for (int c=0;c<k;++c) {
        const int x=count[size_t(v)*k+c];
        if (x<value) { chosen=c;value=x;ties=1; }
        else if (x==value && rng()%uint64_t(++ties)==0) chosen=c;
      }
      color[v]=chosen;
      for (int j=offsets[v];j<offsets[v+1];++j) ++count[size_t(adjacency[j])*k+chosen];
    }
    int conflicts=0;
    for (const auto &e : edges) conflicts += color[e[0]]==color[e[1]];
    const int initial_conflicts=conflicts;
    int best=conflicts;
    std::vector<int> best_color(color), active, position(n,-1);
    auto refresh = [&](int v) {
      const bool bad=count[size_t(v)*k+color[v]]>0;
      if (bad && position[v]<0) { position[v]=int(active.size());active.push_back(v); }
      else if (!bad && position[v]>=0) {
        const int p=position[v], last=active.back(); active[p]=last;position[last]=p;
        active.pop_back();position[v]=-1;
      }
    };
    for (int v=0;v<n;++v) refresh(v);
    std::vector<uint64_t> tabu(size_t(n)*k,0);
    const auto start=std::chrono::steady_clock::now();
    auto elapsed = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); };
    auto save_model = [&]() {
      std::ofstream out(argv[3]);
      if (!out) throw std::runtime_error("cannot write best model");
      for (int c : best_color) out << c << '\n';
    };
    save_model();
    std::cout << "{\"stage\":\"initialized\",\"initial_conflicts\":" << conflicts
              << ",\"conflicted_vertices\":" << active.size() << ",\"free_vertices\":" << n << "}" << std::endl;
    uint64_t moves=0, best_move=0, forced_tabu_moves=0;
    double checkpoint=0;
    while (conflicts>0 && moves<max_moves) {
      if ((moves&255)==0 && elapsed()>=seconds) break;
      int move_vertex=-1, move_color=-1, move_delta=std::numeric_limits<int>::max(), ties=0;
      const int samples=std::min(sample_size,int(active.size()));
      for (int sample=0;sample<samples;++sample) {
        const int v=active[rng()%active.size()], previous=count[size_t(v)*k+color[v]];
        for (int c=0;c<k;++c) if (c!=color[v]) {
          const int delta=count[size_t(v)*k+c]-previous;
          if (tabu[size_t(v)*k+c]>moves && conflicts+delta>=best) continue;
          if (delta<move_delta) { move_vertex=v;move_color=c;move_delta=delta;ties=1; }
          else if (delta==move_delta && rng()%uint64_t(++ties)==0) { move_vertex=v;move_color=c; }
        }
      }
      if (move_vertex<0) {
        ++forced_tabu_moves;
        move_vertex=active[rng()%active.size()]; move_color=(color[move_vertex]+1+int(rng()%uint64_t(k-1)))%k;
        move_delta=count[size_t(move_vertex)*k+move_color]-count[size_t(move_vertex)*k+color[move_vertex]];
      }
      const int v=move_vertex, previous=color[v];
      const uint64_t tenure=5+rng()%10+uint64_t(active.size()*0.6);
      tabu[size_t(v)*k+previous]=moves+tenure;
      color[v]=move_color; conflicts+=move_delta; ++moves;
      for (int j=offsets[v];j<offsets[v+1];++j) {
        const int w=adjacency[j]; --count[size_t(w)*k+previous];++count[size_t(w)*k+move_color];refresh(w);
      }
      refresh(v);
      if (conflicts<best) { best=conflicts;best_color=color;best_move=moves; }
      if ((moves&4095)==0 && elapsed()-checkpoint>=5) {
        checkpoint=elapsed();save_model();
        std::cout << "{\"stage\":\"progress\",\"moves\":" << moves << ",\"conflicts\":" << conflicts
                  << ",\"best_conflicts\":" << best << ",\"conflicted_vertices\":" << active.size()
                  << ",\"elapsed_seconds\":" << checkpoint << "}" << std::endl;
      }
    }
    save_model();
    int checked=0;for (const auto &e : edges) checked+=best_color[e[0]]==best_color[e[1]];
    if (checked!=best) throw std::runtime_error("incremental conflict count mismatch");
    std::ofstream stats(argv[4]);
    stats << "{\"status\":\"" << (best==0?"PROPER_MODEL_CANDIDATE":"UNKNOWN_IMPROPER_BEST")
          << "\",\"seed\":" << seed << ",\"vertices\":" << n << ",\"edges\":" << m << ",\"k\":" << k
          << ",\"moves\":" << moves << ",\"best_move\":" << best_move
          << ",\"initial_conflicts\":" << initial_conflicts << ",\"best_conflicts\":" << best
          << ",\"final_conflicts\":" << conflicts << ",\"forced_tabu_moves\":" << forced_tabu_moves
          << ",\"elapsed_seconds\":" << elapsed() << ",\"all_vertices_free\":true,\"goal_satisfied\":false}\n";
    return best==0?10:0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';return 1;
  }
}
