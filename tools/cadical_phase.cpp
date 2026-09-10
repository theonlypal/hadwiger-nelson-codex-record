// Thin CaDiCaL API driver. Phase hints change search preferences only.
// The sole clause input is the original DIMACS file; no assumptions,
// symmetry-breaking clauses, or constraints are added by this driver.
#include "cadical.hpp"
#include <csignal>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

static volatile std::sig_atomic_t stopped = 0;
static void stop(int) { stopped = 1; }
struct Stopper : CaDiCaL::Terminator {
  bool terminate() override { return stopped != 0; }
};

int main(int argc, char **argv) {
  std::string solution, phases, configuration;
  std::vector<std::string> positional;
  int seed = 0;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--version") { std::cout << CaDiCaL::Solver::version() << "\n"; return 0; }
    if (arg == "--no-colors") continue;
    if (arg.rfind("--seed=", 0) == 0) { seed = std::stoi(arg.substr(7)); continue; }
    if (arg == "--sat" || arg == "--unsat") { configuration = arg.substr(2); continue; }
    if ((arg == "-w" || arg == "--phase-file") && i+1 < argc) {
      (arg == "-w" ? solution : phases) = argv[++i]; continue;
    }
    if (arg.size() && arg[0] == '-') { std::cerr << "Unsupported option\n"; return 1; }
    positional.push_back(arg);
  }
  if (positional.size() != 2 || solution.empty()) {
    std::cerr << "Expected -w SOLUTION CNF PROOF [--phase-file FILE]\n"; return 1;
  }
  CaDiCaL::Solver solver;
  solver.set("seed", seed);
  solver.set("quiet", 0);
  solver.set("binary", 1);
  if (!configuration.empty() && !solver.configure(configuration.c_str())) return 1;
  if (!solver.trace_proof(positional[1].c_str())) return 1;
  int variables = 0;
  if (const char *error = solver.read_dimacs(positional[0].c_str(), variables)) {
    std::cerr << error << "\n"; return 1;
  }
  int phase_count = 0;
  if (!phases.empty()) {
    std::ifstream input(phases);
    if (!input) return 1;
    std::vector<bool> seen(variables+1, false);
    int literal = 0;
    bool terminated = false;
    while (input >> literal) {
      if (terminated) { std::cerr << "Trailing phase data\n"; return 1; }
      if (!literal) { terminated = true; continue; }
      if (literal < -variables || literal > variables || seen[std::abs(literal)]) {
        std::cerr << "Invalid or duplicate phase variable\n"; return 1;
      }
      seen[std::abs(literal)] = true;
      solver.phase(literal);
      ++phase_count;
    }
    if (!terminated || !input.eof()) { std::cerr << "Malformed phase file\n"; return 1; }
  }
  std::cout << "c phase hints " << phase_count << " (preferences only)\n";
  Stopper stopper;
  solver.connect_terminator(&stopper);
  std::signal(SIGTERM, stop); std::signal(SIGINT, stop);
  int result = solver.solve();
  std::ofstream output(solution);
  if (!output) return 1;
  if (result == 10) {
    output << "s SATISFIABLE\nv";
    for (int v = 1; v <= variables; ++v) output << ' ' << (solver.val(v) > 0 ? v : -v);
    output << " 0\n";
  } else output << (result == 20 ? "s UNSATISFIABLE\n" : "s UNKNOWN\n");
  solver.statistics();
  solver.close_proof_trace();
  return result;
}
