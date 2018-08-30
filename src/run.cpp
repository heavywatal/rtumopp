// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <tumopp/simulation.hpp>

//' Run C++ simulation
//' @param args command line arguments as a string vector
//' @rdname rcpprun
// [[Rcpp::export]]
Rcpp::CharacterVector cpp_tumopp(const std::vector<std::string>& args) {
    try {
        tumopp::Simulation simulation(args);
        simulation.run();
        return Rcpp::CharacterVector::create(
            Rcpp::Named("config", simulation.config_string()),
            Rcpp::Named("population", simulation.history()),
            Rcpp::Named("snapshots", simulation.snapshots()),
            Rcpp::Named("drivers", simulation.drivers()),
            Rcpp::Named("distances", simulation.pairwise_distance()),
            Rcpp::Named("ms", simulation.ms())
        );
    } catch (const std::runtime_error& e) {
        Rcpp::Rcerr << e.what() << "\n";
    }
    return {};
}

//' @param nsam number of samples for ms-like output
//' @rdname rcpprun
// [[Rcpp::export]]
Rcpp::CharacterVector cpp_tumopp_ms(unsigned int nsam, const std::vector<std::string>& args) {
    try {
        std::vector<std::string> rargs{{std::to_string(nsam), "1"}};
        rargs.insert(rargs.end(), args.begin(), args.end());
        tumopp::Simulation simulation(rargs);
        simulation.run();
        return simulation.ms();
    } catch (const std::runtime_error& e) {
        Rcpp::Rcerr << e.what() << "\n";
    }
    return {};
}
