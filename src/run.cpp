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
            Rcpp::Named("config", simulation.config()),
            Rcpp::Named("population", simulation.history()),
            Rcpp::Named("snapshots", simulation.snapshots()),
            Rcpp::Named("drivers", simulation.drivers())
        );
    } catch (const std::runtime_error& e) {
        Rcpp::Rcerr << e.what() << "\n";
    }
    return {};
}
