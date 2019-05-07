#include <Rcpp.h>
#include <tumopp/simulation.hpp>

// [[Rcpp::export]]
Rcpp::CharacterVector cpp_tumopp(const std::vector<std::string>& args) {
    try {
        std::streambuf* obuf = tumopp::std_cout_rdbuf(Rcpp::Rcout.rdbuf());
        std::streambuf* ebuf = tumopp::std_cerr_rdbuf(Rcpp::Rcerr.rdbuf());
        tumopp::Simulation simulation(args);
        simulation.run();
        tumopp::std_cout_rdbuf(obuf);
        tumopp::std_cout_rdbuf(ebuf);
        return Rcpp::CharacterVector::create(
            Rcpp::Named("config", simulation.config()),
            Rcpp::Named("population", simulation.history()),
            Rcpp::Named("snapshots", simulation.snapshots()),
            Rcpp::Named("drivers", simulation.drivers()),
            Rcpp::Named("benchmark", simulation.benchmark())
        );
    } catch (const std::runtime_error& e) {
        Rcpp::Rcerr << e.what() << "\n";
    }
    return {};
}
