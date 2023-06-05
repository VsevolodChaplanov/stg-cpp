#include <stg/kriging.hpp>
#include <stg/spectral_method.hpp>

using namespace stg::kriging;
using namespace stg::spectral;
using namespace std::string_literals;


int main() {
    const std::string work_dir = "./kriging_result/";
    const std::string cov_filename = "r_11.vtk";
    const std::string fert_filename = "phi_11.vtk";
    IntegrateFert1D<double> fert_integrator{stg::spectral::DataLoader{work_dir}, cov_filename, fert_filename};

    fert_integrator.integrate_covariations();
    fert_integrator.save_calculated_fert(work_dir + "fert_11.vtk"s);
    return EXIT_SUCCESS;
}
