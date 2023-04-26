#include <stg/kriging.hpp>

using namespace stg::kriging;


int main() {
  DataLoader data_loader{"test_resources/"};
  KrigingAnalysis<double> kriging_test{data_loader};

  const auto& calc_cov_test = kriging_test.calculate_covariations();
  const auto& cov_data = kriging_test.covariations_from_data();
  return EXIT_SUCCESS;
}
