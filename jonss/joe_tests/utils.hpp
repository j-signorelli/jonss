#include <mfem/mfem.hpp>

namespace jonss
{
namespace joe_tests
{

/**
 * @brief Generate a 2D rectangular grid with x-length set to ensure
 * exactly \p num_waves in the domain for a single acoustic wave in
 * a freestream, both in the x-direction.
 * 
 * @warning The grid order is kept at 1. For isoparametric elements,
 * set using mfem::Mesh::SetCurvature().
 * 
 * @param k                  Acoustic wave wavenumber.
 * @param num_waves          Number of acoustic waves to include.
 * @param n_elems_per_wave   Number of elements to include per
 *                           wavelength.
 * @param n_elems_y          Number of elements to include in y-direction.
 * @param grid_ar            Aspect ratio (length_x/length_y) of grid.
 *                           This sets the y-spacing.
 */
mfem::Mesh Create2DWaveDomain(const double k,
                              const int num_waves=1,
                              const int n_elems_per_wave=60,
                              const int n_elems_y=10,
                              const double grid_ar=6);

} // namespace joe_tests

} // namespace jonss
