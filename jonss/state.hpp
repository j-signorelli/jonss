#ifndef JONSS_STATE
#define JONSS_STATE

#include "options.hpp"

#include <mfem/mfem.hpp>

namespace jonss
{

/**
 * @brief Struct to hold model-specific conservative variables.
 * 
 * @details Default template contains the standard conserved variables. 
 * This may be explicitly specialized to contain more for a specific 
 * \ref FluidOption, such as for species transport.
 * 
 * @tparam TFluid    Fluid model.
 * @tparam TDim      Spatial dimension.
 */
template<FluidOption TFluid, int TDim>
struct State
{
   /// Density.
   mfem::real_t rho;

   /// Momentum.
   mfem::real_t rhoV[TDim];
   
   /// Energy.
   mfem::real_t rhoE;
};

} // namespace jonss

#endif // JONSS_STATE
