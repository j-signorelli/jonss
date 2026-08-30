#ifndef JONSS_STATE
#define JONSS_STATE

#include "options.hpp"

namespace jonss
{

/**
 * @brief Struct to hold model-specific conservative variables.
 * 
 * @details This should be explicitly specialized for each 
 * \ref FluidModelOption.
 * 
 * @tparam TModel    Fluid model.
 * @tparam TDim      Spatial dimension.
 */
template<FluidModelOption TModel, int TDim>
struct State
{
   MFEM_HOST_DEVICE State() = delete;
};

/**
 * @brief State specialization for FluidModelOption::AirCPG.
 */
template<int TDim>
struct State<FluidModelOption::AirCPG, TDim>
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
