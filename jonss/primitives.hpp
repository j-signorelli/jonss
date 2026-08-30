#ifndef JONSS_PRIMITIVES
#define JONSS_PRIMITIVES

#include "options.hpp"
#include "constants.hpp"

namespace jonss
{

/**
 * @brief Struct to hold relevant, model-specific primitive and auxiliary
 * variables.
 * 
 * @details For variables, dependent on the state, that are to be computed 
 * once to avoid expensive re-computations on-the-fly and support storage 
 * in smem.  This should be explicitly specialized for each 
 * \ref FluidOption.
 * 
 * @tparam TFluid    Fluid model.
 * @tparam TDim      Spatial dimension.
 */
template<FluidOption TFluid, int TDim>
struct Primitives
{
   MFEM_HOST_DEVICE Primitives() = delete;
};

/**
 * @brief Primitives specialization for FluidOption::AirCPG.
 * 
 * @tparam TDim      Spatial dimension.
 */
template<int TDim>
struct Primitives<FluidOption::EulerCPG, TDim>
{
   /// Velocity.
   mfem::real_t vel[TDim];

   /// Velocity magnitude, squared.
   mfem::real_t vel_sq;

   /// Pressure.
   mfem::real_t p;

   /// Specific total enthalpy.
   mfem::real_t H;
};

/**
 * @brief Function for initializing a Primitives struct.
 * 
 * @tparam TFluid    Fluid model.
 * @tparam TDim      Spatial dimension. Templated to allow compiler to unroll
 *                   loops.
 * 
 * @param[in] state     State struct.
 * @param[out] prim     Primitives struct.
 */
template<FluidOption TFluid, int TDim>
MFEM_HOST_DEVICE inline
void ComputePrimitives(const FluidConstants<TFluid> &constants,
                       const State<TFluid,TDim> &state,
                       Primitives<TFluid,TDim> &prim)
{
   using enum FluidOption;

   if constexpr (TFluid == EulerCPG)
   {
      // See "I do like CFD" for details.

      prim.vel_sq = 0.0;
      for (int d = 0; d < TDim; d++)
      {
         prim.vel[d] = state.rhoV[d]/state.rho;
         prim.vel_sq += prim.vel[d]*prim.vel[d];
      }
      prim.p = (constants.gamma-1.0)*(state.rhoE-0.5*state.rho*prim.vel_sq);
      prim.H = (state.rhoE + prim.p)/state.rho;
   }
   else
   {
      static_assert(TDim != TDim, "Unimplemented model.");
   }
}

} // namespace jonss

#endif // JONSS_PRIMITIVES
