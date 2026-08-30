#ifndef JONSS_FLUXES
#define JONSS_FLUXES

#include "options.hpp"
#include "state.hpp"
#include "primitives.hpp"

#include <cstdint>

#include <mfem/mfem.hpp>

namespace jonss
{

/**
 * @brief Function for computing the inviscid fluxes.
 * 
 * @tparam TModel    Fluid model.
 * @tparam TDim      Spatial dimension. Templated to allow compiler to unroll
 *                   loops.
 * 
 * @param state      State struct.
 * @param prim       Primitives struct.
 * 
 * @return fluxes    Computed fluxes.
 */
template<FluidModelOption TModel, int TDim>
MFEM_HOST_DEVICE inline
void ComputeInviscidFluxes(const State<TModel,TDim> &state,
                           const Primitives<TModel,TDim> &prim,
                           State<TModel,TDim> (&fluxes)[TDim]) 
{
   using enum FluidModelOption;

   if constexpr (TModel == AirCPG)
   {
      // See "I do like CFD" for details.

      for (int di = 0; di < TDim; di++)
      {
         // Set density fluxes
         fluxes[di].rho = state.rhoV[di];

         // Set momentum fluxes
         for (int dj = 0; dj < TDim; dj++)
         {
            fluxes[di].rhoV[dj] = state.rhoV[di]*prim.vel[dj];
         }
         fluxes[di].rhoV[di] += prim.p;

         // Set energy fluxes
         fluxes[di].rhoE = state.rhoV[di]*prim.H;
      }
   }
   else
   {
      static_assert(TDim != TDim, "Unimplemented model.");
   }
}

/**
 * @brief Compute the numerical fluxes.
 * 
 * @tparam TModel    Fluid model.
 * @tparam TFlux     Numerical flux.
 * @tparam TDim      Spatial dimension. Templated to allow compiler to unroll
 *                   loops.
 * @tparam TStab     If true, include stabilization (for face fluxes).
 *                   If false, do not include stabilization (for volume fluxes).
 * 
 * @param state1     First state.
 * @param state2     Second state.
 * 
 * @return fluxes    Computed numerical fluxes.
 */
template<FluidModelOption TModel, NumericalFluxOption TFlux, int TDim, 
         bool TStab>
MFEM_HOST_DEVICE inline void ComputeNumericalFluxes(
      const State<TModel,TDim> &state1,
      const State<TModel,TDim> &state2,
      State<TModel,TDim> (&fluxes)[TDim])
{
   using enum NumericalFluxOption;

   if constexpr (TFlux == LocalLaxFriedrichs)
   { 
      State<TModel,TDim> F_1[TDim];
      State<TModel,TDim> F_2[TDim];

      //  TODO: Will need to think about viscous soon. Hurray.
      ComputeInviscidFluxes<TModel,TDim>(state1, F_1);
      ComputeInviscidFluxes<TModel,TDim>(state2, F_2);

      if constexpr (TModel == FluidModelOption::AirCPG)
      {
         for (int di = 0; di < TDim; di++)
         {
            fluxes[di].rho = 0.5*(F_1[di].rho + F_2[di].rho);
            for (int dj = 0; dj < TDim; dj++)
            {
               fluxes[di].rhoV[dj] = 0.5*(F_1[di].rhoV[dj] + F_2[di].rhoV[dj]);
            }
            fluxes[di].rhoE = 0.5*(F_1[di].rhoE + F_2[di].rhoE);

            if constexpr (TStab)
            {
               // mfem::real_t lambda_max;
               // TODO
               // fluxes[d][c] -= 0.5*lambda_max*(state2[c] - state1[c]);
            }
         }
      }
      else
      {
         static_assert(TDim != TDim, "Unimplemented numerical flux.");
      }

   }
   else if constexpr (TFlux == Chandrashekar)
   {

   }
   else
   {
      static_assert(TDim != TDim, "Unimplemented numerical flux.");
   }
}

} // namespace jonss

#endif  // JONSS_FLUXES
