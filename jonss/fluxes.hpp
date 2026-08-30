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
 * @brief Function for computing the exact, analytical fluxes.
 * 
 * @tparam TFluid    Fluid model.
 * @tparam TVisc     Viscosity model.
 * @tparam TDim      Spatial dimension. Templated to allow compiler to unroll
 *                   loops.
 * 
 * @param[in] state      State struct.
 * @param[in] prim       Primitives struct.
 * 
 * @param[out] fluxes    Computed fluxes.
 */
template<FluidOption TFluid, ViscosityOption TVisc, int TDim>
MFEM_HOST_DEVICE inline
void ComputeExactFluxes(const State<TFluid,TDim> &state,
                        const Primitives<TFluid,TVisc,TDim> &prim,
                        State<TFluid,TDim> (&fluxes)[TDim]) 
{
   using enum FluidOption;

   if constexpr (TFluid == CPG)
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
 * @tparam TFluid    Fluid model.
 * @tparam TFlux     Numerical flux.
 * @tparam TDim      Spatial dimension. Templated to allow compiler to unroll
 *                   loops.
 * @tparam TStab     If true, include stabilization (for face fluxes).
 *                   If false, do not include stabilization (for volume
 *                   fluxes).
 * 
 * @param[in] state1     First state.
 * @param[in] state2     Second state.
 * 
 * @param[out] fluxes    Computed numerical fluxes.
 */
template<FluidOption TFluid, NumericalFluxOption TFlux, int TDim, bool TStab>
MFEM_HOST_DEVICE inline
void ComputeNumericalFluxes(const State<TFluid,TDim> &state1,
                            const State<TFluid,TDim> &state2,
                            State<TFluid,TDim> (&fluxes)[TDim])
{
   using enum NumericalFluxOption;
   using enum FluidOption;
   
   if constexpr (TFlux == LocalLaxFriedrichs)
   { 
      State<TFluid,TDim> F_1[TDim];
      State<TFluid,TDim> F_2[TDim];

      //  TODO: Will need to think about viscous soon. Hurray.
      ComputeExactFluxes<TFluid,TDim>(state1, F_1);
      ComputeExactFluxes<TFluid,TDim>(state2, F_2);

      if constexpr (TFluid == CPG)
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
