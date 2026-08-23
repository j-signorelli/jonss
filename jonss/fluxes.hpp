#ifndef JONSS_FLUXES
#define JONSS_FLUXES

#include <cstdint>

#include <mfem/mfem.hpp>

namespace jonss
{

/**
 * @brief Modeling options.
 * 
 * @details These are used to assist in determining the 
 * number of equations/conserved variables and models
 * at compile-time to prevent any vtable lookups.
 * 
 */
enum class ModelOption : std::uint8_t
{
   /// Air, assuming calorically-perfect gas.
   AirCPG,

   /// Number of ModelOption enumerators.
   Size
};

/**
 * @brief Compile-time variable for the number of equations
 * for a given ModelOption.
 * 
 * @tparam TDim    Spatial dimension.
 * @tparam TModel  Fluid model.
 */
template<int TDim, ModelOption TModel>
constexpr static int kNumEq = 
[]()
{
   using enum ModelOption;

   if constexpr (TModel == AirCPG)
   {
      // Density, Momentum (TDim), and Energy
      return TDim + 2;
   }
   else
   {
      static_assert(TDim != TDim, "Unimplemented model.");
   }
}();

/**
 * @brief Function for computing the inviscid fluxes.
 * 
 * @details See "I do like CFD" for details.
 * 
 * @warning Only ModelOption::AirCPG is implemented.
 * 
 * @tparam TDim      Spatial dimension. Templated to allow compiler to unroll
 *                   loops.
 * @tparam TModel    Fluid model.
 * 
 * @param state      State vector.
 * @param fluxes     Flux vectors.
 */
template<int TDim, ModelOption TModel>
MFEM_HOST_DEVICE inline
void ComputeInviscidFluxes(const mfem::real_t (&state)[kNumEq<TDim,TModel>],
                           mfem::real_t (&fluxes)[TDim][kNumEq<TDim,TModel>]) 
{
   using enum ModelOption;

   if constexpr (TModel == AirCPG)
   {
      mfem::real_t vel[TDim];
      mfem::real_t vel_sq = 0.0;
      for (int d = 0; d < TDim; d++)
      {
         vel[d] = state[1+d]/state[0];
         vel_sq += vel[d]*vel[d];
      }
      const mfem::real_t p = (1.4-1.0)*(state[TDim+1]-state[0]*vel_sq/2.0);
      const mfem::real_t H = state[TDim+1]/state[0] + p/state[0];

      for (int di = 0; di < TDim; di++)
      {
         // Set density fluxes
         fluxes[di][0] = state[1+di];

         // Set momentum fluxes
         for (int dj = 0; dj < TDim; dj++)
         {
            fluxes[di][dj+1] = state[1+di]*vel[dj];
         }
         fluxes[di][di+1] += p;

         // Set energy fluxes
         fluxes[di][TDim+1] = state[1+di]*H;
      }
   }
   else
   {
      static_assert(TDim != TDim, "Unimplemented model.");
   }
}

/// Numerical flux options.
enum class FluxOption : std::uint8_t
{
   LocalLaxFriedrichs,
   Chandrashekar,
   Size
};

/**
 * @brief Compute the numerical fluxes.
 * 
 * @tparam TDim      Spatial dimension. Templated to allow compiler to unroll
 *                   loops.
 * @tparam TFlux     Numerical flux.
 * @tparam TStab     If true, include stabilization (for face fluxes).
 *                   If false, do not include stabilization (for volume fluxes).
 * @tparam TModel    Fluid model.
 * 
 * @param state1     First state.
 * @param state2     Second state.
 * @return fluxes    Computed numerical fluxes.
 */
template<int TDim, bool TStab, FluxOption TFlux, ModelOption TModel>
MFEM_HOST_DEVICE void ComputeNumericalFluxes(
      const mfem::real_t (&state1)[kNumEq<TDim,TModel>],
      const mfem::real_t (&state2)[kNumEq<TDim,TModel>],
      mfem::real_t (&fluxes)[TDim][kNumEq<TDim,TModel>])
{
   using enum FluxOption;

   if constexpr (TFlux == LocalLaxFriedrichs)
   { 
      mfem::real_t F_1[TDim][kNumEq<TDim,TModel>];
      mfem::real_t F_2[TDim][kNumEq<TDim,TModel>];
      ComputeInviscidFluxes<TDim,TModel>(state1, F_1);
      ComputeInviscidFluxes<TDim,TModel>(state2, F_2);

      for (int d = 0; d < TDim; d++)
      {
         for (int c = 0; c < kNumEq<TDim,TModel>; c++)
         {
            fluxes[d][c] = 0.5*(F_1[d][c] + F_2[d][c]);
            if constexpr (TStab)
            {
               // Compute lambda_max on-the-fly
               mfem::real_t lambda_max;
               // TODO: Figure out exactly what these should be...
               // I think they are just the eigenvalues in the I do like CFD.
               // But Theseus confusing me.
               fluxes[d][c] -= 0.5*lambda_max*(state2[c] - state1[c]);
            }
         }
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
