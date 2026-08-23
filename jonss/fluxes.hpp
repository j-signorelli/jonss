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
 * @brief Compile-time variable for getting the number of equations
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
 * @brief Generic, "compile-time-abstract" struct representing a numerical flux.
 * 
 * 
 * @details This should be AT LEAST partially specialized for every FluxOption!
 * 
 * A struct is used to support member variables in derived types, if one
 * wanted to have runtime customization of a numerical flux scheme. Additionally, partial
 * specialization is supported, such that an explicit specialization is not needed for
 * every single combination of \p TFlux and \p TModel.
 * 
 * @tparam TFlux     Numerical flux.
 * @tparam TModel    Fluid model.
 */
template<FluxOption TFlux, ModelOption TModel>
struct NumericalFlux
{
   /**
    * @brief Compute the numerical flux.
    * 
    * @tparam TDim      Spatial dimension. Templated to allow compiler to unroll
    *                   loops.
    * @tparam TStab     If true, enable stabilization term. \p nor must be 
    *                   set! If false, \p nor is not used. In other words, if
    *                   true, compute face fluxes in DGSEM. If false, compute 
    *                   volume fluxes in DGSEM.
    * @param state1     First state.
    * @param state2     Second state.
    * @param nor        Pointer to normal vector - only used if \p TStab is
    *                   true.
    * @return fluxes    Computed numerical fluxes.
    */
   template<int TDim, bool TStab>
   MFEM_HOST_DEVICE void ComputeFlux(
         const mfem::real_t (&state1)[kNumEq<TDim,TModel>],
         const mfem::real_t (&state2)[kNumEq<TDim,TModel>],
         const mfem::real_t *nor[TDim],
         mfem::real_t (&fluxes)[TDim][kNumEq<TDim,TModel>])
   {
      static_assert(TDim != TDim, "Unimplemented flux.");
   }
};

/**
 * @brief THE function for computing volume fluxes.
 * 
 */
template<int TDim, FluxOption TFlux, ModelOption TModel>
void MFEM_HOST_DEVICE ComputeVolumeFlux(
         const NumericalFlux<TFlux,TModel> &numerical_flux,
         const mfem::real_t (&state1)[kNumEq<TDim,TModel>],
         const mfem::real_t (&state2)[kNumEq<TDim,TModel>],
         mfem::real_t (&fluxes)[TDim][kNumEq<TDim,TModel>])
{
   numerical_flux.ComputeFlux<false>(state1, state2, nullptr, fluxes);
}

/**
 * @brief THE function for computing face fluxes.
 */
template<int TDim, FluxOption TFlux, ModelOption TModel>
void MFEM_HOST_DEVICE ComputeFaceFlux(
         const NumericalFlux<TFlux,TModel> &numerical_flux,
         const mfem::real_t (&state1)[kNumEq<TDim,TModel>],
         const mfem::real_t (&state2)[kNumEq<TDim,TModel>],
         const mfem::real_t (&nor)[TDim],
         mfem::real_t (&fluxes)[TDim][kNumEq<TDim,TModel>])
{
   numerical_flux.ComputeFlux<true>(state1, state2, &nor, fluxes);
}

/**
 * @brief Local Lax-Friedrichs numerical flux specialization.
 * 
 */
template<ModelOption TModel>
struct NumericalFlux<FluxOption::LocalLaxFriedrichs,TModel>
{
   template<int TDim, bool TStab>
   MFEM_HOST_DEVICE void ComputeFlux(
         const mfem::real_t (&state1)[kNumEq<TDim,TModel>],
         const mfem::real_t (&state2)[kNumEq<TDim,TModel>],
         const mfem::real_t *nor[TDim],
         mfem::real_t (&fluxes)[TDim][kNumEq<TDim,TModel>])
   {
      using enum ModelOption;

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
               // 
            }
         }
      }
   }
   }
};

/**
 * @brief Chandrashekar numerical flux specialization.
 */
template<ModelOption TModel>
struct NumericalFlux<FluxOption::Chandrashekar,TModel>
{
   template<int TDim, bool TStab>
   MFEM_HOST_DEVICE void ComputeFlux(
         const mfem::real_t (&state1)[kNumEq<TDim,TModel>],
         const mfem::real_t (&state2)[kNumEq<TDim,TModel>],
         const mfem::real_t *nor[TDim],
         mfem::real_t (&fluxes)[TDim][kNumEq<TDim,TModel>])
   {
      using enum ModelOption;

      
   }
};


} // namespace jonss

#endif  // JONSS_FLUXES
