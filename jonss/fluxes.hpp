#ifndef JONSS_FLUXES
#define JONSS_FLUXES

#include <cstdint>

#include <mfem/mfem.hpp>

namespace jonss
{

/// Numerical flux options.
enum class FluxOption : std::uint8_t
{
   LocalLaxFriedrichs,
   Chandrashekar,
   Size
};

/**
 * @brief Function for computing the inviscid fluxes for CPG.
 * 
 * @details See "I do like CFD" for details. For now, it is
 * assumed that \f$\gamma=1.4\f$.
 * 
 * @tparam TDim      Spatial dimension.
 * 
 * @param state      State vector.
 * @param fluxes     Flux vectors.
 */
template<int TDim>
MFEM_HOST_DEVICE inline
void ComputeInviscidFluxes(const mfem::real_t (&state)[TDim+2],
                           mfem::real_t (&fluxes)[TDim][TDim+2]) 
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

template<int TDim, bool TStab>
MFEM_HOST_DEVICE inline
void LocalLaxFriedrichs(const mfem::real_t (&state1)[TDim+2],
                        const mfem::real_t (&state2)[TDim+2],
                        mfem::real_t (&flux)[TDim][TDim+2])
{
   mfem::real_t F_1[TDim][TDim+2];
   mfem::real_t F_2[TDim][TDim+2];
   ComputeInviscidFluxes(state1, F_1);
   ComputeInviscidFluxes(state2, F_2);

   for (int d = 0; d < TDim; d++)
   {
      for (int c = 0; c < TDim+2; c++)
      {
         flux[d][c] = 0.5*(F_1[d][c] + F_2[d][c]);
         if constexpr (TStab)
         {
            // TODO.
         }
      }
   }
}


} // namespace jonss

#endif  // JONSS_FLUXES
