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
 * in smem.  This should be specialized for every valid
 * \ref FluidOption and \ref ViscosityOption combination.
 *
 * @tparam TFluid    Fluid model.
 * @tparam TDim      Spatial dimension.
 */
template<FluidOption TFluid, ViscosityOption TVisc, int TDim>
struct Primitives
{
   MFEM_HOST_DEVICE Primitives() = delete;
};

/**
 * @brief Primitives explicit specialization for an inviscid,
 * calorically-perfect gas.
 *
 * @tparam TDim      Spatial dimension.
 */
template<int TDim>
struct Primitives<FluidOption::CPG, ViscosityOption::Inviscid, TDim>
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
 * @brief Primitives partial specialization for a viscous,
 * calorically-perfect gas.
 *
 * @tparam TDim      Spatial dimension.
 * @tparam TVisc     Viscosity model.
 */
template<ViscosityOption TVisc, int TDim>
struct Primitives<FluidOption::CPG, TVisc, TDim>
   : public Primitives<FluidOption::CPG, ViscosityOption::Inviscid, TDim>
{
   /// Viscous stress tensor.
   mfem::real_t tau[TDim][TDim];

   /// Heat flux.
   mfem::real_t q[TDim];
};

/**
 * @brief Function for initializing a Primitives struct.
 *
 * @tparam TFluid    Fluid model.
 * @tparam TVisc     Viscosity model.
 * @tparam TDim      Spatial dimension. Templated to allow compiler to unroll
 *                   loops.
 *
 * @param[in] constants   Constants struct.
 * @param[in] state       State struct.
 * @param[out] prim       Primitives struct.
 */
template<FluidOption TFluid, ViscosityOption TVisc, int TDim>
MFEM_HOST_DEVICE inline
void ComputePrimitives(const FluidConstants<TFluid,TVisc> &constants,
                       const State<TFluid,TDim> &state,
                       Primitives<TFluid,TVisc,TDim> &prim)
{
   using enum FluidOption;
   using enum ViscosityOption;

   if constexpr (TFluid == CPG)
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

      if constexpr (TVisc != Inviscid)
      {
         // TODO...
      }
   }
   else
   {
      static_assert(TDim != TDim, "Unimplemented model.");
   }
}



} // namespace jonss

#endif // JONSS_PRIMITIVES
