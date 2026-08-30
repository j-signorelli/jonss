#ifndef JONSS_PRIMITIVES
#define JONSS_PRIMITIVES

#include "options.hpp"


namespace jonss
{

/**
 * @brief Struct to hold relevant, model-specific primitive and auxiliary
 * variables.
 * 
 * @details For variables, dependent on the state, that are to be computed 
 * once to avoid expensive re-computations on-the-fly and support storage 
 * in smem.  This should be explicitly specialized for each 
 * \ref FluidModelOption.
 * 
 * @tparam TModel    Fluid model.
 * @tparam TDim      Spatial dimension.
 */
template<FluidModelOption TModel, int TDim>
struct Primitives
{
   MFEM_HOST_DEVICE Primitives() = delete;
};

/**
 * @brief Primitives specialization for FluidModelOption::AirCPG.
 * 
 * @tparam TDim      Spatial dimension.
 */
template<int TDim>
struct Primitives<FluidModelOption::AirCPG, TDim>
{
   /// Specific heat ratio.
   static constexpr double kGamma = 1.4;

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
 * @tparam TModel    Fluid model.
 * @tparam TDim      Spatial dimension. Templated to allow compiler to unroll
 *                   loops.
 * 
 * @param state      State struct.
 * @return prim      Primitives struct.
 */
template<FluidModelOption TModel, int TDim>
MFEM_HOST_DEVICE inline
void ComputePrimitives(const State<TModel,TDim> &state,
                       Primitives<TModel,TDim> &prim)
{
   using enum FluidModelOption;

   if constexpr (TModel == AirCPG)
   {
      // See "I do like CFD" for details.
      
      prim.vel_sq = 0.0;
      for (int d = 0; d < TDim; d++)
      {
         prim.vel[d] = state.rhoV[d]/state.rho;
         prim.vel_sq += prim.vel[d]*prim.vel[d];
      }
      prim.p = (prim.kGamma-1.0)*(state.rhoE-0.5*state.rho*prim.vel_sq);
      prim.H = (state.rhoE + prim.p)/state.rho;
   }
   else
   {
      static_assert(TDim != TDim, "Unimplemented model.");
   }
}

} // namespace jonss

#endif // JONSS_PRIMITIVES
