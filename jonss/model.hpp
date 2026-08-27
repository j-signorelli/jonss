#ifndef JONSS_MODEL
#define JONSS_MODEL

#include <cstdint>

#include <mfem/mfem.hpp>

namespace jonss
{

/**
 * @brief Modeling options.
 * 
 * @details These are used to assist in determining the 
 * number of equations/conserved variables and models
 * at compile-time to prevent any vtable lookups at 
 * runtime.
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
 * @brief Struct to hold relevant, model-specific
 * primitive and auxiliary variables to be computed once, 
 * to avoid expensive re-computations on-the-fly and 
 * support storage in smem. Templated to support 
 * explicit specializations for a given ModelOption.
 * 
 * @tparam TDim      Spatial dimension.
 * @tparam TModel    Fluid model.
 */
template<int TDim, ModelOption TModel>
struct Primitives
{

};

/// Primitives explicit-specialization for ModelOption::AirCPG.
template<int TDim>
struct Primitives<TDim, ModelOption::AirCPG>
{
   /// @brief  Velocity.
   mfem::real_t vel[TDim];

   /// @brief Velocity magnitude, squared.
   mfem::real_t vel_sq;

   /// @brief Pressure.
   mfem::real_t p;

   /// @brief @todo.
   mfem::real_t H;
};

template<int TDim, ModelOption TModel>
MFEM_HOST_DEVICE inline
void ComputePrimitives(const mfem::real_t (&state)[kNumEq<TDim,TModel>]
                        Primitives<TDim,TModel> &prim)
{
   if constexpr (TModel == AirCPG)
   {
      mfem::real_t vel[TDim];
      prim.vel_sq = 0.0;
      for (int d = 0; d < TDim; d++)
      {
         prim.vel[d] = state[1+d]/state[0];
         prim.vel_sq += prim.vel[d]*prim.vel[d];
      }
      prim.p = (1.4-1.0)*(state[TDim+1]-state[0]*0.5*prim.vel_sq);
      prim.H = state[TDim+1]/state[0] + prim.p/state[0];
   }
}

} // namespace jonss

#endif // JONSS_MODEL