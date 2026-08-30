#ifndef JONSS_CONSTANTS
#define JONSS_CONSTANTS

#include "options.hpp"

#include <mfem/mfem.hpp>

namespace jonss
{

/**
 * @brief Constants associated with a particular ViscosityOption.
 *
 * @details This should be explicitly specialized for every
 * \ref ViscosityOption.
 *
 * @tparam TFluid    Fluid model.
 */
template<ViscosityOption TVisc>
struct ViscosityConstants
{
   ViscosityConstants() = delete;
};

/// ViscosityConstants explicit specialization for Sutherland's Law.
template<>
struct ViscosityConstants<ViscosityOption::Sutherland>
{
   /// Reference temperature
   const mfem::real_t T_ref;

   /// Viscosity at \ref T_ref.
   const mfem::real_t mu_ref;

   /// Sutherland temperature.
   const mfem::real_t S;
};

/**
 * @brief Constants associated with a particular \ref FluidOption and
 * \ref ViscosityOption.
 *
 * @details This should be specialized for every valid
 * \ref FluidOption and \ref ViscosityOption combination.
 *
 * @tparam TFluid    Fluid model.
 * @tparam TVisc     Viscosity model.
 */
template<FluidOption TFluid, ViscosityOption TVisc>
struct FluidConstants
{
   FluidConstants() = delete;
};

/**
 * @brief FluidConstants explicit specialization for an inviscid,
 * calorically-perfect gas.
 */
template<>
struct FluidConstants<FluidOption::CPG, ViscosityOption::Inviscid>
{
   const mfem::real_t gamma;

   FluidConstants<FluidOption::CPG, ViscosityOption::Inviscid>(
      const mfem::real_t &gamma_)
      : gamma(gamma_) {}
};

/**
 * @brief FluidConstants partial specialization for a viscous,
 * calorically-perfect gas.
 *
 */
template<ViscosityOption TVisc>
struct FluidConstants<FluidOption::CPG, TVisc>
   : public FluidConstants<FluidOption::CPG, ViscosityOption::Inviscid>
{
   /// Prandtl number.
   const mfem::real_t Pr;

   /// Viscosity model.
   const ViscosityConstants<TVisc> visc_constants;

   /// Bulk viscosity factor.
   const mfem::real_t bulk_visc_fac;


   FluidConstants<FluidOption::CPG,TVisc>(
      const mfem::real_t &gamma_,
      const mfem::real_t &Pr_,
      const ViscosityConstants<TVisc> &visc_constants_,
      const mfem::real_t &bulk_visc_fac_)
      : FluidConstants<FluidOption::CPG, ViscosityOption::Inviscid>(gamma_),
        Pr(Pr_),
        visc_constants(visc_constants_),
        bulk_visc_fac(bulk_visc_fac_) {}
};

} // namespace jonss

#endif // JONSS_CONSTANTS
