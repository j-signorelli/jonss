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
 * @tparam TVisc  Viscosity model.
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
   /// Reference temperature. 273.15 for standard air.
   const mfem::real_t T_ref;

   /// Reference dynamic viscosity at \ref T_ref. 1.716e-5 for standard air.
   const mfem::real_t mu_ref;

   /// Sutherland temperature. 110.4 for standard air.
   const mfem::real_t S;

   /// Evaluated, simplified prefactor 
   const mfem::real_t C1;

   ViscosityConstants<ViscosityOption::Sutherland>(
      const mfem::real_t &T_ref_,
      const mfem::real_t &mu_ref_,
      const mfem::real_t &S_)
   : T_ref(T_ref_),
     mu_ref(mu_ref_),
     S(S_),
     C1(mu_ref*(T_ref+S)/std::pow(T_ref, 2.0/3.0)) {}
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
   /// Specific heat ratio. 1.4 for standard air.
   const mfem::real_t gamma;

   /// Specific gas constant. 287.05 for standard air.
   const mfem::real_t R;

   FluidConstants<FluidOption::CPG, ViscosityOption::Inviscid>(
      const mfem::real_t &gamma_,
      const mfem::real_t &R_)
      : gamma(gamma_), R(R_) {}
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
   /// Prandtl number. 0.72 for standard air.
   const mfem::real_t Pr;

   /// Bulk viscosity factor. 0.6 for standard air ( @todo cite ). 
   const mfem::real_t bulk_visc_fac;

   /// Viscosity model.
   const ViscosityConstants<TVisc> visc_constants;

   FluidConstants<FluidOption::CPG,TVisc>(
      const mfem::real_t &gamma_,
      const mfem::real_t &R_,
      const mfem::real_t &Pr_,
      const ViscosityConstants<TVisc> &visc_constants_,
      const mfem::real_t &bulk_visc_fac_)
      : FluidConstants<FluidOption::CPG, ViscosityOption::Inviscid>(gamma_, R_),
        Pr(Pr_),
        visc_constants(visc_constants_),
        bulk_visc_fac(bulk_visc_fac_) {}
};

} // namespace jonss

#endif // JONSS_CONSTANTS
