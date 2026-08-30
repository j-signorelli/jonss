#ifndef JONSS_CONSTANTS
#define JONSS_CONSTANTS

#include "options.hpp"

#include <mfem/mfem.hpp>

namespace jonss
{

/**
 * @brief Constants associated with a particular FluidOption.
 * 
 * @details This should be explicitly specialized for every
 * \ref FluidOption.
 * 
 * @tparam TFluid 
 */
template<FluidOption TFluid>
struct FluidConstants
{
   FluidConstants() = delete;
};

/**
 * @brief FluidConstants explicit specialization for \ref
 * FluidOption::EulerCPG.
 */
template<>
struct FluidConstants<FluidOption::EulerCPG>
{
   const mfem::real_t gamma;

   FluidConstants<FluidOption::EulerCPG>(const mfem::real_t &gamma_)
   : gamma(gamma_) {}
};

/**
 * @brief FluidConstants explicit specialization for \ref
 * FluidOption::NavierStokesCPG.
 */
template<>
struct FluidConstants<FluidOption::NavierStokesCPG> 
   : public FluidConstants<FluidOption::EulerCPG>
{
   /// Prandtl number.
   const mfem::real_t Pr;

   /// Viscosity model.
   const ViscosityOption visc_model;

   /// Bulk viscosity factor.
   const mfem::real_t bulk_visc_fac;


   FluidConstants<FluidOption::NavierStokesCPG>(
         const mfem::real_t &gamma_, 
         const mfem::real_t &Pr_, 
         const ViscosityOption &visc_model_, 
         const mfem::real_t &bulk_visc_fac_)
   : FluidConstants<FluidOption::EulerCPG>(gamma_),
     Pr(Pr_),
     visc_model(visc_model_),
     bulk_visc_fac(bulk_visc_fac_) {}
};

} // namespace jonss

#endif // JONSS_CONSTANTS
