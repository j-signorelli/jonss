#ifndef JONSS_INTEGRATORS
#define JONSS_INTEGRATORS

#include "fluxes.hpp"

#include <mfem/mfem.hpp>

namespace jonss
{

/**
 * @brief Integrator for the volume term in split-form DGSEM @cite gassner2016,
 * given by
 *
 * \f[
 * \left(\displaystyle\int_E\phi_{ijk}\dfrac{\partial \phi_{\lambda\mu\nu}}
 *       {\partial\xi^l}d\tilde{V}\right)\tilde{F}^l_{\lambda\mu\nu}.
 *
 * \f]
 *
 * @details This integrator requires tensor-product elements using a nodal
 * basis, where collocation of solution nodes and quadrature points
 * (Gauss-Lobatto) is performed. No other `mfem::IntegrationRule` is supported.
 *
 * For `mfem::AssemblyLevel::LEGACY`, the legacy, CPU-only implementation using
 * `mfem::NonlinearFormIntegrator::AssembleElementVector()` is employed.
 *
 */
class SplitFormVolumeFluxIntegrator : public mfem::NonlinearFormIntegrator
{
protected:

public:
   SplitFormVolumeFluxIntegrator() {}

   void AssembleElementVector(const mfem::FiniteElement &el,
                              mfem::ElementTransformation &trans,
                              const mfem::Vector &elfun,
                              mfem::Vector &elvect) override;

// const mfem::IntegrationRule* GetDefaultIntegrationRule(
//                            const mfem::FiniteElement &trial_fe,
//                            const mfem::FiniteElement &test_fe,
//                            const mfem::ElementTransformation &trans) const
// {
//    // TODO.
//    return &(IntRules.Get(test_fe.GetGeomType(), test_fe.GetOrder()));
// }


};

} // namespace jonss

#endif // JONSS_INTEGRATORS
