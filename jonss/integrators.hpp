#include <mfem/mfem.hpp>

namespace jonss
{

class SplitFormVolumeFluxIntegrator : public mfem::NonlinearFormIntegrator
{
private:
   const mfem::NumericalFlux &num_flux_;
public:
   SplitFormVolumeFluxIntegrator(const mfem::NumericalFlux &num_flux)
   : num_flux_(num_flux) {}

   void AssembleElementVector(const mfem::FiniteElement &el,
                              mfem::ElementTransformation &trans,
                              const mfem::Vector &elfun,
                              mfem::Vector &elvect) override; 
   // Partial assembly incoming...
};

} // namespace jonss
