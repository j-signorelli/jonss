#include "integrators.hpp"

using namespace mfem;

namespace jonss
{

void SplitFormVolumeFluxIntegrator::AssembleElementVector(
   const FiniteElement &el, 
   ElementTransformation &trans,
   const Vector &elfun,
   Vector &elvect)
{
   MFEM_ASSERT(dynamic_cast<const TensorBasisElement*>(&fe) != nullptr,
               "Only tensor product elements are currently supported.");
   
   el.GetDofToQuad
   
}

} // namespace jonss
