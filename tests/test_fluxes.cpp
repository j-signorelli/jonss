#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <jonss/jonss.hpp>

using namespace jonss;

TEMPLATE_TEST_CASE_SIG("Compute inviscid fluxes" ,"[Fluxes]",
                        ((FluidOption TModel), TModel),
                        FluidOption::EulerCPG,
                        FluidOption::NavierStokesCPG)
{
   using enum FluidOption;
   if constexpr (TModel == EulerCPG)
   {
      // TODO: Initialize each and verify.
      const FluidConstants<EulerCPG> constants(1.4);
      State<EulerCPG,3> state;
      Primitives<EulerCPG,3> prim;
      ComputePrimitives(constants,state,prim);    
      
      State<EulerCPG,3> fluxes[3];
      ComputeInviscidFluxes(state,prim,fluxes);
   }
   else
   {
      FAIL("Unimplemented test for given FluidOption");
   }
}
