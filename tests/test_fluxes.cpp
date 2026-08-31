#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <jonss/jonss.hpp>

using namespace jonss;

using enum FluidOption;
using enum ViscosityOption;

TEMPLATE_TEST_CASE_SIG("Compute inviscid fluxes" ,"[Fluxes]",
                        ((FluidOption TFluid, ViscosityOption TVisc), TFluid, TVisc),
                        (CPG, Inviscid),
                        (CPG, Sutherland))
{
   using enum FluidOption;
   if constexpr (TFluid == CPG)
   {
      if constexpr (TVisc == Inviscid)
      {
         // TODO: Initialize each and verify.
         const FluidConstants<CPG, Inviscid> constants(1.4, 287.05);
         State<CPG,3> state;
         Primitives<CPG,Inviscid,3> prim;
         ComputePrimitives(constants,state,prim);    
         
         State<CPG,3> fluxes[3];
         //ComputeInviscidFluxes(state,prim,fluxes);
      }
      else
      {
         FAIL("Unimplemented test for given ViscosityOption");
      }
   }
   else
   {
      FAIL("Unimplemented test for given FluidOption");
   }
}
