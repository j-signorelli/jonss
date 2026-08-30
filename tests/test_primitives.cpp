#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <jonss/jonss.hpp>

using namespace jonss;

TEMPLATE_TEST_CASE_SIG("Compute primitives" ,"[Primitives]",
                        ((FluidModelOption TModel), TModel),
                        FluidModelOption::AirCPG)
{
   using enum FluidModelOption;
   if constexpr (TModel == AirCPG)
   {
      // TODO: Initialize each and verify.
      State<AirCPG,3> state;
      Primitives<AirCPG,3> prim;
      ComputePrimitives(state,prim);
   }
   else
   {
      FAIL("Unimplemented test for given FluidModelOption");
   }
}