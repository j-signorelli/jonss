#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <jonss/jonss.hpp>

using namespace jonss;
using namespace Catch::Matchers;

TEMPLATE_TEST_CASE_SIG("Compute primitives" ,"[Primitives]",
                        ((FluidModelOption TModel), TModel),
                        FluidModelOption::AirCPG)
{
   using enum FluidModelOption;
   if constexpr (TModel == AirCPG)
   {
      // Initialize test state
      State<AirCPG,3> state;
      state.rho = 0.031762238707;
      state.rhoV[0] = state.rho*471.82963046368343;
      state.rhoV[1] = 0.0;
      state.rhoV[2] = state.rho*817.232892479548;
      state.rhoE = 0.0; // TODO

      // Compute primitives from test state
      Primitives<AirCPG,3> prim;
      ComputePrimitives(state,prim);

      // Initialize exact primitives
      Primitives<AirCPG,3> exact_prim;
      exact_prim.p = 568.749414924658;
      exact_prim.vel[0] = 471.82963046368343;
      exact_prim.vel[1] = 0.0;
      exact_prim.vel[2] = 817.232892479548;
      exact_prim.vel_sq = 890492.8007339844;
      exact_prim.H = 0; // TODO;

      // Check equality
      constexpr static double kTol = 1e-12;
      CHECK_THAT(prim.p, WithinRel(exact_prim.p, kTol));
      CHECK_THAT(prim.vel[0], WithinRel(exact_prim.vel[0], kTol));
      CHECK_THAT(prim.vel[1], WithinRel(exact_prim.vel[1], kTol));
      CHECK_THAT(prim.vel[2], WithinRel(exact_prim.vel[2], kTol));
      CHECK_THAT(prim.vel_sq, WithinRel(exact_prim.vel_sq, kTol));
      CHECK_THAT(prim.H, WithinRel(exact_prim.H, kTol));
   }
   else
   {
      FAIL("Unimplemented test for given FluidModelOption");
   }
}
