#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <jonss/jonss.hpp>

#include <type_traits>

using namespace jonss;
using namespace Catch::Matchers;

using enum FluidOption;
using enum ViscosityOption;

TEMPLATE_TEST_CASE_SIG("Compute primitives" ,"[Primitives]",
                        ((FluidOption TFluid, 
                          ViscosityOption TVisc),
                         TFluid, TVisc),
                        (CPG, Inviscid),
                        (CPG, Sutherland))
{
   if constexpr (TFluid == CPG)
   {
      if constexpr (TVisc == Inviscid)
      {
         // Initialize constants
         const FluidConstants<CPG,Inviscid> constants(1.4, 287.05);

         // Initialize test state
         State<CPG,3> state;
         state.rho = 0.031762238707;
         state.rhoV[0] = state.rho*471.82963046368343;
         state.rhoV[1] = 0.0;
         state.rhoV[2] = state.rho*817.232892479548;
         state.rhoE = 0.0; // TODO

         // Compute primitives from test state
         Primitives<CPG,Inviscid,3> prim;
         ComputePrimitives(constants,state,prim);

         // Initialize exact primitives
         Primitives<CPG,Inviscid,3> exact_prim;
         exact_prim.p = 568.749414924658;
         exact_prim.vel[0] = 471.82963046368343;
         exact_prim.vel[1] = 0.0;
         exact_prim.vel[2] = 817.232892479548;
         exact_prim.vel_sq = 890492.8007339844;
         exact_prim.H = 0; // TODO;

         // Check equality
         constexpr static mfem::real_t kTol = 
         [&]()
         {
            if constexpr (std::is_same_v<mfem::real_t,double>)
            {
               return 1e-12;
            }
            else
            {
               return 1e-9;
            }
         }();

         CHECK_THAT(prim.p, WithinRel(exact_prim.p, kTol));
         CHECK_THAT(prim.vel[0], WithinRel(exact_prim.vel[0], kTol));
         CHECK_THAT(prim.vel[1], WithinRel(exact_prim.vel[1], kTol));
         CHECK_THAT(prim.vel[2], WithinRel(exact_prim.vel[2], kTol));
         CHECK_THAT(prim.vel_sq, WithinRel(exact_prim.vel_sq, kTol));
         CHECK_THAT(prim.H, WithinRel(exact_prim.H, kTol));
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
