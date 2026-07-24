/**
 * @brief Single acoustic wave, periodic-domain 2D test.
 */
#include <jonss/jonss.hpp>
#include <jabber/jabber.hpp>
#include <mfem/mfem.hpp>

using namespace mfem;
using namespace jonss::joe_tests;

struct TestContext
{
   /// Freestream x-velocity.
   double u_infty =943.659260927367;

   /// Freestream density.
   double rho_infty = 0.031762238707;

   /// Freestream pressure.
   double p_infty = 568.749414924658;

   /// Specific heat ratio.
   double gamma = 1.4;

   /// Number of waves to include in domain.
   int num_waves = 1;

   /// Number of elements / wave in x-direction
   int num_elems_per_wave = 10;

   /// Pressure amplitude of single acoustic wave.
   double amp = 10;

   /// Frequency of planar wave.
   double freq = 136e3;

   /// Phase of planar wave.
   double phase = 0.0;

   /// Planar wave speed. true = slow, false = fast.
   bool slow = true;
   
   /// Enable/disable visualization w/ GLVis.
   bool visualization = true;

} ctx;

int main(int argc, char* argv[])
{
   // Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int rank = Mpi::WorldRank();
   Hypre::Init();

   // Parse command line arguments
   OptionsParser args(argc, argv);
   args.AddOption(&ctx.u_infty, "-u", "--u-infty", 
                  "Freestream velocity vector.");
   args.AddOption(&ctx.rho_infty, "-rho", "--rho-infty", 
                  "Freestream density.");
   args.AddOption(&ctx.p_infty, "-p", "--p-infty",
                  "Freestream pressure.");
   args.AddOption(&ctx.gamma, "-g", "--gamma",
                  "Specific heat ratio.");
   args.AddOption(&ctx.num_waves, "-nw", "--num-waves",
                  "Number of waves to include in test.");
   args.AddOption(&ctx.num_elems_per_wave, "-ne", "--num-elems",
                  "Number of elements **per wave** in x.");
   args.AddOption(&ctx.amp, "-a", "--amplitude", 
                  "Acoustic wave amplitude.");
   args.AddOption(&ctx.freq, "-f", "--freq",
                  "Acoustic wave frequency.");
   args.AddOption(&ctx.phase, "-psi", "--phase",
                  "Acoustic wave phase.");
   args.AddOption(&ctx.slow, "-slow", "--slow", "-fast", "--fast",
                  "Acoustic wave speed flag.");
   args.AddOption(&ctx.visualization, "-vis", "--visualization",
                  "-no-vis", "--no-visualization",
                  "Enable or disable GLVis visualization.");                  
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(std::cout);
      return 1;
   }
   if (rank == 0)
   {
      args.PrintOptions(std::cout);
   }
   
   const double k = jabber::ComputeWavenumber()
   Mesh mesh = Create2DWaveDomain()

   return 0;
}