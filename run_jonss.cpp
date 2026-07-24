#include <mfem/mfem.hpp>

#include <jabber/jabber.hpp>

using namespace mfem;

/**
 * @brief Generate periodic rectangular grid with x-length set to ensure
 * exactly \p num_waves in the domain.
 *
 * @details Default x-spacing matches PC2 Low Re x-spacing.
 * 
 * @note The grid order is not set to \p order. For isoparametric elements,
 * set mfem::Mesh::SetCurvature() using \p order.
 */
mfem::Mesh CreateTestDomain(const mfem::real_t x_spacing=0.5e-4,
                              const mfem::real_t y_spacing=1e-4,
                              const int num_waves=1,
                              const mfem::real_t length_y=1e-3);

struct
{
   /// Freestream x-velocity.
   Vector u_infty({943.659260927367, 0.0});

   /// Freestream density.
   double rho_infty = 0.031762238707;

   /// Freestream pressure.
   double p_infty = 568.749414924658;

   /// Specific heat ratio.
   double gamma = 1.4;

   /// Pressure amplitude of single acoustic wave.
   double amp = 10;

   /// Frequency of planar wave.
   double freq = 136e3;

   /// Phase of planar wave.
   double phase = 0.0;

   /// Planar wave speed. true = slow, false = fast.
   bool slow = true;

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
   args.AddOption(&ctx.amp, "-a", "--amplitude", 
                  "Acoustic wave amplitude.");
   args.AddOption(&ctx.freq, "-f", "--freq",
                  "Acoustic wave frequency.");
   args.AddOption(&ctx.phase, "-psi", "--phase",
                  "Acoustic wave phase.");
   args.AddOption(&ctx.slow, "-slow", "--slow", "-fast", "--fast",
                  "Acoustic wave speed flag.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   if (rank == 0)
   {
      args.PrintOptions(cout);
   }

   // Initialize the mesh



   // Initialize acoustic field
   jabber::AcousticField field(2, coords, kPInfty, kRhoInfty, 
                                 kUInfty, kGamma);

   // Add wave
   jabber::Wave wave;
   wave.amplitude = kAmp;
   wave.frequency = kFreq;
   wave.k_hat = std::vector<double>({1.0,0.0});
   wave.phase = kPhase;
   wave.speed = kSpeed;
   field.AddWave(wave);

   field.Finalize();
   field.Compute(0.0);

   return 0;
}
