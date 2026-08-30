#ifndef JONSS_MODEL
#define JONSS_MODEL

#include <cstdint>

#include <mfem/mfem.hpp>

namespace jonss
{

/// Fluid modeling options.
enum class FluidOption : std::uint8_t
{
   /// Standard, CPG air.
   AirCPG,

   /// Number of ModelOption enumerators.
   Size
};

/// Numerical flux options.
enum class NumericalFluxOption : std::uint8_t
{
   LocalLaxFriedrichs,
   Chandrashekar,
   Size
};

} // namespace jonss

#endif // JONSS_MODEL
