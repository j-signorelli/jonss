#ifndef JONSS_MODEL
#define JONSS_MODEL

#include <cstdint>

namespace jonss
{

/// Fluid model options.
enum class FluidOption : std::uint8_t
{
   /// Calorically-perfect gas.
   CPG
};

/// Viscosity model options.
enum class ViscosityOption : std::uint8_t
{
   /// Inviscid (Euler).
   Inviscid,

   /// Sutherland's law for viscosity.
   Sutherland
};

/// Numerical flux options.
enum class NumericalFluxOption : std::uint8_t
{
   LocalLaxFriedrichs,
   Chandrashekar
};

} // namespace jonss

#endif // JONSS_MODEL
