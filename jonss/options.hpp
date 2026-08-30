#ifndef JONSS_MODEL
#define JONSS_MODEL

#include <cstdint>

namespace jonss
{

/// Fluid model options.
enum class FluidOption : std::uint8_t
{
   /// Inviscid, calorically perfect gas.
   EulerCPG,

   /// Viscous, calorically perfect gas.
   NavierStokesCPG,

   /// Number of FluidOption enumerators.
   Size
};

/// Viscosity model options.
enum class ViscosityOption : std::uint8_t
{
   /// Sutherland's law for viscosity.
   Sutherland,

   /// Number of ViscosityOption enumerators.
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
