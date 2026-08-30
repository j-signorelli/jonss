#include "utils.hpp"

using namespace mfem;


namespace joe_tests
{

Mesh Create2DWaveDomain(const double k,
                        const int num_waves,
                        const int n_elems_per_wave,
                        const int n_elems_y,
                        const double grid_ar)
{
   const double length_x = 2*M_PI*num_waves/k;
   const double x_spacing = length_x/n_elems_per_wave;
   const double length_y = length_x/grid_ar;
   const double y_spacing = length_y/n_elems_y;

   std::vector<double> xs(n_elems_per_wave+1), ys(n_elems_y+1);

   for (std::size_t i = 0; i < xs.size(); i++)
   {
      xs[i] = i*x_spacing;
   }

   for (std::size_t i = 0; i < ys.size(); i++)
   {
      ys[i] = i*y_spacing;
   }

   const int num_vert = (n_elems_per_wave+1)*(n_elems_y+1);
   const int num_elem = n_elems_per_wave*n_elems_y;
   const int num_bdr_elem = 2*n_elems_per_wave + 2*n_elems_y;
   Mesh mesh(2, num_vert, num_elem, num_bdr_elem);

   // Set vertices
   for (int j = 0; j < ys.size(); j++)
   {
      for (int i = 0; i < xs.size(); i++)
      {
         mesh.AddVertex(xs[i], ys[j]);
      }
   }

   // Set elements
   for (int j = 0; j < ys.size() - 1; j++)
   {
      for (int i = 0; i < xs.size() - 1; i++)
      {
         mesh.AddQuad(i   + j*xs.size(),
                      i+1 + j*xs.size(),
                      i+1 + (j+1)*xs.size(),
                      i   + (j+1)*xs.size());
      }
   }

   // Set boundary elements
   for (int j = 0; j < ys.size() - 1; j++)
   {
      // Left: bdr_attr 1
      mesh.AddBdrSegment(j*xs.size(), (j+1)*xs.size(), 1);

      // Right: bdr_attr 3
      mesh.AddBdrSegment((j+1)*xs.size() - 1, (j+2)*xs.size() - 1, 3);
   }

   // Bottom + Top
   for (int i = 0; i < xs.size() - 1; i++)
   {
      // Bottom: bdr_attr 2
      mesh.AddBdrSegment(i, i+1, 2);

      // Top: bdr_attr 4
      mesh.AddBdrSegment(i   + (ys.size() - 1)*xs.size(),
                         i+1 + (ys.size() - 1)*xs.size(),
                         4);
   }

   // Finalize quad mesh, generating boundary edge elements.
   mesh.FinalizeQuadMesh(true, 0, true);

   // Create periodic mesh
   // (https://mfem.org/howto/periodic-boundaries/#square)
   Vector x_translations({xs.back(), 0.0});
   Vector y_translations({0.0, ys.back()});
   std::vector<Vector> translations({x_translations, y_translations});
   Mesh periodic_mesh =
      Mesh::MakePeriodic(mesh,mesh.CreatePeriodicVertexMapping(translations));

   return periodic_mesh;
}


} // namespace joe_tests
