/*
  Copyright (C) 2018 - 2023 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
 */


#include <aspect/mesh_deformation/phase_boundary.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>

#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    PhaseBoundary<dim>::PhaseBoundary()
    {}



    template <int dim>
    void
    PhaseBoundary<dim>::update ()
    {

    }



    template <int dim>
    void PhaseBoundary<dim>::phase_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                            const IndexSet &mesh_locally_owned,
                                            const IndexSet &mesh_locally_relevant,
                                            LinearAlgebra::Vector &output,
                                            const std::set<types::boundary_id> &boundary_ids) const
    {
      // Do nothing at time zero
      if (this->get_timestep_number() < 1)
        return;

      LinearAlgebra::Vector phase_boundary;
      phase_boundary.reinit(mesh_locally_owned, this->get_mpi_communicator());

      // Initialize quadrature on the support points of the system.
      // This is where we want to compute the mesh deformation.
      std::vector<Point<dim-1>> face_support_points = mesh_deformation_dof_handler.get_fe().base_element(0).get_unit_face_support_points();
      Quadrature<dim-1> face_support_quadrature(face_support_points);
      FEFaceValues<dim> fe_support_values (this->get_mapping(),
                                           mesh_deformation_dof_handler.get_fe(),
                                           face_support_quadrature,
                                           update_values | update_normal_vectors
                                           | update_gradients | update_quadrature_points);

      // We need another fevalues object using the finite element
      // that contains the temperature and pressure.
      FEFaceValues<dim> fe_values_p_T (this->get_mapping(),
                                       this->get_fe(),
                                       face_support_quadrature,
                                       update_values | update_gradients | update_quadrature_points);

      // The number of quadrature points on a mesh deformation surface face
      const unsigned int n_face_q_points = fe_support_values.n_quadrature_points;

      // What we need to build our system on the mesh deformation element

      // The nr of shape functions per face of the mesh deformation element
      const unsigned int dofs_per_face = mesh_deformation_dof_handler.get_fe().dofs_per_face;

      // Map of local to global cell dof indices
      std::vector<types::global_dof_index> face_dof_indices (dofs_per_face);

      // Vector for getting the local dim displacement and initial topography values
      std::vector<Tensor<1, dim>> displacement_values(n_face_q_points);
      std::vector<Tensor<1, dim>> initial_topography_values(n_face_q_points);

      // The global displacements/initial topography on the MeshDeformation FE
      LinearAlgebra::Vector displacements = this->get_mesh_deformation_handler().get_mesh_displacements();
      LinearAlgebra::Vector initial_topography = this->get_mesh_deformation_handler().get_initial_topography();

      // An extractor for the dim-valued displacement and initial topography vectors.
      FEValuesExtractors::Vector extract_vertical_displacements(0);
      FEValuesExtractors::Vector extract_initial_topography(0);

      // Cell iterator over the FE that contains the solution.
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active();

      // Iterate over all cells to find those at the mesh deformation boundary.
      for (const auto &phase_boundary_cell : mesh_deformation_dof_handler.active_cell_iterators())
        {
          if (phase_boundary_cell->at_boundary() && phase_boundary_cell->is_locally_owned())
            for (const unsigned int face_no : phase_boundary_cell->face_indices())
              if (phase_boundary_cell->face(face_no)->at_boundary())
                {
                  // Boundary indicator of current cell face
                  const types::boundary_id boundary_indicator
                    = phase_boundary_cell->face(face_no)->boundary_id();

                  // Only apply phase boundary to the requested boundaries.
                  if (boundary_ids.find(boundary_indicator) == boundary_ids.end())
                    continue;

                  // Recompute values, gradients, etc on the faces
                  fe_support_values.reinit (phase_boundary_cell, face_no);
                  fe_values_p_T.reinit (cell, face_no);

                  // Evaluate the material model on the cell face.
                  // For this, we need the finite element that contains the solution.
                  MaterialModel::MaterialModelInputs<dim> in(fe_values_p_T, phase_boundary_cell, this->introspection(), this->get_solution());
                  MaterialModel::MaterialModelOutputs<dim> out(fe_values_p_T.n_quadrature_points, this->n_compositional_fields());
                  in.requested_properties = MaterialModel::MaterialProperties::density;
                  this->get_material_model().evaluate(in, out);

                  // Get the global numbers of the local DoFs of the mesh deformation cell
                  phase_boundary_cell->face(face_no)->get_dof_indices (face_dof_indices);

                  // Extract the displacement and initial topography values
                  fe_support_values[extract_vertical_displacements].get_function_values (displacements, displacement_values);
                  fe_support_values[extract_initial_topography].get_function_values (initial_topography, initial_topography_values);

                  // Loop over the quadrature points of the current face
                  for (unsigned int i=0; i<face_dof_indices.size(); ++i)
                    {
                      // Given the face dof, we get the component and overall cell dof index.
                      const std::pair<unsigned int, unsigned int> component_index = mesh_deformation_dof_handler.get_fe().face_system_to_component_index(i);
                      const unsigned int component = component_index.first;
                      const unsigned int support_index = component_index.second;

                      // Given the face dof, we get the component and overall cell dof index.
                      const std::pair<unsigned int, unsigned int> component_index_p_T = this->get_fe().face_system_to_component_index(i);
                      const unsigned int support_index_p_T = component_index_p_T.second;

                      // Get the gravity vector to compute the outward direction of displacement
                      const Point<dim> point = fe_support_values.quadrature_point(support_index);
                      Tensor<1,dim> direction = -(this->get_gravity_model().gravity_vector(point));

                      // Normalize direction vector
                      if (direction.norm() > 0.0)
                        direction *= 1./direction.norm();
                      else
                        AssertThrow(direction.norm() > 0.0,
                                    ExcMessage("Gravity must be non-zero in models with phase boundary."));

                      // Compute the total displacement in the gravity direction,
                      // i.e. the initial topography + any additional mesh displacement.
                      //const double delta_p = (in.pressure[support_index_p_T] - phase_transition_pressure
                      const double delta_p = (this->get_adiabatic_conditions().pressure(point) - phase_transition_pressure
                                              + clapeyron_slope * (in.temperature[support_index_p_T] - phase_transition_temperature));
                      const double delta_r = delta_p / (out.densities[support_index_p_T] * this->get_gravity_model().gravity_vector(point).norm());

                      const unsigned int dimension_of_dof = mesh_deformation_dof_handler.get_fe().system_to_component_index(support_index).first;

                      const double new_surface = displacement_values[support_index][dimension_of_dof]
                                                 + initial_topography_values[support_index][dimension_of_dof]
                                                 + delta_r * direction[component];

                      phase_boundary[face_dof_indices[i]] = new_surface;
                    }
                }
          // Make sure we have the right cell for the FE that contains temperature and pressure.
          ++cell;
        }

      phase_boundary.compress(VectorOperation::insert);

      // The phase_boundary vector contains the new displacements, but we need to return a velocity.
      // Therefore, we compute v=d_displacement/d_t.
      // d_displacement are the new mesh node locations
      // minus the old locations, which are initial_topography + displacements.
      LinearAlgebra::Vector velocity(mesh_locally_owned, mesh_locally_relevant, this->get_mpi_communicator());
      velocity = phase_boundary;
      velocity -= initial_topography;
      velocity -= displacements;

      // The velocity
      if (this->get_timestep() > 0.)
        velocity /= this->get_timestep();
      else
        AssertThrow(false, ExcZero());

      output = velocity;
    }



    template <int dim>
    Tensor<1,dim>
    PhaseBoundary<dim>::compute_initial_deformation_on_boundary(const types::boundary_id /*boundary_indicator*/,
                                                                const Point<dim> &position) const
    {
      // Compute initial phase boundary
      const double topo = 0.0;
      const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(position);

      Tensor<1,dim> topography_direction;
      if (gravity.norm() > 0.0)
        topography_direction = -gravity / gravity.norm();

      return topo * topography_direction;
    }



    /**
     * A function that creates constraints for the velocity of certain mesh
     * vertices (e.g. the surface vertices) for a specific boundary.
     * The calling class will respect
     * these constraints when computing the new vertex positions.
     */
    template <int dim>
    void
    PhaseBoundary<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                                 AffineConstraints<double> &mesh_velocity_constraints,
                                                                 const std::set<types::boundary_id> &boundary_id) const
    {
      LinearAlgebra::Vector boundary_velocity;

      const IndexSet &mesh_locally_owned = mesh_deformation_dof_handler.locally_owned_dofs();
      const IndexSet mesh_locally_relevant = DoFTools::extract_locally_relevant_dofs (mesh_deformation_dof_handler);
      boundary_velocity.reinit(mesh_locally_owned, mesh_locally_relevant,
                               this->get_mpi_communicator());

      // Determine the phase boundary surface based on the current temperature
      // and pressure
      phase_boundary(mesh_deformation_dof_handler, mesh_locally_owned,
                     mesh_locally_relevant, boundary_velocity, boundary_id);

      // now insert the relevant part of the solution into the mesh constraints
      const IndexSet constrained_dofs =
        DoFTools::extract_boundary_dofs(mesh_deformation_dof_handler,
                                        ComponentMask(dim, true),
                                        boundary_id);

      for (unsigned int i = 0; i < constrained_dofs.n_elements();  ++i)
        {
          types::global_dof_index index = constrained_dofs.nth_index_in_set(i);
          if (mesh_velocity_constraints.can_store_line(index))
            if (mesh_velocity_constraints.is_constrained(index)==false)
              {
                mesh_velocity_constraints.add_line(index);
                mesh_velocity_constraints.set_inhomogeneity(index, boundary_velocity[index]);
              }
        }
    }



    template <int dim>
    bool
    PhaseBoundary<dim>::
    needs_surface_stabilization () const
    {
      return false;
    }



    template <int dim>
    void PhaseBoundary<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Phase boundary");
        {
          prm.declare_entry ("Phase transition temperature", "0",
                             Patterns::Double (0.),
                             "The value of the temperature at the phase boundary. "
                             "Units: \\si{\\kelvin}.");
          prm.declare_entry ("Phase transition pressure", "0",
                             Patterns::Double (0.),
                             "The value of the pressure at the phase boundary. "
                             "Units: \\si{\\pascal}.");
          prm.declare_entry ("Phase transition Clapeyron slope", "0",
                             Patterns::Double (),
                             "The value of the Clapeyron slope of the phase boundary. "
                             "Units: \\si{\\pascal\\per\\kelvin}.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void PhaseBoundary<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("Phase boundary");
        {
          phase_transition_temperature = prm.get_double ("Phase transition temperature");
          phase_transition_pressure = prm.get_double ("Phase transition pressure");
          clapeyron_slope = prm.get_double ("Phase transition Clapeyron slope");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(PhaseBoundary,
                                           "phase boundary",
                                           "TODO.")
  }
}
