/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


#include <aspect/mesh_deformation/free_surface.h>
#include <aspect/simulator.h>
#include <aspect/global.h>
#include <aspect/simulator/assemblers/interface.h>
#include <aspect/melt.h>


namespace aspect
{
  namespace Assemblers
  {
    template <int dim>
    ApplyStabilization<dim>::ApplyStabilization(const double stabilization_theta)
      :
      free_surface_theta(stabilization_theta)
    {}

    template <int dim>
    void
    ApplyStabilization<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>       &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim>      &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      if (!this->get_parameters().free_surface_enabled)
        return;

      if (this->get_parameters().include_melt_transport)
        {
          this->get_melt_handler().apply_free_surface_stabilization_with_melt (this->get_free_surface_handler().get_stabilization_term(),
                                                                               scratch.cell,
                                                                               scratch,
                                                                               data);
          return;
        }

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const typename DoFHandler<dim>::active_cell_iterator cell (&this->get_triangulation(),
                                                                 scratch.finite_element_values.get_cell()->level(),
                                                                 scratch.finite_element_values.get_cell()->index(),
                                                                 &this->get_dof_handler());

      const unsigned int n_face_q_points = scratch.face_finite_element_values.n_quadrature_points;
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();

      // only apply on free surface faces
      if (cell->at_boundary() && cell->is_locally_owned())
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
          if (cell->face(face_no)->at_boundary())
            {
              const types::boundary_id boundary_indicator
                = cell->face(face_no)->boundary_id();

              if (this->get_parameters().free_surface_boundary_indicators.find(boundary_indicator)
                  == this->get_parameters().free_surface_boundary_indicators.end())
                continue;

              scratch.face_finite_element_values.reinit(cell, face_no);

              this->compute_material_model_input_values (this->get_solution(),
                                                         scratch.face_finite_element_values,
                                                         cell,
                                                         false,
                                                         scratch.face_material_model_inputs);

              this->get_material_model().evaluate(scratch.face_material_model_inputs, scratch.face_material_model_outputs);

              for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
                {
                  for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
                    {
                      if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                        {
                          scratch.phi_u[i_stokes] = scratch.face_finite_element_values[introspection.extractors.velocities].value(i, q_point);
                          ++i_stokes;
                        }
                      ++i;
                    }

                  const Tensor<1,dim>
                  gravity = this->get_gravity_model().gravity_vector(scratch.face_finite_element_values.quadrature_point(q_point));
                  double g_norm = gravity.norm();

                  // construct the relevant vectors
                  const Tensor<1,dim> n_hat = scratch.face_finite_element_values.normal_vector(q_point);
                  const Tensor<1,dim> g_hat = (g_norm == 0.0 ? Tensor<1,dim>() : gravity/g_norm);

                  const double pressure_perturbation = scratch.face_material_model_outputs.densities[q_point] *
                                                       this->get_timestep() *
                                                       free_surface_theta *
                                                       g_norm;

                  // see Kaus et al 2010 for details of the stabilization term
                  for (unsigned int i=0; i< stokes_dofs_per_cell; ++i)
                    for (unsigned int j=0; j< stokes_dofs_per_cell; ++j)
                      {
                        // The fictive stabilization stress is (phi_u[i].g)*(phi_u[j].n)
                        const double stress_value = -pressure_perturbation*
                                                    (scratch.phi_u[i]*g_hat) * (scratch.phi_u[j]*n_hat)
                                                    *scratch.face_finite_element_values.JxW(q_point);

                        data.local_matrix(i,j) += stress_value;
                      }
                }
            }
    }
  }

  namespace MeshDeformation
  {
    template <int dim>
    void
    FreeSurface<dim>::initialize ()
    {
      this->get_signals().set_assemblers.connect(std_cxx11::bind(&FreeSurface<dim>::set_assemblers,
                                                                 std_cxx11::cref(*this),
                                                                 std_cxx11::_1,
                                                                 std_cxx11::_2));
    }



    template <int dim>
    void FreeSurface<dim>::set_assemblers(const SimulatorAccess<dim> &,
                                          aspect::Assemblers::Manager<dim> &assemblers) const
    {
      aspect::Assemblers::ApplyStabilization<dim> *surface_stabilization
        = new aspect::Assemblers::ApplyStabilization<dim>(free_surface_theta);

      assemblers.stokes_system.push_back(
        std_cxx11::unique_ptr<aspect::Assemblers::ApplyStabilization<dim> > (surface_stabilization));

      // Note that we do not want face_material_model_data, because we do not
      // connect to a face assembler. We instead connect to a normal assembler,
      // and compute our own material_model_inputs in apply_stabilization
      // (because we want to use the solution instead of the current_linearization_point
      // to compute the material properties).
      assemblers.stokes_system_assembler_on_boundary_face_properties.needed_update_flags |= (update_values  |
          update_gradients |
          update_quadrature_points |
          update_normal_vectors |
          update_JxW_values);
    }



    template <int dim>
    void FreeSurface<dim>::project_velocity_onto_boundary(const DoFHandler<dim> &free_surface_dof_handler,
                                                          const IndexSet &mesh_locally_owned,
                                                          const IndexSet &mesh_locally_relevant,
                                                          LinearAlgebra::Vector &output) const
    {
      // TODO: should we use the extrapolated solution?

      // stuff for iterating over the mesh
      QGauss<dim-1> face_quadrature(free_surface_dof_handler.get_fe().degree+1);
      UpdateFlags update_flags = UpdateFlags(update_values | update_quadrature_points
                                             | update_normal_vectors | update_JxW_values);
      FEFaceValues<dim> fs_fe_face_values (this->get_mapping(), free_surface_dof_handler.get_fe(), face_quadrature, update_flags);
      FEFaceValues<dim> fe_face_values (this->get_mapping(), this->get_fe(), face_quadrature, update_flags);
      const unsigned int n_face_q_points = fe_face_values.n_quadrature_points,
                         dofs_per_cell = fs_fe_face_values.dofs_per_cell;

      // stuff for assembling system
      std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
      Vector<double> cell_vector (dofs_per_cell);
      FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

      // stuff for getting the velocity values
      std::vector<Tensor<1,dim> > velocity_values(n_face_q_points);

      // set up constraints
      ConstraintMatrix mass_matrix_constraints(mesh_locally_relevant);
      DoFTools::make_hanging_node_constraints(free_surface_dof_handler, mass_matrix_constraints);

      typedef std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundary_pairs;
      periodic_boundary_pairs pbp = this->get_geometry_model().get_periodic_boundary_pairs();
      for (periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
        DoFTools::make_periodicity_constraints(free_surface_dof_handler,
                                               (*p).first.first, (*p).first.second, (*p).second, mass_matrix_constraints);

      mass_matrix_constraints.close();

      // set up the matrix
      LinearAlgebra::SparseMatrix mass_matrix;
#ifdef ASPECT_USE_PETSC
      LinearAlgebra::DynamicSparsityPattern sp(mesh_locally_relevant);

#else
      TrilinosWrappers::SparsityPattern sp (mesh_locally_owned,
                                            mesh_locally_owned,
                                            mesh_locally_relevant,
                                            this->get_mpi_communicator());
#endif
      DoFTools::make_sparsity_pattern (free_surface_dof_handler, sp, mass_matrix_constraints, false,
                                       Utilities::MPI::this_mpi_process(this->get_mpi_communicator()));
#ifdef ASPECT_USE_PETSC
      SparsityTools::distribute_sparsity_pattern(sp,
                                                 free_surface_dof_handler.n_locally_owned_dofs_per_processor(),
                                                 this->get_mpi_communicator(), mesh_locally_relevant);

      sp.compress();
      mass_matrix.reinit (mesh_locally_owned, mesh_locally_owned, sp, this->get_mpi_communicator());
#else
      sp.compress();
      mass_matrix.reinit (sp);
#endif

      FEValuesExtractors::Vector extract_vel(0);

      // make distributed vectors.
      LinearAlgebra::Vector rhs, dist_solution;
      rhs.reinit(mesh_locally_owned, this->get_mpi_communicator());
      dist_solution.reinit(mesh_locally_owned, this->get_mpi_communicator());

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(), endc= this->get_dof_handler().end();
      typename DoFHandler<dim>::active_cell_iterator
      fscell = free_surface_dof_handler.begin_active();

      for (; cell!=endc; ++cell, ++fscell)
        if (cell->at_boundary() && cell->is_locally_owned())
          for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
            if (cell->face(face_no)->at_boundary())
              {
                const types::boundary_id boundary_indicator
                  = cell->face(face_no)->boundary_id();
                if (this->get_parameters().free_surface_boundary_indicators.find(boundary_indicator)
                    == this->get_parameters().free_surface_boundary_indicators.end())
                  continue;

                fscell->get_dof_indices (cell_dof_indices);
                fs_fe_face_values.reinit (fscell, face_no);
                fe_face_values.reinit (cell, face_no);
                fe_face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), velocity_values);

                cell_vector = 0;
                cell_matrix = 0;
                for (unsigned int point=0; point<n_face_q_points; ++point)
                  {
                    // Select the direction onto which to project the velocity solution
                    Tensor<1,dim> direction;
                    if ( advection_direction == SurfaceAdvection::normal ) // project onto normal vector
                      direction = fs_fe_face_values.normal_vector(point);
                    else if ( advection_direction == SurfaceAdvection::vertical ) // project onto local gravity
                      direction = this->get_gravity_model().gravity_vector(fs_fe_face_values.quadrature_point(point));
                    else
                      AssertThrow(false, ExcInternalError());

                    direction *= ( direction.norm() > 0.0 ? 1./direction.norm() : 0.0 );

                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                      {
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                          {
                            cell_matrix(i,j) += (fs_fe_face_values[extract_vel].value(j,point) *
                                                 fs_fe_face_values[extract_vel].value(i,point) ) *
                                                fs_fe_face_values.JxW(point);
                          }

                        cell_vector(i) += (fs_fe_face_values[extract_vel].value(i,point) * direction)
                                          * (velocity_values[point] * direction)
                                          * fs_fe_face_values.JxW(point);
                      }
                  }

                mass_matrix_constraints.distribute_local_to_global (cell_matrix, cell_vector,
                                                                    cell_dof_indices, mass_matrix, rhs, false);
              }

      rhs.compress (VectorOperation::add);
      mass_matrix.compress(VectorOperation::add);

      // Jacobi seems to be fine here.  Other preconditioners (ILU, IC) run into troubles
      // because the matrix is mostly empty, since we don't touch internal vertices.
      LinearAlgebra::PreconditionJacobi preconditioner_mass;
      preconditioner_mass.initialize(mass_matrix);

      SolverControl solver_control(5*rhs.size(), this->get_parameters().linear_stokes_solver_tolerance*rhs.l2_norm());
      SolverCG<LinearAlgebra::Vector> cg(solver_control);
      cg.solve (mass_matrix, dist_solution, rhs, preconditioner_mass);

      mass_matrix_constraints.distribute (dist_solution);
      output = dist_solution;
    }


    template <int dim>
    void
    FreeSurface<dim>::deformation_constraints(const DoFHandler<dim> &free_surface_dof_handler,
                                              ConstraintMatrix &mesh_constraints) const
    {
      // For the free surface indicators we constrain the displacement to be v.n
      LinearAlgebra::Vector boundary_velocity;

      const IndexSet mesh_locally_owned = free_surface_dof_handler.locally_owned_dofs();
      IndexSet mesh_locally_relevant;
      DoFTools::extract_locally_relevant_dofs (free_surface_dof_handler,
                                               mesh_locally_relevant);
      boundary_velocity.reinit(mesh_locally_owned, mesh_locally_relevant, this->get_mpi_communicator());
      project_velocity_onto_boundary(free_surface_dof_handler,mesh_locally_owned,mesh_locally_relevant,boundary_velocity);

      // now insert the relevant part of the solution into the mesh constraints
      IndexSet constrained_dofs;
      DoFTools::extract_boundary_dofs(free_surface_dof_handler, ComponentMask(dim, true),
                                      constrained_dofs, this->get_parameters().free_surface_boundary_indicators);
      for (unsigned int i = 0; i < constrained_dofs.n_elements();  ++i)
        {
          types::global_dof_index index = constrained_dofs.nth_index_in_set(i);
          if (mesh_constraints.can_store_line(index))
            if (mesh_constraints.is_constrained(index)==false)
              {
                mesh_constraints.add_line(index);
                mesh_constraints.set_inhomogeneity(index, boundary_velocity[index]);
              }
        }
    }



    template <int dim>
    void FreeSurface<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Free surface");
      {
        prm.declare_entry("Free surface stabilization theta", "0.5",
                          Patterns::Double(0,1),
                          "Theta parameter described in Kaus et. al. 2010. "
                          "An unstabilized free surface can overshoot its "
                          "equilibrium position quite easily and generate "
                          "unphysical results.  One solution is to use a "
                          "quasi-implicit correction term to the forces near the "
                          "free surface.  This parameter describes how much "
                          "the free surface is stabilized with this term, "
                          "where zero is no stabilization, and one is fully "
                          "implicit.");
        prm.declare_entry("Surface velocity projection", "normal",
                          Patterns::Selection("normal|vertical"),
                          "After each time step the free surface must be "
                          "advected in the direction of the velocity field. "
                          "Mass conservation requires that the mesh velocity "
                          "is in the normal direction of the surface. However, "
                          "for steep topography or large curvature, advection "
                          "in the normal direction can become ill-conditioned, "
                          "and instabilities in the mesh can form. Projection "
                          "of the mesh velocity onto the local vertical direction "
                          "can preserve the mesh quality better, but at the "
                          "cost of slightly poorer mass conservation of the "
                          "domain.");
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void FreeSurface<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Free surface");
      {
        free_surface_theta = prm.get_double("Free surface stabilization theta");
        std::string advection_dir = prm.get("Surface velocity projection");

        if ( advection_dir == "normal")
          advection_direction = SurfaceAdvection::normal;
        else if ( advection_dir == "vertical")
          advection_direction = SurfaceAdvection::vertical;
        else
          AssertThrow(false, ExcMessage("The surface velocity projection must be ``normal'' or ``vertical''."));

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
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(FreeSurface,
                                           "free surface",
                                           "")
  }
}
