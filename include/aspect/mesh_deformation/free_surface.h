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


#ifndef _aspect_mesh_deformation_free_surface_h
#define _aspect_mesh_deformation_free_surface_h

#include <aspect/mesh_deformation/interface.h>

#include <aspect/simulator_access.h>
#include <aspect/simulator/assemblers/interface.h>


namespace aspect
{
  using namespace dealii;

  namespace Assemblers
  {
    /**
     * Apply stabilization to a cell of the system matrix.  The
     * stabilization is only added to cells on a free surface.  The
     * scheme is based on that of Kaus et. al., 2010.  Called during
     * assembly of the system matrix.
     */
    template <int dim>
    class ApplyStabilization: public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        ApplyStabilization(const double stabilization_theta);

        virtual ~ApplyStabilization () {};

        virtual
        void
        execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                 internal::Assembly::CopyData::CopyDataBase<dim> &data) const;

      private:
        /**
         * Stabilization parameter for the free surface.  Should be between
         * zero and one. A value of zero means no stabilization.  See Kaus
         * et. al. 2010 for more details.
         */
        double free_surface_theta;
    };
  }

  namespace MeshDeformation
  {
    template<int dim>
    class FreeSurface : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        virtual void initialize();

        /**
         * Called by Simulator::set_assemblers() to allow the FreeSurfaceHandler
         * to register its assembler.
         */
        void set_assemblers(const SimulatorAccess<dim> &simulator_access,
                            aspect::Assemblers::Manager<dim> &assemblers) const;

        virtual
        void
        deformation_constraints(const DoFHandler<dim> &free_surface_dof_handler,
                                ConstraintMatrix &mesh_velocity_constraints) const;

        /**
         * Declare parameters for the free surface handling.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters for the free surface handling.
         */
        void parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Project the Stokes velocity solution onto the
         * free surface. Called by make_constraints()
         */
        void project_velocity_onto_boundary (const DoFHandler<dim> &free_surface_dof_handler,
                                             const IndexSet &mesh_locally_owned,
                                             const IndexSet &mesh_locally_relevant,
                                             LinearAlgebra::Vector &output) const;

        /**
         * Stabilization parameter for the free surface.  Should be between
         * zero and one. A value of zero means no stabilization.  See Kaus
         * et. al. 2010 for more details.
         */
        double free_surface_theta;

        /**
         * A struct for holding information about how to advect the free surface.
         */
        struct SurfaceAdvection
        {
          enum Direction { normal, vertical };
        };

        /**
         * Stores whether to advect the free surface in the normal direction
         * or the direction of the local vertical.
         */
        typename SurfaceAdvection::Direction advection_direction;
    };
  }
}


#endif
