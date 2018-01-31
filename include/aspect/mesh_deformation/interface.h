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


#ifndef _aspect_mesh_deformation_interface_h
#define _aspect_mesh_deformation_interface_h

#include <aspect/plugins.h>
#include <aspect/simulator.h>

namespace aspect
{
  using namespace dealii;

  /**
   * A namespace that contains everything that is related to the deformation
   * of the mesh vertices over time.
   */
  namespace MeshDeformation
  {
    /**
     * A base class for mesh deformation plugins. Each of these plugins should
     * implement a function that determines the displacement for mesh vertices
     * given the current position of the mesh vertex.
     */
    template<int dim>
    class Interface
    {
      public:
        /**
         * Destructor. Made virtual to enforce that derived classes also have
         * virtual destructors.
         */
        virtual ~Interface();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        virtual void initialize ();

        /**
         * A function that is called at the beginning of each time step and
         * that allows the implementation to update internal data structures.
         * This is useful, for example, if you have mesh deformation that
         * depends on time, or on the solution of the previous step.
         *
         * The default implementation of this function does nothing.
         */
        virtual void update();

        virtual
        void
        deformation_constraints(const DoFHandler<dim> &free_surface_dof_handler,
                                ConstraintMatrix &mesh_constraints) const = 0;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };

    template<int dim>
    class FreeSurfaceHandler: public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize the free surface handler, allowing it to read in
         * relevant parameters as well as giving it a reference to the
         * Simulator that owns it, since it needs to make fairly extensive
         * changes to the internals of the simulator.
         */
        FreeSurfaceHandler(Simulator<dim> &simulator);

        /**
         * Destructor for the free surface handler.
         */
        ~FreeSurfaceHandler();

        void initialize();

        void update();

        /**
         * The main execution step for the free surface implementation. This
         * computes the motion of the free surface, moves the boundary nodes
         * accordingly, redistributes the internal nodes in order to
         * preserve mesh regularity, and calculates the Arbitrary-
         * Lagrangian-Eulerian correction terms for advected quantities.
         */
        void execute();

        /**
         * Allocates and sets up the members of the FreeSurfaceHandler. This
         * is called by Simulator<dim>::setup_dofs()
         */
        void setup_dofs();

        /**
         * Declare parameters for the free surface handling.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters for the free surface handling.
         */
        void parse_parameters (ParameterHandler &prm);

        /**
         * A function that is used to register mesh deformation objects in such
         * a way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * plugins are implemented to register these plugins, rather than also
         * having to modify the Manager class by adding the new initial
         * temperature plugin class.
         *
         * @param name A string that identifies the mesh deformation model
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this mesh deformation model
         * wants to read from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this mesh deformation model.
         */
        static
        void
        register_mesh_deformation
        (const std::string &name,
         const std::string &description,
         void (*declare_parameters_function) (ParameterHandler &),
         Interface<dim> *(*factory_function) ());


        /**
         * Return a list of names of all mesh deformation models currently
         * used in the computation, as specified in the input file.
         */
        const std::vector<std::string> &
        get_active_mesh_deformation_names () const;

        /**
         * Return a list of pointers to all mesh deformation models
         * currently used in the computation, as specified in the input file.
         */
        const std::vector<std_cxx11::shared_ptr<Interface<dim> > > &
        get_active_mesh_deformation_models () const;

        /**
         * Go through the list of all mesh deformation models that have been selected in
         * the input file (and are consequently currently active) and see if one
         * of them has the desired type specified by the template argument. If so,
         * return a pointer to it. If no mesh deformation model is active
         * that matches the given type, return a NULL pointer.
         */
        template <typename MeshDeformationType>
        MeshDeformationType *
        find_mesh_deformation_model () const;

        /**
         * For the current plugin subsystem, write a connection graph of all of the
         * plugins we know about, in the format that the
         * programs dot and neato understand. This allows for a visualization of
         * how all of the plugins that ASPECT knows about are interconnected, and
         * connect to other parts of the ASPECT code.
         *
         * @param output_stream The stream to write the output to.
         */
        static
        void
        write_plugin_graph (std::ostream &output_stream);

        /**
         * Exception.
         */
        DeclException1 (ExcMeshDeformationNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered mesh deformation objects.");

      private:
        /**
         * A list of mesh deformation objects that have been requested in the
         * parameter file.
         */
        std::vector<std_cxx11::shared_ptr<Interface<dim> > > mesh_deformation_objects;

        /**
         * A list of names of mesh deformation objects that have been requested
         * in the parameter file.
         */
        std::vector<std::string> model_names;

        /**
         * Set the boundary conditions for the solution of the elliptic
         * problem, which computes the displacements of the internal
         * vertices so that the mesh does not become too distorted due to
         * motion of the free surface.  Velocities of vertices on the free
         * surface are set to be the normal of the Stokes velocity solution
         * projected onto that surface.  Velocities of vertices on free-slip
         * boundaries are constrained to be tangential to those boundaries.
         * Velocities of vertices on no-slip boundaries are set to be zero.
         */
        void make_constraints ();

        /**
         * Solve vector Laplacian equation for internal mesh displacements.
         */
        void compute_mesh_displacements ();

        /**
         * Calculate the velocity of the mesh for ALE corrections.
         */
        void interpolate_mesh_velocity ();

        /**
         * Reference to the Simulator object to which a FreeSurfaceHandler
         * instance belongs.
         */
        Simulator<dim> &sim;

        /**
        * Finite element for the free surface implementation, which is
        * used for tracking mesh deformation.
        */
        const FESystem<dim> free_surface_fe;

        /**
         * DoFHandler for the free surface implementation.
         */
        DoFHandler<dim> free_surface_dof_handler;

        /**
         * BlockVector which stores the mesh velocity.
         * This is used for ALE corrections.
         */
        LinearAlgebra::BlockVector mesh_velocity;

        /**
         * Vector for storing the positions of the mesh vertices. This
         * is used for calculating the mapping from the reference cell to
         * the position of the cell in the deformed mesh. This must be
         * redistributed upon mesh refinement.
         */
        LinearAlgebra::Vector mesh_displacements;

        /**
         * Vector for storing the mesh velocity in the free surface finite
         * element space, which is, in general, not the same finite element
         * space as the Stokes system. This is used for interpolating
         * the mesh velocity in the free surface finite element space onto
         * the velocity in the Stokes finite element space, which is then
         * used for making the ALE correction in the advection equations.
         */
        LinearAlgebra::Vector fs_mesh_velocity;

        /**
         * IndexSet for the locally owned DoFs for the mesh system
         */
        IndexSet mesh_locally_owned;

        /**
         * IndexSet for the locally relevant DoFs for the mesh system
         */
        IndexSet mesh_locally_relevant;

        /**
         * Storage for the mesh displacement constraints for solving the
         * elliptic problem
         */
        ConstraintMatrix mesh_displacement_constraints;

        /**
         * Storage for the mesh vertex constraints for keeping the mesh conforming
         * upon redistribution.
         */
        ConstraintMatrix mesh_vertex_constraints;

        /**
         * A struct for holding information about how to advect the free surface.
         */
        struct SurfaceVelocity
        {
          enum Type { free_surface, function };
        };

        /**
         * Stores whether to advect the free surface in the normal direction
         * or the direction of the local vertical.
         */
        typename SurfaceVelocity::Type surface_velocity;

        /**
         * A set of boundary indicators that denote those boundaries that are
         * allowed to move their mesh tangential to the boundary. All
         * boundaries that have tangential material velocity boundary
         * conditions are in this set by default, but it can be extended by
         * open boundaries, boundaries with traction boundary conditions, or
         * boundaries with prescribed material velocities if requested in
         * the parameter file.
         */
        std::set<types::boundary_id> tangential_mesh_boundary_indicators;

        friend class Simulator<dim>;
        friend class SimulatorAccess<dim>;
    };



    template <int dim>
    template <typename MeshDeformationType>
    inline
    MeshDeformationType *
    FreeSurfaceHandler<dim>::find_mesh_deformation_model () const
    {
      for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
           p = mesh_deformation_objects.begin();
           p != mesh_deformation_objects.end(); ++p)
        if (MeshDeformationType *x = dynamic_cast<MeshDeformationType *> ( (*p).get()) )
          return x;
      return NULL;
    }



    /**
     * Return a string that consists of the names of mesh deformation models that can
     * be selected. These names are separated by a vertical line '|' so
     * that the string can be an input to the deal.II classes
     * Patterns::Selection or Patterns::MultipleSelection.
     */
    template <int dim>
    std::string
    get_valid_model_names_pattern ();



    /**
     * Given a class name, a name, and a description for the parameter file
     * for a mesh deformation model, register it with the functions that can
     * declare their parameters and create these objects.
     *
     * @ingroup MeshDeformations
     */
#define ASPECT_REGISTER_MESH_DEFORMATION_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_MESH_DEFORMATION_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::MeshDeformation::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::MeshDeformation::FreeSurfaceHandler<2>::register_mesh_deformation, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::MeshDeformation::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::MeshDeformation::FreeSurfaceHandler<3>::register_mesh_deformation, \
                                name, description); \
  }
  }
}

#endif
