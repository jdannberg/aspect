/*
  Copyright (C) 2012 - 2017 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_temperature_adiabatic_h
#define _aspect_initial_temperature_adiabatic_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    /**
     * A class that implements adiabatic initial conditions for the
     * temperature field and, optional, upper and lower thermal boundary
     * layers calculated using the half-space cooling model. The age of the
     * boundary layers are input parameters.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class Adiabatic : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Age of the upper thermal boundary layer at the surface of the
         * model. If set to zero, no boundary layer will be present in the
         * model.
         */
        double age_top_boundary_layer;
        /* Age of the lower thermal boundary layer. */
        double age_bottom_boundary_layer;

        /**
         * Radius (in m) of the initial temperature perturbation at the bottom
         * of the model domain.
         */
        double radius;
        /**
         * Amplitude (in K) of the initial temperature perturbation at the
         * bottom of the model domain.
         */
        double amplitude;
        /*
         * Position of the initial temperature perturbation (in the
         * center or at the boundary of the model domain).
         */
        std::string perturbation_position;

        /*
         * Deviation from adiabaticity in a subadiabatic mantle
         * temperature profile. 0 for an adiabatic temperature
         * profile.
         */
        double subadiabaticity;

        /*
         * Whether the initial temperature profile should take into account
         * the effect of latent heat of melting, which increases the temperature
         * gradient.
         */
        bool include_latent_heat_of_melt;

        /*
         * Whether to use a function (instead of the constant value given in
         * age_top_boundary_layer to prescribe the age of the top boundary
         * layer.
         */
        bool use_age_function_for_top_boundary;

        /**
         * A function object representing the compositional fields that will
         * be used as a reference profile for calculating the thermal
         * diffusivity. The function depends only on depth.
         */
        std_cxx11::unique_ptr<Functions::ParsedFunction<1> > function;

        /**
        * A function object representing the age of the lithosphere. The
        * function depends only on the horizontal coordinate(s), and can be
        * used to specify different ages of the lithosphere for different places
        * in the model.
        */
        std_cxx11::unique_ptr<Functions::ParsedFunction<dim-1> > age_function;
    };
  }
}


#endif
