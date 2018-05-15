/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_composition_porosity_h
#define _aspect_initial_composition_porosity_h

#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements initial conditions for the porosity field
     * by computing the equilibrium melt fraction for the given initial
     * condition and reference pressure profile. Note that this plugin only
     * works if there is a compositional field called 'porosity', and the
     * used material model implements the 'MeltFractionModel' interface.
     * All compositional fields except porosity are not changed by this plugin.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class Porosity : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        virtual
        double initial_composition (const Point<dim> &position,
                                    const unsigned int compositional_index) const;

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

      private:
        /**
         * A function object representing the compositional fields that this
         * initial composition should apply to.
         */
        std::vector<std::string> names_of_melt_fields;
    };
  }
}


#endif
