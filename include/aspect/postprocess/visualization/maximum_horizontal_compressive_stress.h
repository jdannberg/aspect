/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef __aspect__postprocess_visualization_maximum_horizontal_compressive_stress_h
#define __aspect__postprocess_visualization_maximum_horizontal_compressive_stress_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_postprocessor.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class that computes a field of horizontal vectors that
       * represent the direction of maximal horizontal compressive
       * stress. For an exact definition, see the documentation of
       * this plugin in the manual.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class MaximumHorizontalCompressiveStress
        : public DataPostprocessor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          virtual
          void
          compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                             const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                             const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                             const std::vector<Point<dim> >                  &normals,
                                             const std::vector<Point<dim> >                  &evaluation_points,
                                             std::vector<Vector<double> >                    &computed_quantities) const;

          /**
           * Return the vector of strings describing the names of the computed
           * quantities. Given the purpose of this class, this is a vector
           * with entries all equal to the name of the plugin.
           */
          virtual std::vector<std::string> get_names () const;

          /**
           * This functions returns information about how the individual
           * components of output files that consist of more than one data set
           * are to be interpreted. The returned value is
           * DataComponentInterpretation::component_is_scalar repeated
           * SymmetricTensor::n_independent_components times. (These
           * components should really be part of a symmetric tensor, but
           * deal.II does not allow marking components as such.)
           */
          virtual
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const;

          /**
           * Return which data has to be provided to compute the derived
           * quantities. The flags returned here are the ones passed to the
           * constructor of this class.
           */
          virtual UpdateFlags get_needed_update_flags () const;
      };
    }
  }
}

#endif
