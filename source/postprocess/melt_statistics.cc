/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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



#include <aspect/postprocess/melt_statistics.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/melt_global.h>
#include <aspect/material_model/melt_simple.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    MeltStatistics<dim>::execute (TableHandler &statistics)
    {
      Assert (dynamic_cast<const MaterialModel::MeltGlobal<dim>*> (&this->get_material_model()) != 0 ||
              dynamic_cast<const MaterialModel::MeltSimple<dim>*> (&this->get_material_model()) != 0,
              ExcMessage ("This postprocessor can only be used with the melt simple "
                          "or melt global material model."));
      // TODO: this could easily be extended to also include the latent heat melt material model

      // create a quadrature formula based on the temperature element alone.
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());

      std::ostringstream output;
      output.precision(4);

      double local_melt_integral = 0.0;
      double local_min_melt = std::numeric_limits<double>::max();
      double local_max_melt = -std::numeric_limits<double>::max();

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      // compute the integral quantities by quadrature
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.temperature]
            .get_function_values (this->get_solution(),
                                  in.temperature);
            in.position = fe_values.get_quadrature_points();

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                double melt_fraction = 0.0;
                if (const MaterialModel::MeltGlobal<dim> *
                    melt_material_model = dynamic_cast <const MaterialModel::MeltGlobal<dim>*> (&this->get_material_model()))
                  melt_fraction = melt_material_model->melt_fraction(in.temperature[q],
                                                                     this->get_adiabatic_conditions().pressure(in.position[q]),
                                                                     0.0);
                else if (const MaterialModel::MeltSimple<dim> *
                         melt_material_model = dynamic_cast <const MaterialModel::MeltSimple<dim>*> (&this->get_material_model()))
                  melt_fraction = melt_material_model->melt_fraction(in.temperature[q],
                                                                     this->get_adiabatic_conditions().pressure(in.position[q]));

                local_melt_integral += melt_fraction * fe_values.JxW(q);
                local_min_melt       = std::min(local_min_melt, melt_fraction);
                local_max_melt       = std::max(local_max_melt, melt_fraction);
              }

          }

      const double global_melt_integral
        = Utilities::MPI::sum (local_melt_integral, this->get_mpi_communicator());
      double global_min_melt = 0;
      double global_max_melt = 0;

      // now do the reductions that are
      // min/max operations. do them in
      // one communication by multiplying
      // one value by -1
      {
        double local_values[2] = { -local_min_melt, local_max_melt };
        double global_values[2];

        Utilities::MPI::max (local_values, this->get_mpi_communicator(), global_values);

        global_min_melt = -global_values[0];
        global_max_melt = global_values[1];
      }


      // finally produce something for the statistics file
      statistics.add_value ("Minimal melt fraction",
                            global_min_melt);
      statistics.add_value ("Total melt fraction",
                            global_melt_integral);
      statistics.add_value ("Maximal melt fraction",
                            global_max_melt);

      // also make sure that the other columns filled by the this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Minimal melt fraction",
                                  "Total melt fraction",
                                  "Maximal melt fraction"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }

      output << global_min_melt << ", "
             << global_melt_integral << ", "
             << global_max_melt;

      return std::pair<std::string, std::string> ("Melt fraction min/total/max:",
                                                  output.str());

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MeltStatistics,
                                  "melt statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the melt fraction, averaged by volume. ")
  }
}
