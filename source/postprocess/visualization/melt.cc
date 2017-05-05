/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/melt.h>
#include <aspect/melt.h>
#include <aspect/utilities.h>
#include <aspect/simulator.h>
#include <aspect/material_model/interface.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/data_out.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {

      template <int dim>
      MeltMaterialProperties<dim>::
      MeltMaterialProperties ()
        :
        DataPostprocessor<dim> ()
      {}

      template <int dim>
      std::vector<std::string>
      MeltMaterialProperties<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;

        for (unsigned int i=0; i<property_names.size(); ++i)
          if (property_names[i] == "fluid density gradient")
            for (unsigned int i=0; i<dim; ++i)
              solution_names.push_back ("fluid_density_gradient");
          else
            {
              solution_names.push_back(property_names[i]);
              std::replace(solution_names.back().begin(),solution_names.back().end(),' ', '_');
            }

        return solution_names;
      }

      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      MeltMaterialProperties<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
        for (unsigned int i=0; i<property_names.size(); ++i)
          {
            if (property_names[i] == "fluid density gradient")
              {
                for (unsigned int c=0; c<dim; ++c)
                  interpretation.push_back (DataComponentInterpretation::component_is_part_of_vector);
              }
            else
              interpretation.push_back (DataComponentInterpretation::component_is_scalar);
          }

        return interpretation;
      }

      template <int dim>
      UpdateFlags
      MeltMaterialProperties<dim>::
      get_needed_update_flags () const
      {
        return update_gradients | update_values  | update_q_points;
      }

      template <int dim>
      void
      MeltMaterialProperties<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        AssertThrow(this->include_melt_transport()==true,
                    ExcMessage("'Include melt transport' has to be on when using melt transport postprocessors."));

        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points, this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points, this->n_compositional_fields());
        MeltHandler<dim>::create_material_model_outputs(out);

        in.position = input_data.evaluation_points;
        Point<dim> mid_point;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            in.pressure[q] = input_data.solution_values[q][this->introspection().component_indices.pressure];
            in.temperature[q] = input_data.solution_values[q][this->introspection().component_indices.temperature];
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = solution_gradients[q][d];
            in.strain_rate[q] = symmetrize (grad_u);


            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = input_data.solution_values[q][this->introspection().component_indices.compositional_fields[c]];

            mid_point += evaluation_points[q]/n_quadrature_points;
          }

        typename DoFHandler<dim>::active_cell_iterator cell;
        cell = (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_dof_handler(), mid_point)).first;
        in.cell = &cell;

        this->get_material_model().evaluate(in, out);
        MaterialModel::MeltOutputs<dim> *melt_outputs = out.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
        AssertThrow(melt_outputs != NULL,
                    ExcMessage("Need MeltOutputs from the material model for computing the melt properties."));


        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            unsigned output_index = 0;
            for (unsigned int i=0; i<property_names.size(); ++i, ++output_index)
              {
                if (property_names[i] == "compaction viscosity")
                  computed_quantities[q][output_index] = melt_outputs->compaction_viscosities[q];
                else if (property_names[i] == "fluid viscosity")
                  computed_quantities[q][output_index] = melt_outputs->fluid_viscosities[q];
                else if (property_names[i] == "permeability")
                  computed_quantities[q][output_index] = melt_outputs->permeabilities[q];
                else if (property_names[i] == "fluid density")
                  computed_quantities[q][output_index] = melt_outputs->fluid_densities[q];
                else if (property_names[i] == "fluid density gradient")
                  {
                    for (unsigned int k=0; k<dim; ++k, ++output_index)
                      {
                        computed_quantities[q][output_index] = melt_outputs->fluid_density_gradients[q][k];
                      }
                    --output_index;
                  }
                else
                  AssertThrow(false, ExcNotImplemented());
              }
          }
      }

      template <int dim>
      void
      MeltMaterialProperties<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Melt material properties");
            {
              const std::string pattern_of_names
                = "compaction viscosity|fluid viscosity|permeability|"
                  "fluid density|fluid density gradient";

              prm.declare_entry("List of properties",
                                "compaction viscosity,permeability",
                                Patterns::MultipleSelection(pattern_of_names),
                                "A comma separated list of melt properties that should be "
                                "written whenever writing graphical output. "
                                "The following material properties are available:\n\n"
                                +
                                pattern_of_names);
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      MeltMaterialProperties<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Melt material properties");
            {
              property_names = Utilities::split_string_list(prm.get ("List of properties"));
              AssertThrow(Utilities::has_unique_entries(property_names),
                          ExcMessage("The list of strings for the parameter "
                                     "'Postprocess/Visualization/Melt material properties/List of properties' contains entries more than once. "
                                     "This is not allowed. Please check your parameter file."));
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MeltMaterialProperties,
                                                  "melt material properties",
                                                  "A visualization output object that generates output "
                                                  "for melt related properties of the material model.")

    }
  }
}
