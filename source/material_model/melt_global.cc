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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/material_model/melt_global.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>


namespace aspect
{
  namespace MaterialModel
  {
    namespace internal
    {

      class MeltFractionLookup
      {
        public:
          MeltFractionLookup(const std::string &filename,
                             const bool interpol,
                             const std::string &p_unit,
                             const std::string &T_unit,
                             const MPI_Comm &comm)
          {
            /* Initializing variables */
            interpolation = interpol;
            delta_press=-1.0;
            min_press=-1.0;
            delta_temp=-1.0;
            min_temp=-1.0;
            numtemp=0;
            numpress=0;

            std::string temp;
            // Read data from disk and distribute among processes
            std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

            getline(in, temp); // eat first line
            getline(in, temp); // eat next line
            getline(in, temp); // eat next line
            getline(in, temp); // eat next line

            // we may have to convert temperature and pressure to SI units
            double pressure_scaling_factor = numbers::signaling_nan<double>();
            double temperature_offset = numbers::signaling_nan<double>();

            if (p_unit == "GPa")
              pressure_scaling_factor = 1.e9;
            else if (p_unit == "Pa")
              pressure_scaling_factor = 1.0;
            else if (p_unit == "bar")
              pressure_scaling_factor = 1.e5;
            else if (p_unit == "kbar")
              pressure_scaling_factor = 1.e8;
            else
              AssertThrow (false,
                           ExcMessage ("The value <" + p_unit + "> for a pressure unit "
                                       "is not one of the valid values."));

            if (T_unit == "Kelvin")
              temperature_offset = 0.0;
            else if (p_unit == "Celsius")
              temperature_offset = -273.15;
            else
              AssertThrow (false,
                           ExcMessage ("The value <" + T_unit + "> for a temperature unit "
                                       "is not one of the valid values."));

            in >> min_temp;
            min_temp += temperature_offset;
            getline(in, temp);
            in >> delta_temp;
            getline(in, temp);
            in >> numtemp;
            getline(in, temp);
            getline(in, temp);
            in >> min_press;
            min_press *= pressure_scaling_factor;  // conversion from [GPa] to [Pa]
            getline(in, temp);
            in >> delta_press;
            delta_press *= pressure_scaling_factor; // conversion from [GPa] to [Pa]
            getline(in, temp);
            in >> numpress;
            getline(in, temp);
            getline(in, temp);
            getline(in, temp);

            Assert(min_temp >= 0.0, ExcMessage("Read in of Material header failed (mintemp)."));
            Assert(delta_temp > 0, ExcMessage("Read in of Material header failed (delta_temp)."));
            Assert(numtemp > 0, ExcMessage("Read in of Material header failed (numtemp)."));
            Assert(min_press >= 0, ExcMessage("Read in of Material header failed (min_press)."));
            Assert(delta_press > 0, ExcMessage("Read in of Material header failed (delta_press)."));
            Assert(numpress > 0, ExcMessage("Read in of Material header failed (numpress)."));


            max_temp = min_temp + (numtemp-1) * delta_temp;
            max_press = min_press + (numpress-1) * delta_press;

            peridotite_melt_fractions.reinit(numtemp,numpress);
            basalt_melt_fractions.reinit(numtemp,numpress);

            unsigned int i = 0;
            while (!in.eof())
              {
                double temp1,temp2;
                double peridotite,basalt;
                in >> temp1 >> temp2;
                in >> peridotite;
                if (in.fail())
                  {
                    in.clear();
                    peridotite = peridotite_melt_fractions[(i-1)%numtemp][(i-1)/numtemp];
                  }
                in >> basalt;
                if (in.fail())
                  {
                    in.clear();
                    basalt = basalt_melt_fractions[(i-1)%numtemp][(i-1)/numtemp];
                  }

                getline(in, temp);
                if (in.eof())
                  break;

                peridotite_melt_fractions[i%numtemp][i/numtemp]=peridotite;
                basalt_melt_fractions[i%numtemp][i/numtemp]=basalt;

                i++;
              }
            Assert(i==numtemp*numpress, ExcMessage("Melt fraction table size not consistent with header."));

          }


          double
          peridotite_melt_fraction(double temperature,
                                   double pressure) const
          {
            return value(temperature,pressure,peridotite_melt_fractions,interpolation);
          }


          double
          basalt_melt_fraction(double temperature,
                               double pressure) const
          {
            return value(temperature,pressure,basalt_melt_fractions,interpolation);
          }


          double
          value (const double temperature,
                 const double pressure,
                 const dealii::Table<2,
                 double> &values,
                 bool interpol) const
          {
            const double nT = get_nT(temperature);
            const unsigned int inT = static_cast<unsigned int>(nT);

            const double np = get_np(pressure);
            const unsigned int inp = static_cast<unsigned int>(np);

            Assert(inT<values.n_rows(), ExcMessage("Attempting to look up a temperature value with index greater than the number of rows."));
            Assert(inp<values.n_cols(), ExcMessage("Attempting to look up a pressure value with index greater than the number of columns."));

            if (!interpol)
              return values[inT][inp];
            else
              {
                // compute the coordinates of this point in the
                // reference cell between the data points
                const double xi = nT-inT;
                const double eta = np-inp;

                Assert ((0 <= xi) && (xi <= 1), ExcInternalError());
                Assert ((0 <= eta) && (eta <= 1), ExcInternalError());

                // use these coordinates for a bilinear interpolation
                return ((1-xi)*(1-eta)*values[inT][inp] +
                        xi    *(1-eta)*values[inT+1][inp] +
                        (1-xi)*eta    *values[inT][inp+1] +
                        xi    *eta    *values[inT+1][inp+1]);
              }
          }

        private:
          double get_nT(double temperature) const
          {
            temperature=std::max(min_temp, temperature);
            temperature=std::min(temperature, max_temp-delta_temp);
            Assert(temperature>=min_temp, ExcMessage("ASPECT found a temperature less than min_T."));
            Assert(temperature<=max_temp, ExcMessage("ASPECT found a temperature greater than max_T."));
            return (temperature-min_temp)/delta_temp;
          }

          double get_np(double pressure) const
          {
            pressure=std::max(min_press, pressure);
            pressure=std::min(pressure, max_press-delta_press);
            Assert(pressure>=min_press, ExcMessage("ASPECT found a pressure less than min_p."));
            Assert(pressure<=max_press, ExcMessage("ASPECT found a pressure greater than max_p."));
            return (pressure-min_press)/delta_press;
          }

          dealii::Table<2,double> peridotite_melt_fractions;
          dealii::Table<2,double> basalt_melt_fractions;

          double delta_press;
          double min_press;
          double max_press;
          double delta_temp;
          double min_temp;
          double max_temp;
          unsigned int numtemp;
          unsigned int numpress;
          bool interpolation;
      };
    }

    template <int dim>
    double
    MeltGlobal<dim>::
    reference_viscosity () const
    {
      return eta_0;
    }

    template <int dim>
    double
    MeltGlobal<dim>::
    reference_darcy_coefficient () const
    {
      // 0.01 = 1% melt
      return reference_permeability * std::pow(0.01,3.0) / eta_f;
    }

    template <int dim>
    bool
    MeltGlobal<dim>::
    is_compressible () const
    {
      return false;
    }


    template <int dim>
    void
    MeltGlobal<dim>::
    initialize()
    {
      melt_fraction_lookup.reset(new internal::MeltFractionLookup(data_directory+melt_fraction_file_name,
                                                                  interpolation,
                                                                  pressure_unit,
                                                                  temperature_unit,
                                                                  this->get_mpi_communicator()));
    }


    template <int dim>
    double
    MeltGlobal<dim>::
    melt_fraction (const double temperature,
                   const double pressure,
                   const double depletion) const
    {
      const double T_solidus  = surface_solidus
                                + pressure_solidus_change * pressure
                                + std::max(depletion_solidus_change * depletion, -200.0);
      const double T_liquidus = T_solidus + 500.0;

      double melt_fraction;
      if (temperature < T_solidus)
        melt_fraction = 0.0;
      else if (temperature > T_liquidus)
        melt_fraction = 1.0;
      else
        melt_fraction = (temperature - T_solidus) / (T_liquidus - T_solidus);

      return melt_fraction;
    }


    template <int dim>
    void
    MeltGlobal<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      if (read_melt_from_file)
        {
          // if we read the melt from a file, the melt fraction depends on the path we are on:
          // the solid-->melt phase transition uses a different diagram than the melt-->solid
          // phase trasition
          // As we do not know if material is melting or freezing at the moment, we will here
          // only output the melt fraction on the solid-->melt side of the path
          for (unsigned int q=0; q<in.temperature.size(); ++q)
            {
              double peridotite_melt = melt_fraction_lookup->peridotite_melt_fraction(in.temperature[q],in.pressure[q]);

              if (this->introspection().compositional_name_exists("peridotite"))
                {
                  const double basalt_melt = melt_fraction_lookup->basalt_melt_fraction(in.temperature[q],in.pressure[q]);
                  const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
                  const double peridotite_fraction = (1.0 - std::max(-in.composition[q][peridotite_idx],0.0));
                  const double basalt_fraction = std::min(std::max(-in.composition[q][peridotite_idx],0.0),1.0);

                  melt_fractions[q] = peridotite_fraction * peridotite_melt
                                      + basalt_fraction * basalt_melt;
                }
              else
                melt_fractions[q] = peridotite_melt;
            }
        }
      else
        {
          double depletion = 0.0;

          for (unsigned int q=0; q<in.temperature.size(); ++q)
            {
              if (this->include_melt_transport())
                {
                  const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
                  const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
                  depletion = in.composition[q][peridotite_idx] - in.composition[q][porosity_idx];
                }
              melt_fractions[q] = this->melt_fraction(in.temperature[q],
                                                      std::max(0.0, in.pressure[q]),
                                                      depletion);
            }
        }
      return;
    }


    template <int dim>
    void
    MeltGlobal<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      std::vector<double> old_porosity(in.position.size());
      std::vector<double> old_depletion(in.position.size());
      std::vector<double> old_melt_composition(in.position.size());

      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim> >();

      // make sure the compositional fields we want to use exist
      if (this->include_melt_transport())
        AssertThrow(this->introspection().compositional_name_exists("porosity"),
                    ExcMessage("Material model Melt simple with melt transport only "
                               "works if there is a compositional field called porosity."));

      if (this->include_melt_transport() && include_melting_and_freezing)
        AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                    ExcMessage("Material model Melt global only works if there is a "
                               "compositional field called peridotite."));

      if (this->include_melt_transport() && include_melting_and_freezing && read_melt_from_file)
        AssertThrow(this->introspection().compositional_name_exists("crystallized_fraction"),
                    ExcMessage("Reading in melt fractions from a file only works if there "
                               "is a field called crystallized_fraction."));

      // we want to get the porosity field from the old solution here,
      // because we need a field that is not updated in the nonlinear iterations
      if (this->include_melt_transport() && in.current_cell.state() == IteratorState::valid
          && this->get_timestep_number() > 0 && !this->get_parameters().use_operator_splitting)
        {
          // Prepare the field function
          Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector>
          fe_value(this->get_dof_handler(), this->get_old_solution(), this->get_mapping());

          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          fe_value.set_active_cell(in.current_cell);
          fe_value.value_list(in.position,
                              old_porosity,
                              this->introspection().component_indices.compositional_fields[porosity_idx]);

          if(read_melt_from_file && include_melting_and_freezing)
            {
              const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
              const unsigned int crystallization_idx = this->introspection().compositional_index_for_name("crystallized_fraction");
              fe_value.value_list(in.position,
                                  old_depletion,
                                  this->introspection().component_indices.compositional_fields[peridotite_idx]);
              fe_value.value_list(in.position,
                                  old_melt_composition,
                                  this->introspection().component_indices.compositional_fields[crystallization_idx]);
            }
        }
      else if (this->get_parameters().use_operator_splitting)
        for (unsigned int i=0; i<in.position.size(); ++i)
          {
            const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
            old_porosity[i] = in.composition[i][porosity_idx];

            if(read_melt_from_file)
              {
                old_depletion[i] = in.composition[i][this->introspection().compositional_index_for_name("peridotite")];
                old_melt_composition[i] = in.composition[i][this->introspection().compositional_index_for_name("crystallized_fraction")];
              }
          }

      for (unsigned int i=0; i<in.position.size(); ++i)
        {
          // calculate density first, we need it for the reaction term
          // temperature dependence of density is 1 - alpha * (T - T(adiabatic))
          double temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
            temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                      * thermal_expansivity;
          else
            temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;

          // calculate composition dependence of density
          const double delta_rho = this->introspection().compositional_name_exists("peridotite")
                                   ?
                                   depletion_density_change * in.composition[i][this->introspection().compositional_index_for_name("peridotite")]
                                   :
                                   0.0;
          out.densities[i] = (reference_rho_s + delta_rho) * temperature_dependence
                             * std::exp(compressibility * (in.pressure[i] - this->get_surface_pressure()));

          // now compute melting and crystallization
          if (this->include_melt_transport() && include_melting_and_freezing)
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
              const unsigned int crystallization_idx = (read_melt_from_file
                                                        ?
                                                        this->introspection().compositional_index_for_name("crystallized_fraction")
                                                        :
                                                        numbers::invalid_unsigned_int);
              double eq_melt_composition = 0.0;

              // calculate the melting rate as difference between the equilibrium melt fraction
              // and the solution of the previous time step
              // solidus is lowered by previous melting events (fractional melting)
              // we can either use a simplified, linear parametrization, or read the melt fraction
              // from a file
              double melting_rate = 0.0;
              if(read_melt_from_file)
                {
                  const double eq_melt_fraction = melt_fraction_lookup->peridotite_melt_fraction(in.temperature[i],in.pressure[i]);
                  // TODO: interpolate between basalt and peridotite if depletion < 0!

                  eq_melt_composition = melt_fraction_lookup->basalt_melt_fraction(in.temperature[i],in.pressure[i]);
                  if (eq_melt_fraction >= std::max(old_depletion[i],0.0))
                    melting_rate = eq_melt_fraction - std::max(old_depletion[i],0.0);
                  else
                    melting_rate = (old_melt_composition[i] > std::max(eq_melt_composition, 0.0)
                                   ?
                                   (eq_melt_composition - old_melt_composition[i]) / old_melt_composition[i]
                                   :
                                   0.0);
                }
              else
                {
                  const double eq_melt_fraction = melt_fraction(in.temperature[i],
                                                                this->get_adiabatic_conditions().pressure(in.position[i]),
                                                                in.composition[i][peridotite_idx] - in.composition[i][porosity_idx]);
                  melting_rate = eq_melt_fraction - old_porosity[i];
                }

              // do not allow negative porosity or porosity > 1
              if (old_porosity[i] + melting_rate < 0)
                melting_rate = -old_porosity[i];
              if (old_porosity[i] + melting_rate > 1)
                melting_rate = 1.0 - old_porosity[i];

              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                {
                  if (c == peridotite_idx && this->get_timestep_number() > 1 && (in.strain_rate.size()))
                    out.reaction_terms[i][c] = melting_rate
                                               - in.composition[i][peridotite_idx] * trace(in.strain_rate[i]) * this->get_timestep();
                  else if (c == porosity_idx && this->get_timestep_number() > 1 && (in.strain_rate.size()))
                    out.reaction_terms[i][c] = melting_rate
                                               * out.densities[i]  / this->get_timestep();
                  else if (c == crystallization_idx && this->get_timestep_number() > 1 && (in.strain_rate.size()))
                    out.reaction_terms[i][c] = (melting_rate > 0
                                                ?
                                                0.0 /*std::min(std::max(old_porosity[i], 0.0), 1.0) * old_melt_composition[i]
                                                + melting_rate * eq_melt_composition*/
                                                :
                                                std::max(melting_rate/old_melt_composition[i], -old_melt_composition[i]));
                  else
                    out.reaction_terms[i][c] = 0.0;

                  // fill reaction rate outputs if the model uses operator splitting
                  if (this->get_parameters().use_operator_splitting)
                    {
                      if (reaction_rate_out != NULL)
                        {
                          if ((c == peridotite_idx || c == crystallization_idx) && this->get_timestep_number() > 0)
                            reaction_rate_out->reaction_rates[i][c] = out.reaction_terms[i][c] / this->get_timestep() ;
                          else if (c == porosity_idx && this->get_timestep_number() > 0)
                            reaction_rate_out->reaction_rates[i][c] = melting_rate / this->get_timestep();
                          else
                            reaction_rate_out->reaction_rates[i][c] = 0.0;
                        }
                      out.reaction_terms[i][c] = 0.0;
                    }
                }

              const double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));
              out.viscosities[i] = eta_0 * exp(- alpha_phi * porosity);
            }
          else
            {
              out.viscosities[i] = eta_0;
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                out.reaction_terms[i][c] = 0.0;
            }

          out.entropy_derivative_pressure[i]    = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = thermal_conductivity;
          out.compressibilities[i] = 0.0;

          double visc_temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
            {
              const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
              visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
            }
          else
            {
              const double delta_temp = in.temperature[i]-reference_T;
              visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/reference_T),1e4),1e-4);
            }
          out.viscosities[i] *= visc_temperature_dependence;
        }

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim> >();

      if (melt_out != NULL)
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              double porosity = std::max(in.composition[i][porosity_idx],0.0);

              melt_out->fluid_viscosities[i] = eta_f;
              melt_out->permeabilities[i] = reference_permeability * std::pow(porosity,3) * std::pow(1.0-porosity,2);
              melt_out->fluid_density_gradients[i] = Tensor<1,dim>();

              // temperature dependence of density is 1 - alpha * (T - T(adiabatic))
              double temperature_dependence = 1.0;
              if (this->include_adiabatic_heating ())
                temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                          * thermal_expansivity;
              else
                temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
              melt_out->fluid_densities[i] = reference_rho_f * temperature_dependence
                                             * std::exp(melt_compressibility * (in.pressure[i] - this->get_surface_pressure()));

              melt_out->compaction_viscosities[i] = xi_0 * exp(- alpha_phi * porosity);

              double visc_temperature_dependence = 1.0;
              if (this->include_adiabatic_heating ())
                {
                  const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
                  visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
                }
              else
                {
                  const double delta_temp = in.temperature[i]-reference_T;
                  visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/reference_T),1e4),1e-4);
                }
              melt_out->compaction_viscosities[i] *= visc_temperature_dependence;
            }
        }
    }



    template <int dim>
    void
    MeltGlobal<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt global");
        {
          prm.declare_entry ("Reference solid density", "3000",
                             Patterns::Double (0),
                             "Reference density of the solid $\\rho_{s,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference melt density", "2500",
                             Patterns::Double (0),
                             "Reference density of the melt/fluid$\\rho_{f,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: $K$.");
          prm.declare_entry ("Reference shear viscosity", "5e20",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa s$.");
          prm.declare_entry ("Reference bulk viscosity", "1e22",
                             Patterns::Double (0),
                             "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa s$.");
          prm.declare_entry ("Reference melt viscosity", "10",
                             Patterns::Double (0),
                             "The value of the constant melt viscosity $\\eta_f$. Units: $Pa s$.");
          prm.declare_entry ("Exponential melt weakening factor", "27",
                             Patterns::Double (0),
                             "The porosity dependence of the viscosity. Units: dimensionless.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of the shear viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal bulk viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of the bulk viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $C_p$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Reference permeability", "1e-8",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: $m^2$.");
          prm.declare_entry ("Depletion density change", "0.0",
                             Patterns::Double (),
                             "The density contrast between material with a depletion of 1 and a "
                             "depletion of zero. Negative values indicate lower densities of "
                             "depleted material. Depletion is indicated by the compositional "
                             "field with the name peridotite. Not used if this field does not "
                             "exist in the model. "
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Surface solidus", "1300",
                             Patterns::Double (0),
                             "Solidus for a pressure of zero. "
                             "Units: $K$.");
          prm.declare_entry ("Depletion solidus change", "200.0",
                             Patterns::Double (),
                             "The solidus temperature change for a depletion of 100\\%. For positive "
                             "values, the solidus gets increased for a positive peridotite field "
                             "(depletion) and lowered for a negative peridotite field (enrichment). "
                             "Scaling with depletion is linear. Only active when fractional melting "
                             "is used. "
                             "Units: $K$.");
          prm.declare_entry ("Pressure solidus change", "6e-8",
                             Patterns::Double (),
                             "The linear solidus temperature change with pressure. For positive "
                             "values, the solidus gets increased for positive pressures. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Solid compressibility", "0.0",
                             Patterns::Double (0),
                             "The value of the compressibility of the solid matrix. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Melt compressibility", "0.0",
                             Patterns::Double (0),
                             "The value of the compressibility of the melt. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Melt bulk modulus derivative", "0.0",
                             Patterns::Double (0),
                             "The value of the pressure derivative of the melt bulk "
                             "modulus. "
                             "Units: None.");
          prm.declare_entry ("Include melting and freezing", "true",
                             Patterns::Bool (),
                             "Whether to include melting and freezing (according to a simplified "
                             "linear melting approximation in the model (if true), or not (if "
                             "false).");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/melt-fraction-model/melt_global/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Melt fraction file name", "peridotite_melt_contour_output.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file names of the melt fraction data (melt fraction "
                             "data is assumed to be in order with the ordering "
                             "of the compositional fields). Note that there are "
                             "three options on how many files need to be listed "
                             "here: 1. If only one file is provided, it is used "
                             "for the whole model domain, and compositional fields "
                             "are ignored. 2. If there is one more file name than the "
                             "number of compositional fields, then the first file is "
                             "assumed to define a `background composition' that is "
                             "modified by the compositional fields. If there are "
                             "exactly as many files as compositional fields, the fields are "
                             "assumed to represent the fractions of different materials "
                             "and the average property is computed as a sum of "
                             "the value of the compositional field times the "
                             "material property of that field.");
          prm.declare_entry ("Read melt fraction from file", "false",
                             Patterns::Bool (),
                             "Whether to read the melt fraction from a data file (if true) "
                             "or to use a simple linearized, analytical melting model.");
          prm.declare_entry ("Pressure unit in melt fraction file", "Pa",
                             Patterns::Selection("Pa|GPa|bar|kbar"),
                             "What unit the pressure should have in the data file that "
                             "determines the melt fraction."
                             "\n\n"
                             "Possible choices: Pa|GPa|bar|kbar"
                             "\n\n"
                             "This option is ignored if no such data file is used in the "
                             "computation..");
          prm.declare_entry ("Temperature unit in melt fraction file", "Kelvin",
                             Patterns::Selection("Kelvin|Celsius"),
                             "What unit the temperature should have in the data file that "
                             "determines the melt fraction."
                             "\n\n"
                             "Possible choices: Kelvin|Celsius"
                             "\n\n"
                             "This option is ignored if no such data file is used in the "
                             "computation..");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MeltGlobal<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt global");
        {
          reference_rho_s                   = prm.get_double ("Reference solid density");
          reference_rho_f                   = prm.get_double ("Reference melt density");
          reference_T                       = prm.get_double ("Reference temperature");
          eta_0                             = prm.get_double ("Reference shear viscosity");
          xi_0                              = prm.get_double ("Reference bulk viscosity");
          eta_f                             = prm.get_double ("Reference melt viscosity");
          reference_permeability            = prm.get_double ("Reference permeability");
          thermal_viscosity_exponent        = prm.get_double ("Thermal viscosity exponent");
          thermal_bulk_viscosity_exponent   = prm.get_double ("Thermal bulk viscosity exponent");
          thermal_conductivity              = prm.get_double ("Thermal conductivity");
          reference_specific_heat           = prm.get_double ("Reference specific heat");
          thermal_expansivity               = prm.get_double ("Thermal expansion coefficient");
          alpha_phi                         = prm.get_double ("Exponential melt weakening factor");
          depletion_density_change          = prm.get_double ("Depletion density change");
          surface_solidus                   = prm.get_double ("Surface solidus");
          depletion_solidus_change          = prm.get_double ("Depletion solidus change");
          pressure_solidus_change           = prm.get_double ("Pressure solidus change");
          compressibility                   = prm.get_double ("Solid compressibility");
          melt_compressibility              = prm.get_double ("Melt compressibility");
          include_melting_and_freezing      = prm.get_bool ("Include melting and freezing");

          data_directory                    = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          melt_fraction_file_name           = prm.get ("Melt fraction file name");
          read_melt_from_file               = prm.get_bool ("Read melt fraction from file");

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model Melt simple with Thermal viscosity exponent can not have reference_T=0."));

          pressure_unit                     = prm.get ("Pressure unit in melt fraction file");
          temperature_unit                  = prm.get ("Temperature unit in melt fraction file");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    MeltGlobal<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (this->get_parameters().use_operator_splitting
          && out.template get_additional_output<ReactionRateOutputs<dim> >() == NULL)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std_cxx11::shared_ptr<MaterialModel::AdditionalMaterialOutputs<dim> >
            (new MaterialModel::ReactionRateOutputs<dim> (n_points, this->n_compositional_fields())));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MeltGlobal,
                                   "melt global",
                                   "A material model that implements a simple formulation of the "
                                   "material parameters required for the modelling of melt transport, "
                                   "including a source term for the porosity according to a simplified "
                                   "linear melting model similar to \\cite{schmeling2006}:\n"
                                   "$\\phi_\\text{equilibrium} = \\frac{T-T_\\text{sol}}{T_\\text{liq}-T_\\text{sol}}$\n"
                                   "with "
                                   "$T_\\text{sol} = T_\\text{sol,0} + \\Delta T_p \\, p + \\Delta T_c \\, C$ \n"
                                   "$T_\\text{liq} = T_\\text{sol}  + \\Delta T_\\text{sol-liq}$.")
  }
}
