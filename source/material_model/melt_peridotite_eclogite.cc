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


#include <aspect/material_model/melt_peridotite_eclogite.h>
#include <aspect/simulator.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    MeltPeridotiteEclogite<dim>::initialize ()
    {
      // check if the applicable compositional fields exist
      AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                  ExcMessage("Material model Melt peridotite eclogite only works if there is a "
                             "compositional field called peridotite."));

      if (this->include_melt_transport())
        {
          AssertThrow(this->introspection().compositional_name_exists("porosity"),
                      ExcMessage("Material model Melt peridotite eclogite with melt transport only "
                                 "works if there is a compositional field called porosity."));
        }
    }


    template <int dim>
    double
    MeltPeridotiteEclogite<dim>::
    reference_viscosity () const
    {
      return eta_0;
    }


    template <int dim>
    bool
    MeltPeridotiteEclogite<dim>::
    is_compressible () const
    {
      return model_is_compressible;
    }

    template <int dim>
    double
    MeltPeridotiteEclogite<dim>::
    melt_fraction (const double temperature,
                   const double pressure) const
    {
      // anhydrous melting of peridotite after Katz, 2003
      const double T_solidus  = A1 + 273.15
                                + A2 * pressure
                                + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus = C1 + 273.15
                                + C2 * pressure
                                + C3 * pressure * pressure;

      // melt fraction for peridotite with clinopyroxene
      double peridotite_melt_fraction;
      if (temperature < T_solidus || pressure > 1.3e10)
        peridotite_melt_fraction = 0.0;
      else if (temperature > T_lherz_liquidus)
        peridotite_melt_fraction = 1.0;
      else
        peridotite_melt_fraction = std::pow((temperature - T_solidus) / (T_lherz_liquidus - T_solidus),beta);

      // melt fraction after melting of all clinopyroxene
      const double R_cpx = r1 + r2 * std::max(0.0, pressure);
      const double F_max = M_cpx / R_cpx;

      if (peridotite_melt_fraction > F_max && temperature < T_liquidus)
        {
          const double T_max = std::pow(F_max,1/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
          peridotite_melt_fraction = F_max + (1 - F_max) * pow((temperature - T_max) / (T_liquidus - T_max),beta);
        }
      return peridotite_melt_fraction;
    }


    template <int dim>
    double
    MeltPeridotiteEclogite<dim>::
    entropy_change (const double temperature,
                    const double pressure,
                    const double maximum_melt_fraction,
                    const NonlinearDependence::Dependence dependence) const
    {
      double melt_fraction_derivative = 0;

      // calculate latent heat of melting
      // we need the change of melt fraction in dependence of pressure and temperature

      // for peridotite after Katz, 2003
      const double T_solidus        = A1 + 273.15
                                      + A2 * pressure
                                      + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus       = C1 + 273.15
                                      + C2 * pressure
                                      + C3 * pressure * pressure;

      const double dT_solidus_dp        = A2 + 2 * A3 * pressure;
      const double dT_lherz_liquidus_dp = B2 + 2 * B3 * pressure;
      const double dT_liquidus_dp       = C2 + 2 * C3 * pressure;

      if (temperature > T_solidus && temperature < T_liquidus && pressure < 1.3e10)
        {
          // melt fraction when clinopyroxene is still present
          double melt_fraction_derivative_temperature
            = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
              / (T_lherz_liquidus - T_solidus);

          double melt_fraction_derivative_pressure
            = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
              * (dT_solidus_dp * (temperature - T_lherz_liquidus)
                 + dT_lherz_liquidus_dp * (T_solidus - temperature))
              / pow(T_lherz_liquidus - T_solidus,2);

          // melt fraction after melting of all clinopyroxene
          const double R_cpx = r1 + r2 * std::max(0.0, pressure);
          const double F_max = M_cpx / R_cpx;

          if (melt_fraction(temperature, pressure) > F_max)
            {
              const double T_max = std::pow(F_max,1.0/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
              const double dF_max_dp = - M_cpx * std::pow(r1 + r2 * pressure,-2) * r2;
              const double dT_max_dp = dT_solidus_dp
                                       + 1.0/beta * std::pow(F_max,1.0/beta - 1.0) * dF_max_dp * (T_lherz_liquidus - T_solidus)
                                       + std::pow(F_max,1.0/beta) * (dT_lherz_liquidus_dp - dT_solidus_dp);

              melt_fraction_derivative_temperature
                = (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                  / (T_liquidus - T_max);

              melt_fraction_derivative_pressure
                = dF_max_dp
                  - dF_max_dp * std::pow((temperature - T_max)/(T_liquidus - T_max),beta)
                  + (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                  * (dT_max_dp * (T_max - T_liquidus) - (dT_liquidus_dp - dT_max_dp) * (temperature - T_max)) / std::pow(T_liquidus - T_max, 2);
            }

          if (dependence == NonlinearDependence::temperature)
            melt_fraction_derivative = melt_fraction_derivative_temperature;
          else if (dependence == NonlinearDependence::pressure)
            melt_fraction_derivative = melt_fraction_derivative_pressure;
          else
            AssertThrow(false, ExcMessage("not implemented"));
        }
      return melt_fraction_derivative;
    }


    template <int dim>
    void
    MeltPeridotiteEclogite<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      for (unsigned int q=0; q<in.temperature.size(); ++q)
        melt_fractions[q] = this->melt_fraction(in.temperature[q],
                                                std::max(0.0, in.pressure[q]));
      return;
    }


    template <int dim>
    void
    MeltPeridotiteEclogite<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      DiffusionDislocation<dim>::evaluate(in, out);

      std::vector<double> maximum_melt_fractions(in.position.size());
      std::vector<double> old_porosity(in.position.size());
      std::vector<double> old_temperature(in.position.size());

      // we want to get the peridotite field from the old solution here,
      // because it tells us how much of the material was already molten
      if (this->include_melt_transport() && in.cell
          && this->get_timestep_number() > 0)
        {
          // Prepare the field function
          Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector>
          fe_value(this->get_dof_handler(), this->get_old_solution(), this->get_mapping());

          // get peridotite and porosity field from the old the solution
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
          const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");

          fe_value.set_active_cell(*in.cell);
          fe_value.value_list(in.position,
                              maximum_melt_fractions,
                              this->introspection().component_indices.compositional_fields[peridotite_idx]);

          fe_value.value_list(in.position,
                              old_porosity,
                              this->introspection().component_indices.compositional_fields[porosity_idx]);

          fe_value.value_list(in.position,
                              old_temperature,
                              this->introspection().component_indices.temperature);
        }

      for (unsigned int i=0; i<in.position.size(); ++i)
        {

          // we can not use the densities from the diffusion dislocation model
          // as they assume compositional fields between 0 and 1, so we overwrite them
          // calculate density first, we need it for the reaction term
          // first, calculate temperature dependence of density
          double temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
            {
              // temperature dependence is 1 - alpha * (T - T(adiabatic))
              temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                        * thermal_expansivity;
            }
          else
            temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;

          // calculate composition dependence of density
          const double delta_rho = this->introspection().compositional_name_exists("peridotite")
                                   ?
                                   depletion_density_change * in.composition[i][this->introspection().compositional_index_for_name("peridotite")]
                                   :
                                   0.0;
          const double delta_rho_lithosphere = this->introspection().compositional_name_exists("lithosphere")
                                               ?
                                               lithosphere_density_change * in.composition[i][this->introspection().compositional_index_for_name("lithosphere")]
                                               :
                                               0.0;
          out.densities[i] = (reference_rho_s + delta_rho + delta_rho_lithosphere)
                             * temperature_dependence * std::exp(compressibility * (in.pressure[i] - this->get_surface_pressure()));


          if (this->include_melt_transport())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");

              // calculate the melting rate as difference between the equilibrium melt fraction
              // and the solution of the previous time step
              double melting = 0.0;
              if (fractional_melting)
                {
                  // solidus is lowered by previous melting events (fractional melting)
                  const double solidus_change = (in.composition[i][peridotite_idx] - in.composition[i][porosity_idx]) * depletion_solidus_change;
                  const double eq_melt_fraction = melt_fraction(in.temperature[i] - solidus_change, this->get_adiabatic_conditions().pressure(in.position[i]));
                  melting = eq_melt_fraction - old_porosity[i];
                }
              else if (this->get_timestep_number() > 0)
                {
                  // batch melting
                  // melting rate = df/dT dT/dt + df/dp grad p u_s
                  const double df_dp = entropy_change (in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]), 0.0, NonlinearDependence::pressure);
                  const double df_dT = entropy_change (in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]), 0.0, NonlinearDependence::temperature);
                  // melting = df_dT * (in.temperature[i] - old_temperature[i]) + df_dp * (in.pressure_gradient[i] * in.velocity[i]) * this->get_timestep();
                  melting = melt_fraction(in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]))
                            - std::max(maximum_melt_fractions[i], 0.0);
                }
              // freezing of melt below the solidus
              {
                const double freezing_potential = melt_fraction(in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i])) - old_porosity[i];
                const double freezing = freezing_rate * this->get_timestep() / year_in_seconds * 0.5 * (freezing_potential - std::abs(freezing_potential));
                melting += freezing;
              }

              // do not allow negative porosity
              if (old_porosity[i] + melting < 0)
                melting = -old_porosity[i];

              // because depletion is a volume-based, and not a mass-based property that is advected,
              // additional scaling factors on the right hand side apply
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                {
                  if (c == peridotite_idx && this->get_timestep_number() > 0 && (in.strain_rate.size()))
                    out.reaction_terms[i][c] = melting * (1 - maximum_melt_fractions[i])
                                               / (1 - maximum_melt_fractions[i]);
                  else if (c == porosity_idx && this->get_timestep_number() > 0 && (in.strain_rate.size()))
                    out.reaction_terms[i][c] = melting
                                               * out.densities[i] / this->get_timestep();
                  else
                    out.reaction_terms[i][c] = 0.0;
                }

              // reduce viscosity if there is melt present
              if (in.strain_rate.size())
                {
                  const double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));
                  out.viscosities[i] *= exp(- alpha_phi * porosity);
                }
            }

          out.compressibilities[i] = compressibility;
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
              melt_out->permeabilities[i] = (old_porosity[i] > this->get_melt_handler().melt_transport_threshold
                                             ?
                                             std::max(reference_permeability * std::pow(porosity,3) * std::pow(1.0-porosity,2),0.0)
                                             :
                                             0.0);

              // first, calculate temperature dependence of density
              double temperature_dependence = 1.0;
              if (this->include_adiabatic_heating ())
                {
                  // temperature dependence is 1 - alpha * (T - T(adiabatic))
                  temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                            * thermal_expansivity;
                }
              else
                temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;

              // the fluid compressibility includes two parts, a constant compressibility, and a pressure-dependent one
              // this is a simplified formulation, experimental data are often fit to the Birch-Murnaghan equation of state
              const double fluid_compressibility = melt_compressibility / (1.0 + in.pressure[i] * melt_bulk_modulus_derivative * melt_compressibility);

              melt_out->fluid_densities[i] = reference_rho_f * std::exp(fluid_compressibility * (in.pressure[i] - this->get_surface_pressure()))
                                             * temperature_dependence;

              melt_out->fluid_density_gradients[i] = melt_out->fluid_densities[i] * melt_out->fluid_densities[i]
                                                     * fluid_compressibility
                                                     * this->get_gravity_model().gravity_vector(in.position[i]);

              const double phi_0 = 0.05;
              porosity = std::max(std::min(porosity,0.995),1.e-3);
              melt_out->compaction_viscosities[i] = xi_0 * phi_0 / porosity;
              if (in.strain_rate.size())
                {
                  melt_out->compaction_viscosities[i] *= out.viscosities[i]/(eta_0 * exp(- alpha_phi * porosity));
                }
            }
        }
    }


    template <int dim>
    void
    MeltPeridotiteEclogite<dim>::declare_parameters (ParameterHandler &prm)
    {
      DiffusionDislocation<dim>::declare_parameters(prm);

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt peridotite eclogite");
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
                             "The value of the pressure derivative of the melt bulk"
                             "modulus. "
                             "Units: None.");
          prm.declare_entry ("Use full compressibility", "false",
                             Patterns::Bool (),
                             "If the compressibility should be used everywhere in the code"
                             "(if true), changing the volume of material when the density changes, "
                             "or only in the momentum conservation and advection equations "
                             "(if false).");
          prm.declare_entry ("Use fractional melting", "false",
                             Patterns::Bool (),
                             "If fractional melting should be used (if true), including a solidus "
                             "change based on depletion (in this case, the amount of melt that has "
                             "migrated away from its origin), and freezing of melt when it has moved "
                             "to a region with temperatures lower than the solidus; or if batch "
                             "melting should be used (if false), assuming that the melt fraction only "
                             "depends on temperature and pressure, and how much melt has already been "
                             "generated at a given point, but not considering movement of melt in "
                             "the melting parameterization.");
          prm.declare_entry ("Freezing rate", "0.0",
                             Patterns::Double (0),
                             "Freezing rate of melt when in subsolidus regions."
                             "Units: $1/yr$.");
          prm.declare_entry ("Depletion density change", "0.0",
                             Patterns::Double (),
                             "The density contrast between material with a depletion of 1 and a "
                             "depletion of zero. Negative values indicate lower densities of"
                             "depleted material. Depletion is indicated by the compositional"
                             "field with the name peridotite. Not used if this field does not "
                             "exist in the model."
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Lithosphere density change", "0.0",
                             Patterns::Double (),
                             "The density contrast between the lithosphere and the background mantle. "
                             "Not used if this field does not exist in the model."
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Depletion solidus change", "200.0",
                             Patterns::Double (0),
                             "The solidus temperature change for a depletion of 100\\%. For positive "
                             "values, the solidus gets increased for a positive peridotite field "
                             "(depletion) and lowered for a negative peridotite field (enrichment)."
                             "Scaling with depletion is linear. Only active when fractional melting "
                             "is used. "
                             "Units: $K$.");
          prm.declare_entry ("A1", "1085.7",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of peridotite. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("A2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("A3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("B1", "1475.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the lherzolite "
                             "liquidus used for calculating the fraction "
                             "of peridotite-derived melt. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("B2", "8.0e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("B3", "-3.2e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("C1", "1780.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the liquidus "
                             "of peridotite. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("C2", "4.50e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("C3", "-2.0e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("r1", "0.5",
                             Patterns::Double (),
                             "Constant in the linear function that "
                             "approximates the clinopyroxene reaction "
                             "coefficient. "
                             "Units: non-dimensional.");
          prm.declare_entry ("r2", "8e-11",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the linear function that approximates "
                             "the clinopyroxene reaction coefficient. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("beta", "1.5",
                             Patterns::Double (),
                             "Exponent of the melting temperature in "
                             "the melt fraction calculation. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Mass fraction cpx", "0.15",
                             Patterns::Double (),
                             "Mass fraction of clinopyroxene in the "
                             "peridotite to be molten. "
                             "Units: non-dimensional.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MeltPeridotiteEclogite<dim>::parse_parameters (ParameterHandler &prm)
    {
      DiffusionDislocation<dim>::parse_parameters(prm);

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt peridotite eclogite");
        {
          reference_rho_s            = prm.get_double ("Reference solid density");
          reference_rho_f            = prm.get_double ("Reference melt density");
          reference_T                = prm.get_double ("Reference temperature");
          eta_0                      = prm.get_double ("Reference shear viscosity");
          xi_0                       = prm.get_double ("Reference bulk viscosity");
          eta_f                      = prm.get_double ("Reference melt viscosity");
          reference_permeability     = prm.get_double ("Reference permeability");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          thermal_bulk_viscosity_exponent = prm.get_double ("Thermal bulk viscosity exponent");
          thermal_conductivity       = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_expansivity        = prm.get_double ("Thermal expansion coefficient");
          alpha_phi                  = prm.get_double ("Exponential melt weakening factor");
          compressibility            = prm.get_double ("Solid compressibility");
          melt_compressibility       = prm.get_double ("Melt compressibility");
          model_is_compressible      = prm.get_bool ("Use full compressibility");
          fractional_melting         = prm.get_bool ("Use fractional melting");
          freezing_rate              = prm.get_double ("Freezing rate");
          melt_bulk_modulus_derivative = prm.get_double ("Melt bulk modulus derivative");
          depletion_density_change   = prm.get_double ("Depletion density change");
          depletion_solidus_change   = prm.get_double ("Depletion solidus change");
          lithosphere_density_change = prm.get_double ("Lithosphere density change");

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model Melt simple with Thermal viscosity exponent can not have reference_T=0."));

          A1              = prm.get_double ("A1");
          A2              = prm.get_double ("A2");
          A3              = prm.get_double ("A3");
          B1              = prm.get_double ("B1");
          B2              = prm.get_double ("B2");
          B3              = prm.get_double ("B3");
          C1              = prm.get_double ("C1");
          C2              = prm.get_double ("C2");
          C3              = prm.get_double ("C3");
          r1              = prm.get_double ("r1");
          r2              = prm.get_double ("r2");
          beta            = prm.get_double ("beta");
          M_cpx           = prm.get_double ("Mass fraction cpx");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MeltPeridotiteEclogite,
                                   "melt peridotite eclogite",
                                   "A material model that implements a simple formulation of the "
                                   "material parameters required for the modelling of melt transport, "
                                   "including a source term for the porosity according to the melting "
                                   "model for dry peridotite of \\cite{KSL2003}. This also includes a "
                                   "computation of the latent heat of melting (if the 'latent heat' "
                                   "heating model is active)."
                                   "\n\n"
                                   "Most of the material properties are constant, except for the shear, "
                                   "viscosity $\\eta$, the compaction viscosity $\\xi$, and the "
                                   "permeability $k$, which depend on the porosity; and the solid and melt "
                                   "densities, which depend on temperature and pressure:\n "
                                   "$\\eta(\\phi,T) = \\eta_0 e^{\\alpha(\\phi-\\phi_0)} e^{-\\beta(T-T_0)/T_0}$, "
                                   "$\\xi(\\phi,T) = \\xi_0 \\frac{\\phi_0}{\\phi} e^{-\\beta(T-T_0)/T_0}$, "
                                   "$k=k_0 \\phi^n (1-\\phi)^m$, "
                                   "$\\rho=\\rho_0 (1 - \\alpha (T - T_\\text{adi})) e^{\\kappa p}$."
                                   "\n\n"
                                   "The model is compressible only if this is specified in the input file, "
                                   "and contains compressibility for both solid and melt.")
  }
}
