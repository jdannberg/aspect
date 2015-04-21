/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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


#ifndef __aspect__model_damage_rheology_h
#define __aspect__model_damage_rheology_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/std_cxx1x/array.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Lookup
    {
      class MaterialLookup;
    }

    /**
     * A material model that consists of globally constant values for all
     * material parameters except that the density decays linearly with the
     * temperature and the viscosity, which depends on the temperature,
     * pressure, strain rate and grain size.
     *
     * The grain size evolves in time, dependent on strain rate, temperature,
     * creep regime, and phase transitions.
     *
     * The model is considered compressible.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class DamageRheology : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. Loads the material data and sets up
         * pointers.
         */
        virtual
        void
        initialize ();

        /**
         * Called at the beginning of each time step and allows the material
         * model to update internal data structures.
         */
        virtual void update();

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return true if the viscosity() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the density() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        density_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the compressibility() function returns something
         * that may depend on the variable identifies by the argument.
         *
         * This function must return false for all possible arguments if the
         * is_compressible() function returns false.
         */
        virtual bool
        compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the specific_heat() function returns something that
         * may depend on the variable identifies by the argument.
         */
        virtual bool
        specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the thermal_conductivity() function returns
         * something that may depend on the variable identifies by the
         * argument.
         */
        virtual bool
        thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the contuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        virtual unsigned int thermodynamic_phase (const double temperature,
                                                  const double pressure,
                                                  const std::vector<double> &compositional_fields) const;

        virtual void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                              typename Interface<dim>::MaterialModelOutputs &out) const;
        /**
         * @}
         */

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
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
        /**
         * @}
         */

        /**
         * Returns the ratio of dislocation to diffusion viscosity. Useful
         * for postprocessing purposes to determine the regime of deformation
         * in the viscosity ratio postprocessor.
         */
        double
        viscosity_ratio (const double temperature,
                         const double pressure,
                         const std::vector<double> &composition,
                         const SymmetricTensor<2,dim> &strain_rate,
                         const Point<dim> &position) const;

        /**
         * Returns the enthalpy as calculated by HeFESTo.
         */
        virtual double enthalpy (const double      temperature,
                                 const double      pressure,
                                 const std::vector<double> &compositional_fields,
                                 const Point<dim> &position) const;

        /**
         * Returns the enthalpy derivatives for the evaluate function and
         * postprocessors.
         */
        std_cxx1x::array<std::pair<double, unsigned int>,2>
        enthalpy_derivative (const typename Interface<dim>::MaterialModelInputs &in) const;

        /**
         * Returns the p-wave velocity as calculated by HeFESTo.
         */
        virtual double seismic_Vp (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &position) const;

        /**
         * Returns the s-wave velocity as calculated by HeFESTo.
         */
        virtual double seismic_Vs (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &position) const;

      protected:
        double reference_rho;
        double reference_T;
        double eta;
        double composition_viscosity_prefactor_1;
        double composition_viscosity_prefactor_2;
        double compositional_delta_rho;
        double thermal_alpha;
        double reference_specific_heat;

        /**
         * The constant compressibility.
         */
        double reference_compressibility;

        /**
         * The thermal conductivity.
         */
        double k_value;

        // grain evolution parameters
        double gas_constant; // J/(K*mol)
        std::vector<double> grain_growth_activation_energy;
        std::vector<double> grain_growth_activation_volume;
        std::vector<double> grain_growth_rate_constant;
        std::vector<double> grain_growth_exponent;
        std::vector<double> reciprocal_required_strain;
        std::vector<double> recrystallized_grain_size;

        // for paleowattmeter
        bool use_paleowattmeter;
        std::vector<double> grain_boundary_energy;
        std::vector<double> boundary_area_change_work_fraction;
        std::vector<double> geometric_constant;

        // rheology parameters
        double dislocation_viscosity_iteration_threshold;
        unsigned int dislocation_viscosity_iteration_number;
        std::vector<double> dislocation_creep_exponent;
        std::vector<double> dislocation_activation_energy;
        std::vector<double> dislocation_activation_volume;
        std::vector<double> dislocation_creep_prefactor;
        std::vector<double> diffusion_creep_exponent;
        std::vector<double> diffusion_activation_energy;
        std::vector<double> diffusion_activation_volume;
        std::vector<double> diffusion_creep_prefactor;
        std::vector<double> diffusion_creep_grain_size_exponent;

        // Because of the nonlinear nature of this material model many
        // parameters need to be kept within bounds to ensure stability of the
        // solution. These bounds can be adjusted as input parameters.
        double max_temperature_dependence_of_eta;
        double min_eta;
        double max_eta;
        double min_specific_heat;
        double max_specific_heat;
        double min_thermal_expansivity;
        double max_thermal_expansivity;
        unsigned int max_latent_heat_substeps;
        double min_grain_size;
        double pv_grain_size_scaling;

        bool advect_log_gransize;


        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;

        virtual double diffusion_viscosity (const double      temperature,
                                            const double      pressure,
                                            const std::vector<double>    &compositional_fields,
                                            const SymmetricTensor<2,dim> &,
                                            const Point<dim> &position) const;

        /**
         * This function calculates the dislocation viscosity. For this purpose
         * we need the dislocation component of the strain rate, which we can
         * only compute by knowing the dislocation viscosity. Therefore, we
         * iteratively solve for the dislocation viscosity and update the
         * dislocation strain rate in each iteration using the new value
         * obtained for the dislocation viscosity. The iteration is started
         * with a dislocation viscosity calculated for the whole strain rate
         * unless a guess for the viscosity is provided, which can reduce the
         * number of iterations significantly.
         */
        virtual double dislocation_viscosity (const double      temperature,
                                              const double      pressure,
                                              const std::vector<double>    &compositional_fields,
                                              const SymmetricTensor<2,dim> &strain_rate,
                                              const Point<dim> &position,
                                              const double viscosity_guess = 0) const;

        /**
         * This function calculates the dislocation viscosity for a given
         * dislocation strain rate.
         */
        double dislocation_viscosity_fixed_strain_rate (const double      temperature,
                                                        const double      pressure,
                                                        const std::vector<double> &,
                                                        const SymmetricTensor<2,dim> &dislocation_strain_rate,
                                                        const Point<dim> &position) const;

        virtual double density (const double temperature,
                                const double pressure,
                                const std::vector<double> &compositional_fields,
                                const Point<dim> &position) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position) const;

        virtual double thermal_expansion_coefficient (const double      temperature,
                                                      const double      pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const;
        /**
         * Rate of grain size growth (Ostwald ripening) or reduction
         * (due to phase transformations) in dependence on temperature
         * pressure, strain rate, mineral phase and creep regime.
         * We use the grain size evolution laws described in Solomatov
         * and Reese, 2008. Grain size variations in the Earth’s mantle
         * and the evolution of primordial chemical heterogeneities,
         * J. Geophys. Res., 113, B07408.
         */
        virtual
        double
        grain_size_growth_rate (const double                  temperature,
                                const double                  pressure,
                                const std::vector<double>    &compositional_fields,
                                const SymmetricTensor<2,dim> &strain_rate,
                                const Tensor<1,dim>          &velocity,
                                const Point<dim>             &position,
                                const unsigned int            phase_index,
                                const int                     crossed_transition) const;

        /**
         * Function that defines the phase transition interface
         * (0 above, 1 below the phase transition).This is done
         * individually for each transition and summed up in the end.
         */
        virtual
        double
        phase_function (const Point<dim> &position,
                        const double temperature,
                        const double pressure,
                        const unsigned int phase) const;

        /**
         * Function that returns the phase for a given
         * position, temperature, pressure and compositional
         * field index.
         */
        virtual
        unsigned int
        get_phase_index (const Point<dim> &position,
                         const double temperature,
                         const double pressure) const;

        /**
         * Function that takes an object in the same format
         * as in.composition as argument and converts the
         * vector that corresponds to the grain size to its
         * logarithms and back and limits the grain size to
         * a global minimum.
         * @in normal_to_log: if true, convert from the grain
         * size to its logarithm, otherwise from log to grain
         * size
         */
        virtual
        void
        convert_log_grain_size (const bool normal_to_log,
                                std::vector<double> &compositional_fields) const;

        // list of depth, width and Clapeyron slopes for the different phase
        // transitions and in which phase they occur
        std::vector<double> transition_depths;
        std::vector<double> transition_temperatures;
        std::vector<double> transition_slopes;
        std::vector<std::string> transition_phases;
        std::vector<double> transition_widths;


        /* The following variables are properties of the material files
         * we read in.
         */
        std::string datadirectory;
        std::vector<std::string> material_file_names;
        std::vector<std::string> derivatives_file_names;
        unsigned int n_material_data;
        bool use_table_properties;
        bool use_enthalpy;
        bool use_bilinear_interpolation;


        /**
         * The format of the provided material files. Currently we support
         * the PERPLEX and HeFESTo data formats.
         */
        enum formats
        {
          perplex,
          hefesto
        } material_file_format;

        /**
         * List of pointers to objects that read and process data we get from
         * Perplex files. There is one pointer/object per compositional field
         * data provided.
         */
        std::vector<std_cxx1x::shared_ptr<MaterialModel::Lookup::MaterialLookup> > material_lookup;
    };

  }
}

#endif
