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


#ifndef __aspect__model_melt_simple_h
#define __aspect__model_melt_simple_h

#include <aspect/material_model/melt_interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that implements a simple formulation of the
     * material parameters required for the modelling of melt transport,
     * including a source term for the porosity according to the melting
     * model for dry peridotite of Katz, 2003.
     *
     * Most of the material properties are constant, except for the shear,
     * compaction and melt viscosities and the permeability, which depend on
     * the porosity.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MeltSimple : public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
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

        virtual void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                              typename Interface<dim>::MaterialModelOutputs &out) const;

        virtual void evaluate_with_melt(const typename MeltInterface<dim>::MaterialModelInputs &in,
                                        typename MeltInterface<dim>::MaterialModelOutputs &out) const;

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

      private:
        double reference_rho_s;
        double reference_rho_f;
        double reference_T;
        double eta_0;
        double eta_f;
        double thermal_viscosity_exponent;
        double thermal_expansivity;
        double reference_specific_heat;
        double thermal_conductivity;
        double reference_permeability;
        double alpha_phi;

        /**
         * Parameters for anhydrous melting of peridotite after Katz, 2003
         */

        // for the solidus temperature
        double A1;   // °C
        double A2; // °C/Pa
        double A3; // °C/(Pa^2)

        // for the lherzolite liquidus temperature
        double B1;   // °C
        double B2;   // °C/Pa
        double B3; // °C/(Pa^2)

        // for the liquidus temperature
        double C1;   // °C
        double C2;  // °C/Pa
        double C3; // °C/(Pa^2)

        // for the reaction coefficient of pyroxene
        double r1;     // cpx/melt
        double r2;     // cpx/melt/GPa
        double M_cpx;  // mass fraction of pyroxene

        // melt fraction exponent
        double beta;

        // entropy change upon melting
        double peridotite_melting_entropy_change;

        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure) const;

        virtual
        double
        entropy_change (const double temperature,
                        const double pressure,
                        const double maximum_melt_fraction,
                        const NonlinearDependence::Dependence dependence) const;
    };

  }
}

#endif
