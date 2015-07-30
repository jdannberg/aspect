#include <aspect/material_model/melt_interface.h>

/** 
 * This material model extends the simple material model
 * to include melt transport.
 */

namespace aspect
{
  template <int dim>
  class SimpleWithMelt:
      public MaterialModel::MeltInterface<dim>
  {
      public:
      virtual bool
      viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        if ((dependence & MaterialModel::NonlinearDependence::compositional_fields) != MaterialModel::NonlinearDependence::none)
          return true;
        return false;
      }

      virtual bool
      density_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      compressibility_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      specific_heat_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      thermal_conductivity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }

      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_viscosity () const
      {
        return 0.2;
      }

      virtual double reference_density () const
      {
        return 1.0;
      }
      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                 typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        for (unsigned int i=0;i<in.position.size();++i)
          {
            out.viscosities[i] = 1.0;
            out.thermal_expansion_coefficients[i] = 0.01;
            out.specific_heat[i] = 1250.0;
            out.thermal_conductivities[i] = 1e-6;
            out.compressibilities[i] = 0.0;
            out.densities[i] = 1.0 * (1 - out.thermal_expansion_coefficients[i] * in.temperature[i]) + 100*in.composition[i][0];

            // Pressure derivative of entropy at the given positions.
            out.entropy_derivative_pressure[i] = 0.0;
            // Temperature derivative of entropy at the given positions.
            out.entropy_derivative_temperature[i] = 0.0;
            for (unsigned int c=0;c<in.composition[i].size();++c)
              out.reaction_terms[i][c] = 0.0;
          }
      }

      virtual void evaluate_with_melt(const typename MaterialModel::MeltInterface<dim>::MaterialModelInputs &in,
          typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs &out) const
      {
        evaluate(in, out);

        for (unsigned int i=0;i<in.position.size();++i)
          {
            double porosity = in.composition[i][0];
            out.compaction_viscosities[i] = exp(porosity);
            out.fluid_viscosities[i] = 1.0;
            out.permeabilities[i] = porosity * porosity;
            out.fluid_compressibilities[i] = 0.0;
            out.fluid_densities[i] = 0.5;
          }

      }

  };
}
  

// explicit instantiations
namespace aspect
{
    ASPECT_REGISTER_MATERIAL_MODEL(SimpleWithMelt,
                                   "simple with melt",
                                   "A simple material model that is like the "
				   "'simple' model, but with melt migration.")
}
