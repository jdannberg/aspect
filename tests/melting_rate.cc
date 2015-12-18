#include <aspect/material_model/melt_interface.h>
#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;


namespace aspect
{
  template <int dim>
  class MeltingRate:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
      virtual bool
      viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        if ((dependence & MaterialModel::NonlinearDependence::temperature) != MaterialModel::NonlinearDependence::none)
          return true;
        else
          return false;
      }

      virtual bool
      density_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        if ((dependence & MaterialModel::NonlinearDependence::temperature) != MaterialModel::NonlinearDependence::none)
          return true;
        else if ((dependence & MaterialModel::NonlinearDependence::compositional_fields) != MaterialModel::NonlinearDependence::none)
          return true;
        else
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
        return 5e20;
      }

      virtual double reference_density () const
      {
        return 3000;
      }
      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        for (unsigned int i=0; i<in.position.size(); ++i)
          {
            const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
            const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
            const double porosity = std::max(in.composition[i][porosity_idx],0.0);

            for (unsigned int c=0; c<in.composition[i].size(); ++c)
              {
                const double x = in.position[i][0] - 100000.0;
                const double y = in.position[i][1] - 50000.0;
                const double c0 = 10000.0;
                const double deltaT = this->get_timestep();
                if (c == peridotite_idx && deltaT>0.0)
                  out.reaction_terms[i][c] = 0.0001 * std::exp(-(x*x+y*y)/(2*c0*c0));
                else if (c == porosity_idx && deltaT>0.0)
                  out.reaction_terms[i][c] = 0.0001 * std::exp(-(x*x+y*y)/(2*c0*c0)) * 3000.0 / this->get_timestep();
                else
                  out.reaction_terms[i][c] = 0.0;
              }

            out.viscosities[i] = 5e20;
            out.densities[i] = 3000.0 * (1 - 2e-5 * (in.temperature[i] - 293.0));
            out.thermal_expansion_coefficients[i] = 2e-5;
            out.specific_heat[i] = 1250.0;
            out.thermal_conductivities[i] = 4.7;
            out.compressibilities[i] = 0.0;
          }
      }

      virtual void evaluate_with_melt(const typename MaterialModel::MeltInterface<dim>::MaterialModelInputs &in,
                                      typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs &out) const
      {
        evaluate(in, out);
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        for (unsigned int i=0; i<in.position.size(); ++i)
          {
            const double porosity = in.composition[i][porosity_idx];

            out.compaction_viscosities[i] = 5e20;
            out.fluid_viscosities[i]= 10.0;
            out.permeabilities[i]= 1e-8 * std::pow(porosity,3) * std::pow(1.0-porosity,2);
            out.fluid_densities[i]= 2500.0;
            out.fluid_compressibilities[i] = 0.0;
          }

      }

  };


}



// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(MeltingRate,
                                 "melting rate",
                                 "")
}
