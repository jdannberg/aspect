#include <aspect/material_model/melt_interface.h>
#include <aspect/compositional_initial_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>



namespace aspect
{
  /**
   * This is the "Solitary wave" benchmark defined in the following paper:
   * @code
   *  @Article{DMGT11,
   *    author =       {T. Keller and D. A. May and B. J. P. Kaus},
   *    title =        {Numerical modelling of magma dynamics coupled
   *                    to tectonic deformation of lithosphere and crust},
   *    journal =      {Geophysical Journal International},
   *    year =         2013,
   *    volume =       195(3),
   *    pages =        {1406-1442}
   * @endcode
   *
   * To calculate the initial condition, which is a solitary wave solution of
   * the magma dynamics equations, we use the equation for the one-dinemsional
   * case and the non-dimensionalization as it is described in
   * @code
   *  @Article{SS11,
   *    author =       {G. Simpson and M. Spiegelman},
   *    title =        {Solitary Wave Benchmarks in Magma Dynamics},
   *    journal =      {Journal of Scientific Computing},
   *    year =         2011,
   *    volume =       49(3),
   *    pages =        {268-290}
   * @endcode
   *
   * Specifically, this means that we scale the porosity with the background
   * porosity, and the coordinates with the compaction length $\delta_0$, which is
   * defined as $\sqrt \frac{k(\phi_0) \xi^{*}+4/3 \eta^{*}}{\eta_f}$.  $k(\phi_0)$ is the
   * permeability at background porosity, $\xi^{*}$ is the compaction viscosity,
   * $\eta^{*}$ is the shear viscosity of the fluid and $\eta_f$ is the shear viscosity
   * of the melt.
   */
  namespace SolitaryWaveBenchmark
  {
    using namespace dealii;

    namespace AnalyticSolutions
    {
      // vectors to store the porosity field and the corresponding coordinate in
      const unsigned int max_points = 2e8;
      std::vector<double> porosity, coordinate;

      /**
       * @note The solitary wave solution only exists as a function x = func(phi)
       * and not phi = func(x), which is what we would like to have for describing
       * the shape of the wave. Thus, we calculate x = func(phi) for a range of phis
       * between the background porosity and the amplitude of the wave. In a next
       * step, we interpolate these values to the grid.
       *
       * @param phi The characteristic shape of the wave, with phi --> 1
       * for x --> +- infinity
       *
       * @param amplitude The amplitude of the solitary wave, which is always
       * greater than 1.
       */
       double solitary_wave_solution (const double phi, const double amplitude)
       {
         AssertThrow(phi > 1.0 && phi <= amplitude,
                     ExcMessage("The solitary wave solution can only be computed "
                         "for porosities larger than the background porosity of 1"
                         "and smaller than or equal to the amplitude of the wave."));
         AssertThrow(amplitude > 1,
                     ExcMessage("Amplitude of the solitary wave must be larger than 1!"));
         const double A_1   = std::sqrt(amplitude - 1.0);
         const double A_phi = std::sqrt(amplitude - phi);
         return std::sqrt(amplitude + 0.5)
                * (2 * A_phi - 1.0/A_1 * std::log((A_1 - A_phi)/(A_1 + A_phi)));
       }


       /**
        * This function reads the coordinate and the porosity of the solitary wave
        * from an input file.
        *
        * @param filename Name of the input file.
        */
       void read_solitary_wave_solution (const std::string &filename)
       {
         std::string temp;
         std::ifstream in(filename.c_str(), std::ios::in);
         AssertThrow (in,
                      ExcMessage (std::string("Couldn't open file <") + filename));

         while (!in.eof())
           {
             double x, f;
             in >> x >> f;
             if (in.eof())
               break;
             getline(in, temp);

             coordinate.insert(coordinate.begin(),x);
             porosity.insert(porosity.begin(),f);
           }
       }


       /**
        * This function gets the coordinate as an input parameters and gives
        * back the porosity of the solitary wave. As this function is only defined
        * implicitly, we have to interpolate from the coordinates where we have the
        * porosity to our mesh.
        *
        * @param amplitude The amplitude of the solitary wave, which is always
        * greater than 1.
        * @param offset The offset of the center of the solitary wave from the
        * boundary of the domain.
        */
       void compute_porosity (const double amplitude,
                              const double background_porosity,
                              const double offset,
                              const double compaction_length,
                              const bool read_solution,
                              const std::string file_name)
       {
         // non-dimensionalize the amplitude
         const double non_dim_amplitude = amplitude / background_porosity;

         if(read_solution)
           read_solitary_wave_solution(file_name);
         else
           {
             porosity.resize(max_points);
             coordinate.resize(max_points);

             // get the coordinates where we have the solution
             for (unsigned int i=0;i<max_points;++i)
               {
                 porosity[i] = 1.0 + 1e-10*non_dim_amplitude
                               + double(i)/double(max_points-1) * (non_dim_amplitude * (1.0 - 1e-10) - 1.0);
                 coordinate[i] = solitary_wave_solution(porosity[i], non_dim_amplitude);
               }
           }

             for (unsigned int i=0;i<coordinate.size();++i)
               {
                 // re-scale porosity and position
                 porosity[i] *= background_porosity;
                 coordinate[i] *= compaction_length;
               }
       }


       double interpolate (const double position,
                           const double offset)
       {
         // interpolate from the solution grid to the mesh used in the simulation
         // solitary wave is a monotonically decreasing function, so the coordinates
         // should be in descending order

         // we only have the solution of the solitary wave for
         // coordinates larger than 0 (one half of the wave)
         const double x = (position > offset
                          ?
                          position - offset
                          :
                          offset - position);

         if(x > coordinate[0])
           return porosity[0];

         unsigned int j= coordinate.size()-2;
         unsigned int i = j/2;
         while (!(x < coordinate[j] && x >= coordinate[j+1]))
         {
           if(x < coordinate[j])
        	 j += i;
           else
        	 j -= i;
           if (i>1)
             i /= 2;
         }

         const double distance = (x - coordinate[j+1])
                                 /(coordinate[j] - coordinate[j+1]);
         return porosity[j+1] + distance * (porosity[j] - porosity[j+1]);
       }

       /**
        * The exact solution for the Solitary wave benchmark.
        */
       template <int dim>
       class FunctionSolitaryWave : public Function<dim>
       {
         public:
    	   FunctionSolitaryWave (const double offset, const double delta, const std::vector<double> &initial_pressure, const double max_z)
    	     :
    		 Function<dim>(dim+4),
    		 offset_(offset),
    	     delta_(delta),
    	     initial_pressure_(initial_pressure),
    	     max_z_(max_z)
    	   {}

    	   virtual void set_delta(const double delta)
    	   {
    		   delta_ = delta;
    	   }

           virtual void vector_value (const Point< dim > &p,
                                      Vector< double >   &values) const
           {
        	 const double index = static_cast<int>((p[dim-1]-delta_)/max_z_ * (initial_pressure_.size()-1));
        	 AssertThrow(index < initial_pressure_.size(), ExcMessage("not in range"));
        	 const double z_coordinate1 = static_cast<double>(index)/static_cast<double>(initial_pressure_.size()-1) * max_z_;
        	 const double z_coordinate2 = static_cast<double>(index+1)/static_cast<double>(initial_pressure_.size()-1) * max_z_;
        	 const double interpolated_pressure = (index == initial_pressure_.size()-1)
                    		                      ?
                    		                      initial_pressure_[index]
                    		                      :
                    		                      initial_pressure_[index] + (initial_pressure_[index+1] - initial_pressure_[index])
                    		                                  * (p[dim-1]-delta_ - z_coordinate1) / (z_coordinate2 - z_coordinate1);

             values[dim+3] = AnalyticSolutions::interpolate(p[dim-1]-delta_,offset_); //porosity
             values[dim+1] = interpolated_pressure;                                   //compaction pressure
           }

         private:
           const double offset_;
           double delta_;
           const std::vector<double> initial_pressure_;
           const double max_z_;
       };
     }


    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class SolitaryWaveMaterial : public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
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
        if ((dependence & MaterialModel::NonlinearDependence::compositional_fields) != MaterialModel::NonlinearDependence::none)
          return true;
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
        return eta_0;
      }

      virtual double reference_density () const
      {
        return reference_rho_s;
      }

      double length_scaling (const double porosity) const
      {
        return std::sqrt(reference_permeability * std::pow(porosity,3) * (xi_0 + 4.0/3.0 * eta_0) / eta_f);
      }

      double velocity_scaling (const double porosity) const
      {
    	const Point<dim> surface_point = this->get_geometry_model().representative_point(0.0);
        return reference_permeability * std::pow(porosity,2) * (reference_rho_s - reference_rho_f)
        	   * this->get_gravity_model().gravity_vector(surface_point).norm() / eta_f;
      }

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

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

        for (unsigned int i=0;i<in.position.size();++i)
          {
            double porosity = in.composition[i][porosity_idx];

            out.viscosities[i] = eta_0 * (1.0 - porosity);
            out.densities[i] = reference_rho_s;
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 0.0;
            out.compressibilities[i] = 0.0;
          }
      }

      virtual void evaluate_with_melt(const typename MaterialModel::MeltInterface<dim>::MaterialModelInputs &in,
                                      typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs &out) const
      {
        evaluate(in, out);
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

        for (unsigned int i=0;i<in.position.size();++i)
          {
            double porosity = in.composition[i][porosity_idx];

            out.compaction_viscosities[i] = xi_0 * (1.0 - porosity);
            out.fluid_viscosities[i]= eta_f;
            out.permeabilities[i]= reference_permeability * std::pow(porosity,3);
            out.fluid_densities[i]= reference_rho_f;
            out.fluid_compressibilities[i] = 0.0;
          }

      }

      private:
        double reference_rho_s;
        double reference_rho_f;
        double eta_0;
        double xi_0;
        double eta_f;
        double reference_permeability;
    };

    template <int dim>
    void
    SolitaryWaveMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Solitary wave");
        {
          prm.declare_entry ("Reference solid density", "3000",
              Patterns::Double (0),
              "Reference density of the solid $\\rho_{s,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference melt density", "2500",
              Patterns::Double (0),
              "Reference density of the melt/fluid$\\rho_{f,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference shear viscosity", "1e20",
              Patterns::Double (0),
              "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
              "Units: $Pa s$.");
          prm.declare_entry ("Reference compaction viscosity", "1e20",
              Patterns::Double (0),
              "The value of the constant volumetric viscosity $\\xi_0$ of the solid matrix. "
              "Units: $Pa s$.");
          prm.declare_entry ("Reference melt viscosity", "100.0",
              Patterns::Double (0),
              "The value of the constant melt viscosity $\\eta_f$. Units: $Pa s$.");

          prm.declare_entry ("Reference permeability", "5e-9",
              Patterns::Double(),
              "Reference permeability of the solid host rock."
              "Units: $m^2$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SolitaryWaveMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Solitary wave");
        {
          reference_rho_s            = prm.get_double ("Reference solid density");
          reference_rho_f            = prm.get_double ("Reference melt density");
          eta_0                      = prm.get_double ("Reference shear viscosity");
          xi_0                       = prm.get_double ("Reference compaction viscosity");
          eta_f                      = prm.get_double ("Reference melt viscosity");
          reference_permeability     = prm.get_double ("Reference permeability");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    /**
     * An initial conditions model for the solitary waves benchmark.
     */
    template <int dim>
    class SolitaryWaveInitialCondition : public CompositionalInitialConditions::Interface<dim>,
                                         public ::aspect::SimulatorAccess<dim>
    {
      public:

      /**
       * Initialization function. Take references to the material model and
       * get the compaction length, so that it can be used subsequently to
       * compute the analytical solution for the shape of the solitary wave.
       */
      void
      initialize ();

      /**
       * Return the boundary velocity as a function of position.
       */
      virtual
      double
      initial_composition (const Point<dim> &position, const unsigned int n_comp) const;

      static
      void
      declare_parameters (ParameterHandler &prm);

      virtual
      void
      parse_parameters (ParameterHandler &prm);

      double
      get_amplitude () const;

      double
      get_background_porosity () const;

      double
      get_offset () const;

      private:
        double amplitude;
        double background_porosity;
        double offset;
        double compaction_length;
        bool read_solution;
        std::string file_name;
    };

    template <int dim>
    double
    SolitaryWaveInitialCondition<dim>::get_amplitude () const
    {
      return amplitude;
    }

    template <int dim>
    double
    SolitaryWaveInitialCondition<dim>::get_background_porosity () const
    {
      return background_porosity;
    }

    template <int dim>
    double
    SolitaryWaveInitialCondition<dim>::get_offset () const
    {
      return offset;
    }

    template <int dim>
    void
    SolitaryWaveInitialCondition<dim>::initialize ()
    {
      std::cout << "Initialize solitary wave solution"
                << std::endl;

      if (dynamic_cast<const SolitaryWaveMaterial<dim> *>(&this->get_material_model()) != NULL)
        {
          const SolitaryWaveMaterial<dim> *
          material_model
            = dynamic_cast<const SolitaryWaveMaterial<dim> *>(&this->get_material_model());

          compaction_length = material_model->length_scaling(background_porosity);
        }
      else
        {
          AssertThrow(false,
                      ExcMessage("Initial condition Solitary Wave only works with the material model Solitary wave."));
        }

      AnalyticSolutions::compute_porosity(amplitude,
                                          background_porosity,
                                          offset,
                                          compaction_length,
                                          read_solution,
                                          file_name);
    }


    template <int dim>
    double
    SolitaryWaveInitialCondition<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      return AnalyticSolutions::interpolate(position[dim-1],
                                            offset);
    }


    template <int dim>
    void
    SolitaryWaveInitialCondition<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Compositional initial conditions");
      {
        prm.enter_subsection("Solitary wave initial condition");
        {
          prm.declare_entry ("Amplitude", "0.01",
              Patterns::Double (0),
              "Amplitude of the solitary wave. Units: none.");
          prm.declare_entry ("Background porosity", "0.001",
              Patterns::Double (0),
              "Background porosity of the solitary wave. Units: none.");
          prm.declare_entry ("Offset", "150",
              Patterns::Double (0),
              "Offset of the center of the solitary wave from the boundary"
              "of the domain. "
              "Units: $m$.");
          prm.declare_entry ("Read solution from file", "false",
              Patterns::Bool (),
              "Whether to read the porosity initial condition from "
              "a file or to compute it.");
          prm.declare_entry ("File name", "solitary_wave.txt",
              Patterns::Anything (),
              "The file name of the porosity initial condition data. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SolitaryWaveInitialCondition<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Compositional initial conditions");
      {
        prm.enter_subsection("Solitary wave initial condition");
        {
          amplitude            = prm.get_double ("Amplitude");
          background_porosity  = prm.get_double ("Background porosity");
          offset               = prm.get_double ("Offset");
          read_solution        = prm.get_bool ("Read solution from file");
          file_name            = prm.get ("File name");

          AssertThrow(amplitude > background_porosity,
                      ExcMessage("Amplitude of the solitary wave must be larger "
                          "than the background porosity."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    /**
      * A postprocessor that evaluates the accuracy of the solution.
      *
      * The implementation of error evaluators that correspond to the
      * benchmarks defined in the paper Keller et al. reference above.
      */
    template <int dim>
    class SolitaryWavePostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        /**
         * Initialization function. Take references to the material model and
         * initial conditions model to get the parameters necessary for computing
         * the analytical solution for the shape of the solitary wave and store them.
         */
        void
        initialize ();

        void
        store_initial_pressure ();

        double
        compute_phase_shift ();

      private:
        double amplitude;
        double background_porosity;
        double offset;
        double compaction_length;
        double velocity_scaling;
        double boundary_velocity;
        unsigned int max_points;
        std::vector<double> initial_pressure;
        double maximum_pressure;
        std_cxx1x::shared_ptr<Function<dim> > ref_func;

    };

    template <int dim>
    void
    SolitaryWavePostprocessor<dim>::initialize ()
    {
        // verify that we are using the "Solitary wave" initial conditions and material model,
    	// then get the parameters we need

        if (dynamic_cast<const SolitaryWaveInitialCondition<dim> *>(&this->get_compositional_initial_conditions()) != NULL)
          {
            const SolitaryWaveInitialCondition<dim> *
            initial_conditions
              = dynamic_cast<const SolitaryWaveInitialCondition<dim> *>(&this->get_compositional_initial_conditions());

            amplitude           = initial_conditions->get_amplitude();
            background_porosity = initial_conditions->get_background_porosity();
            offset              = initial_conditions->get_offset();
          }
        else
          {
            AssertThrow(false,
                        ExcMessage("Postprocessor Solitary Wave only works with the initial conditions model Solitary wave."));
          }

        if (dynamic_cast<const SolitaryWaveMaterial<dim> *>(&this->get_material_model()) != NULL)
          {
            const SolitaryWaveMaterial<dim> *
            material_model
              = dynamic_cast<const SolitaryWaveMaterial<dim> *>(&this->get_material_model());

            compaction_length = material_model->length_scaling(background_porosity);
            velocity_scaling = material_model->velocity_scaling(background_porosity);
          }
        else
          {
            AssertThrow(false,
                        ExcMessage("Postprocessor Solitary Wave only works with the material model Solitary wave."));
          }

        // we also need the boundary velocity, but we can not get it from simulator access
        // TODO: write solitary wave boundary condition where the phase speed is calculated!

        max_points = 1e6;
        initial_pressure.resize(max_points);
        maximum_pressure = 0.0;
    }

    template <int dim>
    void
    SolitaryWavePostprocessor<dim>::store_initial_pressure ()
    {
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.pressure).degree);
      const unsigned int n_q_points = quadrature_formula.size();
      const double max_depth = this->get_geometry_model().maximal_depth();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      // do the same stuff we do in depth average
      std::vector<double> volume(max_points,0.0);
      std::vector<double> p_s(n_q_points);
      std::vector<double> p_f(n_q_points);
      std::vector<double> phi(n_q_points);
      double local_max_pressure = 0.0;

      const unsigned int porosity_index = this->introspection().compositional_index_for_name("porosity");

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
            		                                                                  p_s);
            fe_values[this->introspection().extractors.compaction_pressure].get_function_values (this->get_solution(),
            		                                                                             p_f);
            fe_values[this->introspection().extractors.compositional_fields[porosity_index]].get_function_values (this->get_solution(),
            		                                                                             phi);
            for (unsigned int q=0; q<n_q_points; ++q)
              {
            	double z = fe_values.quadrature_point(q)[dim-1];
                const unsigned int idx = static_cast<unsigned int>((z*(max_points-1))/max_depth);
                AssertThrow(idx < max_points, ExcInternalError());

                initial_pressure[idx] += (1.0-phi[q]) * (p_s[q] - p_f[q]) * fe_values.JxW(q);
                volume[idx] += fe_values.JxW(q);

                local_max_pressure = std::max (local_max_pressure, std::abs((1.0-phi[q]) * (p_s[q] - p_f[q])));
              }
          }

      std::vector<double> pressure_all(max_points, 0.0);
      std::vector<double> volume_all(max_points, 0.0);
      Utilities::MPI::sum(volume, this->get_mpi_communicator(), volume_all);
      Utilities::MPI::sum(initial_pressure, this->get_mpi_communicator(), pressure_all);
      maximum_pressure = Utilities::MPI::max (local_max_pressure, this->get_mpi_communicator());

      // fill the first and last element of the initial_pressure vector if they are empty
      // this makes sure they can be used for the interpolation later on
      if(pressure_all[0] == 0.0)
      {
    	unsigned int j = 1;
    	while (pressure_all[j] == 0.0)
    	  j++;
    	pressure_all[0] = pressure_all[j];
    	volume_all[0]   = volume_all[j];
    	initial_pressure[0] = pressure_all[j] / (static_cast<double>(volume_all[j])+1e-20);
      }

      if(pressure_all[max_points-1] == 0.0)
      {
    	unsigned int k = max_points-2;
    	while (pressure_all[k] == 0.0)
    	  k--;
    	pressure_all[max_points-1] = pressure_all[k];
    	volume_all[max_points-1]   = volume_all[k];
    	initial_pressure[max_points-1] = pressure_all[k] / (static_cast<double>(volume_all[k])+1e-20);
      }

      // interpolate between the non-zero elements to fill the elements that are 0
      for (unsigned int i=1; i<max_points-1; ++i)
      {
    	if (pressure_all[i] != 0.0)
    	  initial_pressure[i] = pressure_all[i] / (static_cast<double>(volume_all[i])+1e-20);
    	else
    	{
    	  // interpolate between the values we have
    	  unsigned int k = i-1;
    	  while (pressure_all[k] == 0.0)
    		k--;
    	  Assert(k >= 0, ExcInternalError());
    	  unsigned int j = i+1;
    	  while (pressure_all[j] == 0.0)
    		j++;
    	  Assert(j < max_points, ExcInternalError());
    	  const double pressure_k = pressure_all[k] / (static_cast<double>(volume_all[k])+1e-20);
    	  const double pressure_j = pressure_all[j] / (static_cast<double>(volume_all[j])+1e-20);
    	  initial_pressure[i] = pressure_k + (pressure_j - pressure_k) * static_cast<double>(i-k)/static_cast<double>(j-k);
    	}
      }
    }

    template <int dim>
    double
    SolitaryWavePostprocessor<dim>::compute_phase_shift ()
    {
      AssertThrow(Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) == 1,
                  ExcNotImplemented());
      // TODO: how do we calculate the coordinate of the wave peak with more than 1 processor?

      AssertThrow(this->introspection().compositional_name_exists("porosity"),
                  ExcMessage("Postprocessor Solitary Wave only works if there is a compositional field called porosity."));
      const unsigned int porosity_index = this->introspection().compositional_index_for_name("porosity");

      // create a quadrature formula based on the compositional element alone.
      // be defensive about determining that a compositional field actually exists
      AssertThrow (this->introspection().base_elements.compositional_fields
                   != numbers::invalid_unsigned_int,
                   ExcMessage("This postprocessor cannot be used without compositional fields."));
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.compositional_fields).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      std::vector<double> compositional_values(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      // The idea here is to first find the maximum, and then use the analytical solution of the
      // solitary wave to calculate a phase shift for every point.
      // This has to be done separately for points left and right of the maximum.
      // In the end, these values for the phase shift are averaged.

      // compute the maximum composition by quadrature (because we also need the coordinate)
      double local_max_composition = -std::numeric_limits<double>::max();
      double local_z_coordinate = -std::numeric_limits<double>::max();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.compositional_fields[porosity_index]].get_function_values (this->get_solution(),
            		                                                                                 compositional_values);
            for (unsigned int q=0; q<n_q_points; ++q)
              {
                const double composition = compositional_values[q];

                if(composition > local_max_composition)
                {
                  local_max_composition = composition;
                  local_z_coordinate = fe_values.quadrature_point(q)[dim-1];
                }
              }
          }

      cell = this->get_dof_handler().begin_active();
      double phase_shift_integral = 0.0;
      unsigned int number_of_points = 0;

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.compositional_fields[porosity_index]].get_function_values (this->get_solution(),
            		                                                                                 compositional_values);

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                const double composition = compositional_values[q];

                if(composition > background_porosity + (amplitude - background_porosity)*0.05  && composition <= amplitude)
                {
                  double z = fe_values.quadrature_point(q)[dim-1];

                  if(z > local_z_coordinate)
                    z -= offset;
                  else
                    z = offset - z;

				  double z_analytical = compaction_length
									    * AnalyticSolutions::solitary_wave_solution(composition/background_porosity,
																				    amplitude/background_porosity);

				  phase_shift_integral += std::abs(z - z_analytical);
				  number_of_points += 1;
                }
              }
          }

      phase_shift_integral /= static_cast<double>(number_of_points);

      // TODO: different case for moving wave (with zero boundary velocity)
      // const double phase_speed = velocity_scaling * (2.0 * amplitude / background_porosity + 1);
      return phase_shift_integral; // + phase_speed * this->get_time();
    }

    template <int dim>
    std::pair<std::string,std::string>
    SolitaryWavePostprocessor<dim>::execute (TableHandler &statistics)
    {
      AssertThrow(Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) == 1,
                  ExcNotImplemented());

      if(this->get_timestep_number()==0)
      {
    	store_initial_pressure();
        ref_func.reset (new AnalyticSolutions::FunctionSolitaryWave<dim>(offset,0.0,initial_pressure,
      		                                                           this->get_geometry_model().maximal_depth()));
      }

      double delta=0;

      if (dynamic_cast<const SolitaryWaveMaterial<dim> *>(&this->get_material_model()) != NULL)
        {
          const SolitaryWaveMaterial<dim> *
          material_model
            = dynamic_cast<const SolitaryWaveMaterial<dim> *>(&this->get_material_model());

          delta = compute_phase_shift();

          // TODO: we should use delta as an argument, but it does not work yet
          // ref_func.set_delta(delta);

        }
      else
        {
          AssertThrow(false,
                      ExcMessage("Postprocessor Solitary Wave only works with the material model Solitary wave."));
        }

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.pressure).degree);

      // we need the compaction pressure, but we only have the solid and the fluid pressure stored in the solution vector.
      // Hence, we create a new vector only with the compaction pressure
      LinearAlgebra::BlockVector compaction_pressure(this->get_solution());

      const unsigned int por_idx = this->introspection().compositional_index_for_name("porosity");
      const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.pressure).get_unit_support_points());
      std::vector<double> porosity_values(quadrature.size());
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values);

      std::vector<types::global_dof_index> local_dof_indices (this->get_fe().dofs_per_cell);
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            cell->get_dof_indices (local_dof_indices);
            fe_values[this->introspection().extractors.compositional_fields[por_idx]].get_function_values (
            		this->get_solution(), porosity_values);
            for (unsigned int j=0; j<this->get_fe().base_element(this->introspection().base_elements.pressure).dofs_per_cell; ++j)
              {
                unsigned int pressure_idx
                = this->get_fe().component_to_system_index(this->introspection().component_indices.pressure,
                    /*dof index within component=*/ j);

                // skip entries that are not locally owned:
                if (!this->get_dof_handler().locally_owned_dofs().is_element(pressure_idx))
                  continue;

                unsigned int p_f_idx
                = this->get_fe().component_to_system_index(this->introspection().component_indices.compaction_pressure,
                    /*dof index within component=*/ j);

                double p_s = this->get_solution()(local_dof_indices[pressure_idx]);
                double p_f = this->get_solution()(local_dof_indices[p_f_idx]);
                double phi = porosity_values[j];
                double p_c;
                p_c = (1.0-phi) * (p_s - p_f);

                compaction_pressure(local_dof_indices[p_f_idx]) = p_c;
              }
          }


      // what we want to compare:
      // (1) error of the numerical phase speed c:
      // c_numerical = c_analytical - Delta / time;
      const double c_analytical = velocity_scaling * (2.0 * amplitude / background_porosity + 1);
      const double c_numerical = c_analytical + delta / this->get_time();
      const double error_c = std::abs (c_numerical / c_analytical - 1);

      // (3) preservation of shape of melt fraction
      // (4) preservation of the shape of compaction pressure

      Vector<float> cellwise_errors_f (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());

      // get correct components for porosity and compaction pressure
      ComponentSelectFunction<dim> comp_f(dim+3, dim+4);
      ComponentSelectFunction<dim> comp_p(dim+1, dim+4);

      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_f,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_f);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
    		                             compaction_pressure,
                                         *ref_func,
                                         cellwise_errors_p,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_p);

      std::ostringstream os;
      os << std::scientific << cellwise_errors_f.l2_norm() / (amplitude * sqrt(cellwise_errors_f.size()))
                    << ", " << cellwise_errors_p.l2_norm() / (maximum_pressure * sqrt(cellwise_errors_p.size()))
                    << ", " << error_c
                    << ", " << delta;

      return std::make_pair("Errors e_f, e_p, e_c, delta:", os.str());
    }

  }
}



// explicit instantiations
namespace aspect
{
  namespace SolitaryWaveBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SolitaryWaveMaterial,
                                   "Solitary Wave",
                                   "A material model that corresponds to the 'SolitaryWave' benchmark "
                                   "defined in Keller et al., JGI, 2013.")

    ASPECT_REGISTER_POSTPROCESSOR(SolitaryWavePostprocessor,
                                  "solitary wave statistics",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "the Keller et al., JGI, 2013, paper with the one computed by ASPECT "
                                  "and reports the error.")

   ASPECT_REGISTER_COMPOSITIONAL_INITIAL_CONDITIONS(SolitaryWaveInitialCondition,
                                                    "Solitary wave initial condition",
                                                    "Composition is set to a solitary wave function.")
  }
}
