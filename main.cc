//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Basic framework for Phase Field (Cahn-Hilliard mixed formulation)
//Created May 2018
//authors: rudraa (2018)
//

//deal.II headers
#include "include/headers.h"
//input parameter headers
#include "parameters.h"
//physics headers
#include "include/chemo.h"

//Namespace
namespace phaseField1
{
  using namespace dealii;

  //Initial conditions
  template <int dim>
  class InitalConditions: public Function<dim> {
  public:
    InitalConditions (): Function<dim>(TotalDOF){std::srand(5);}
    /*
    double value (const Point<dim>   &p) const {
      return (0.63 + 0.02*(0.5 -(double)(std::rand() % 100 )/100.0)); //c
      }*/
    void vector_value (const Point<dim>   &p, Vector<double>   &values) const {
      Assert (values.size() == TotalDOF, ExcDimensionMismatch (values.size(),TotalDOF));
     
      values(0)=0.00;
      if (std::sqrt(p.square())<0.22) {values(0)=1.00;}
      // values(1)=((double)(std::rand()%100))/100.;
      values(1)=0.01 + ((double)(std::rand()%25)/100.)*0.02;
    }
  };
  
  template <int dim>
  class phaseField{
  public:
    phaseField ();
    ~phaseField ();
    void run ();

  private:
    void applyBoundaryConditions(const unsigned int increment);
    void setup_system ();
    void assemble_system ();
    void solveIteration (bool isProject=false);
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int increment);
    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;
    IndexSet                                  locally_owned_dofs;
    IndexSet                                  locally_relevant_dofs;
    ConstraintMatrix                          constraints;
    LA::MPI::SparseMatrix                     system_matrix;
    LA::MPI::Vector                           locally_relevant_solution, U, Un, UGhost, UnGhost, dU;
    LA::MPI::Vector                           system_rhs;
    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;
    double                                    bulkenergy, volumephase, interfaceenergy;
    
    //solution variables
    unsigned int currentIncrement, currentIteration;
    double totalTime, currentTime, dt;
    Sacado::Fad::DFad<double> free_energy;
    std::vector<std::string> nodal_solution_names; std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
    std::ofstream Bulk, Face, Volume;
  };

  template <int dim>
  phaseField<dim>::phaseField ():
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening)),
    fe(FE_Q<dim>(1),TotalDOF),
    dof_handler (triangulation),
    pcout (std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator)== 0)),
    computing_timer (mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times){
    //solution variables
    dt=TimeStep; totalTime=TotalTime;
    currentIncrement=0; currentTime=0;
    free_energy=0.;
    //nodal Solution names
    nodal_solution_names.push_back("c"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    nodal_solution_names.push_back("solute"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    // nodal_solution_names.push_back("mu"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
      Bulk.open("Bulk.txt");
      Volume.open("Volume.txt");
      Face.open("Face.txt");
    }
  }
  
  template <int dim>
  phaseField<dim>::~phaseField () {
    dof_handler.clear ();
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
      Bulk.close();
      Volume.close();
      Face.close();
    }
  }

  //Apply boundary conditions
  template <int dim>
  void phaseField<dim>::applyBoundaryConditions(const unsigned int increment){
    constraints.clear (); 
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    
    //Setup boundary conditions
    //No Dirchlet BC are necessary for the parabolic problem
    
    constraints.close ();
  }
  
  //Setup
  template <int dim>
  void phaseField<dim>::setup_system (){
    TimerOutput::Scope t(computing_timer, "setup");
    dof_handler.distribute_dofs (fe);
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);
    
    locally_relevant_solution.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    //Non-ghost vectors
    system_rhs.reinit (locally_owned_dofs, mpi_communicator);
    U.reinit (locally_owned_dofs, mpi_communicator);
    Un.reinit (locally_owned_dofs, mpi_communicator);
    dU.reinit (locally_owned_dofs, mpi_communicator);
    //Ghost vectors
    UGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    UnGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    //call applyBoundaryConditions to setup constraints matrix needed for generating the sparsity pattern
    applyBoundaryConditions(0);
    
    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
  }

  //Assembly
  template <int dim>
  void phaseField<dim>::assemble_system (){
    TimerOutput::Scope t(computing_timer, "assembly");
    system_rhs=0.0; system_matrix=0.0;
    const QGauss<dim>  quadrature_formula(3);
    const QGauss<dim-1>	face_quadrature_formula (2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);
    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors);
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    // std::cout<<"dof per cell"<<dofs_per_cell;
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;
  
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    Sacado::Fad::DFad<double> BulkEnergy=0., PhaseVolume=0., FaceEnergy=0.;
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix = 0; local_rhs = 0; 
	cell->get_dof_indices (local_dof_indices);
	 //AD variables
	
	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); Table<1, double > ULocalConv(dofs_per_cell);
	for (unsigned int i=0; i<dofs_per_cell; ++i){
	  if (std::abs(UGhost(local_dof_indices[i]))<1.0e-16){ULocal[i]=0.0;}
	  else{ULocal[i]=UGhost(local_dof_indices[i]);}
	  ULocal[i].diff (i, dofs_per_cell);
	  ULocalConv[i]= UnGhost(local_dof_indices[i]);
	}
	//
	//	std::cout<<dofs_per_cell<<"     " ;
	//for(unsigned int i=0;i<dofs_per_cell;i++)std::cout<<ULocal[i]<<" ";exit(-1);
	dealii::Table<1,double> c_conv(n_q_points);
	//setup residual vector
	free_energy=0.;
	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
	for (unsigned int i=0; i<dofs_per_cell; ++i) {R[i]=0.0;}
	
	//populate residual vector 
	residualForChemo(fe_values, 0, fe_face_values, cell, dt, ULocal, ULocalConv, R,/* currentTime, totalTime,*/ c_conv,local_matrix,currentIncrement,free_energy,BulkEnergy, PhaseVolume, FaceEnergy);
	
	//evaluate Residual(R) and Jacobian(R')

	for(unsigned int i=0;i<dofs_per_cell;i++){
	  local_rhs(i)=-R[i].val();
	}
	for(unsigned int i=0;i<dofs_per_cell;i++){
	  for(unsigned int j=0;j<dofs_per_cell;j++){
	    local_matrix(i,j)=R[i].fastAccessDx(j);
	  }
	}
	/*for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  for (unsigned int j=0; j<dofs_per_cell; ++j){
	    // R' by AD
	    local_matrix(i,j)= R[i].fastAccessDx(j);
	  }
	  //R
	  local_rhs(i) = -R[i].val();
	  }*/
	/*for(int i=0;i<dofs_per_cell;i++){
	  for(int j=0;j<dofs_per_cell;j++){
	    std::cout<<local_matrix(i,j)<<" ";
	  }std::cout<<"\n";
	}
	std::cout<<"\n\n\n";
	for(unsigned int i=0;i<dofs_per_cell;i++)std::cout<<R[i]<<" ";
	exit(-1);*/

	constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
      }
    //bulkenergy, volumephase, interfaceenergy;
    bulkenergy=BulkEnergy.val(); Utilities::MPI::sum(bulkenergy, mpi_communicator);
    volumephase=PhaseVolume.val(); Utilities::MPI::sum(volumephase, mpi_communicator);
    interfaceenergy=FaceEnergy.val(); Utilities::MPI::sum(interfaceenergy, mpi_communicator);
    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
  }
  
template <int dim>
  void phaseField<dim>::solveIteration(bool isProject){
    TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);
      
    //check for convergence of iterative solver, and in case of slow convergence for smaller problem switch to Direct Solver.  
    //try
    {
      //Iterative solvers from Petsc and Trilinos
      SolverControl solver_control (dof_handler.n_dofs(), 1e-12);
#ifdef USE_PETSC_LA
      LA::SolverGMRES solver(solver_control, mpi_communicator);
#else
      LA::SolverGMRES solver(solver_control);
#endif
      LA::MPI::PreconditionJacobi preconditioner;
      LA::MPI::PreconditionJacobi::AdditionalData data;
#ifdef USE_PETSC_LA
      //data.symmetric_operator = true;
#else
      // Trilinos defaults are good 
#endif
      if(!isProject){
	preconditioner.initialize(system_matrix, data);
	solver.solve (system_matrix, completely_distributed_solution, system_rhs, preconditioner);
	constraints.distribute (completely_distributed_solution);
	locally_relevant_solution = completely_distributed_solution;
	dU = completely_distributed_solution;
      }
      /*else{
	preconditioner.initialize(mass_matrix, data);
	solver.solve (mass_matrix, completely_distributed_solution, system_rhs, preconditioner);
	constraints_L2.distribute(completely_distributed_solution);
	locally_relevant_solution_L2=completely_distributed_solution;
	U_L2=completely_distributed_solution;
	UGhost_L2=U_L2;
      }*/
      pcout << "   Solved in " << solver_control.last_step()
	    << " iterations." << std::endl;
    }
    /*
    catch(...){
      //Direct solver MUMPS
      SolverControl cn;
      PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
      if(!isProject){
	solver.set_symmetric_mode(false);
	solver.solve(system_matrix, completely_distributed_solution, system_rhs);
	constraints.distribute (completely_distributed_solution);
	locally_relevant_solution = completely_distributed_solution;
	dU = completely_distributed_solution;
      }
      else{
	solver.set_symmetric_mode(true);
	solver.solve(mass_matrix, completely_distributed_solution, system_rhs);
	constraints_L2.distribute(completely_distributed_solution);
	locally_relevant_solution_L2=completely_distributed_solution;
	U_L2=completely_distributed_solution;
	UGhost_L2=U_L2;
      }
    */
  }

 

  //Solve
  template <int dim>
  void phaseField<dim>::solve(){
    double res=1, tol=1.0e-6, abs_tol=1.0e-9, initial_norm=0, current_norm=0;
    double machineEPS=1.0e-15;
    currentIteration=0;
    char buffer[200];
    while (true){
      if (currentIteration>=50){sprintf(buffer, "maximum number of iterations reached without convergence. \n"); pcout<<buffer; break; exit (1);}
      if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); pcout<<buffer; break; exit (1);}
      assemble_system();
      current_norm=system_rhs.l2_norm();
      initial_norm=std::max(initial_norm, current_norm);
      res=current_norm/initial_norm;
      sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); pcout<<buffer; 
      if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); pcout<<buffer; pcout<<"\n"; /*pcout<<currentIncrement<<"\t "<<free_energy.val();/*myfile<<free_energy<<" ";*/ break;}
      solveIteration();
      U+=dU; UGhost=U; 
      ++currentIteration;
    }
    Un=U; UnGhost=Un;
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
      Bulk <<currentIncrement << "\t" << bulkenergy << "\n" << std::flush;
      Volume <<currentIncrement << "\t" << volumephase << "\n" << std::flush;;
      Face << currentIncrement << "\t" << interfaceenergy << "\n" << std::flush;
    }
  }

  //Output
  template <int dim>
  void phaseField<dim>::output_results (const unsigned int cycle) {
    TimerOutput::Scope t(computing_timer, "output");
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (UnGhost, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);    

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");
    
    data_out.build_patches ();
    const std::string filename = ("solution-" +
                                  Utilities::int_to_string (cycle, 2) +
                                  "." +
                                  Utilities::int_to_string
                                  (triangulation.locally_owned_subdomain(), 4));
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
      std::vector<std::string> filenames;
      for (unsigned int i=0;
	   i<Utilities::MPI::n_mpi_processes(mpi_communicator);
	   ++i)
	filenames.push_back ("solution-" +
			     Utilities::int_to_string (cycle, 2) +
			     "." +
			     Utilities::int_to_string (i, 4) +
			     ".vtu");
      
      std::ofstream master_output (("solution-" +
				    Utilities::int_to_string (cycle, 2) +
				    ".pvtu").c_str());
      data_out.write_pvtu_record (master_output, filenames);
    }
  }

  //Solve problem
  template <int dim>
  void phaseField<dim>::run (){
    //setup problem geometry and mesh
    GridGenerator::hyper_cube (triangulation, -problemWidth/2.0, problemWidth/2.0, true);
    triangulation.refine_global (refinementFactor);
    setup_system ();
    pcout << "   Number of active cells:       "
	  << triangulation.n_global_active_cells()
	  << std::endl
	  << "   Number of degrees of freedom: "
	  << dof_handler.n_dofs()
	  << std::endl;
    
    //setup initial conditions
    VectorTools::interpolate(dof_handler, InitalConditions<dim>(), U); Un=U;
    
    //sync ghost vectors to non-ghost vectors
    UGhost=U;  UnGhost=Un;
    output_results (0);
    //myfile.open("freeEnergy.txt",ios::out);
    //Time stepping
    currentIncrement=0;
    for (currentTime=0; currentTime<totalTime; currentTime+=dt){
      currentIncrement++;
      solve();
      /*if(currentIncrement<20 || (currentIncrement%100==0)){
	output_results(currentIncrement);
      }*/
      output_results(currentIncrement);
      pcout << std::endl;
    }
   
    //computing_timer.print_summary ();
  }
}


int main(int argc, char *argv[]){
  try
    {
      using namespace dealii;
      using namespace phaseField1;
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      phaseField<2> problem;
      //fstream myfile;
      problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
