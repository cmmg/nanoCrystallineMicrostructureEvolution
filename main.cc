
//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Basic framework for Finite Strain Elasticity
//Created May 2018
//authors: rudraa (2018)
//
//main
//deal.II headers
#include "include/headers.h"
//input parameter headers
#include "parameters.h"
//physics headers
#include "include/mechanics.h"
#include "include/chemo.h"

//Namespace
namespace elasticity1
{
  using namespace dealii;

  template<int dim>
  class InitialCondition: public Function<dim>{
  public:
    std::vector<Point<dim> >* grainPoints;
    std::vector<unsigned int> * grainID;
    InitialCondition(std::vector<Point<dim> >*_grainPoints, std::vector<unsigned int>*_grainId):Function<dim>(totalDOF), grainPoints(_grainPoints),grainID(_grainId){}
    void vector_value(const Point<dim> &p, Vector<double>&values)const{
      Assert(values.size()==totalDOF, ExcDimensionMismatch(values.size(),totalDOF));
      // initial 3 components of values denote displacement, rest of the components denote phase field parameter
      values(0)=0.;values(1)=0.; values(2)=0.;
      /*for(int i=0;i<n_diff_grains;i++){
	values(i)=((double)(std::rand()%100))/100;
      }*/
      Table<1, double> distance(n_seed_points);
      for(unsigned int i=0;i<n_seed_points;i++){
	distance[i]=p.distance((*grainPoints)[i]);
      }
      
      int min=0;
      for(unsigned int i=0;i<n_seed_points;i++){
	if(distance[i]<distance[min])min=i;
      }
      unsigned int g_id=(*grainID)[min];
      for(unsigned int i=0;i<n_diff_grains;i++){
	if(i==g_id) values(dim+i)=(double)(std::rand()%100)/100;
	else values(dim+i)=0.;
	}
      
      
    }
  };
  
  
  template <int dim>
  class elasticity{
  public:
    elasticity ();
    ~elasticity ();
    void run ();

  private:
    void applyBoundaryConditions(const unsigned int increment);
    void setup_system ();
    void assemble_system ();
    void solveIteration (bool isProject=false);
    void grain_generation();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int increment, bool isProject=false);
    void l2_projection();
    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;
    IndexSet                                  locally_owned_dofs;
    IndexSet                                  locally_relevant_dofs;
    ConstraintMatrix                          constraints, constraints2,constraints_L2;
    LA::MPI::SparseMatrix                     system_matrix,mass_matrix;
    LA::MPI::Vector                           locally_relevant_solution, U, Un, UGhost, UnGhost, dU;
    LA::MPI::Vector                           locally_relevant_solution_L2, U_L2, UGhost_L2;
    LA::MPI::Vector                           system_rhs;
    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;
    std::vector<Point<dim> >                  grain_seeds;
    std::vector<unsigned int>                 grain_ID;
    //solution variables
    unsigned int currentIncrement, currentIteration;
    double totalTime, currentTime, dt;
    std::vector<std::string> nodal_solution_names; std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
 std::vector<std::string> nodal_solution_names_L2; std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation_L2;
    //history variables
    std::map<typename DoFHandler<dim>::active_cell_iterator, std::vector< historyVariables<dim>* > > history;
  };

  template <int dim>
  elasticity<dim>::elasticity ():
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening)),
    fe(FE_Q<dim>(1), totalDOF),// change the degree of freedoms per node from displacement to displacement+grain_parameter
    dof_handler (triangulation),
    pcout (std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator)== 0)),
    computing_timer (mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times){
    //solution variables
    dt=TimeStep; totalTime=TotalTime;
    currentIncrement=0; currentTime=0;
    
    //nodal Solution names
    for (unsigned int i=0; i<dim; ++i){
      nodal_solution_names.push_back("u"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    }
    // nodal_solution_names.push_back("phase"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);


    char buffer[100];
    for(char i=0;i<n_diff_grains;i++){
      sprintf(buffer, "eta%u",i);
      nodal_solution_names.push_back(buffer);nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    }
    
    
    nodal_solution_names_L2.push_back("Grain_Id"); nodal_data_component_interpretation_L2.push_back(DataComponentInterpretation::component_is_scalar);
    nodal_solution_names_L2.push_back("VonMissesStrain"); nodal_data_component_interpretation_L2.push_back(DataComponentInterpretation::component_is_scalar);
    if (dim==3) {nodal_solution_names_L2.push_back("Gamma3"); nodal_data_component_interpretation_L2.push_back(DataComponentInterpretation::component_is_scalar);}
    nodal_solution_names_L2.push_back("random"); nodal_data_component_interpretation_L2.push_back(DataComponentInterpretation::component_is_scalar);
  }
  
  template <int dim>
  elasticity<dim>::~elasticity (){
    dof_handler.clear ();
  }
  
  //Apply boundary conditions
  template <int dim>
  void elasticity<dim>::applyBoundaryConditions(const unsigned int increment){
    constraints.clear (); constraints2.clear (); constraints_L2.clear ();
    constraints.reinit (locally_relevant_dofs);
    constraints2.reinit (locally_relevant_dofs);
    constraints_L2.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints2);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints_L2);
    
    //Setup boundary conditions
    /*
    std::vector<bool> uBCX0 (dim, true); 
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(dim), constraints, uBCX0);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(dim), constraints2, uBCX0);
    std::vector<bool> uBCX1 (dim, false); uBCX1[0]=true; 
    VectorTools::interpolate_boundary_values (dof_handler, 1, ConstantFunction<dim>(0.001, dim), constraints, uBCX1);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(dim), constraints2, uBCX1);
    */
    
    
    //for shear uncomment this
    
    /* std::vector<bool> uBCX0 (dim, true); uBCX0[0]=true;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(dim), constraints, uBCX0);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(dim), constraints2, uBCX0);
    std::vector<bool> uBCX1 (dim, true); uBCX1[2]=false;
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(dim), constraints, uBCX1);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(dim), constraints2, uBCX1);
    std::vector<bool> uBCX1Z (dim, false); uBCX1Z[2]=true;
    VectorTools::interpolate_boundary_values (dof_handler, 1, ConstantFunction<dim>(0.0001, dim), constraints, uBCX1Z);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(dim), constraints2, uBCX1Z);*/



     
    
    // for simple shear uncomment this
    /*
    std::vector<bool> uBCX1Z (dim, false); uBCX1Z[2]=true;
    VectorTools::interpolate_boundary_values (dof_handler, 1, ConstantFunction<dim>(0.0001, dim), constraints, uBCX1Z);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(dim), constraints2, uBCX1Z);*/
     




    // for pure shear uncomment this    
    /* std::vector<bool> uBCX1Y (dim, true); uBCX1Y[2]=false;
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(dim), constraints, uBCX1Y);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(dim), constraints2, uBCX1Y);*/

    
   //simple tension boundary condition
      
   std::vector<bool> uBCX0 (totalDOF, false); uBCX0[0]=true;
   VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(totalDOF), constraints, uBCX0);
   VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(totalDOF), constraints2, uBCX0);
   
   std::vector<bool> uBCY0 (totalDOF, false); uBCY0[1]=true;
   VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(totalDOF), constraints, uBCY0);
   VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(totalDOF), constraints2, uBCY0);
   if (dim==3) {
     std::vector<bool> uBCZ0 (totalDOF, false); uBCZ0[2]=true;
     VectorTools::interpolate_boundary_values (dof_handler, 4, ZeroFunction<dim>(totalDOF), constraints, uBCZ0);
     VectorTools::interpolate_boundary_values (dof_handler, 4, ZeroFunction<dim>(totalDOF), constraints2, uBCZ0);
   }
   std::vector<bool> uBCX1 (totalDOF, false); uBCX1[0]=true;
   if (currentIncrement<10){
     VectorTools::interpolate_boundary_values (dof_handler, 1, ConstantFunction<dim>(0.00, totalDOF), constraints, uBCX1);
   }
   else{
     VectorTools::interpolate_boundary_values (dof_handler, 1, ConstantFunction<dim>(0.00001, totalDOF), constraints, uBCX1);
   }
   VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(totalDOF), constraints2, uBCX1);
    
   constraints.close ();
   constraints2.close ();    
   constraints_L2.close ();

   }
  //grain id generation
  template<int dim>
  void elasticity<dim>::grain_generation(){

    for(unsigned int I=0;I<n_seed_points;I++){
      Point<dim> grain_id;
      std::srand(5.2*I);
      grain_id[0]=((double)((std::rand())%100)/100)-0.5;
      grain_id[1]=((double)((std::rand())%100)/100)-0.5;
      grain_id[2]=((double)((std::rand())%100)/100)-0.5;
      
      grain_seeds.push_back(grain_id);
    }
    for(unsigned int i=0;i<n_seed_points;i++){
      if(i<n_diff_grains)grain_ID.push_back(i);
      else{
	Table<1, double> distance(i),temp_d(i);unsigned int var,findid_j, findid_k;
	for(unsigned int j=0;j<i;j++){
	  distance[j]=grain_seeds[i].distance(grain_seeds[j]);
	  temp_d[j]=distance[j];
	}
	for(unsigned int j=0;j<i;j++){
	  for(unsigned int k=0;k<i;k++){
	    if(temp_d[k]>temp_d[j]){double t=temp_d[k];temp_d[k]=temp_d[j];temp_d[j]=t;}
	  }
	}
	for(unsigned int j=0;j<i;j++){ 
	  var=0;
	  for(unsigned int l=0;l<i;l++){if(temp_d[i-1-j]==distance[l])findid_j=grain_ID[l]; }
	  for(unsigned int k=0;k<n_diff_grains-1;k++){
	    for(unsigned int l=0;l<i;l++)if(temp_d[k]==distance[l])findid_k=grain_ID[l]; 
	    if(findid_j==findid_k){var=1;break;}
	  }
	  if(var==1)continue;
	  else{grain_ID.push_back(findid_j);break;}
	}
	//else ends
      }
    }
    
  }
  

  //Setup
  template <int dim>
  void elasticity<dim>::setup_system (){
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
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints2, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

    //create data structure for L2 projection
     DynamicSparsityPattern dsp_L2 (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp_L2, constraints_L2, false);
    SparsityTools::distribute_sparsity_pattern (dsp_L2, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    mass_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp_L2, mpi_communicator);
    locally_relevant_solution_L2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    U_L2.reinit (locally_owned_dofs, mpi_communicator);
    UGhost_L2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    
    //setup history variables
    const QGauss<dim>  quadrature_formula(3);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell){
      if (cell->is_locally_owned()){
	for (unsigned int q=0; q<fe_values.n_quadrature_points; q++){
	  history[cell].push_back(new historyVariables<dim>); //create histroy variables object at each quad point of the cell.
	  //history[cell].back()->Grain_Id=0;
	  //history[cell].back()->orientationAngle=0;
	  for(unsigned int i=0;i<n_slip_systems;i++){
	    history[cell].back()->CRSS[i]=0.0016;
	    history[cell].back()->CRSS_iteration[i]=0.0016;
	    history[cell].back()->Gamma_iteration[i]=0.;
	    history[cell].back()->Gamma[i]=0.;
	    history[cell].back()->shear_Stress_iteration[i]=0.;
	    history[cell].back()->shear_Stress[i]=0.;
	    history[cell].back()->Gamma_total[i]=0.;
	  }
	  for(unsigned int i=0;i<dim;i++)
	    {   
	      for(unsigned int j=0;j<dim;j++)
		{	 
		  history[cell].back()->Fp_previous(i,j)=(i==j);
		  history[cell].back()->Fp_previousIteration(i,j)=(i==j);
		  
		}
	    }
	}
      }
    }
  }
  
  //Assembly
  template <int dim>
  void elasticity<dim>::assemble_system (){
    TimerOutput::Scope t(computing_timer, "assembly");
    system_rhs=0.0; system_matrix=0.0;
    const QGauss<dim>  quadrature_formula(3);
    const QGauss<dim-1>	face_quadrature_formula (3);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);
    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors);
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell); 
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;
  
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix = 0; local_rhs = 0; 
	cell->get_dof_indices (local_dof_indices);
	 //AD variables
	Table<1, double> ULocal(dofs_per_cell); Table<1, double > ULocalConv(dofs_per_cell);
	for (unsigned int i=0; i<dofs_per_cell; ++i){
	  ULocal[i]=UGhost(local_dof_indices[i]);
	  ULocalConv[i]= UnGhost(local_dof_indices[i]);
	}
	//get defomration map
	deformationMap<double, dim> defMap(n_q_points); 
	getDeformationMap<double, dim>(fe_values, 0, ULocal, defMap, currentIteration);
	
	//pasing a reference of map to residual function in order to store all stress, strain, back_stress, slip_rate at current cell
	residualForMechanics(fe_values, 0, ULocal, ULocalConv, defMap, currentIteration, currentIncrement, history[cell], local_rhs, local_matrix,  grain_seeds,grain_ID);
	
	residualForChemo(fe_values, dim, fe_face_values, cell,  dt, ULocal, ULocalConv,local_rhs, local_matrix);//DOF component 0 if only chemo, dim for coupling---for mechanics switch residualChemo off
	
	for(int i=0;i<dofs_per_cell;i++){
	  local_rhs[i]=-local_rhs[i];
	}
	
	if ((currentIteration==0)){
	  constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
	}
	else{
	  constraints2.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
	}	
      }
   
    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
  }

   template <int dim>
   void elasticity<dim>::l2_projection (){
     TimerOutput::Scope t(computing_timer, "projection");
     system_rhs=0.0; mass_matrix=0.0;
     const QGauss<dim>  quadrature_formula(3);
     FEValues<dim> fe_values (fe, quadrature_formula,
			      update_values   |
			      update_quadrature_points |
			      update_JxW_values);
     const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell); 
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;
  
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix = 0.0; local_rhs = 0.0; 
	cell->get_dof_indices (local_dof_indices);
	
	//implement L2 projection
	for (unsigned int q=0; q<n_q_points; q++){
	  for(unsigned int i=0;i<dofs_per_cell;i++){
	    const unsigned int ci = fe_values.get_fe().system_to_component_index(i).first;
	    if (ci==0){
	      local_rhs(i)+=fe_values.shape_value(i,q)*(history[cell][q]->Grain_Id)*fe_values.JxW(q);
	    }
	    else if (ci==1){
	      local_rhs(i)+= fe_values.shape_value(i,q)*history[cell][q]->vonMissesStrain*fe_values.JxW(q);
	    }
	    for(unsigned int j=0;j<dofs_per_cell;j++){
	      const unsigned int cj = fe_values.get_fe().system_to_component_index(j).first;
	      if (ci==cj){
		local_matrix(i,j)+= fe_values.shape_value(i,q)*fe_values.shape_value(j,q)*fe_values.JxW(q);
	      }
	    }
	  }
	}
	
	//assemble
	constraints_L2.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, mass_matrix, system_rhs);
      }
    mass_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);

    //solve
    //solveIteration(true);
   
    //
    
  }

  //Solve
  template <int dim>
  void elasticity<dim>::solveIteration(bool isProject){
    TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);
    /*    
    //Iterative solvers from Petsc and Trilinos
    SolverControl solver_control (dof_handler.n_dofs(), 1e-12);
#ifdef USE_PETSC_LA
    LA::SolverGMRES solver(solver_control, mpi_communicator);
#else
    LA::SolverGMRES solver(solver_control);
#endif
    LA::MPI::PreconditionAMG preconditioner;
    LA::MPI::PreconditionAMG::AdditionalData data;
#ifdef USE_PETSC_LA
    //data.symmetric_operator = true;
#else
    // Trilinos defaults are good 
#endif
    preconditioner.initialize(system_matrix, data);
    solver.solve (system_matrix, completely_distributed_solution, system_rhs, preconditioner);
    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;
    */
    //Direct solver MUMPS
    SolverControl cn;
    PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
    if(!isProject){
      //solver.set_symmetric_mode(false);
      //system_matrix.print(std::cout);
      solver.solve(system_matrix, completely_distributed_solution, system_rhs);
      //
      if ((currentIteration==0)){
	constraints.distribute (completely_distributed_solution);
      }
      else{
	constraints2.distribute (completely_distributed_solution);
      }
      locally_relevant_solution = completely_distributed_solution;
      dU = completely_distributed_solution;
    }
    else
      {
	//solver.set_symmetric_mode(true);
	solver.solve(mass_matrix, completely_distributed_solution, system_rhs);
	constraints_L2.distribute (completely_distributed_solution);
	locally_relevant_solution_L2 = completely_distributed_solution;
	U_L2 = completely_distributed_solution;
	UGhost_L2 = U_L2;
      }
  }

  //Solve
  template <int dim>
  void elasticity<dim>::solve(){
    double res=1, tol=1.0e-8, abs_tol=1.0e-14, initial_norm=0, current_norm=0;
    double machineEPS=1.0e-15;
    currentIteration=0;
    char buffer[200];
    while (true){
      if (currentIteration>=200){sprintf(buffer, "maximum number of iterations reached without convergence. \n"); pcout<<buffer; break; exit (1);}
      if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); pcout<<buffer; break; exit (1);}
      assemble_system();
      current_norm=system_rhs.l2_norm();
      initial_norm=std::max(initial_norm, current_norm);
      res=current_norm/initial_norm;
      sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); pcout<<buffer; 
      if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); pcout<<buffer; break;}
      solveIteration();
      U+=dU; UGhost=U; 
      ++currentIteration;
      //break;
    }
    Un=U; UnGhost=Un;
  }

  //Output
  template <int dim>
  void elasticity<dim>::output_results (const unsigned int cycle, bool isProject) {
    TimerOutput::Scope t(computing_timer, "output");
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
   
    data_out.add_data_vector (UnGhost, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
    //data_out.add_data_vector (UGhost_L2, nodal_solution_names_L2, DataOut<dim>::type_dof_data, nodal_data_component_interpretation_L2);

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
  void elasticity<dim>::run (){
    //setup problem geometry and mesh
    GridGenerator::hyper_cube (triangulation, -problemWidth/2.0, problemWidth/2.0, true);
    triangulation.refine_global (refinementFactor);
    grain_generation();
    setup_system ();
    pcout << "   Number of active cells:       "
	  << triangulation.n_global_active_cells()
	  << std::endl
	  << "   Number of degrees of freedom: "
	  << dof_handler.n_dofs()
	  << std::endl;
   
    //setup initial conditions
    
    VectorTools::interpolate(dof_handler, InitialCondition<dim>(&grain_seeds, &grain_ID), U); Un=U;
    //U=0.0;
    // 
    //sync ghost vectors to non-ghost vectors
    UGhost=U;  UnGhost=Un;
    output_results (0);
    
    //Time stepping
    currentIncrement=0;
    for (currentTime=0; currentTime<totalTime; currentTime+=dt){
      currentIncrement++;
     
      applyBoundaryConditions(currentIncrement);

      solve();
      //l2_projection();
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
      using namespace elasticity1;
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      elasticity<DIMS> problem;
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
