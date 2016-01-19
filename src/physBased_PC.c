#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"
#include "physBased_PC.h"
#include <unistd.h>
#include "pcSW.h"
#include "pcWave.h"

void GlobalInit_PBPC(PB_PC* pb_pc, LinearSolver* lsol, Simulation* simu, int* mat2assemble){

  switch(lsol->pc_type){
    case PHY_BASED_P1: 
      if(lsol->problem_type == WAVE){
        Init_PBPC_Wave_SchurPressure_BCVelocity(simu, pb_pc, mat2assemble);}
      else if(lsol->problem_type == SW){
        printf("The Shallow Water problem has no implementation of the Schur on the pressure. Please, choose another preconditioner!\n");
        exit(1);}//Init_PBPC_SW_SchurPressure_BCVelocity(simu, pb_pc, mat2assemble);}
      else{
        printf("Please, provide the implementation of the specified problem, or choose a different one.\n");
        exit(1);}
      break;

    case PHY_BASED_P2:
      if(lsol->problem_type == WAVE){
        Init_PBPC_Wave_SchurPressure_BCPressure(simu, pb_pc, mat2assemble);}
      else if(lsol->problem_type == SW){
        printf("The Shallow Water problem has no implementation of the Schur on the pressure. Please, choose another preconditioner!\n");
        exit(1);}//Init_PBPC_SW_SchurPressure_BCVelocity(simu, pb_pc, mat2assemble);}
      else{
        printf("Please, provide the implementation of the specified problem, or choose a different one.\n");
        exit(1);}
      break;
     
    case PHY_BASED_U1:
      if(lsol->problem_type == WAVE){
        Init_PBPC_Wave_SchurVelocity_BCVelocity(simu, pb_pc, mat2assemble);}
      else if(lsol->problem_type == SW){
        Init_PBPC_SW_SchurVelocity_BCVelocity(simu, pb_pc, mat2assemble);}
      else{
        printf("Please, provide the implementation of the specified problem, or choose a different one.\n");
        exit(1);}
      break;

    case PHY_BASED_U2:
      if(lsol->problem_type == WAVE){
        Init_PBPC_Wave_SchurVelocity_BCPressure(simu, pb_pc, mat2assemble);}
      else if(lsol->problem_type == SW){
        Init_PBPC_SW_SchurVelocity_BCPressure(simu, pb_pc, mat2assemble);}
      else{
        printf("Please, provide the implementation of the specified problem, or choose a different one.\n");
        exit(1);}
      break;

    default: printf("No preconditioner was specified. Aborting...\n");
             exit(1);
  }

  Init_Parameters_PhyBasedPC(pb_pc);
}

void Init_Parameters_PhyBasedPC(PB_PC* pb_pc){

  pb_pc->solver_prediction=LU;//PAR_CG;
  pb_pc->solver_propagation=LU;//PAR_CG;
  pb_pc->solver_correction=LU;//PAR_CG;

  pb_pc->pc_prediction=PAR_JACOBI;
  pb_pc->pc_propagation=PAR_JACOBI;//MULTICOLOREDGS;//PAR_ILU;//JACOBI;
  pb_pc->pc_correction=PAR_JACOBI;

  pb_pc->tol_prediction=1.e-9;
  pb_pc->tol_propagation=1e-10;
  pb_pc->tol_correction=1.e-9;

  pb_pc->itermax_prediction=1000;
  pb_pc->itermax_propagation=1000;
  pb_pc->itermax_correction=1000;
  
  pb_pc->restart_prediction=30;
  pb_pc->restart_propagation=30;
  pb_pc->restart_correction=30;
}


void PhyBased_PC_CG(void* pb_pc, Simulation *simu, real* globalSol, real*globalRHS){
  
  PB_PC* pc = (PB_PC*) pb_pc;

  // 1) PREDICTION STEP

  pc->D.lsol.solver_type=pc->solver_prediction;
  pc->D.lsol.tol=pc->tol_prediction;
  pc->D.lsol.pc_type=pc->pc_prediction;
  pc->D.lsol.iter_max=pc->itermax_prediction;
  pc->D.lsol.restart_gmres=pc->restart_prediction;

  //printf("RHS assembly.....\n");
  for (int i=0;i<pc->D.nb_fe_nodes;i++){
    pc->D.lsol.rhs[i] = globalRHS[i*3];
  }
  //printf("Solution...\n");
  //DisplayLinearSolver(&pc->D.lsol);
  Advanced_SolveLinearSolver(&pc->D.lsol,simu);

 
  // 2) PROPAGATION STEP

  pc->Schur.lsol.solver_type=pc->solver_propagation;
  pc->Schur.lsol.tol=pc->tol_propagation;
  pc->Schur.lsol.pc_type=pc->pc_propagation;
  pc->Schur.lsol.iter_max=pc->itermax_propagation;
  pc->Schur.lsol.restart_gmres=pc->restart_propagation;
  
  // Parsing L1P, L2P into the "sol" of L1 and L2 (since it is unused).
  pc->L1.lsol.MatVecProduct(&pc->L1.lsol,pc->D.lsol.sol,pc->L1.lsol.sol);
  pc->L2.lsol.MatVecProduct(&pc->L2.lsol,pc->D.lsol.sol,pc->L2.lsol.sol);

  //printf("RHS assembly.....\n");
  for (int i=0;i<pc->D.nb_fe_nodes;i++){
    pc->Schur.lsol.rhs[i*2]   = globalRHS[i*3+1] - pc->L1.lsol.sol[i];
    pc->Schur.lsol.rhs[i*2+1] = globalRHS[i*3+2] - pc->L2.lsol.sol[i];
  }
  
  //DisplayLinearSolver(&pc->Schur.lsol);
  Advanced_SolveLinearSolver(&pc->Schur.lsol,simu);

  // 3) CORRECTION STEP
  // Extracting both U1 and U2 from the previous solution of the propagation step
  real *solU1=calloc(pc->U1.nb_fe_nodes, sizeof(real));
  real *solU2=calloc(pc->U2.nb_fe_nodes, sizeof(real));
  extract2CGVectors(&pc->U1,&pc->U2,pc->Schur.lsol.sol,solU1,solU2);

  pc->D.lsol.solver_type=pc->solver_correction;
  pc->D.lsol.tol=pc->tol_correction;
  pc->D.lsol.pc_type=pc->pc_correction;
  pc->D.lsol.iter_max=pc->itermax_correction;
  pc->D.lsol.restart_gmres=pc->restart_correction;
  
  pc->U1.lsol.MatVecProduct(&pc->U1.lsol,solU1,pc->U1.lsol.sol);
  pc->U2.lsol.MatVecProduct(&pc->U2.lsol,solU2,pc->U2.lsol.sol);
  pc->D.lsol.MatVecProduct(&pc->D.lsol,pc->D.lsol.sol,pc->D.lsol.rhs);

  //printf("RHS assembly.....\n");
  for (int i=0;i<pc->D.nb_fe_nodes;i++){
    pc->D.lsol.rhs[i] += - pc->U1.lsol.sol[i] - pc->U2.lsol.sol[i];
  }
  //printf("Solution...\n");
  //DisplayLinearSolver(&pc->D.lsol);
  Advanced_SolveLinearSolver(&pc->D.lsol,simu);
  
 
  // 4) OUTPUT STEP Final concatenation
  cat2CGVectors(&pc->D,&pc->Schur,pc->D.lsol.sol,pc->Schur.lsol.sol,globalSol);

  free(solU1);
  free(solU2);
}





void PhyBased_PC_InvertSchur_CG(void* pb_pc, Simulation *simu, real* globalSol, real* globalRHS){
  
  PB_PC* pc = (PB_PC*) pb_pc;

  // 1) PREDICTION STEP
  printf(" prediction\n");
  pc->D.lsol.solver_type=pc->solver_prediction;
  pc->D.lsol.tol=pc->tol_prediction;
  pc->D.lsol.pc_type=pc->pc_prediction;
  pc->D.lsol.iter_max=pc->itermax_prediction;
  pc->D.lsol.restart_gmres=pc->restart_prediction;

  //printf("RHS assembly.....\n");
  for (int i=0;i<pc->D.nb_fe_nodes;i++){
    pc->D.lsol.rhs[i*2]   = globalRHS[i*3+1];
    pc->D.lsol.rhs[i*2+1] = globalRHS[i*3+2];
  }
  //printf("Solution prediction\n");
  Advanced_SolveLinearSolver(&pc->D.lsol,simu);
  /*for (int i=0;i<pc->D.nb_fe_nodes;i++){
    printf(" 1) U[%d] = %8e, V[%d] = %8e\n", i, pc->D.lsol.sol[i*2], i, pc->D.lsol.sol[i*2+1]);
    }*/

  // 2) PROPAGATION STEP
   printf(" propagation\n");
  pc->Schur.lsol.solver_type=pc->solver_propagation;
  pc->Schur.lsol.tol=pc->tol_propagation;
  pc->Schur.lsol.pc_type=pc->pc_propagation;
  pc->Schur.lsol.iter_max=pc->itermax_propagation;
  pc->Schur.lsol.restart_gmres=pc->restart_propagation;

  real *solU1=calloc(pc->U1.nb_fe_nodes, sizeof(real));
  real *solU2=calloc(pc->U2.nb_fe_nodes, sizeof(real));
  extract2CGVectors(&pc->U1,&pc->U2,pc->D.lsol.sol,solU1,solU2);
  
  // Parsing L1P, L2P into the "sol" of L1 and L2 (since it is unused).
  pc->L1.lsol.MatVecProduct(&pc->L1.lsol,solU1,pc->L1.lsol.sol);
  pc->L2.lsol.MatVecProduct(&pc->L2.lsol,solU2,pc->L2.lsol.sol);

  //printf("RHS assembly.....\n");
  for (int i=0;i<pc->D.nb_fe_nodes;i++){
    pc->Schur.lsol.rhs[i]   = globalRHS[i*3] - pc->L1.lsol.sol[i] - pc->L2.lsol.sol[i];
  }
  
  //printf("Solution propagation\n");
  Advanced_SolveLinearSolver(&pc->Schur.lsol,simu);
  /*for (int i=0;i<pc->D.nb_fe_nodes;i++){
    printf(" 2) P[%d] = %8e\n", i, pc->Schur.lsol.sol[i]);
    }*/

  
  // 3) CORRECTION STEP
 printf(" correction\n");
  pc->D.lsol.solver_type=pc->solver_correction;
  pc->D.lsol.tol=pc->tol_correction;
  pc->D.lsol.pc_type=pc->pc_correction;
  pc->D.lsol.iter_max=pc->itermax_correction;
  pc->D.lsol.restart_gmres=pc->restart_correction;
  
  pc->U1.lsol.MatVecProduct(&pc->U1.lsol,pc->Schur.lsol.sol,pc->U1.lsol.sol);
  pc->U2.lsol.MatVecProduct(&pc->U2.lsol,pc->Schur.lsol.sol,pc->U2.lsol.sol);
  pc->D.lsol.MatVecProduct(&pc->D.lsol,pc->D.lsol.sol,pc->D.lsol.rhs);

  //printf("RHS assembly.....\n");
  for (int i=0;i<pc->D.nb_fe_nodes;i++){
    pc->D.lsol.rhs[i*2] = pc->D.lsol.rhs[i*2]- pc->U1.lsol.sol[i];
    pc->D.lsol.rhs[i*2+1] =  pc->D.lsol.rhs[i*2+1] - pc->U2.lsol.sol[i];
  }

  //printf("Solution correction\n");
  Advanced_SolveLinearSolver(&pc->D.lsol,simu);
  /*for (int i=0;i<pc->D.nb_fe_nodes;i++){
    printf("Prout Schur P : P[%d] = %8e, U[%d] = %8e, V[%d] = %8e\n", i, pc->Schur.lsol.sol[i], i, pc->D.lsol.sol[2*i], i, pc->D.lsol.sol[2*i+1]);
    }*/

  
  // 4) OUTPUT STEP Final concatenation
  cat2CGVectors(&pc->D,&pc->Schur,pc->D.lsol.sol,pc->Schur.lsol.sol,globalSol);

  free(solU1);
  free(solU2);
}

void GenericOperator_PBPC_NonLinear(void* pb_pc, ContinuousSolver* cs){

  ContinuousSolver* ps = cs;
  PB_PC* pc = (PB_PC*) pb_pc;
  int nnodes = pc->D.nnodes;
  field* f0 = &pc->D.simu->fd[0];
  // TODO This boolean is currently useless, since we rebuild all the matrices
  // at each time-step. It is even unclear if in every stage of the process
  // these booleans are not overwritten...
  bool assembled = pc->D.lsol.mat_is_assembly && 
                   pc->L1.lsol.mat_is_assembly &&
                   pc->L2.lsol.mat_is_assembly &&
                   pc->U1.lsol.mat_is_assembly &&
                   pc->U2.lsol.mat_is_assembly &&
                   pc->Schur.lsol.mat_is_assembly;
  int nbel = pc->D.nbel;
  int dg_to_fe_index[pc->D.nb_dg_nodes];
  for (int i=0; i<pc->D.nb_dg_nodes; i++){
    dg_to_fe_index[i]=pc->D.dg_to_fe_index[i];
  }

  if(!assembled){
    for(int ie = 0; ie < nbel; ie++){  

        
      int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
      int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);

      // Build Derivatives
      // value = 0, dx = 1, dy = 2, dz = 3
      real var[ps->nb_phy_vars][4][ps->nnodes]; 
      for(int iv = 0; iv<ps->nb_phy_vars; iv++){
        for(int i_der = 0; i_der<4; i_der++){
          for(int iloc = 0; iloc < ps->nnodes; iloc++){
            var[iv][i_der][iloc] = 0.0;
          }
        }
      }
      // Computation of the variable and its derivatives for each Gauss point
      for(int iloc = 0; iloc < ps->nnodes; iloc++){
        real wpg;
        real xref[3];
        real dtau[3][3],codtau[3][3];
        real dphiref_j[3];
        real dphi_j[3];
        real basisPhi_j[4];

        int ipgmacro = iloc + isubcell * ps->nnodes;
        ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);
        
        // Loop on variables
        for(int iv=0; iv<ps->nb_phy_vars; iv++){
          int ino_dg = iloc + ie * ps->nnodes;
          int ino_fe = ps->dg_to_fe_index[ino_dg];
          var[iv][0][iloc] = ps->lsol.sol[iv+ino_fe*ps->nb_phy_vars];

          for(int jGauss = 0; jGauss < ps->nnodes; jGauss++){
            int jGaussmacro = jGauss + isubcell * ps->nnodes;
            int jGauss_dg = jGauss + ie * ps->nnodes;
            int jGauss_fe = ps->dg_to_fe_index[jGauss_dg];
            grad_psi_pg(f0->deg,f0->raf,jGaussmacro,ipgmacro,dphiref_j);
            Ref2Phy(ps->simu->fd[iemacro].physnode,
                    xref,dphiref_j,0,NULL,
                    dtau,codtau,dphi_j,NULL);
            real det = dot_product(dtau[0], codtau[0]);
            basisPhi_j[1]=dphi_j[0]/det;
            basisPhi_j[2]=dphi_j[1]/det;
            basisPhi_j[3]=dphi_j[2]/det;	    
            for (int i=1;i<4;i++){
              var[iv][i][iloc] += basisPhi_j[i]*ps->lsol.sol[iv+jGauss_fe*ps->nb_phy_vars];
            } // end for i (derivative)
          } // end for jGauss
        } // end for iv
      } // end for iloc

      // Init: local matrix and rhs
      real aloc[ps->nnodes*ps->nb_phy_vars][ps->nnodes*ps->nb_phy_vars];
      real rhsloc[ps->nnodes*ps->nb_phy_vars];
      for(int iloc = 0; iloc < ps->nnodes*ps->nb_phy_vars; iloc++){
        rhsloc[iloc] = 0.0;
        for(int jloc = 0; jloc < ps->nnodes*ps->nb_phy_vars; jloc++){
          aloc[iloc][jloc] = 0.0;
        }
      }

      /////////////////////////////////////////////////
      /////////// Now compute matrix and RHS //////////
      /////////////////////////////////////////////////
      for(int ipg = 0;ipg < nnodes; ipg++){
        real wpg;
        real xref[3];
        int ipgmacro = ipg + isubcell * nnodes;
        
        ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);
        real loc_rhs[ps->nb_phy_vars];
        for (int i=0; i<ps->nb_phy_vars; i++){
          loc_rhs[i]=0;
        }
        
        for(int iloc = 0; iloc < nnodes; iloc++){
          real dtau[3][3],codtau[3][3];
          real dphiref_i[3],dphiref_j[3];
          real dphi_i[3],dphi_j[3];
          real basisPhi_i[4], basisPhi_j[4];
          int ilocmacro = iloc + isubcell * nnodes;
          int ino_dg = iloc + ie * nnodes;
          int ino_fe = dg_to_fe_index[ino_dg];
          grad_psi_pg(f0->deg,f0->raf,ilocmacro,ipgmacro,dphiref_i);
          Ref2Phy(pc->D.simu->fd[iemacro].physnode,
                  xref,dphiref_i,0,NULL,
                  dtau,codtau,dphi_i,NULL);
        
          real det = dot_product(dtau[0], codtau[0]);
          if (ilocmacro==ipgmacro){
            basisPhi_i[0]=1;
          }
          else
          {
            basisPhi_i[0]=0;
          }
          basisPhi_i[1]=dphi_i[0]/det;
          basisPhi_i[2]=dphi_i[1]/det;
          basisPhi_i[3]=dphi_i[2]/det;

          ////////////////////////////////////
          ///// BUILDING LOCAL MATRICES //////
          ////////////////////////////////////
          real h  = var[0][0][iloc];
          real hx = var[0][1][iloc];
          real hy = var[0][2][iloc];
          real hz = var[0][3][iloc];
          real u  = var[1][0][iloc];
          real ux = var[1][1][iloc];
          real uy = var[1][2][iloc];
          real uz = var[1][3][iloc];
          real v  = var[2][0][iloc];
          real vx = var[2][1][iloc];
          real vy = var[2][2][iloc];
          real vz = var[2][3][iloc];
          
          real var2[12] = {h, hx, hy, hz, u, ux, uy, uz, v, vx, vy ,vz};
   
          pc->loc_mat_assembly(pc,var2);

          for(int jloc = 0; jloc < nnodes; jloc++){
            int jlocmacro = jloc + isubcell * nnodes;
            int jno_dg = jloc + ie * nnodes;
            int jno_fe = dg_to_fe_index[jno_dg];
            grad_psi_pg(f0->deg,f0->raf,jlocmacro,ipgmacro,dphiref_j);
            Ref2Phy(pc->D.simu->fd[iemacro].physnode,
                    xref,dphiref_j,0,NULL,
                    dtau,codtau,dphi_j,NULL);
            if (jlocmacro==ipgmacro){
              basisPhi_j[0]=1;
            }
            else
            {
              basisPhi_j[0]=0;
            }
            basisPhi_j[1]=dphi_j[0]/det;
            basisPhi_j[2]=dphi_j[1]/det;
            basisPhi_j[3]=dphi_j[2]/det;
            real val;
            real res[4];
            ContinuousSolver * cs;

            /////////////////////////////////////
            ///// BUILDING GLOBAL MATRICES //////
            /////////////////////////////////////
            // Building D Matrix
            cs = &pc->D;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2
            

            // Building L1 Matrix
            cs = &pc->L1;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

            // Building L2 Matrix
            cs = &pc->L2;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

            // Building U1 Matrix
            cs = &pc->U1;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

            // Building U2 Matrix
            cs = &pc->U2;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

            // Building Schur Matrix
            cs = &pc->Schur;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

          } // end for jloc
        } // end for iloc
      } // end for ipg
    } // end for ie
  } // end if assembled

  // Applying boundaries right after building the preconditioner matrices, because there is no 
  // need to differ them (RHS is built otherwise, and bcs are applied for them too) 
  // -> Check penalization parameters to be corresponding one to another between preconditioner
  // and full problem.
  boundary_assembly(pc);
}


void GenericOperator_PBPC(void* pb_pc, ContinuousSolver* cs){

  PB_PC* pc = (PB_PC*) pb_pc;
  int nnodes = pc->D.nnodes;
  field* f0 = &pc->D.simu->fd[0];
  bool assembled = pc->D.lsol.mat_is_assembly && 
                   pc->L1.lsol.mat_is_assembly &&
                   pc->L2.lsol.mat_is_assembly &&
                   pc->U1.lsol.mat_is_assembly &&
                   pc->U2.lsol.mat_is_assembly &&
                   pc->Schur.lsol.mat_is_assembly;
  int nbel = pc->D.nbel;
  int dg_to_fe_index[pc->D.nb_dg_nodes];
  for (int i=0; i<pc->D.nb_dg_nodes; i++){
    dg_to_fe_index[i]=pc->D.dg_to_fe_index[i];
  }

  if(!assembled){
    for(int ie = 0; ie < nbel; ie++){  

        
      int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
      int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);
        
      for(int ipg = 0;ipg < nnodes; ipg++){
        real wpg;
        real xref[3];
        int ipgmacro = ipg + isubcell * nnodes;
        
        ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);
        
        for(int iloc = 0; iloc < nnodes; iloc++){
          real dtau[3][3],codtau[3][3];
          real dphiref_i[3],dphiref_j[3];
          real dphi_i[3],dphi_j[3];
          real basisPhi_i[4], basisPhi_j[4];
          int ilocmacro = iloc + isubcell * nnodes;
          int ino_dg = iloc + ie * nnodes;
          int ino_fe = dg_to_fe_index[ino_dg];
          grad_psi_pg(f0->deg,f0->raf,ilocmacro,ipgmacro,dphiref_i);
          Ref2Phy(pc->D.simu->fd[iemacro].physnode,
                  xref,dphiref_i,0,NULL,
                  dtau,codtau,dphi_i,NULL);
        
          real det = dot_product(dtau[0], codtau[0]);
          if (ilocmacro==ipgmacro){
            basisPhi_i[0]=1;
          }
          else
          {
            basisPhi_i[0]=0;
          }
          basisPhi_i[1]=dphi_i[0]/det;
          basisPhi_i[2]=dphi_i[1]/det;
          basisPhi_i[3]=dphi_i[2]/det;
          for(int jloc = 0; jloc < nnodes; jloc++){
            int jlocmacro = jloc + isubcell * nnodes;
            int jno_dg = jloc + ie * nnodes;
            int jno_fe = dg_to_fe_index[jno_dg];
            grad_psi_pg(f0->deg,f0->raf,jlocmacro,ipgmacro,dphiref_j);
            Ref2Phy(pc->D.simu->fd[iemacro].physnode,
                    xref,dphiref_j,0,NULL,
                    dtau,codtau,dphi_j,NULL);
            if (jlocmacro==ipgmacro){
              basisPhi_j[0]=1;
            }
            else
            {
              basisPhi_j[0]=0;
            }
            basisPhi_j[1]=dphi_j[0]/det;
            basisPhi_j[2]=dphi_j[1]/det;
            basisPhi_j[3]=dphi_j[2]/det;
            real val;
            real res[4];
            ContinuousSolver * cs ;

            // Building D Matrix
            cs = &pc->D;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

            // Building L1 Matrix
            cs = &pc->L1;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

            // Building L2 Matrix
            cs = &pc->L2;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

            // Building U1 Matrix
            cs = &pc->U1;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

            // Building U2 Matrix
            cs = &pc->U2;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

            // Building Schur Matrix
            cs = &pc->Schur;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

          } // end for jloc
        } // end for iloc
      } // end for ipg
    } // end for ie
  } // end if assembled
}

void GenericRHS_PBPC_NonLinear(void* pb_pc, ContinuousSolver* cs){

  ContinuousSolver* ps = cs;
  PB_PC* pc = (PB_PC*) pb_pc;
  int nnodes = pc->D.nnodes;
  field* f0 = &pc->D.simu->fd[0];
  int nbel = pc->D.nbel;
  int dg_to_fe_index[pc->D.nb_dg_nodes];
  for (int i=0; i<pc->D.nb_dg_nodes; i++){
    dg_to_fe_index[i]=pc->D.dg_to_fe_index[i];
  }

  for(int ie = 0; ie < nbel; ie++){  

      
    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);


    // Build Derivatives
    // value = 0, dx = 1, dy = 2, dz = 3
    real var[ps->nb_phy_vars][4][ps->nnodes]; 
    for(int iv = 0; iv<ps->nb_phy_vars; iv++){
      for(int i_der = 0; i_der<4; i_der++){
        for(int iloc = 0; iloc < ps->nnodes; iloc++){
          var[iv][i_der][iloc] = 0.0;
        }
      }
    }
    // Computation of the variable and its derivatives for each Gauss point
    for(int iloc = 0; iloc < ps->nnodes; iloc++){
      real wpg;
      real xref[3];
      real dtau[3][3],codtau[3][3];
      real dphiref_j[3];
      real dphi_j[3];
      real basisPhi_j[4];

      int ipgmacro = iloc + isubcell * ps->nnodes;
      ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);
      
      // Loop on variables
      for(int iv=0; iv<ps->nb_phy_vars; iv++){
        int ino_dg = iloc + ie * ps->nnodes;
        int ino_fe = ps->dg_to_fe_index[ino_dg];
        var[iv][0][iloc] = ps->lsol.sol[iv+ino_fe*ps->nb_phy_vars];

        for(int jGauss = 0; jGauss < ps->nnodes; jGauss++){
          int jGaussmacro = jGauss + isubcell * ps->nnodes;
          int jGauss_dg = jGauss + ie * ps->nnodes;
          int jGauss_fe = ps->dg_to_fe_index[jGauss_dg];
          grad_psi_pg(f0->deg,f0->raf,jGaussmacro,ipgmacro,dphiref_j);
          Ref2Phy(ps->simu->fd[iemacro].physnode,
                  xref,dphiref_j,0,NULL,
                  dtau,codtau,dphi_j,NULL);
          real det = dot_product(dtau[0], codtau[0]);
          basisPhi_j[1]=dphi_j[0]/det;
          basisPhi_j[2]=dphi_j[1]/det;
          basisPhi_j[3]=dphi_j[2]/det;
          for (int i=1;i<4;i++){
            var[iv][i][iloc] += basisPhi_j[i]*ps->lsol.sol[iv+jGauss_fe*ps->nb_phy_vars];
          } // end for i (derivative)
        } // end for jGauss
      } // end for iv
    } // end for iloc

    // Init: local rhs
    real rhsloc[ps->nnodes*ps->nb_phy_vars];
    for(int iloc = 0; iloc < ps->nnodes*ps->nb_phy_vars; iloc++){
      rhsloc[iloc] = 0.0;
    }

    //////////////////////////////////////
    /////////// Now compute RHS //////////
    //////////////////////////////////////
    for(int ipg = 0;ipg < nnodes; ipg++){
      real wpg;
      real xref[3];
      int ipgmacro = ipg + isubcell * nnodes;
      
      ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);
      real loc_rhs[ps->nb_phy_vars];
      for (int i=0; i<ps->nb_phy_vars; i++){
        loc_rhs[i]=0;
      }
      
      for(int iloc = 0; iloc < nnodes; iloc++){
        real dtau[3][3],codtau[3][3];
        real dphiref_i[3],dphiref_j[3];
        real dphi_i[3],dphi_j[3];
        real basisPhi_i[4], basisPhi_j[4];
        int ilocmacro = iloc + isubcell * nnodes;
        int ino_dg = iloc + ie * nnodes;
        int ino_fe = dg_to_fe_index[ino_dg];
        grad_psi_pg(f0->deg,f0->raf,ilocmacro,ipgmacro,dphiref_i);
        Ref2Phy(pc->D.simu->fd[iemacro].physnode,
                xref,dphiref_i,0,NULL,
                dtau,codtau,dphi_i,NULL);
      
        real det = dot_product(dtau[0], codtau[0]);
        if (ilocmacro==ipgmacro){
          basisPhi_i[0]=1;
        }
        else
        {
          basisPhi_i[0]=0;
        }
        basisPhi_i[1]=dphi_i[0]/det;
        basisPhi_i[2]=dphi_i[1]/det;
        basisPhi_i[3]=dphi_i[2]/det;

        //////////////////////////////
        ///// BUILDING MATRICES //////
        //////////////////////////////
        real h  = var[0][0][iloc];
        real hx = var[0][1][iloc];
        real hy = var[0][2][iloc];
        real hz = var[0][3][iloc];
        real u  = var[1][0][iloc];
        real ux = var[1][1][iloc];
        real uy = var[1][2][iloc];
        real uz = var[1][3][iloc];
        real v  = var[2][0][iloc];
        real vx = var[2][1][iloc];
        real vy = var[2][2][iloc];
        real vz = var[2][3][iloc];
        
        real var2[12] = {h, hx, hy, hz, u, ux, uy, uz, v, vx, vy ,vz};
  
        pc->loc_rhs_assembly(pc,var2,ps->simu->dt*wpg*det*basisPhi_i[0],loc_rhs);

        for (int iv=0; iv<ps->nb_phy_vars; iv++){
          rhsloc[iloc*ps->nb_phy_vars+iv] += loc_rhs[iv];
        }
      } // end for iloc
    } // end for ipg

    // Building global rhs 
    for(int iloc = 0; iloc < ps->nnodes; iloc++){
      int ino_dg = iloc + ie * ps->nnodes;
      int ino_fe = ps->dg_to_fe_index[ino_dg];
      for (int iv=0; iv<ps->nb_phy_vars; iv++){
        int ipot_fe = ino_fe*ps->nb_phy_vars+iv;
        ps->lsol.rhs[ino_fe*ps->nb_phy_vars+iv] += rhsloc[iloc*ps->nb_phy_vars+iv];
      }
    }

  } // end for ie
}

void freePB_PC(PB_PC* pb_pc){
  free(pb_pc->list_mat2assemble);
  freeContinuousSolver(&pb_pc->D);
  freeContinuousSolver(&pb_pc->L1);
  freeContinuousSolver(&pb_pc->L2);
  freeContinuousSolver(&pb_pc->U1);
  freeContinuousSolver(&pb_pc->U2);
  freeContinuousSolver(&pb_pc->Schur);
  freeContinuousSolver(&pb_pc->fullSolver);
  free(pb_pc->rhs_prediction);
  free(pb_pc->rhs_propagation);
  free(pb_pc->rhs_correction);
  
}

void reset(PB_PC* pb_pc){
  
  for (int i=0;i<pb_pc->D.nb_fe_dof; i++){
    pb_pc->D.lsol.sol[i] = 0.0;
    pb_pc->D.lsol.rhs[i] = 0.0;
  }
  for (int i=0;i<pb_pc->L1.nb_fe_dof; i++){
    pb_pc->L1.lsol.sol[i] = 0.0;
    pb_pc->L1.lsol.rhs[i] = 0.0;
  }
  for (int i=0;i<pb_pc->L2.nb_fe_dof; i++){
    pb_pc->L2.lsol.sol[i] = 0.0;
    pb_pc->L2.lsol.rhs[i] = 0.0;
  }
  for (int i=0;i<pb_pc->U1.nb_fe_dof; i++){
    pb_pc->U1.lsol.sol[i] = 0.0;
    pb_pc->U1.lsol.rhs[i] = 0.0;
  }
  for (int i=0;i<pb_pc->U2.nb_fe_dof; i++){
    pb_pc->U2.lsol.sol[i] = 0.0;
    pb_pc->U2.lsol.rhs[i] = 0.0;
  }
  for (int i=0;i<pb_pc->Schur.nb_fe_dof; i++){
    pb_pc->Schur.lsol.sol[i] = 0.0;
    pb_pc->Schur.lsol.rhs[i] = 0.0;
  }
  //for (int i=0;i<pb_pc->D.nb_fe_nodes; i++){
  //  // resetting sols
  //  for(int ivar=0;ivar<pb_pc->D.nb_phy_vars; ivar++){
  //    pb_pc->D.lsol.sol[i*pb_pc->D.nb_phy_vars+ivar]=0;
  //    pb_pc->D.lsol.rhs[i*pb_pc->D.nb_phy_vars+ivar]=0;
  //  }
  //   for(int ivar=0;ivar<pb_pc->L1.nb_phy_vars; ivar++){
  //    pb_pc->L1.lsol.sol[i*pb_pc->L1.nb_phy_vars+ivar]=0;
  //    pb_pc->L1.lsol.rhs[i*pb_pc->L1.nb_phy_vars+ivar]=0;
  //  }
  //   for(int ivar=0;ivar<pb_pc->L2.nb_phy_vars; ivar++){
  //    pb_pc->L2.lsol.sol[i*pb_pc->L2.nb_phy_vars+ivar]=0;
  //    pb_pc->L2.lsol.rhs[i*pb_pc->L2.nb_phy_vars+ivar]=0;
  //  }
  //   for(int ivar=0;ivar<pb_pc->U1.nb_phy_vars; ivar++){
  //    pb_pc->U1.lsol.sol[i*pb_pc->U1.nb_phy_vars+ivar]=0;
  //    pb_pc->U1.lsol.rhs[i*pb_pc->U1.nb_phy_vars+ivar]=0;
  //  }
  //   for(int ivar=0;ivar<pb_pc->U2.nb_phy_vars; ivar++){
  //    pb_pc->U2.lsol.sol[i*pb_pc->U2.nb_phy_vars+ivar]=0;
  //    pb_pc->U2.lsol.rhs[i*pb_pc->U2.nb_phy_vars+ivar]=0;
  //  }
  //   for(int ivar=0;ivar<pb_pc->Schur.nb_phy_vars; ivar++){
  //    pb_pc->Schur.lsol.sol[i*pb_pc->Schur.nb_phy_vars+ivar]=0;
  //    pb_pc->Schur.lsol.rhs[i*pb_pc->Schur.nb_phy_vars+ivar]=0;
  //  }
  //}
  
  if (pb_pc->nonlinear){
    pb_pc->D.reset_dt(&pb_pc->D.lsol);
    pb_pc->L1.reset_dt(&pb_pc->L1.lsol);
    pb_pc->L2.reset_dt(&pb_pc->L2.lsol);
    pb_pc->U1.reset_dt(&pb_pc->U1.lsol);
    pb_pc->U2.reset_dt(&pb_pc->U2.lsol);
    pb_pc->Schur.reset_dt(&pb_pc->Schur.lsol);
  }
}



void PiDgToCg(ContinuousSolver * cs, int nbVarIn, real * rhsIn, real * rhsOut){
  
 field* f = &cs->simu->fd[0];

 real ** Restriction;
 real **tabvar;
 real **tabvar_transform;
 real* coeff = calloc(cs->nb_fe_nodes, sizeof(real));

 Restriction = calloc(cs->nb_fe_nodes,sizeof(real*));
 for (int i=0;i<cs->nb_fe_nodes;i++){
   Restriction[i] = calloc(cs->nb_dg_nodes,sizeof(real));
 }
 
 tabvar = calloc(cs->nb_phy_vars,sizeof(real*));
 for (int i=0;i<cs->nb_phy_vars;i++){
   tabvar[i] = calloc(cs->nb_dg_nodes,sizeof(real));
 }

 tabvar_transform = calloc(cs->nb_phy_vars,sizeof(real*));
 for (int i=0;i<cs->nb_phy_vars;i++){
   tabvar_transform[i] = calloc(cs->nb_fe_nodes,sizeof(real));
 }


 for (int i=0;i<cs->nb_dg_nodes;i++){
   for (int iv=0;iv<cs->nb_phy_vars;iv++){
     tabvar[iv][i] = rhsIn[cs->list_of_var[iv]+i*nbVarIn];
   }
 }

for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
   for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
     Restriction[ino_fe][ino_dg]=0.;
   }
 }

 for(int ie = 0; ie < cs->nbel; ie++){  

    int iemacro  = ie / (f->raf[0] * f->raf[1] * f->raf[2]);
    int isubcell = ie % (f->raf[0] * f->raf[1] * f->raf[2]);
    
    for(int iloc = 0; iloc < cs->nnodes; iloc++){
      real wpg;
      real xref[3];
      int ilocmacro = iloc + isubcell * cs->nnodes;
      ref_pg_vol(f->deg,f->raf,ilocmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      Ref2Phy(cs->simu->fd[iemacro].physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      real det = dot_product(dtau[0], codtau[0]);	
      int ino_dg = iloc + ie * cs->nnodes;
      int ino_fe = cs->dg_to_fe_index[ino_dg];
      Restriction[ino_fe][ino_dg]= wpg * det; // We need to multiply by the member of dg node for 1 fe node.
      coeff[ino_fe] += wpg * det; // We need to multiply by the member of dg node for 1 fe node.
    }
 }

 for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
   for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
     Restriction[ino_fe][ino_dg]/=coeff[ino_fe];
   }
 }

 for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
   for(int iv=0;iv<cs->nb_phy_vars;iv++){
	 tabvar_transform[iv][ino_fe]=0.0;
   }
   for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
     for(int iv=0;iv<cs->nb_phy_vars;iv++){
	 tabvar_transform[iv][ino_fe]+=Restriction[ino_fe][ino_dg]*tabvar[iv][ino_dg];
     }
   }
 }

 for (int i=0;i<cs->nb_fe_nodes;i++){
   for (int iv=0;iv<cs->nb_phy_vars;iv++){
     rhsOut[iv+i*cs->nb_phy_vars]=tabvar_transform[iv][i];
   }
 }

 free(coeff);

}

void PiInvertCgToDg(ContinuousSolver * cs, int nbVarOut, real * rhsIn, real * rhsOut){
  
 field* f = &cs->simu->fd[0];
 real ** Reconstruction;
 real **tabvar;
 real **tabvar_transform;
 real* coeff = calloc(cs->nb_dg_nodes, sizeof(real));

 Reconstruction = calloc(cs->nb_dg_nodes,sizeof(real*));
 for (int i=0;i<cs->nb_dg_nodes;i++){
   Reconstruction[i] = calloc(cs->nb_fe_nodes,sizeof(real));
 }
 
 tabvar = calloc(cs->nb_phy_vars,sizeof(real*));
 for (int i=0;i<cs->nb_phy_vars;i++){
   tabvar[i] = calloc(cs->nb_fe_nodes,sizeof(real));
 }

 tabvar_transform = calloc(cs->nb_phy_vars,sizeof(real*));
 for (int i=0;i<cs->nb_phy_vars;i++){
   tabvar_transform[i] = calloc(cs->nb_dg_nodes,sizeof(real));
 }

 for (int i=0;i<cs->nb_fe_nodes;i++){
   for (int iv=0;iv<cs->nb_phy_vars;iv++){
     tabvar[iv][i] = rhsIn[iv+i*cs->nb_phy_vars];
   }
 }

 for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
   for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
     Reconstruction[ino_dg][ino_fe]=0.;
   }
 }

for(int ie = 0; ie < cs->nbel; ie++){  

    int iemacro  = ie / (f->raf[0] * f->raf[1] * f->raf[2]);
    int isubcell = ie % (f->raf[0] * f->raf[1] * f->raf[2]);
    
    for(int iloc = 0; iloc < cs->nnodes; iloc++){
      real wpg;
      real xref[3];
      int ilocmacro = iloc + isubcell * cs->nnodes;
      ref_pg_vol(f->deg,f->raf,ilocmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      Ref2Phy(cs->simu->fd[iemacro].physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      real det = dot_product(dtau[0], codtau[0]);	
      int ino_dg = iloc + ie * cs->nnodes;
      int ino_fe = cs->dg_to_fe_index[ino_dg];
      Reconstruction[ino_dg][ino_fe]= wpg * det; // We need to multiply by the member of dg node for 1 fe node.
      coeff[ino_dg] = wpg * det; // We need to multiply by the member of dg node for 1 fe node.
    }
 }
 for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
   for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
     Reconstruction[ino_dg][ino_fe]/=coeff[ino_dg];
   }
 }

 for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
    for(int iv=0;iv<cs->nb_phy_vars;iv++){
	 tabvar_transform[iv][ino_dg]=0.0;
     }
     for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
       for(int iv=0;iv<cs->nb_phy_vars;iv++){
	 tabvar_transform[iv][ino_dg]+=Reconstruction[ino_dg][ino_fe]*tabvar[iv][ino_fe];
     }
   }
 }

for (int i=0;i<cs->nb_dg_nodes;i++){
   for (int iv=0;iv<cs->nb_phy_vars;iv++){
     rhsOut[cs->list_of_var[iv]+i*nbVarOut]=tabvar_transform[iv][i];
   }
 }
 free(coeff);
}


void RobinFlux_SchurPressure(void * cs, real * xpg, real * w, real *vnorm, real * flux){
  ContinuousSolver * ps=cs;
  real lambda=-1e7;
  real Coef_diff=0;
  real p0=0;
  real Win[ps->simu->fd[0].model.m];
    
  ps->simu->fd[0].model.ImposedData(xpg,ps->simu->dt,Win);
  p0=Win[0];

  Coef_diff=ps->simu->dt*ps->simu->vmax*ps->simu->theta;

  // u=- c dt theta nabla p
  // boudnary condition for wave lambda p-1/c (u,n)=g ==> * Coef_diff * (nabla p,n)+lambda p=g
  // ==> -(nabla p,n)= c* lambda / (Coeff*diff)-g *c /(Coeff*diff)
  // ==> -(Coeff*diff)**2 (nabla p,n)= c* lambda * (Coeff*diff)-g *c *(Coeff*diff)

  flux[0]=ps->simu->vmax * lambda* Coef_diff * w[0];//- ps->simu->vmax * Coef_diff * p0;
}

void Dirichlet_Velocity(void * cs, real * xpg, real * w, real *vnorm, real * flux){
  ContinuousSolver * ps=cs;
  real mu=1.0e20;//15;
  real Coef_diff=0;
  real u0_1 = 0;
  real u0_2 = 0;
  real Win[3];
  int angle=0;

  if(xpg[0]== 0.0 && xpg[1]== 0.0) { angle=1; }
  if(xpg[0]== 0.0 && xpg[1]== 1.0) { angle=1; }
  if(xpg[0]== 1.0 && xpg[1]== 0.0) { angle=1; }
  if(xpg[0]== 1.0 && xpg[1]== 1.0) { angle=1; }

  if(angle ==0){
    Coef_diff=ps->simu->dt*ps->simu->vmax*ps->simu->theta;
    flux[0]= mu * (w[0]*vnorm[0]*vnorm[0] + w[1]*vnorm[0]*vnorm[1])
      - mu * (u0_1*vnorm[0]*vnorm[0] + u0_2*vnorm[0]*vnorm[1]);
    flux[1]=  mu * (w[0]*vnorm[1]*vnorm[0] + w[1]*vnorm[1]*vnorm[1])
      -  mu * (u0_1*vnorm[1]*vnorm[0] + u0_2*vnorm[1]*vnorm[1]);
   }
  else {
    flux[0]= mu * w[0] - mu * u0_1;
    flux[1]= mu * w[1] - mu * u0_2;
  }
 
}

void Dirichlet_vorticity(void * cs, real * xpg, real * w, real *vnorm, real * flux){
  ContinuousSolver * ps=cs;
  real mu=1.0e21;//15;
  real Coef_diff=0;
  real w_0 = 0;

  int angle=0;

  flux[0]= mu * w[0] - mu * w_0;

 
}


void BoundaryTerm_Xderivative(void * cs, real * xpg, real * w, real *vnorm, real * flux){
  ContinuousSolver * ps=cs;
  real lambda=0;
  real Coef_diff=0;
  real p0=0;

  Coef_diff=ps->simu->dt*ps->simu->vmax*ps->simu->theta;
  // term c dt theta n1
  flux[0]= Coef_diff * w[0] * vnorm[0];
}

void BoundaryTerm_Yderivative(void * cs, real * xpg, real * w, real *vnorm, real * flux){
  ContinuousSolver * ps=cs; 
  real lambda=0;
  real Coef_diff=0;
  real p0=0;

  Coef_diff=ps->simu->dt*ps->simu->vmax*ps->simu->theta;
  // term c dt theta n2
  flux[0]= Coef_diff * w[0] * vnorm[1];
}

// This function is there only to ease the writing (since some of the 
// boundaries are not applied in some cases.
void boundary_assembly(PB_PC* pb_pc){
  if (    pb_pc->D.bc_assembly != NULL)     pb_pc->D.bc_assembly(&pb_pc->D);
  if (   pb_pc->L1.bc_assembly != NULL)    pb_pc->L1.bc_assembly(&pb_pc->L1);
  if (   pb_pc->L2.bc_assembly != NULL)    pb_pc->L2.bc_assembly(&pb_pc->L2);
  if (   pb_pc->U1.bc_assembly != NULL)    pb_pc->U1.bc_assembly(&pb_pc->U1);
  if (   pb_pc->U2.bc_assembly != NULL)    pb_pc->U2.bc_assembly(&pb_pc->U2);
  if (pb_pc->Schur.bc_assembly != NULL) pb_pc->Schur.bc_assembly(&pb_pc->Schur);
}
