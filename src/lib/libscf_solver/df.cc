/***************************************************************************
*
*       df.cc in psi4/src/lib/libscf_solver
*       By Rob Parrish, CCMST Georgia Tech
*       robparrish@gmail.com
*       14 June 2010
*
*       Canonical (delocalized) density fitting routines for SCF
*
*
*****************************************************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include "hf.h"
#include "rhf.h"
#include "uhf.h"
#include "rohf.h"

#include <libmints/mints.h>

using namespace std;
using namespace psi;

namespace psi { namespace scf {

void HF::form_B()
{
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                        FORM B
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    //Welcome to form_B, responsible for creation of the three-index tensor with embedded fitting metric
    //on core or disk.    

    //Make sure we're in the right spot
    if (print_)
        fprintf(outfile, "\n  Computing Integrals using Density Fitting\n");
    //TODO: Add support for molecular symmetry
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n"); fflush(outfile);
        abort();
    } 

    //Grab norbs and ndocc and get the ri basis up to the class scope     
    int norbs = basisset_->nbf(); 
    int ndocc = doccpi_[0];
    ribasis_ =shared_ptr<BasisSet>(new BasisSet(chkpt_, "DF_BASIS_SCF"));
    ri_nbf_ = ribasis_->nbf(); 
    
    if (print_>5) {
        basisset_->print(outfile); 
        ribasis_->print(outfile);
        fflush(outfile);
    }

    
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                        RESTART?
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (options_.get_bool("RI_SCF_RESTART")) 
    {
        //restart!! Use existing 3-index tensor on disk
        if (print_)
            fprintf(outfile,"\n  Attempting to restart existing DF-SCF computation\n"); fflush(outfile); 
        
        //First read in tensor sizes to set the bookkeeping up
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->read(PSIF_DFSCF_BJ,"N_TRI",(char *) &(ntri_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->read(PSIF_DFSCF_BJ,"N_TRI_NAIVE",(char *) &(ntri_naive_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        
        //Use ri_pair_mu_ and ri_pair_nu_ to keep track of things
        //Across schwarz sieve and unfortunate shell indexing
        ri_pair_nu_ = init_int_array(ntri_naive_);
        ri_pair_mu_ = init_int_array(ntri_naive_);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->read(PSIF_DFSCF_BJ,"RI_PAIR_MU",(char *) &(ri_pair_mu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->read(PSIF_DFSCF_BJ,"RI_PAIR_NU",(char *) &(ri_pair_nu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        
        //Now determine the storage type. It might change if you switch machines
        string storage_type;
        storage_type = options_.get_str("RI_SCF_STORAGE");

        //Detrmine storage algorithm based on user input
        //Size of the three-index tensor
        unsigned long memA = ntri_*(long)ri_nbf_;
        //Size of the fitting metric 
        unsigned long memJ = ri_nbf_*(long)ri_nbf_;
        if (storage_type == "CORE")
            df_storage_ = core;
        else if (storage_type == "DISK")
            df_storage_ = disk;
        else if (storage_type == "DEFAULT")
        {
    	    //set df_storage_ semi-heuristically based on available memory
    	    if (((long)((memA+memJ)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
                df_storage_ = core; //Full in-core
    	    else
                df_storage_ = disk; //Disk
        }	

        if (df_storage_ == core)
            fprintf(outfile,"\n  Density Fitting Algorithm proceeding on Core.\n"); 
        else if (df_storage_ == disk)
            fprintf(outfile,"\n  Density Fitting Algorithm proceeding on Disk\n"); 
        fflush(outfile);
        
        if (df_storage_ == core) {
            //We need the three-index tensor in the core
            //fprintf(outfile,"  n_tri_ %d, n_tri_naive_ %d, ri_nbf_ %d\n",ntri_,ntri_naive_,ri_nbf_); fflush(outfile);
            B_ia_P_ = block_matrix(ri_nbf_,ntri_naive_);
            next_PSIF_DFSCF_BJ = PSIO_ZERO;
            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(B_ia_P_[0][0]),sizeof(double)*ntri_naive_*ri_nbf_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        }
        psio_->close(PSIF_DFSCF_BJ,1); //we'll need to reuse this guy (probably)
        return;
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                    SCHWARZ SIEVE
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //Form the schwarz sieve
    timer_on("Schwarz Sieve");

    int sig_fun_pairs = 0;
    int sig_shell_pairs = 0;

    int *schwarz_shell_pairs;
    int *schwarz_fun_pairs;
    if (schwarz_ > 0.0) {
        
        schwarz_shell_pairs = init_int_array(basisset_->nshell()*(basisset_->nshell()+1)/2);
        schwarz_fun_pairs = init_int_array(norbs*(norbs+1)/2);
        double* max_shell_val = init_array(basisset_->nshell()*(basisset_->nshell()+1)/2);;
        double* max_fun_val = init_array(norbs*(norbs+1)/2);
        double max_global_val = 0.0;

        IntegralFactory schwarzfactory(basisset_,basisset_,basisset_,basisset_);
        shared_ptr<TwoBodyInt> eri = shared_ptr<TwoBodyInt>(schwarzfactory.eri());
        const double *buffer = eri->buffer();

        int MU, NU, mu, nu,omu,onu, nummu, numnu, index;
        int MUNU = 0;
        int munu = 0;
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU, ++MUNU) {
                numnu = basisset_->shell(NU)->nfunction();
                eri->compute_shell(MU,NU,MU,NU);
                for (mu=0; mu < nummu; ++mu) {
                    omu = basisset_->shell(MU)->function_index() + mu;
                    for (nu=0; nu < numnu; ++nu) {
                        onu = basisset_->shell(NU)->function_index() + nu;
                       
                        if (omu>=onu) {
                            index = mu*(numnu*nummu*numnu+numnu)+nu*(nummu*numnu+1);
                            //int check = mu*numnu*nummu*numnu+nu*nummu*numnu+mu*numnu+nu;
                            //fprintf(outfile,"   Index = %d, (%d %d| %d %d) = %20.15f\n",index,omu,onu,omu,onu, buffer[index] );
                            if (max_global_val<abs(buffer[index]))
                                max_global_val = abs(buffer[index]);
                            if (max_shell_val[MUNU]<abs(buffer[index]))
                                max_shell_val[MUNU] = abs(buffer[index]);
                            if (max_fun_val[omu*(omu+1)/2+onu]<abs(buffer[index]))
                                max_fun_val[omu*(omu+1)/2+onu] = abs(buffer[index]);
                        }
                    }
                }
            }
        }       
        for (int ij = 0; ij < norbs*(norbs+1)/2; ij ++)
            if (max_fun_val[ij]*max_global_val>=schwarz_*schwarz_){
                schwarz_fun_pairs[ij] = 1;
                sig_fun_pairs++;
            }
        for (int ij = 0; ij < basisset_->nshell()*(basisset_->nshell()+1)/2; ij ++)
            if (max_shell_val[ij]*max_global_val>=schwarz_*schwarz_){
                schwarz_shell_pairs[ij] = 1;
                sig_shell_pairs++;
            }
        
        //for (int i = 0, ij = 0; i<norbs; i++)
            //for (int j = 0; j<=i; j++, ij++)
                //fprintf(outfile,"   Function pair %d = (%d,%d), Max val %14.10f, Max Integral %14.10f, Significant %s\n",ij,i,j,max_fun_val[ij],max_fun_val[ij]*max_global_val,(schwarz_fun_pairs[ij])?"YES":"NO");
        //fprintf(outfile,"\n  Shell Pair Schwarz Sieve, schwarz_ = %14.10f:\n",schwarz_);
        //for (int i = 0, ij = 0; i<basisset_->nshell(); i++)
            //for (int j = 0; j<=i; j++, ij++)
                //fprintf(outfile,"   Shell pair %d = (%d,%d), Max val %14.10f, Max Integral %14.10f, Significant %s\n",ij,i,j,max_shell_val[ij],max_shell_val[ij]*max_global_val,(schwarz_shell_pairs[ij])?"YES":"NO");
        //fprintf(outfile, "\n");

        free(max_fun_val);
        free(max_shell_val);
    
        ntri_naive_ = sig_fun_pairs; //Matrix size for most of the algorithm
        ntri_ = ntri_naive_; //For now!

    } else {
        ntri_ = norbs*(norbs+1)/2; //Yeah, eat it 
        ntri_naive_ = norbs*(norbs+1)/2; 
        schwarz_shell_pairs = init_int_array(basisset_->nshell()*(basisset_->nshell()+1)/2);
        schwarz_fun_pairs = init_int_array(norbs*(norbs+1)/2);
        for (int ij = 0; ij < basisset_->nshell()*(basisset_->nshell()+1)/2; ij++)
            schwarz_shell_pairs[ij] = 1;
        for (int ij = 0; ij < ntri_; ij++)
            schwarz_fun_pairs[ij] = 1;
    }

    timer_off("Schwarz Sieve");

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                    DETERMINE STORAGE
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //Detrmine storage algorithm based on user input
    //Size of the three-index tensor
    unsigned long memA = ntri_*(long)ri_nbf_;
    //Size of the fitting metric 
    unsigned long memJ = ri_nbf_*(long)ri_nbf_;

    string storage_type;
    storage_type = options_.get_str("RI_SCF_STORAGE");

    if (storage_type == "CORE")
        df_storage_ = core;
    else if (storage_type == "DISK")
        df_storage_ = disk;
    else if (storage_type == "DEFAULT")
    {
    	//set df_storage_ semi-heuristically based on available memory
    	if (((long)((memA+memJ)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
            df_storage_ = core; //Full in-core
    	else
            df_storage_ = disk; //Disk
    }	

    if (df_storage_ == core)
        fprintf(outfile,"\n  Density Fitting Algorithm proceeding on Core.\n"); 
    else if (df_storage_ == disk)
        fprintf(outfile,"\n  Density Fitting Algorithm proceeding on Disk\n"); 
    fflush(outfile);

    
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                      FORM J^-1/2
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //OK, integrals time. start with the fitting matrix J
    //It takes a lot of work to get a null basis with Psi4! 
    shared_ptr<BasisSet> zero = BasisSet::zero_basis_set();
    
    // Create integral factory for J (Fitting Matrix in form_B)
    IntegralFactory rifactory_J(ribasis_, zero, ribasis_, zero);
    shared_ptr<TwoBodyInt> Jint = shared_ptr<TwoBodyInt>(rifactory_J.eri());

    // Integral buffer
    const double *Jbuffer = Jint->buffer();

    // J Matrix
    double **J = block_matrix(ri_nbf_, ri_nbf_);
    // J^{-1/2}
    double **J_mhalf = block_matrix(ri_nbf_, ri_nbf_);
    
    timer_on("Form J Matrix;");
    // J_{MN} = (0M|N0)
    int index = 0;

    for (int MU=0; MU < ribasis_->nshell(); ++MU) {
        int nummu = ribasis_->shell(MU)->nfunction();

        for (int NU=0; NU <= MU; ++NU) {
            int numnu = ribasis_->shell(NU)->nfunction();

            Jint->compute_shell(MU, 0, NU, 0);

            index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = ribasis_->shell(MU)->function_index() + mu;

                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = ribasis_->shell(NU)->function_index() + nu;

                    J[omu][onu] = Jbuffer[index];
                    J[onu][omu] = Jbuffer[index];
                }
            }
        }
    }
    if (print_>5) {
        fprintf(outfile,"\nJ:\n"); fflush(outfile);
        print_mat(J,ri_nbf_,ri_nbf_,outfile);
    }
    timer_off("Form J Matrix;");
    timer_on("Form J^-1/2;");

    // Form J^-1/2
    // First, diagonalize J
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    double* eigval = init_array(ri_nbf_);
    int lwork = ri_nbf_ * 3;
    double* work = init_array(lwork);
    int stat = C_DSYEV('v','u',ri_nbf_,J[0],ri_nbf_,eigval, work,lwork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }
    free(work);

    // Now J contains the eigenvectors of the original J
    // Copy J to J_copy
    double **J_copy = block_matrix(ri_nbf_, ri_nbf_);
    C_DCOPY(ri_nbf_*ri_nbf_,J[0],1,J_copy[0],1); 

    // Now form J^{-1/2} = U(T)*j^{-1/2}*U,
    // where j^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of J
    int linear_dependencies = 0;
    for (int i=0; i<ri_nbf_; i++) {
        if (eigval[i] < options_.get_double("RI_MIN_EIGENVALUE")) {
            eigval[i] = 0.0;
            linear_dependencies++;
        }
        else 
            eigval[i] = 1.0 / sqrt(eigval[i]);

        // scale one set of eigenvectors by the diagonal elements j^{-1/2}
        C_DSCAL(ri_nbf_, eigval[i], J[i], 1);
    }
    free(eigval);
    if (linear_dependencies)
        fprintf(outfile,"  WARNING: %d linear dependencies found in auxiliary basis set.\n",linear_dependencies);

    // J_mhalf = J_copy(T) * J
    C_DGEMM('t','n',ri_nbf_,ri_nbf_,ri_nbf_,1.0,
            J_copy[0],ri_nbf_,J[0],ri_nbf_,0.0,J_mhalf[0],ri_nbf_);

    free_block(J);
    free_block(J_copy);
    timer_off("Form J^-1/2;");

    if (print_>5) {
        fprintf(outfile,"\nJmhalf:\n"); fflush(outfile);
        print_mat(J_mhalf,ri_nbf_,ri_nbf_,outfile);
    }
    
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                        FORM B (FINALLY)
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //Form the AO tensor (A|mn), transfrom to (B|mn) by embedding J^-1/2 
    timer_on("Overall (B|mn)");
    
    //Use ri_pair_mu_ and ri_pair_nu_ to keep track of things
    //Across schwarz sieve and unfortunate shell indexing
    ri_pair_nu_ = init_int_array(ntri_naive_);
    ri_pair_mu_ = init_int_array(ntri_naive_);
  
    double three_index_cutoff = options_.get_double("THREE_INDEX_CUTOFF");
 
    if (df_storage_ == core)
    {	
    	//Build (A|mn) on core, and then transform in place using as large of a buffer as possible
        IntegralFactory rifactory(basisset_, basisset_, ribasis_, zero);
        shared_ptr<TwoBodyInt> eri = shared_ptr<TwoBodyInt>(rifactory.eri());
        const double *buffer = eri->buffer();
        B_ia_P_ = block_matrix(ri_nbf_,ntri_naive_); 
        
        int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
        int start_index, delta_index, l_index;
        start_index = 0;
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                    delta_index = 0;
                    for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
                        numP = ribasis_->shell(Pshell)->nfunction();
                        timer_on("(B|mn) Integrals");
                        eri->compute_shell(MU, NU, Pshell, 0);
                        timer_off("(B|mn) Integrals");
                        l_index = start_index;
                        for (mu=0 ; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                        
                                    for (P=0; P < numP; ++P) {
                                        PHI = ribasis_->shell(Pshell)->function_index() + P;
                                        B_ia_P_[PHI][l_index]= buffer[mu*numnu*numP+nu*numP+P];
                                    }
                                    if (Pshell == 0) {
                                        delta_index++;
                                        ri_pair_mu_[l_index] = omu;
                                        ri_pair_nu_[l_index] = onu;
                                    }
                                    l_index++;
                                } 
                            }
                        }
                    }
                    start_index+=delta_index;
                }
            } 
        }
        //print_mat(B_ia_P_, ri_nbf_,ntri_ ,outfile);

	// Transformation
        int max_cols = (int) (memory_/sizeof(double)-memA-memJ)/(ri_nbf_*(1+MEMORY_SAFETY_FACTOR));
        if (max_cols > ntri_naive_)
            max_cols = ntri_naive_;
        if (max_cols < 1)
            max_cols = 1; //You need to be able to spare that much at least!
        

        double **Temp1 = block_matrix(ri_nbf_,max_cols);

        //fprintf(outfile,"  Max cols %d\n",max_cols);
        //print_mat(B_ia_P_,ri_nbf_,ntri_naive_,outfile);
	
        for (int index = 0; index<ntri_naive_; index+=max_cols)
	{
            int cols = max_cols;
            if (index+cols>=ntri_naive_) 
                cols = ntri_naive_-index;

            for (int r = 0; r<ri_nbf_; r++)
                C_DCOPY(cols,&(B_ia_P_[r][index]),1,&(Temp1[r][0]),1);

            timer_on("(B|mn) Transform");
            C_DGEMM('N','N',ri_nbf_,cols,ri_nbf_,1.0, J_mhalf[0], ri_nbf_, Temp1[0], max_cols,0.0, &B_ia_P_[0][index],ntri_naive_);
            timer_off("(B|mn) Transform");

	}
	free_block(Temp1);

        //print_mat(B_ia_P_,ri_nbf_,ntri_naive_,outfile);
        timer_on("(B|mn) 3-Sieve");

        if (three_index_cutoff>0.0) {
            int left =  0;
            int right = 0;
            bool negligible_col;
            for (right = 0; right<ntri_naive_; right++) {
                negligible_col = true;
                for (int Q = 0; Q<ri_nbf_; Q++) {
                    if (fabs(B_ia_P_[Q][right])>three_index_cutoff) {
                        negligible_col = false;
                        break;
                    }
                }
                if (!negligible_col) {
                    for (int Q = 0; Q<ri_nbf_; Q++)
                        B_ia_P_[Q][left] = B_ia_P_[Q][right];
                    ri_pair_mu_[left] = ri_pair_mu_[right];
                    ri_pair_nu_[left] = ri_pair_nu_[right]; 
                    left++; 
                } else {
                    ntri_--;
                }
            }
        }
        timer_off("(B|mn) 3-Sieve");
        //print_mat(B_ia_P_, ri_nbf_,ntri_naive_ ,outfile);
        //for (int i = 0; i<ntri_naive_; i++)
        //    fprintf(outfile,"  i = %d, mu = %d, nu = %d\n",i,ri_pair_mu_[i],ri_pair_nu_[i]);

        if (options_.get_bool("RI_SCF_SAVE"))
        {
            write_B();
        }
    }
    else if (df_storage_ == disk)
    {

        //Open the BJ file and prestripe it to avoid wrong block errors
        timer_on("(B|mn) Prestriping");
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
	double *Prestripe = init_array(ntri_naive_);
	for (int Q = 0; Q < ri_nbf_; Q++) {
            psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(Prestripe[0]),sizeof(double)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
	}
        free(Prestripe);	
	next_PSIF_DFSCF_BJ = PSIO_ZERO; 
        timer_off("(B|mn) Prestriping");

        int pass = 0;

        //fprintf(outfile, "  Striped"); fflush(outfile);
        
        //Get an ERI object for the AO three-index integrals 
        IntegralFactory rifactory(basisset_, basisset_, ribasis_,zero);
        shared_ptr<TwoBodyInt> eri = shared_ptr<TwoBodyInt>(rifactory.eri());
        const double *buffer = eri->buffer();
        
        //Determine the maximum nubmer of functions and pairs in the AO basis
        int maxfun = 0;
        for (int m = 0; m<basisset_->nshell(); m++)
            if (maxfun<basisset_->shell(m)->nfunction())
                maxfun=basisset_->shell(m)->nfunction();
        int maxpairs = maxfun*maxfun;

        //Find maximum allowed memory block size (we'll need two of the same size for a good multiply)
        int max_cols = (int)((1.0-MEMORY_SAFETY_FACTOR)*memory_/sizeof(double)-memJ)/(2.0*ri_nbf_);         
        if (max_cols > ntri_naive_ + maxpairs - 1)
            max_cols = ntri_naive_ + maxpairs - 1;
        if (max_cols < maxpairs)
            max_cols = maxpairs; //Gotta give me something to work with         

        //Allocate the fitted and unfitted blocks for three-index integrals
        double **Amn = block_matrix(ri_nbf_,max_cols); //Raw integrals
        double **Bmn = block_matrix(ri_nbf_,max_cols); //Fitted integrals

        int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
        int porous_index = 0; //Index within loose tensor block (before three index sieve), block local
        int dense_index = 0; //Index within dense transformed tensor block (after three index sieve), block local
        int global_index = 0; //Dense block start index after schwarz sieve and three index sieve, global
        int start_index = 0; //Loose block start index after schwarz sieve, global
        int aux_index = 0; //Loose block start index before schwarz sieve, global
        int l_index = 0; //Compact pair index after three_index sieve, global 
        int r_index = 0; //Loose pair index before schwarz sieve, global
        int s_index = 0; //Function pair index within shell pair, shell local
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                //fprintf(outfile, "  MU = %d, NU = %d, Sig = %d\n",MU,NU,schwarz_shell_pairs[MU*(MU+1)/2+NU]); fflush(outfile);
                if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                    for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
                        numP = ribasis_->shell(Pshell)->nfunction();
                        timer_on("(B|mn) Integrals");
                        eri->compute_shell(MU, NU, Pshell, 0);
                        timer_off("(B|mn) Integrals");
                        s_index = 0;
                        for (mu=0 ; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                    for (P=0; P < numP; ++P) {
                                        PHI = ribasis_->shell(Pshell)->function_index() + P;
                                        Amn[PHI][start_index+s_index]= buffer[mu*numnu*numP+nu*numP+P];
                                    }
                                    s_index++;
                                    //if (Pshell == 0)
                                    //    fprintf(outfile,"  %MU = %d, NU = %d, mu = %d, nu = %d, omu = %d, onu = %d, l_index = %d, r_index = %d, s_index = %d\n",MU,NU,mu,nu,omu,onu, l_index,r_index,s_index);  
                                    if (Pshell == 0) {
                                        porous_index++;
                                        ri_pair_mu_[r_index] = omu;
                                        ri_pair_nu_[r_index] = onu;
                                        r_index++;
                                    }
                                } 
                            }
                        }
                    }
                    start_index = porous_index;
                    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    //Porous index is the filled length of the Amn tensor block 
                    //Do the transform, three-index sieve, and dump the integrals if either:
                    //   1) The next batch of integrals might not fit
                    //   2) The last shell pair is done
                    if (porous_index + maxpairs >= max_cols || (MU == basisset_->nshell()-1  && NU == MU)) { 
                        //Do the transform
                        pass++;
                        //print_mat(Amn,ri_nbf_,max_cols,outfile);
                        timer_on("(B|mn) Transform");
                        C_DGEMM('N','N',ri_nbf_,porous_index,ri_nbf_,1.0, J_mhalf[0], ri_nbf_, Amn[0], max_cols,0.0, Bmn[0],max_cols);
                        timer_off("(B|mn) Transform");
                        //print_mat(Bmn,ri_nbf_,max_cols,outfile);
        
                        //Use the three index sieve to compact the tensor
                        timer_on("(B|mn) 3-Sieve");
                        if (three_index_cutoff > 0.0) {
                            for (int pair = 0; pair < porous_index; pair++) {
                                bool sig = false;
                                for (int Q = 0; Q<ri_nbf_;Q++) {
                                    if (fabs(Bmn[Q][pair])>=three_index_cutoff) {
                                        sig = true;
                                        break;
                                    }
                                }
                                if (sig) { 
                                    C_DCOPY(ri_nbf_,&Bmn[0][pair],max_cols,&Amn[0][dense_index],max_cols);
                                    dense_index++;
                                    ri_pair_mu_[l_index] = ri_pair_mu_[aux_index+pair];
                                    ri_pair_nu_[l_index] = ri_pair_nu_[aux_index+pair];
                                    l_index++; 
                                } else {
                                    ntri_--;
                                }
                            }
                        } else {
                            dense_index = porous_index;
                        }
                        timer_off("(B|mn) 3-Sieve");

                        //Write the tensor out with the correct striping (mn is the fast index, but we only have blocks)
                        //NOTE: If three_index_cutoff > 0.0, the tensor is in Amn, otherwise it is in Bmn
                        //This allows for threading of the three index sieve
                        
                        //fprintf(outfile,"\n  Pass = %d. Ready for writing:\n",pass); 
                        //fprintf(outfile,"  Max cols = %d, porous_index = %d, dense_index = %d, global_offset = %d, Tensor %s\n",max_cols,porous_index,dense_index,global_index,((three_index_cutoff > 0.0)?"Amn":"Bmn"));
                        //fflush(outfile);

                        //print_mat(Amn,ri_nbf_,max_cols,outfile);
                        //fflush(outfile);
                        //print_mat(Bmn,ri_nbf_,max_cols,outfile);
                        //fflush(outfile);

                        timer_on("(B|mn) Disk Stripe");
                        for (int Q = 0; Q < ri_nbf_; Q++) {
                            next_PSIF_DFSCF_BJ = psio_get_address(PSIO_ZERO,(ULI)(Q*(ULI)ntri_naive_*sizeof(double)+global_index*sizeof(double)));
                            psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *)((three_index_cutoff > 0.0)?&Amn[Q][0]:&Bmn[Q][0]),sizeof(double)*dense_index,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
                        }
                        timer_off("(B|mn) Disk Stripe");
                        aux_index = r_index;
                        global_index += dense_index;
                        start_index = 0;
                        porous_index = 0;
                        dense_index = 0;
                    } //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                }
            }
        } 

        free_block(Amn);
        free_block(Bmn);
        if (print_>1)
            fprintf(outfile,"\n  Through (B|mn) on disk."); fflush(outfile);
        fflush(outfile);
       
        //next_PSIF_DFSCF_BJ = PSIO_ZERO;
        //double** Bhack = block_matrix(ri_nbf_,ntri_naive_);
        //psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *)&Bhack[0][0],sizeof(double)*ntri_naive_*ri_nbf_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        //print_mat(Bhack,ri_nbf_,ntri_naive_,outfile);
        //free(Bhack);

        //for (int i = 0; i<ntri_naive_; i++)
        //    fprintf(outfile,"  i = %d, mu = %d, nu = %d\n",i,ri_pair_mu_[i],ri_pair_nu_[i]);
 
        //Write the restart data, it's cheap
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->write(PSIF_DFSCF_BJ,"RI_PAIR_MU",(char *) &(ri_pair_mu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->write(PSIF_DFSCF_BJ,"RI_PAIR_NU",(char *) &(ri_pair_nu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->write(PSIF_DFSCF_BJ,"N_TRI",(char *) &(ntri_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->write(PSIF_DFSCF_BJ,"N_TRI_NAIVE",(char *) &(ntri_naive_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        psio_->close(PSIF_DFSCF_BJ,1); //We'll reuse this methinks
    }
    timer_off("Overall (B|mn)");
    
    if (print_>1) {
        if (schwarz_) {
            fprintf(outfile,"\n  Function Pair Schwarz Sieve, Cutoff = %14.10E:\n",schwarz_);
            fprintf(outfile,"  %d out of %d basis function pairs removed, %8.5f%% attenuation.\n",norbs*(norbs+1)/2-sig_fun_pairs,norbs*(norbs+1)/2,100.0*(norbs*(norbs+1)/2-sig_fun_pairs)/(1.0*norbs*(norbs+1)/2));
            int pairs = basisset_->nshell()*(basisset_->nshell()+1)/2;
            fprintf(outfile,"  %d out of %d basis shell pairs removed, %8.5f%% attenuation.\n",pairs-sig_shell_pairs,pairs,100.0*(pairs-sig_shell_pairs)/(1.0*pairs));
        }
        if (three_index_cutoff) {
            int attenuation = ntri_naive_-ntri_;
            fprintf(outfile,"  Direct Three-Index Tensor Sieve, Cutoff = %14.10E:\n",three_index_cutoff);
            fprintf(outfile,"  %d of %d (remaining) basis function pairs removed, %8.5f%% attenuation.\n",attenuation,ntri_naive_, 100.0*attenuation/(1.0*ntri_naive_));
        }
        if (schwarz_>0.0 || three_index_cutoff>0.0)
            fprintf(outfile,"  After sieving, %d out of %d basis function pairs remain, %8.5f%% attenuation.\n\n",ntri_,norbs*(norbs+1)/2,100.0*(1.0-ntri_/(1.0*norbs*(norbs+1)/2)));
    }
    fflush(outfile);
}
void HF::write_B()
{
    if (print_)
        fprintf(outfile,"  Saving three-index tensor to file 97 for restart purposes\n");

    psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
    psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
            
    psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(B_ia_P_[0][0]),sizeof(double)*ntri_naive_*ri_nbf_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
                
    next_PSIF_DFSCF_BJ = PSIO_ZERO;
    psio_->write(PSIF_DFSCF_BJ,"RI_PAIR_MU",(char *) &(ri_pair_mu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    next_PSIF_DFSCF_BJ = PSIO_ZERO;
    psio_->write(PSIF_DFSCF_BJ,"RI_PAIR_NU",(char *) &(ri_pair_nu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    next_PSIF_DFSCF_BJ = PSIO_ZERO;
    psio_->write(PSIF_DFSCF_BJ,"N_TRI",(char *) &(ntri_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    next_PSIF_DFSCF_BJ = PSIO_ZERO;
    psio_->write(PSIF_DFSCF_BJ,"N_TRI_NAIVE",(char *) &(ntri_naive_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    psio_->close(PSIF_DFSCF_BJ,1);
}
void HF::free_B()
{
    if (df_storage_ == core)
        free_block(B_ia_P_);
    free(ri_pair_mu_);
    free(ri_pair_nu_);
}
void RHF::form_G_from_RI()
{
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                        FORM G
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    timer_on("Overall G");
    //Get norbs
    int norbs = basisset_->nbf();     
    //Zero the J matrix
    if (J_is_required_)
        J_->zero();
    //Zero the K matrix
    if (K_is_required_)
        K_->zero();

    //D_->print(outfile);
    //C_->print(outfile);
    
    //Rearrange the D matrix as a vector in terms of ri_pair indices
    //Off diagonal elements get 2x weight due to permutational symmetry
    double* DD = init_array(ntri_);
    
    for (int ij = 0; ij<ntri_; ij++) {
        DD[ij] = D_->get(0,ri_pair_mu_[ij],ri_pair_nu_[ij]); 
        if (ri_pair_mu_[ij] != ri_pair_nu_[ij])
            DD[ij] *= 2.0;
            //only A irrep at the moment!!
    }
    //Get the C matrix (exchange messes things up)
    int ndocc = doccpi_[0];
    double** Cocc = block_matrix(ndocc,norbs);
    for (int i=0; i<norbs; i++) {
        for (int j=0; j<ndocc; j++)
            Cocc[j][i] = C_->get(0,i,j);
        //only A irrep at the moment!!
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                        CORE ALGORITHM
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (df_storage_ == core) {
        //B is in core, DGEMM everything
        if (J_is_required_) {
            /* COULOMB PART */
            //Coulomb convolution vector
            double *L = init_array(ri_nbf_);
            //Temporary J matrix
            double *J = init_array(ntri_);
            //DGEMV -> L:
            //L_Q = (Q|ls)*D_{ls}
            timer_on("J DDOT");
            C_DGEMV('N',ri_nbf_,ntri_,1.0,B_ia_P_[0],ntri_naive_,DD,1,0.0,L,1);
            timer_off("J DDOT");
            //DGEMV -> J:
            //J_{mn} = L_Q(Q|mn)
            timer_on("J DAXPY");
            C_DGEMV('T',ri_nbf_,ntri_,1.0,B_ia_P_[0],ntri_naive_,L,1,0.0,J,1);
            timer_off("J DAXPY");
            //Put everything in J_
            for (int ij = 0; ij < ntri_; ij++) {
                J_->set(0,ri_pair_mu_[ij],ri_pair_nu_[ij],J[ij]);
                J_->set(0,ri_pair_nu_[ij],ri_pair_mu_[ij],J[ij]);
            }

            //J_->print(outfile);

            free(L);
            free(J);
            //J_->print(outfile);
        }
        if (K_is_required_) {
            timer_on("Form E");
            /* EXCHANGE PART */
            

            int max_rows = (memory_/sizeof(double)-ri_nbf_*ntri_naive_)/(norbs*ndocc*(1.0+MEMORY_SAFETY_FACTOR));
            if (max_rows > ri_nbf_)
                max_rows = ri_nbf_;
            if (max_rows < 1)
                max_rows = 1;

            //fprintf(outfile,"  max_rows = %d\n",max_rows);

            //E exchange matrix
            double** E = block_matrix(norbs, ndocc*max_rows);
            //Temporary K matrix
            double** K = block_matrix(norbs,norbs);
            //QS temp matrix for DGEMM
            double** QS = block_matrix(max_rows,norbs);
            //Temp matrix for DGEMM
            double** Temp = block_matrix(ndocc,max_rows);
            // Temp matrix for sparse DGEMM if sieve exists
            double** Ctemp = block_matrix(ndocc,norbs);
            // Index array for non-canonical ordering of mn
            int** m_ij_indices = init_int_matrix(norbs,norbs);
            // Index array of n for given m (in order of above)
            int** n_indices = init_int_matrix(norbs,norbs);
            // sizes of above for schwarz sieve
            int* index_sizes = init_int_array(norbs);

            timer_on("Initial E Indexing");
            
            for (int ij = 0; ij<ntri_; ij++) {
                int m = ri_pair_mu_[ij];
                int n = ri_pair_nu_[ij];
                m_ij_indices[m][index_sizes[m]] = ij;
                n_indices[m][index_sizes[m]] = n;
                index_sizes[m]++;
                if (m != n){
                    m_ij_indices[n][index_sizes[n]] = ij;
                    n_indices[n][index_sizes[n]] = m;
                    index_sizes[n]++;
                }
            }

            timer_off("Initial E Indexing");
            
            //CUTOFF STUFF (KINDA LOOSE)
            //What is the smallest overlap element worth computing exchange for
            double cutoff = options_.get_double("OVERLAP_CUTOFF");
            //Contribution to the exchange matrix
            double contribution;
            //Partial contribution K matrix
            double** Ktemp;
            int att_elements = norbs*(norbs+1)/2;
            
            //MASTER LOOP
            int current_rows;
            for (int row = 0; row<ri_nbf_; row += max_rows) {
                current_rows = max_rows;
                if (row+current_rows>=ri_nbf_)
                    current_rows = ri_nbf_-row; 

                timer_on("E SORT");
                for (int m = 0; m<norbs; m++) {
                    //Find out where the m's are!
                    //timer_on("Intermediate Indexing");
                    int n, ij;
                    for (int index = 0; index<index_sizes[m]; index++) {
                        ij = m_ij_indices[m][index];
                        n = n_indices[m][index];
                 
                        //fprintf(outfile,"  index = %d ij = %d \n",index,ij); fflush(outfile);//fprintf(outfile,"(ij, mu) = (%d, %d)\n",ij,mu);
                        C_DCOPY(current_rows,&B_ia_P_[row][ij],ntri_naive_,&QS[0][index],norbs);
                        //for (int Q = 0; Q<ri_nbf_; Q++) {
                        //     QS[Q][index] = B_ia_P_[Q][ij];
                             //fprintf(outfile," ij = %d, mu = %d, Q = %d, val = %14.10f\n",ij,mu,Q,QS[Q][mu]);
                        //}
                        C_DCOPY(ndocc,&Cocc[0][n],norbs,&Ctemp[0][index],norbs);
                        //for (int o = 0; o<ndocc; o++) {
                        //    Ctemp[o][index] = Cocc[o][n];
                        //}
                    }
                    //timer_off("Intermediate Indexing");
                    //timer_on("DGEMM 2");
                
                    C_DGEMM('N','T',ndocc,current_rows,index_sizes[m],1.0,Ctemp[0],norbs,QS[0],norbs, 0.0, Temp[0], max_rows);
                    //timer_off("DGEMM 2");
                    //timer_on("Final Indexing");
                
                    int offset;
                    for (int Q = 0; Q<current_rows; Q++) {
                        offset = Q*ndocc;
                        for (int i = 0; i<ndocc; i++) {
                            E[m][i+offset] = Temp[i][Q];
                        }
                    }
                    //timer_off("Final Indexing");
                }    
                timer_off("E SORT");
                timer_on("E DGEMM");

                //K_{mn} = E_{im}^QE_{in}^Q

                if (cutoff > 0.0) {
                    for (int i = 0; i< norbs; i++)
                        for (int j = 0; j <= i; j++) {
                            //Is S big enough?
                            if (abs(S_->get(0,i,j))>=cutoff) {
                                contribution = C_DDOT(current_rows*ndocc,&E[i][0],1,&E[j][0],1);
                                K[i][j] += contribution;
                                if (i != j)
                                    K[j][i] += contribution;
                                if (row == 0)
                                    att_elements--;
                            }
                        }
                } else {
                    //DGEMM, usually faster regardless
                    //There it is
                    if (row == 0)
                        Ktemp = block_matrix(norbs,norbs);
                    //If the K matrix is sparse due to low overlap or density
                    //Form the matrix elements individually
                    C_DGEMM('N','T',norbs,norbs,current_rows*ndocc,1.0,E[0],max_rows*ndocc,E[0],max_rows*ndocc, 0.0, Ktemp[0], norbs);

                    C_DAXPY(norbs*norbs,1.0,&Ktemp[0][0],1,&K[0][0],1);
                }
                timer_off("E DGEMM");
            } 

            //print_mat(K,norbs,norbs,outfile); fflush(outfile);

            for (int i = 0; i < norbs; i++)
                for (int j = 0; j<=i; j++) {
                    K_->set(0,j,i,K[i][j]);
                    K_->set(0,i,j,K[i][j]);
                }

            if (cutoff > 0.0 && print_>2) {
                fprintf(outfile,"  K matrix was built with with an overlap cutoff of %14.10E\n",cutoff);  
                fprintf(outfile,"  %d out of %d elements were attenuated due to overlap, %8.5f%% attenuation.\n",att_elements,norbs*(norbs+1)/2,100.0*att_elements/(1.0*norbs*(norbs+1)/2)); 
            } else if (cutoff == 0.0) {
                free_block(Ktemp);
            } 
            
            //Frees
            free_block(E);
            free_block(K);
            free_block(Ctemp);
            free(m_ij_indices[0]);
            free(m_ij_indices);
            free(n_indices[0]);
            free(n_indices);
            free(index_sizes);
            free_block(Temp);
            free_block(QS);
            
            timer_off("Form E");
            
            //fprintf(outfile,"\n E: \n");
            //print_mat(E,norbs,ndocc*ri_nbf_,outfile);
            
            //K_->print(outfile);
            fflush(outfile);
        }
    } else {
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //
        //                      DISK ALGORITHM
        //
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        //B is on disk, E will be in disk, Single disk pass
        //B_ia_P_ stores multiple aux basis function rows of 
        //the three index tensor

        //How many rows to read?
        int max_rows = floor(((memory_/sizeof(double)))/((1.0+MEMORY_SAFETY_FACTOR)*(ntri_naive_+ndocc*norbs)));
	if (max_rows > ri_nbf_)
            max_rows = ri_nbf_;
        if (max_rows < 1)
            max_rows = 1; //Without a row, I can't work

        //fprintf(outfile,"  max_rows = %d\n",max_rows); fflush(outfile);

        //FILE STUFF
        //Transformed integrals are stored in PSIF_DFSCF_BJ
        //Which had better exist at this point
        //timer_on("Open B");
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        //timer_off("Open B");
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        //Row height per read, so that we can tune this value
        int rows_per_read = options_.get_int("ROWS_PER_READ"); 
        
        
        //Chunk of three-index tensor
        B_ia_P_ = block_matrix(max_rows,ntri_naive_);
        
        //COULOMB STUFF
        //Coulomb convolution vector, done element by element
        double *L; //Density coefficients
        double *J; //J register
        double *Jtemp; //J contribution
        if (J_is_required_){
            L = init_array(max_rows);
            Jtemp = init_array(ntri_);
            J = init_array(ntri_);
        }

        //EXCHANGE STUFF
        double **E;
        double **K;
        double **QS;
        double **Temp; 
        double **Ctemp;
        int** m_ij_indices;
        int** n_indices;
        int* index_sizes;
        if (K_is_required_) {    
            //E exchange matrix
            E = block_matrix(norbs, ndocc*max_rows);
            //Temporary K matrix
            K = block_matrix(norbs,norbs);
            //QS temp matrix for DGEMM
            QS = block_matrix(max_rows,norbs);
            //Temp matrix for DGEMM
            Temp = block_matrix(ndocc,max_rows);
            // Temp matrix for sparse DGEMM if sieve exists
            Ctemp = block_matrix(ndocc,norbs);
            // Index array for non-canonical ordering of mn
            m_ij_indices = init_int_matrix(norbs,norbs);
            // Index array of n for given m (in order of above)
            n_indices = init_int_matrix(norbs,norbs);
            // sizes of above for schwarz sieve
            index_sizes = init_int_array(norbs);

            timer_on("Initial E Indexing");
            
            for (int ij = 0; ij<ntri_; ij++) {
                int m = ri_pair_mu_[ij];
                int n = ri_pair_nu_[ij];
                m_ij_indices[m][index_sizes[m]] = ij;
                n_indices[m][index_sizes[m]] = n;
                index_sizes[m]++;
                if (m != n){
                    m_ij_indices[n][index_sizes[n]] = ij;
                    n_indices[n][index_sizes[n]] = m;
                    index_sizes[n]++;
                }
            }

            timer_off("Initial E Indexing");
        }

        //CUTOFF STUFF (KINDA LOOSE)
        //What is the smallest overlap element worth computing exchange for
        double cutoff = options_.get_double("OVERLAP_CUTOFF");
        //Contribution to the exchange matrix
        double contribution;
        //Partial contribution K matrix
        double** Ktemp;
        int att_elements = norbs*(norbs+1)/2;
        
        //MASTER LOOP
        int current_rows,offset;
        //int pass = 0;
        for (int row = 0; row <ri_nbf_; row+=max_rows)
        {

            //Setup block size
            //Careful Starfox
	    current_rows = max_rows;
	    if (row+max_rows>=ri_nbf_)
		current_rows = ri_nbf_-row;
            //Read max_rows of the (B|mn) tensor in, place in B_ia_P
            //fprintf(outfile,"\n  Pass %d, current_rows = %d\n",pass,current_rows);
            //pass++;
            timer_on("Read B");
            
            //New method read in a few rows at a time
            int block_height = rows_per_read;
            for (int Q = 0; Q<current_rows; Q+=rows_per_read) {
                if (Q+rows_per_read>=current_rows)
                    block_height = current_rows-Q;
                psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(B_ia_P_[Q][0]),sizeof(double)*ntri_naive_*block_height,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
                
            }

            //double **Bhack = block_matrix(current_rows,ntri_);
            //C_DCOPY(current_rows*ntri_,in_buffer,1,&Bhack[0][0],1);
            //fprintf(outfile,"  B:\n");
            //print_mat(Bhack,current_rows,ntri_,outfile);
            
            timer_off("Read B");

            //Do J here with blocked DGEMV
            if (J_is_required_) {
                //DGEMV -> L:
                //L_Q = (Q|ls)*D_{ls}
                timer_on("J DDOT");
                C_DGEMV('N',current_rows,ntri_,1.0,B_ia_P_[0],ntri_naive_,DD,1,0.0,L,1);
                timer_off("J DDOT");
                //DGEMV -> J:
                //J_{mn} = L_Q(Q|mn)
                timer_on("J DAXPY");
                C_DGEMV('T',current_rows,ntri_,1.0,B_ia_P_[0],ntri_naive_,L,1,0.0,Jtemp,1);
                timer_off("J DAXPY");

                C_DAXPY(ntri_,1.0,Jtemp,1,J,1);
                
            }

            //Do K here just like the core algorithm 
            if (K_is_required_) {
                timer_on("E SORT");
                for (int m = 0; m<norbs; m++) {
                    //Find out where the m's are!
                    //timer_on("Intermediate Indexing");
                    int n, ij;
                    for (int index = 0; index<index_sizes[m]; index++) {
                        ij = m_ij_indices[m][index];
                        n = n_indices[m][index];
                
 
                        //fprintf(outfile,"  index = %d ij = %d \n",index,ij); fflush(outfile);//fprintf(outfile,"(ij, mu) = (%d, %d)\n",ij,mu);
                        C_DCOPY(current_rows,&B_ia_P_[0][ij],ntri_naive_,&QS[0][index],norbs);
                        //for (int Q = 0; Q<ri_nbf_; Q++) {
                        //     QS[Q][index] = B_ia_P_[Q][ij];
                             //fprintf(outfile," ij = %d, mu = %d, Q = %d, val = %14.10f\n",ij,mu,Q,QS[Q][mu]);
                        //}
                        C_DCOPY(ndocc,&Cocc[0][n],norbs,&Ctemp[0][index],norbs);
                        //for (int o = 0; o<ndocc; o++) {
                        //    Ctemp[o][index] = Cocc[o][n];
                        //}
                    }
                    //timer_off("Intermediate Indexing");
                    //timer_on("DGEMM 2");
                
                    C_DGEMM('N','T',ndocc,current_rows,index_sizes[m],1.0,Ctemp[0],norbs,QS[0],norbs, 0.0, Temp[0], max_rows);
                    //timer_off("DGEMM 2");
                    //timer_on("Final Indexing");
                
                    int offset;
                    for (int Q = 0; Q<current_rows; Q++) {
                        offset = Q*ndocc;
                        for (int i = 0; i<ndocc; i++) {
                            E[m][i+offset] = Temp[i][Q];
                        }
                    }
                    //timer_off("Final Indexing");
                }    
                timer_off("E SORT");
                timer_on("E DGEMM");

                //K_{mn} = E_{im}^QE_{in}^Q

                if (cutoff > 0.0) {
                    for (int i = 0; i< norbs; i++)
                        for (int j = 0; j <= i; j++) {
                            //Is S big enough?
                            if (abs(S_->get(0,i,j))>=cutoff) {
                                contribution = C_DDOT(current_rows*ndocc,&E[i][0],1,&E[j][0],1);
                                K[i][j] += contribution;
                                if (i != j)
                                    K[j][i] += contribution;
                                if (row == 0)
                                    att_elements--;
                            }
                        }
                } else {
                    //DGEMM, usually faster regardless
                    //There it is
                    if (row == 0)
                        Ktemp = block_matrix(norbs,norbs);
                    //If the K matrix is sparse due to low overlap or density
                    //Form the matrix elements individually
                    C_DGEMM('N','T',norbs,norbs,current_rows*ndocc,1.0,E[0],max_rows*ndocc,E[0],max_rows*ndocc, 0.0, Ktemp[0], norbs);

                    C_DAXPY(norbs*norbs,1.0,&Ktemp[0][0],1,&K[0][0],1);
                }
                timer_off("E DGEMM");
            }
        }
        psio_->close(PSIF_DFSCF_BJ,1);
        free_block(B_ia_P_);
        /* Form J */
        if (J_is_required_) {
            for (int ij2 = 0; ij2 < ntri_; ij2++) {
                J_->set(0,ri_pair_mu_[ij2],ri_pair_nu_[ij2],J[ij2]);
                J_->set(0,ri_pair_nu_[ij2],ri_pair_mu_[ij2],J[ij2]);
            }
            free(J);
            free(Jtemp);
            free(L);
        }
        /* Form K */
        if (K_is_required_) {
            for (int i = 0; i < norbs; i++)
                for (int j = 0; j<=i; j++) {
                    K_->set(0,j,i,K[i][j]);
                    K_->set(0,i,j,K[i][j]);
            }
            if (cutoff > 0.0 && print_>2) {
                fprintf(outfile,"  K matrix was built with with an overlap cutoff of %14.10E\n",cutoff);  
                fprintf(outfile,"  %d out of %d elements were attenuated due to overlap, %8.5f%% attenuation.\n",att_elements,norbs*(norbs+1)/2,100.0*att_elements/(1.0*norbs*(norbs+1)/2)); 
            } else if (cutoff == 0.0) {
                free_block(Ktemp);
            } 
            
            //Frees
            free_block(E);
            free_block(K);
            free_block(Ctemp);
            free(m_ij_indices[0]);
            free(m_ij_indices);
            free(n_indices[0]);
            free(n_indices);
            free(index_sizes);
            free_block(Temp);
            free_block(QS);
        } 
    }

    free(DD);
    free_block(Cocc);
    
    //J_->print(outfile);
    //K_->print(outfile);

    /* FORM G_ */
    //This method takes one extra O(N^2) scale,
    //but preserves J and K in place
    G_->zero();
    G_->add(J_);
    G_->scale(-2.0);
    if (K_is_required_)
        G_->add(K_);
    G_->scale(-1.0);
    //G_->print(outfile);
    timer_off("Overall G");

}
void RHF::form_J_from_RI()
{
    fprintf(stderr, "RKS RI  Not implemented yet!\n");
    abort();
    /**
    J_->zero();
    int norbs = basisset_->nbf();
    double** D = D_->to_block_matrix();
    //print_mat(D, norbs, norbs,outfile);

    double **J = block_matrix(norbs, norbs);

    double* D2 = init_array(norbs*(norbs+1)/2);
    for (int i = 0, ij = 0; i<norbs; i++) {
        for (int j = 0; j<=i; ij++, j++)
        {
            D2[ij] = (i==j?1.0:2.0)*D[i][j];
        }
    }

    double *L = init_array(ri_nbf_);
    double *Gtemp = init_array(norbs*(norbs+1)/2);

    //B_ia_P_ in core
    if (df_storage_ == core||df_storage_ == double_core)
    {
        for (int i=0; i<ri_nbf_; i++) {
            L[i]=C_DDOT(norbs*(norbs+1)/2,D2,1,B_ia_P_[i],1);
        }

        C_DGEMM('T','N',1,norbs*(norbs+1)/2,ri_nbf_,1.0,L,1,B_ia_P_[0],norbs*(norbs+1)/2, 0.0, Gtemp, norbs*(norbs+1)/2);    
        free(D2);
    } 
    //B_ia_P_ on disk
    else 
    {
        double *DD = init_array(norbs*(norbs+1)/2);
        for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
            DD[ij] = D2[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]];
        free(D2);

        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        double *in_buffer = init_array(norbs*(norbs+1)/2);
        for (int i=0; i<ri_nbf_; i++) {
            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            L[i]=C_DDOT(norbs*(norbs+1)/2,DD,1,in_buffer,1);
        }

        free(DD);
        psio_->close(PSIF_DFSCF_BJ,1);

        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        double *G2 = init_array(norbs*(norbs+1)/2);
        register double LL;
        for (int Q = 0; Q<ri_nbf_; Q++)
        {
            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            LL = L[Q];
            for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
                G2[ij]+=LL*in_buffer[ij];
        }
        free(in_buffer);
        psio_->close(PSIF_DFSCF_BJ,1);

        for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
            Gtemp[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]] = G2[ij];
        free(G2);
    }

    free(L);

    for (int i = 0, ij=0; i<norbs; i++) {
        for (int j = 0; j<=i; ij++,j++)    
        {
            J[i][j] = Gtemp[ij];
            J[j][i] = Gtemp[ij];
        }
    }
    //fprintf(outfile, "\nJ:\n");
    //print_mat(J,norbs,norbs,outfile); fflush(outfile);
    free(Gtemp);
    free_block(D);
    
    for (int i=0; i<norbs; i++) {
        for (int j=0; j<=i; j++) {
            J_->set(0,i,j,J[i][j]);
            if (i!= j)
                J_->set(0,j,i,J[i][j]);
        }
    }
    //fprintf(outfile,"\n");
    //J_->print();
    free_block(J);
    **/
}
void RHF::form_K_from_RI()
{
    fprintf(stderr, "RKS RI  Not implemented yet!\n");
    abort();
    /**
    K_->zero();
    int norbs = basisset_->nbf();
    int ndocc = doccpi_[0];
    double **K = block_matrix(norbs, norbs);
    double** Cocc = block_matrix(ndocc,norbs);
    for (int i=0; i<norbs; i++) {
        for (int j=0; j<ndocc; j++)
            Cocc[j][i] = C_->get(0,i,j);
    }
    //fprintf(outfile,"\nC:\n");
    //print_mat(Cocc,ndocc,norbs,outfile);
    double** B_im_Q;
    //B_ia_P in core, B_im_Q in core
    if (df_storage_ == core||df_storage_ == double_core)
    {
        B_im_Q = block_matrix(norbs, ndocc*ri_nbf_);
        double** QS = block_matrix(ri_nbf_,norbs);
        double** Temp = block_matrix(ndocc,ri_nbf_);

        //print_mat(B_ia_P_,ri_nbf_,norbs*(norbs+1)/2,outfile);
        //fprintf(outfile,"\nYo\n");
        //print_mat(Cocc,ndocc,norbs,outfile);
        for (int m = 0; m<norbs; m++) {
            for (int Q = 0; Q<ri_nbf_; Q++) {
                for (int s = 0; s<norbs; s++) {
                    QS[Q][s] = B_ia_P_[Q][((s>=m)?ioff[s]+m:ioff[m]+s)];
                }
            }
            C_DGEMM('N','T',ndocc,ri_nbf_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], ri_nbf_);
            //print_mat(Temp,ndocc,ri_nbf_,outfile);
            for (int Q = 0; Q<ri_nbf_; Q++) {
                for (int i = 0; i<ndocc; i++) {
                    B_im_Q[m][i+Q*ndocc] = Temp[i][Q];
                }
            }
        }
        //print_mat(B_im_Q,norbs,ndocc*ri_nbf_,outfile);
        free_block(QS);
        free_block(Temp); 
    }
    //B_ia_P in disk, B_im_Q in core
    else if (df_storage_ == flip_B_disk || df_storage_ == k_incore)
    {
        double *in_buffer = init_array(norbs*(norbs+1)/2);
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        B_im_Q = block_matrix(norbs, ndocc*ri_nbf_);

        int mu, nu;

        for (int Q = 0; Q<ri_nbf_; Q++)
        {
            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
            {
                mu = ri_pair_mu_[ij];
                nu = ri_pair_nu_[ij];
                for (int i = 0; i<ndocc; i++)
                {
                    B_im_Q[mu][i+Q*ndocc]+=Cocc[i][nu]*in_buffer[ij];
                    if (mu != nu)
                        B_im_Q[nu][i+Q*ndocc]+=Cocc[i][mu]*in_buffer[ij];
                }
            }
        }

        psio_->close(PSIF_DFSCF_BJ,1);
        free(in_buffer);
        //B_ia_P in disk, B_im_Q in disk
    } else {
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->open(PSIF_DFSCF_K,PSIO_OPEN_NEW);
        psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

        double *in_buffer = init_array(norbs*(norbs+1)/2);
        double *out_buffer = init_array(norbs*ndocc);

        int mu, nu;
        for (int Q = 0; Q<ri_nbf_; Q++)
        {
            for (int im = 0; im<ndocc*norbs; im++)
                out_buffer[im] = 0.0;

            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
            {
                mu = ri_pair_mu_[ij];
                nu = ri_pair_nu_[ij];
                for (int i = 0; i<ndocc; i++)
                {
                    out_buffer[mu*ndocc+i]+=Cocc[i][nu]*in_buffer[ij];
                    if (mu != nu)
                        out_buffer[nu*ndocc+i]+=Cocc[i][mu]*in_buffer[ij];
                }
            }
            psio_->write(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(out_buffer[0]),sizeof(double)*norbs*ndocc,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
        }

        free(in_buffer);
        free(out_buffer);

        psio_->close(PSIF_DFSCF_BJ,1);
        psio_->close(PSIF_DFSCF_K,1);
    }

    free_block(Cocc);

    if (df_storage_ == core ||df_storage_ == double_core|| df_storage_ == flip_B_disk || df_storage_ == k_incore)
    {
        C_DGEMM('N','T',norbs,norbs,ri_nbf_*ndocc,1.0,B_im_Q[0],ri_nbf_*ndocc,B_im_Q[0],ri_nbf_*ndocc, 0.0, K[0], norbs);
        free_block(B_im_Q);
    } else {
        psio_->open(PSIF_DFSCF_K,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

        double *in_buffer = init_array(norbs*ndocc);

        for (int Q = 0; Q<ri_nbf_; Q++)
        {
            psio_->read(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(in_buffer[0]),sizeof(double)*norbs*ndocc,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);

            for (int m = 0; m<norbs; m++)
                for (int n = 0; n<=m; n++)
                    for (int i = 0; i<ndocc; i++)
                    {
                K[m][n]+=in_buffer[m*ndocc+i]*in_buffer[n*ndocc+i];
                K[n][m] = K[m][n];
            }
        }

        free(in_buffer);
        psio_->close(PSIF_DFSCF_K,1);
    }

    //fprintf(outfile, "\nK:\n");
    //print_mat(K,norbs,norbs,outfile);

    for (int i=0; i<norbs; i++) {
        for (int j=0; j<=i; j++) {
            K_->set(0,i,j,K[i][j]);
            if (i!= j)
                K_->set(0,j,i,K[i][j]);
        }
    }
    //fprintf(outfile,"\n");
    //K_->print();
    free_block(K);**/
}
void UHF::form_G_from_RI()
{
    fprintf(stderr, "UHF RI Not implemented yet!\n");
    abort();
    /**
    * This guy needs an update

    int norbs = basisset_->nbf();
    int nalpha = nalphapi_[0];
    //int nbeta = nbetapi_[0];

    double** Da = Da_->to_block_matrix();
    double** Db = Db_->to_block_matrix();

    Ga_->zero();	
    Gb_->zero();
    
    //fprintf(outfile,"  Arrival in form G from RI \n");fflush(outfile);

    double **J = block_matrix(norbs, norbs);
    
    double* D2 = init_array(norbs*(norbs+1)/2);
    for (int i = 0, ij = 0; i<norbs; i++) {
        for (int j = 0; j<=i; ij++, j++)
        {
            D2[ij] = (i==j?1.0:2.0)*Da[i][j];
            D2[ij] += (i==j?1.0:2.0)*Db[i][j];
        }
    }

    double *L = init_array(ri_nbf_);
    double *Gtemp = init_array(norbs*(norbs+1)/2);
    
    //B_ia_P_ in core
    if (df_storage_ == core||df_storage_ == double_core)
    {                

        for (int i=0; i<ri_nbf_; i++) {
            L[i]=C_DDOT(norbs*(norbs+1)/2,D2,1,B_ia_P_[i],1);
        }


        C_DGEMM('T','N',1,norbs*(norbs+1)/2,ri_nbf_,1.0,L,1,B_ia_P_[0],norbs*(norbs+1)/2, 0.0, Gtemp, norbs*(norbs+1)/2);
        free(D2);
    } 
    //B_ia_P_ on disk
    else 
    {
    	double *DD = init_array(norbs*(norbs+1)/2);
    	for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
        {
            DD[ij] = D2[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]];
            //DD[ij] += D2[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]];
        }
    	free(D2);
    	
    	psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    	psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
    	double *in_buffer = init_array(norbs*(norbs+1)/2);
    	for (int i=0; i<ri_nbf_; i++) {
            psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            L[i]=C_DDOT(norbs*(norbs+1)/2,DD,1,in_buffer,1);
      	}

    	free(DD);
    	psio_close(PSIF_DFSCF_BJ,1);
    	
    	psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    	next_PSIF_DFSCF_BJ = PSIO_ZERO;
    	double *G2 = init_array(norbs*(norbs+1)/2);
    	register double LL;
    	for (int Q = 0; Q<ri_nbf_; Q++)
    	{
            psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            LL = L[Q];
            for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
            {
                G2[ij]+=LL*in_buffer[ij];
            }
    	}
    	free(in_buffer);
    	psio_close(PSIF_DFSCF_BJ,1);
    	
    	for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
        {
            Gtemp[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]] = G2[ij];
        }
    	free(G2);
    }
    
    free(L);
    
    for (int i = 0, ij=0; i<norbs; i++) {
        for (int j = 0; j<=i; ij++,j++)    
        {
            J[i][j] = Gtemp[ij];
            J[j][i] = Gtemp[ij];
        }
    }
    //fprintf(outfile, "\nJ:\n");
    //print_mat(J,norbs,norbs,outfile); fflush(outfile);
    free(Gtemp);
    free_block(Da);
    free_block(Db);
    
    double** B_im_Q;
    double** Cocc = block_matrix(nalpha,norbs);
    
    //First Pass: Ca_
    for (int i=0; i<norbs; i++) {
        for (int j=0; j<nalpha; j++)
            Cocc[j][i] = Ca_->get(0,i,j);
    }
    //B_ia_P in core, B_im_Q in core
    if (df_storage_ == core||df_storage_ == double_core)
    {
    	B_im_Q = block_matrix(norbs, nalpha*ri_nbf_);
    	double** QS = block_matrix(ri_nbf_,norbs);
    	double** Temp = block_matrix(nalpha,ri_nbf_);

        //print_mat(B_ia_P_,ri_nbf_,norbs*(norbs+1)/2,outfile);
        //fprintf(outfile,"\nYo\n");
        //print_mat(Cocc,nalpha,norbs,outfile);
        for (int m = 0; m<norbs; m++) {
            for (int Q = 0; Q<ri_nbf_; Q++) {
                for (int s = 0; s<norbs; s++) {
                    QS[Q][s] = B_ia_P_[Q][((s>=m)?ioff[s]+m:ioff[m]+s)];
                }
            }
            //C_DGEMM('N','T',ndocc,ri_nbf_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], ri_nbf_);
            C_DGEMM('N','T',nalpha,ri_nbf_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], ri_nbf_);
            for (int Q = 0; Q<ri_nbf_; Q++) {
                for (int i = 0; i<nalpha; i++) {
                    B_im_Q[m][i+Q*nalpha] = Temp[i][Q];
                }
            }
    	}
    	free_block(QS);
    	free_block(Temp);
    	
    }
    
    //B_ia_P in disk, B_im_Q in core
    else if (df_storage_ == flip_B_disk || df_storage_ == k_incore)
    {
        double *in_buffer = init_array(norbs*(norbs+1)/2);
    	psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    	psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
    	B_im_Q = block_matrix(norbs, nalpha*ri_nbf_);
    	
    	int mu, nu;
    	
    	for (int Q = 0; Q<ri_nbf_; Q++)
    	{
            psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
            {
                mu = ri_pair_mu_[ij];
                nu = ri_pair_nu_[ij];
                for (int i = 0; i<nalpha; i++)
                {
                    B_im_Q[mu][i+Q*nalpha]+=Cocc[i][nu]*in_buffer[ij];
                    if (mu != nu)
                        B_im_Q[nu][i+Q*nalpha]+=Cocc[i][mu]*in_buffer[ij];
                }
            }
    	}

    	psio_close(PSIF_DFSCF_BJ,1);
    	free(in_buffer);
        //B_ia_P in disk, B_im_Q in disk
    } else {
    	psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    	psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
    	psio_open(PSIF_DFSCF_K,PSIO_OPEN_NEW);
    	psio_address next_PSIF_DFSCF_K = PSIO_ZERO;
    	
    	double *in_buffer = init_array(norbs*(norbs+1)/2);
    	double *out_buffer = init_array(norbs*nalpha);
    	
    	int mu, nu;
    	for (int Q = 0; Q<ri_nbf_; Q++)
    	{
            for (int im = 0; im<nalpha*norbs; im++)
                out_buffer[im] = 0.0;

            psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
            {
                mu = ri_pair_mu_[ij];
                nu = ri_pair_nu_[ij];
                for (int i = 0; i<nalpha; i++)
                {
                    out_buffer[mu*nalpha+i]+=Cocc[i][nu]*in_buffer[ij];
                    if (mu != nu)
                        out_buffer[nu*nalpha+i]+=Cocc[i][mu]*in_buffer[ij];
                }
            }
            psio_write(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(out_buffer[0]),sizeof(double)*norbs*nalpha,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
    	}
    	
    	free(in_buffer);
    	free(out_buffer);
    	
    	psio_close(PSIF_DFSCF_BJ,1);
    	psio_close(PSIF_DFSCF_K,1);
    }
    free_block(Cocc);
    //fprintf(outfile,"  3-Index A Formed \n"); fflush(outfile);
    double **Ka = block_matrix(norbs, norbs);

    if (df_storage_ == core ||df_storage_ == double_core|| df_storage_ == flip_B_disk || df_storage_ == k_incore)
    {
    	C_DGEMM('N','T',norbs,norbs,ri_nbf_*nalpha,1.0,B_im_Q[0],ri_nbf_*nalpha,B_im_Q[0],ri_nbf_*nalpha, 0.0, Ka[0], norbs);
    	free_block(B_im_Q);
    } else {
    	psio_open(PSIF_DFSCF_K,PSIO_OPEN_OLD);
    	psio_address next_PSIF_DFSCF_K = PSIO_ZERO;
    	
    	double *in_buffer = init_array(norbs*nalpha);
    	
    	for (int Q = 0; Q<ri_nbf_; Q++)
    	{
            psio_read(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(in_buffer[0]),sizeof(double)*norbs*nalpha,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);

            for (int m = 0; m<norbs; m++)
                for (int n = 0; n<=m; n++)
                    for (int i = 0; i<nalpha; i++)
                    {
                Ka[m][n]+=in_buffer[m*nalpha+i]*in_buffer[n*nalpha+i];
                Ka[n][m] = Ka[m][n];
            }
    	}
    	
    	free(in_buffer);
    	psio_close(PSIF_DFSCF_K,0);
    }
    //fprintf(outfile,"  3-Index A Used \n"); fflush(outfile);

    //Second Pass: Cb_
    double **Kb; //Has to be defined outside the no work group
    if (nbeta_ > 0)
    {//Do work!

	
	Cocc = block_matrix(nbeta_,norbs);
	for (int i=0; i<norbs; i++) {
            for (int j=0; j<nbeta_; j++)
                Cocc[j][i] = Cb_->get(0,i,j);
        }

        //B_ia_P in core, B_im_Q in core
        if (df_storage_ == core||df_storage_ == double_core)
        {
            B_im_Q = block_matrix(norbs, nbeta_*ri_nbf_);
            double** QS = block_matrix(ri_nbf_,norbs);
            double** Temp = block_matrix(nbeta_,ri_nbf_);

            //print_mat(B_ia_P_,ri_nbf_,norbs*(norbs+1)/2,outfile);
            //fprintf(outfile,"\nYo\n");
            //print_mat(Cocc,nbeta_,norbs,outfile);
            for (int m = 0; m<norbs; m++) {
                for (int Q = 0; Q<ri_nbf_; Q++) {
                    for (int s = 0; s<norbs; s++) {
                        QS[Q][s] = B_ia_P_[Q][((s>=m)?ioff[s]+m:ioff[m]+s)];
                    }
                }
                C_DGEMM('N','T',nbeta_,ri_nbf_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], ri_nbf_);
                //print_mat(Temp,nbeta_,ri_nbf_,outfile);
                for (int Q = 0; Q<ri_nbf_; Q++) {
                    for (int i = 0; i<nbeta_; i++) {
                        B_im_Q[m][i+Q*nbeta_] = Temp[i][Q];
                    }
                }
            }
            free_block(QS);
            free_block(Temp);
        }

        //B_ia_P in disk, B_im_Q in core
        else if (df_storage_ == flip_B_disk || df_storage_ == k_incore)
        {
            double *in_buffer = init_array(norbs*(norbs+1)/2);
            psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
            psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
            B_im_Q = block_matrix(norbs, nbeta_*ri_nbf_);

            int mu, nu;

            for (int Q = 0; Q<ri_nbf_; Q++)
            {
    		psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    		for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
    		{
                    mu = ri_pair_mu_[ij];
                    nu = ri_pair_nu_[ij];
                    for (int i = 0; i<nbeta_; i++)
                    {
                        B_im_Q[mu][i+Q*nbeta_]+=Cocc[i][nu]*in_buffer[ij];
                        if (mu != nu)
                            B_im_Q[nu][i+Q*nbeta_]+=Cocc[i][mu]*in_buffer[ij];
                    }
    		}
            }

            psio_close(PSIF_DFSCF_BJ,1);
            free(in_buffer);
            //B_ia_P in disk, B_im_Q in disk
        } else {
            psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
            psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
            psio_open(PSIF_DFSCF_K,PSIO_OPEN_NEW);
            psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

            double *in_buffer = init_array(norbs*(norbs+1)/2);
            double *out_buffer = init_array(norbs*nbeta_);

            int mu, nu;
            for (int Q = 0; Q<ri_nbf_; Q++)
            {
    		for (int im = 0; im<nbeta_*norbs; im++)
                    out_buffer[im] = 0.0;

    		psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    		for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
    		{
                    mu = ri_pair_mu_[ij];
                    nu = ri_pair_nu_[ij];
                    for (int i = 0; i<nbeta_; i++)
                    {
                        out_buffer[mu*nbeta_+i]+=Cocc[i][nu]*in_buffer[ij];
                        if (mu != nu)
                            out_buffer[nu*nbeta_+i]+=Cocc[i][mu]*in_buffer[ij];
                    }
    		}
    		psio_write(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(out_buffer[0]),sizeof(double)*norbs*nbeta_,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
            }

            free(in_buffer);
            free(out_buffer);

            psio_close(PSIF_DFSCF_BJ,1);
            psio_close(PSIF_DFSCF_K,1);
        }

        free_block(Cocc);
        //fprintf(outfile,"  3-Index B Formed \n"); fflush(outfile);
 	Kb = block_matrix(norbs, norbs);

 	if (df_storage_ == core ||df_storage_ == double_core|| df_storage_ == flip_B_disk || df_storage_ == k_incore)
        {
            C_DGEMM('N','T',norbs,norbs,ri_nbf_*nbeta_,1.0,B_im_Q[0],ri_nbf_*nbeta_,B_im_Q[0],ri_nbf_*nbeta_, 0.0, Kb[0], norbs);
            free_block(B_im_Q);
        } else {
            psio_open(PSIF_DFSCF_K,PSIO_OPEN_OLD);
            psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

            double *in_buffer = init_array(norbs*nbeta_);

            for (int Q = 0; Q<ri_nbf_; Q++)
            {
    		psio_read(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(in_buffer[0]),sizeof(double)*norbs*nbeta_,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
    		
    		for (int m = 0; m<norbs; m++)
                    for (int n = 0; n<=m; n++)
                        for (int i = 0; i<nbeta_; i++)
                        {
                    Kb[m][n]+=in_buffer[m*nbeta_+i]*in_buffer[n*nbeta_+i];
                    Kb[n][m] = Kb[m][n];
                }
            }

            free(in_buffer);
            psio_close(PSIF_DFSCF_K,0);
        }
    } //Stop Work!


    for (int i=0; i<norbs; i++) {
        for (int j=0; j<=i; j++) {
            Ga_->set(0,i,j,J[i][j]-Ka[i][j]);
            if (nbeta_ > 0)
                Gb_->set(0,i,j,J[i][j]-Kb[i][j]);
            if (i!= j)
            {
                Ga_->set(0,j,i,J[i][j]-Ka[i][j]);
                if (nbeta_ > 0)
                    Gb_->set(0,j,i,J[i][j]-Kb[i][j]);
            }
        }
    }
    //Ga_->print();
    //Gb_->print();
    //fprintf(outfile,"\n");
    //G_.print();
    free_block(J);
    free_block(Ka);
    if (nbeta_ > 0)
    	free_block(Kb);

    **/

}
void ROHF::form_G_from_RI()
{
    fprintf(stderr, "ROHF RI Not implemented yet!\n");
    abort();
}

}}
