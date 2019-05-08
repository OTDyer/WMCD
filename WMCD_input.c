 /*
 * Functions for passing information into WMCD.
 * This version of WMCD includes passive polymeric systems, smart MC,
 *      and affine deformations.
 *
 * Copyright 2019 Oliver T. Dyer
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "WMCD_header.h"
#include "pcg_basic.h"


void input(System_parameters *params_input, Potentials *pot_input,
           In_out_options *out_input, double arg_option)
{
    /**-----------------------------*/
    /** Box parameters              */
    /**-----------------------------*/

    /** Use periodic or unbounded boundary conditions */
    params_input->periodic = false;

    /* Set wrapping functions by the boundary conditions */
    if(params_input->periodic == false)
    {params_input->pwrap = &dont_wrap;
     params_input->pshort_wrap = &dont_short_wrap;}
    else
    {params_input->pwrap = &wrap;
     params_input->pshort_wrap = &short_wrap;}


    /** Sub-box side length */
    params_input->sub_box_size = 1;
    // NOTE: should be at least as large as the maximum range of the pot() function

    /** Box side length */
    params_input->L = 20;
    // NOTE: must be an integer multiple of bs

    /** Number of particles allowed in a box */
    params_input->n_box_max = 51;
    // NOTE: the value is arbitrary, but balances memory usage with risk of aborting
    //       due to overfull sub-boxes. Higher risk with theta-solvent.


    /**-----------------------------*/
    /** Chain parameters            */
    /**-----------------------------*/

    /** Number of monomers in each polymer */
    //params_input->Nmono = 4;
    params_input->Nmono = (int)arg_option;

    /** Particle concentration = Ntot/L^{3} (periodic system) */
    double concentration = 0.183;
    // NOTE: can be used to pick Npoly.

    /** Number of polymer chains */
    params_input->Npoly = 1;
    //params_input->Npoly = nint(concentration*pow_int(params_input->L,3)/params_input->Nmono);

    /** Total number of particles */
    params_input->Ntot = params_input->Nmono * params_input->Npoly;
    // MUST BE Nmono*Npoly.



    /**-----------------------------*/
    /** Move parameters             */
    /**-----------------------------*/

    /** Minimum wavelet radius */
    params_input->lambda_min = 0.7;

    /** Maximum wavelet radius */
    params_input->lambda_max = 0.5*(params_input->L - params_input->sub_box_size);
    // NOTE: should not exceed 0.5*(L-bs). This will later be reduced further if shear flows are imposed

    /** Determine maximum wavenumber in periodic systems */
    params_input->Lcart = 15;
    // Needs to be high enough for good accuracy but if it is too high you
    // the computation time can increase appreciably. 15 is a good value.

    /** Move amplitude */
    params_input->A0 = 0.25;
//    params_input->A0 = arg_option;

    /** Turn plane wave moves on or off */
    params_input->use_Fourier = true;

    /** Pick total number of moves or number of time units */
    params_input->set_time_max = true;
    // This is typically set to true for real simulations.
    // Setting to false can be very useful for testing and debugging with small numbers of moves.
    // For reference, ~10^5 moves per time unit is quite usual.
    if(params_input->set_time_max == false)
    {params_input->movenum_max = 1E3;}
    else if (params_input->set_time_max == true)
    {params_input->time_max = 1E5;}
    else
    {printf("Error: Nmoves not chosen\n");
     abort();}

    /** Use smart Monte Carlo */
    params_input->smart_MC = true;

    /** Use the code as a Langevin solver */
    /* Basically smart MC without the MC test.
     * Note it will only be treated as true if smart_MC is also true.
     * It is not advisable with diverging potentials such as WCA and FENE unless
     * A0 is sufficiently small (FENE is safer than WCA). */
    params_input->Langevin = false;

    /** Choice of Monte Carlo test: true = Glauber; false = Metropolis */
    /* either can be used with smart_MC == true */
    params_input->use_Glauber = false;

    /** Use move recycling to correct move distribution */
    params_input->recycle = false;
    // This is only important when smart_MC == false.
    // If smart_MC == true, best to set recycle = false.

    /**-----------------------------*/
    /** Physical parameters         */
    /**-----------------------------*/

    /** Thermal energy */
    params_input->kT = 1.0;
    params_input->beta = 1.0/params_input->kT;
    // beta MUST be set to 1/kT.

    /** Ratio between lambda_min and particle hydrodynamic radius */
    params_input->lambda_a = 2.316;
    // NOTE: this is fixed by the choice of mother wavelet
    params_input->hydro_rad = params_input->lambda_min/params_input->lambda_a;
    // MUST be lambda_min / lambda_a.

    /** average time increment after a single move */
    params_input->dt = 6 * pow_int(params_input->lambda_a * params_input->A0, 2) / (105*params_input->Ntot);
    if(params_input->use_Fourier == false)
    {params_input->dt = params_input->dt * (1 - (params_input->lambda_min/params_input->lambda_max) )/(1 -
                                                 pow_int((params_input->lambda_min/params_input->lambda_max), 3) );
    }
    // The details for these can be found in
    //      O T Dyer & R C Ball, J Chem Phys vol. 146, 124111 (2017).
    // Note the numerical values depend on the choice of mother wavelet.

    /** move independent factors in variance and mobility*dt for generation of amplitude*orientation */
    params_input->sMC_var = params_input->A0 * params_input->A0 / 3.0;
    params_input->sMC_mob = params_input->sMC_var / (2 * params_input->kT);
    // These are inter-related by the Fluctuation Dissipation Theorem, so must have sMC_mob = sMC_var/(2kT).
    // How sMC_var depends on WMCD parameters is set out in ??????

    /** -----------------------   Solvent shear rate ---------------------------- */

    params_input->shear_rate = 0.0; /* note 0.01 is a large value in our units */

    /* Reduce maximum wavelet size to fit inside sheared boxes */
    if(params_input->shear_rate != 0.0)
    {params_input->lambda_max = params_input->lambda_max / 1.118;}

    /* Time interval between remapping back onto the initial sheared coordinates */
    if(params_input->shear_rate != 0.0)
    {params_input->reset_time = (params_input->sub_box_size)/params_input->shear_rate;}
    else if(params_input->set_time_max == true)
    {params_input->reset_time = 2.0*params_input->time_max;}
    else
    {params_input->reset_time = 2.0*params_input->movenum_max;}


    /** ----------------------------- solvent type ------------------------------- */

    //char solvent_str[10];
    char solvent_str[10] = "good";
    //char solvent_str[10] = "theta";
    // In good solvent particles repel one another.
    // In theta solvent only bond forces are present.
    // This terminology is standard in polymer literature.

    strcpy(pot_input->solvent, solvent_str);

    /** Estimate of Zimm time, for when to output uncorrelated data */
    if(strcmp(solvent_str,"good") == 0)
    {params_input->Zimm_time = 1000.0 * pow((double)params_input->Nmono / 20, 3 * 0.588);}
    else if(strcmp(solvent_str,"theta") == 0)
    {params_input->Zimm_time = 1000.0 * pow((double)params_input->Nmono / 20, 3 * 0.5);}
    // 1000 time units is approximately the Zimm time for Nmono = 20 chains with bond length approx. 1.
    // The ideal value to set is one measured in simulations.

    char pot_str[10];
    /** potential type */
    if(strcmp(solvent_str,"good") == 0)
    {
        /* comment out lines to pick the desired excluded volume (EV) potential */
        strcpy(pot_str,"WCA");   params_input->EV_fn = &EV_WCA;
        //strcpy(pot_str,"NGauss");   params_input->EV_fn = &EV_NGauss;
        //strcpy(pot_str,"Exp");   params_input->EV_fn = &EV_Exp;
    }
    else
    {
        strcpy(pot_str,"free");
        params_input->EV_fn = &EV_free;
        pot_input->pot_length = 0.0;
        pot_input->pot_energy = 0.0;
    }
    strcpy(pot_input->pot_type, pot_str);

    if(strcmp(pot_input->pot_type,"WCA") == 0)
    {
        pot_input->pot_length = 0.89090;
        pot_input->pot_energy = 1.0;
    }
    else if(strcmp(pot_input->pot_type,"NGauss") == 0)
    {/** PRELIMINARY VALUES!! */
        pot_input->pot_length = 0.3; /* must be less than 1/3 */
        pot_input->pot_energy = params_input->kT * pow_int(1.0/params_input->hydro_rad, 3) * 8;
    }
    else if(strcmp(pot_input->pot_type,"Exp") == 0)
    {
        pot_input->pot_length = 0.25; /* must be less than 1/3 */
        pot_input->pot_energy = 75;
        //pot_input->pot_energy = arg_option;
    }


/* ------------------------------------------------*/


    char bond_str[10];
    /** Bond potential type */
    if(strcmp(solvent_str,"good") == 0)
    {
        /* comment out lines to pick the desired bond potential */
        strcpy(bond_str,"FENE");  params_input->bond_fn = &bond_FENE;
        //strcpy(bond_str,"quadratic");   params_input->bond_fn = &bond_quad;
    }
    else
    {
        strcpy(bond_str,"quadratic");
        params_input->bond_fn = &bond_quad;
    }
    strcpy(pot_input->bond_type, bond_str);

    /** Length scale in bondpot, R0 in FENE, mean separation in Hookian */
    if(strcmp(pot_input->bond_type,"FENE") == 0)
    {pot_input->bond_length = 1.7818;}
    else
    {pot_input->bond_length = 0.9691;}

    /** Spring constant in bondpot*/
    pot_input->spring_const = 8.8194;
    if(pot_input->bond_length == 0.0)
    {pot_input->spring_const = params_input->kT / (4*params_input->hydro_rad*params_input->hydro_rad);}


    /**-----------------------------*/
    /** Output options              */
    /**-----------------------------*/

    /* Note the 'rescale' stuff is now redundant, but I haven't yet unpicked it */
    /** Number of times to output in each output regime */
    /** get N_out outputs in total */
    if(params_input->set_time_max == false)
    {out_input->N_out = 0.1*params_input->movenum_max/params_input->Ntot;}
    else if (params_input->set_time_max == true)
    {out_input->N_out = 1.0 * params_input->time_max;}


    /* In all below, choose true to turn on the output */

    /** Output system parameters for the run */
    out_input->system_params = true;

    /** Whether to make and/or read in an equilibrated config */
    out_input->eq_config_out = false;
    out_input->eq_config_in  = false;

    /** Whether to make and/or read in a config for an unfinished run */
    out_input->cont_config_out = false;
    out_input->cont_config_in  = false;
    /* NOTE: > cont_config_in = true supersedes eq_config_in = true
     *       > it is normally correct for *_out to be true if *_in is so that the
     *             configuration gets updated at the end (and not reset)
     *       > the reverse is not the case as the first run should read in
     *             an equilibrated config
     */

    /** Unwrapped particle locations */
    out_input->unwrapped_config = true;

    /** Wrapped particle locations */
    out_input->wrapped_config = false;

    /** Snapshots of forces on particles */
    out_input->forces_out = false;
    /* NOTE: this will only be meaningful if smart MC is used */

    /** Accept-reject statistics */
    out_input->accept_reject = false;

    /** Evolution of end-to-end and gyration radii squared */
    out_input->polymer_radii = false;

    /** List of centres of mass of polymers at output time **/
    out_input->polymer_CoM = false;
    out_input->off_diag = false;

    /** Evolution of centre of mass position --> chain diffusivity */
    out_input->CoM_displacement = false;

    /** Structure factor */
    out_input->static_structure_factor = false;


}


/* read in equilibrated particle configurations */
void read_equilib_config(double xyz[][3], double xyzn[][3], double xyz0[][3], double C_to_xyz[3][3],
                         System_parameters params, Potentials pot_params, char *directory)
{
    int i,j;
    int Ntot = params.Nmono*params.Npoly;
    int Nconfigs; /* stores number of available configurations to start from */
    double Cart_coords[3];

    char max_name[150], file_name[150]; /* strings to hold file names and directories */
    FILE *config_file;

    /*----------------------------------------*/
    /* obtain the number of available configs */
    /*----------------------------------------*/

    /* build file name from relevant system parameters */
    sprintf(max_name,"%sEquilib_configs/Nconfigs_%d_%d_%s.txt",
            directory, params.Nmono, params.Npoly,pot_params.solvent);

    /* open file and read in the number of available configs */
    config_file = fopen(max_name,"r");
    if(config_file!=NULL)
    {
        fscanf(config_file,"%d",&Nconfigs);

        fclose(config_file);
    }
    else
    {printf("Error: equilibrated configuration not found.\n");
     abort();}


    /*----------------------------------------*/
    /* read in particle configurations        */
    /*----------------------------------------*/

    /* build file name from relevant system parameters */
    sprintf(file_name,
            "%sEquilib_configs/config_%d_%d_%s_%d.txt",
            directory, params.Nmono, params.Npoly,
            pot_params.solvent,
            pcg32_boundedrand(Nconfigs) + 1 );


    config_file = fopen(file_name,"r");

    if(config_file!=NULL)
    {
        /* read in coordinates */
        for(i=0; i<Ntot; i++)
        {for(j=0; j<3; j++)
            {fscanf(config_file,"%lf",&Cart_coords[j]);}

         transform_coords(Cart_coords, xyz[i], C_to_xyz);
        }

        fclose(config_file);
    }
    else
    {printf("Error: equilibrated configuration not found.\n");
     abort();}

    /* copy to 'new' and 'original' locations */
    for(i=0; i<Ntot; i++)
    {for(j=0; j<3; j++)
        {
            xyzn[i][j] = xyz[i][j];
            xyz0[i][j] = xyz[i][j];
        }
    }
}


/* read in configuration from an unfinished run */
void read_continue_config(double xyz[][3], double xyzn[][3], double xyz0[][3], int Ntot,
                    double *time, System_parameters params, char *directory, char *suffix)
{
    int i,j;
    double Cart_coords[3];

    double a1[3], a2[3], a3[3];
    double C_to_xyz[3][3];

    char file_name[150]; /* strings to hold file names and directories */
    FILE *config_file;

    /* build file name from relevant system parameters */
    sprintf(file_name, "%sContinue_configs/config%s", directory, suffix);

    config_file = fopen(file_name,"r");

    if(config_file!=NULL)
    {
        /* read in time */
        fscanf(config_file, "%lf", time);

        lattice_vectors(a1, a2, a3, *time, params);
        Cart_to_xyz_transform(C_to_xyz, a1, a2, a3);

        /* read in coordinates */
        for(i=0; i<Ntot; i++)
        {for(j=0; j<3; j++)
            {fscanf(config_file,"%lf",&Cart_coords[j]);}
        transform_coords(Cart_coords, xyz[i], C_to_xyz);
        }

        fclose(config_file);
    }
    else
    {printf("Error: 'continue' configuration not found.\n");
     abort();}

    /* copy to 'new' and 'original' locations */
    for(i=0; i<Ntot; i++)
    {for(j=0; j<3; j++)
        {
            xyzn[i][j] = xyz[i][j];
            xyz0[i][j] = xyz[i][j];
        }
    }

}
