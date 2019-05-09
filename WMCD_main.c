/*
 * Core functions for WMCD.
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


/* regarding passing arrays as arguments, I use:
 *       *array      when I modify it in the function
 *       array[n]    when I don't modify it and know its size
 *       array[]     when I don't modify it and it size isn't always the same
 *
 * the aim is to clearly see what is being done with array arguments
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "WMCD_header.h"
#include "pcg_basic.h"

// max random numbers to give closed set [0,1] upon multiplication of RNG
#define RAND_MAX_INV 1.0/RAND_MAX
#define POW2_32_INV 1.0/(pow(2,32)-1)
// 'max' random numbers to give open set [0,1) upon multiplication of RNG
#define POW2_32_INV_OPEN 1.0/pow(2,32)

/* All functions that differ between periodic and infinite boxes */

void get_wavelet_parameters(wavelet_parameters *pwavelet, double Rmin, double Rmax,
             bool recycle, Recycle_wavelet *recycled_w, int Ntot, double xyz[][3], double transform[3][3]);

void planewave(Fourier_parameters *pplane_wave, double CDF[][2]);
void planewave_cart(Fourier_parameters *pplane_wave, double k0, double Rmax, int Lcart,
                    double kextent[3], double kmetric[3][3]);

void movew(double *position, System_parameters params, wavelet_parameters wavelet, double xyz_to_C[3][3],
           double C_to_xyz[3][3], double metric[3][3], double images[8][3],double image_dist[8]);

void movepw(double *position, Fourier_parameters plane_wave, double metric[3][3]);

double UF_move(double xyz[][3],System_parameters params,double forces[][3],bool find_forces,
               double metric[3][3],int *anchorxyz,int *spanxyz,
               double (*EV_pot)(double, Potentials, double *), double (*bond_pot)(double, Potentials, double *),
               Potentials pot_params,int Lbs,
               int movers[Lbs][Lbs][Lbs][params.n_box_max],int * non_movers[Lbs][Lbs][Lbs],
               int nmovers,long movenum,long *bondref,const int casesall[13][3],
               double images[8][3],double image_dist[8]);

void evolve(double xyz[][3],double xyzn[][3],double forces[][3],System_parameters params,
            Potentials pot_params,double time, int Lbs,int contents[Lbs][Lbs][Lbs][params.n_box_max],
            int *movers, int moving_pcls[Lbs][Lbs][Lbs][params.n_box_max],
            int non_movers[Lbs][Lbs][Lbs][params.n_box_max], int * pnon_movers[Lbs][Lbs][Lbs],
            long movenum,long *paccptnum,long bondref[][params.Nmono*params.Npoly],const int casesall[13][3],
            double move_accpt[3],Recycle_wavelet *recycled_w,Recycle_Fourier *recycled_F,
            double Ppw, double pwCDF[][2], double R_short[3],
            double dR_acc[3], double dR_acc_reduced[3]);

/* In this program, the main function is structured as:
 *      > declares variables and arrays needed throughout the simulation
 *      > loop calling the evolve function and most outputs
 *      > final outputs etc. and program end.
 */
int main(int argc, char **argv)
{

    srand((long)time(NULL) ^ atoi(argv[1]));
    int seed[2] = {rand(), atoi(argv[2])+1};

    /* seed random numbers using run/job number */
    pcg32_srandom(seed[0], seed[1]);

    long i,j,k; /* looping indices */
    double run_time; /* run timer. uses clock() */

    /* get input parameters from input.c */
    System_parameters params;
    System_parameters *pparams=&params;

    Potentials pot_params;
    Potentials *ppot_params=&pot_params;

    In_out_options out_options;
    In_out_options *pout_options=&out_options;

    input(pparams, ppot_params, pout_options, atof(argv[4]));

    /* change a chosen parameter if want multiple jobs that change this */
    //    if(atof(argv[4]) != 0.0)
    // {params.A0 = atof(argv[4]);}

    /* Get parts of file names common to all outputs */
    char file_directory[50], file_suffix[60];
    get_directory_and_suffix(file_directory, file_suffix, params, pot_params,
                             (atoi(argv[3]) == 0)? seed[0] : atoi(argv[1]),
                             atoi(argv[2]), false);


    /**-----------------------------*/
    /** Box parameters              */
    /**-----------------------------*/

    printf("\n%s box wavelet+Fourier code\n",(params.periodic==true)?"Periodic":"Infinite");

    const int Lbs = params.L/params.sub_box_size; /* L and bs must be comensurate */

    printf("\nINPUTS:\nL = %d, with sub-box length %d\n", params.L, params.sub_box_size);

    /* index for the sub-boxes the particles are found in */
    int (*contents)[Lbs][Lbs][params.n_box_max] = calloc(Lbs, sizeof(*contents));

    /* arrays which split contents into movers and non-movers each move.
     * these are only updated in the 'spans' region, so non_move_ptr efficently directs searches to
     * use either non_move_pcls or contents if inside or outside this region respectively.
     */
    int (*moving_pcls)[Lbs][Lbs][params.n_box_max] = calloc(Lbs, sizeof(*moving_pcls));
    int (*non_move_pcls)[Lbs][Lbs][params.n_box_max] = calloc(Lbs, sizeof(*non_move_pcls));
    int * (*non_move_ptr)[Lbs][Lbs] = malloc(Lbs * sizeof(*non_move_ptr));

    /* Check all arrays are created */
    if(contents == NULL || moving_pcls == NULL || non_move_pcls == NULL || non_move_ptr == NULL)
    {printf("!! Insufficient memory for contents !!\n");
     return 1;}

    for(i=0; i<Lbs; i++)
    {for(j=0; j<Lbs; j++)
        {for(k=0; k<Lbs; k++)
            {non_move_ptr[i][j][k] = contents[i][j][k];}}}

    /* vectors that identify the region containing moving particles */
    int anchormain[3] = {0,0,0}, spansmain[3] = {Lbs-1,Lbs-1,Lbs-1};

    /* list of directions for neighbouring boxes without double counting */
    const int casesall[13][3] = {{1,0,0}, {1,0,-1}, {1,-1,0}, {1,-1,-1}, {0,1,0},
         {0,1,-1},{0,0,1},{1,1,0}, {1,1,-1}, {0,1,1}, {1,0,1}, {1,-1,1}, {1,1,1}};

    /* lattice vectors and coordinate transform */
    double a1[3], a2[3], a3[3];
    double xyz_to_C[3][3], metric[3][3];


    /**-----------------------------*/
    /** Chain parameters            */
    /**-----------------------------*/

    /* Number(s) of particles. */
    const int Ntot = params.Npoly*params.Nmono;

    printf("N = %d, Npoly = %d, Nmono = %d\n",Ntot, params.Npoly, params.Nmono);

    /* list of current (xyz) and 'new' (xyzn) position vectors of particles */
    double xyz[Ntot][3];
    double xyzn[Ntot][3];
    /* list of original positions, for use in measuring net motion */
    double xyz0[Ntot][3];

    /* lists tracking what particles do in a move */
    int movers[Ntot]; /* tracks which particles move in a wavelet */
    long bondref[3][Ntot]; /* tracks which bond energies have been calculated in a move */
    for(i=0; i<Ntot; i++)
    {movers[i] = i+1;
     bondref[0][i] = 0;
     bondref[2][i] = 0;
     bondref[1][i] = 0;}


    /**-----------------------------*/
    /** Wavelet parameters          */
    /**-----------------------------*/

    printf("lambda_min = %g, lambda_max = %g\n",params.lambda_min, params.lambda_max);
    printf("A_0 = %g\n",params.A0);


    /**-----------------------------*/
    /** Fourier parameters          */
    /**-----------------------------*/

    /* Probability of picking a plane wave over a wavelet
     * this is different (and simpler) for the continuous Fourier spectrum
     */
    double Ppw; /* Ppw[0] = prob pick a PW, Ppw[1] = prob pick it via discrete method */
    double Sigmapw; /* sum total (unnormalised) probabilities and those out to edge of Cartesian/discrete region */

    if(params.use_Fourier == false) // if plane waves turned off
    {Ppw = 0.0;}
    else if(params.periodic == false)
    {Ppw = pow_int(params.lambda_min/params.lambda_max, 3);}
    else
    {
        /* sum out a number of shells far enough to have converged. this is L dependent */
        Sigmapw = pwProbSum(2*M_PI*params.lambda_max/params.L, 500);

        Ppw = 1.0/( 1.0 + (4*M_PI*pow_int(params.L/params.lambda_max, 3)*
                    (pow_int(params.lambda_max/params.lambda_min, 3) - 1.0 ))/(15*Sigmapw) );
    }

    printf("prob. of selecting a plane wave move: %G\n",Ppw);


    /* CDF for choosing wavenumber, converts between
     * random no.s with uniform dist (..[][0]) and dist from PDF (..[][1])
     * NOTE: this is only used with an unbounded system
     */
    int NCDFelements=38751; //Number of lines in FourierCDF.txt

    double pwCDF[NCDFelements][2];

    FILE *pwCDFfile;
    pwCDFfile = fopen("FourierCDF.txt","r");

    if(pwCDFfile == NULL)
    {printf("Error: FourierCDF failed to open\n");
     abort();}
    /* Read CDF from file into pwCDF */
    for(i=0; i<NCDFelements; i++)
    {
        fscanf(pwCDFfile, "%lf %lf", &pwCDF[i][0],&pwCDF[i][1]);
        /* Rescale to give k rather than k*Rmax */
        pwCDF[i][1] /= params.lambda_max;
    }
    fclose(pwCDFfile);



    /**-----------------------------*/
    /** Physical parameters         */
    /**-----------------------------*/

    double time, t0, *pt0 = &t0;
    int n_resets=0;

    printf("Time increment after %d accepted %s = %lf\n",Ntot, (Ntot==1)?"move":"moves", params.dt*Ntot);
    printf("Using %s and %s potentials\n", pot_params.pot_type, pot_params.bond_type);
    if(params.Langevin == false)
    {
        if(params.smart_MC == true)
        {printf("Smart MC is used ");
        }
        printf("%sith the %s test\n",(params.smart_MC==true)?"w":"W",(params.use_Glauber==true)?"Glauber":"Metropolis");
    }
    else if(params.smart_MC == true && params.Langevin == true)
    {printf("Run as a Langevin solver\n");}

    /* list of forces on movers.
     * note this does not need (re)initialising as the relevant elements are updated each move
     */
    double forces[Ntot][3];


    /**-----------------------------*/
    /** Move counters               */
    /**-----------------------------*/

    /* Number of accepted moves to make */
    long Nmoves;
    if(params.set_time_max == false)
    {Nmoves = params.movenum_max;}
    else
    {Nmoves = long_nint(params.time_max/params.dt);}

    long movenum = 0; /* current move number */

    /* counters associated with when to output information */
    long output_move = 0, output_moveS = 0, doutput_move = long_nint((double)Nmoves/out_options.N_out);
    long Zimm_count = 1;

    /* used to keep track of failed moves, so that parameters can be reused in the next move */
    Recycle_wavelet recycled_w = {false,1.0,{0.0,0.0,0.0}};
    Recycle_wavelet *precycled_w = &recycled_w;

    Recycle_Fourier recycled_F = {false,1.0,{0.0,1.0,0.0}};
    Recycle_Fourier *precycled_F = &recycled_F;


    /**-----------------------------*/
    /** Accecpt/Reject parameters   */
    /**-----------------------------*/

    long accptnum = 0; /* number of accepted moves */
    long *paccptnum = &accptnum;

    /* variables for accpt/reject statistics */
    int NRbins=20, Nkbins=20; /* number of bins for R and k accept/reject statistics */
    double move_accpt[3]; /* tracks 'size' and success of single move. {w or F, R or k value, success or fail} */
    double R_accpt[NRbins][3], k_accpt[Nkbins][3];
    double wRbins=pow(params.lambda_max/params.lambda_min,1.0/NRbins);/* scale factor of R bins */
    double wkbins=pow(320.0*params.L/(11.0*params.lambda_max),1.0/(Nkbins-1)); /* scale factor of k bins */
    /*  NOTE: 320.0 comes from the upper limit in FourierCDF.txt, while 11.0>2*pi*sqrt(3) */

    /* initialise R_accpt and k_accpt */
    for(i=0; i<NRbins; i++)
    {
        R_accpt[i][0] = pow_int(wRbins,i+1)*params.lambda_min;
        R_accpt[i][1]=0;
        R_accpt[i][2]=0;
    }
    for(i=0; i<Nkbins; i++)
    {
        k_accpt[i][0] = pow_int(wkbins,i+1)*11.0/params.L;
        k_accpt[i][1]=0;
        k_accpt[i][2]=0;
    }


    /**-----------------------------*/
    /** Data analysis objects       */
    /**-----------------------------*/

    /* track CoM position at chosen times */
    int max_power = 17; /* largest power of 2 used for lag separation */
    double R_CoM[3];
    cumulative_data R_CoM_squared[2*max_power+3];
    prev_quantities prev_R_CoM[max_power+1][3];
    /* initialise R_CoM_squared */
    initialise_cumulative_variable(R_CoM_squared, prev_R_CoM, max_power);

    /* track the difference in CoM displacement between rejected and recyled moves */
    double R_short[3] = {0.0,0.0,0.0};
    cumulative_data R_short_squared[2*max_power+3];
    prev_quantities prev_R_short[max_power+1][3];
    /* initialise dacc_rej */
    initialise_cumulative_variable(R_short_squared, prev_R_short, max_power);

    /* adjust physical run time so all lags are commensurate at the end */
    double lag_max = (double)(R_CoM_squared[2*max_power + 2].num_lags)*(out_options.N_out)/(params.time_max);
    if(params.set_time_max == true && out_options.CoM_displacement == true && params.time_max > 2.0*lag_max)
    {
        params.time_max = 2.0 * lag_max * ceil(params.time_max / (2*lag_max));
        Nmoves = long_nint(params.time_max / params.dt) + 1;
        doutput_move = long_nint((double)Nmoves / params.time_max);
    }

    double dR_acc[3] = {0.0,0.0,0.0};
    double dR_acc_red[3] = {0.0,0.0,0.0};


    /**--------------------------------------------------------------------------------*/
    /**--------------------------------------------------------------------------------*/

    printf("------------------------------\n");

    /**---------------------------*/
    /** Construct initial system  */
    /**---------------------------*/

    if(out_options.cont_config_in == true) /* continue from previously generated data   */
    {
        read_continue_config(xyz, xyzn, xyz0, Ntot, pt0, params, file_directory, file_suffix);
        output_move += doutput_move; /* don't output time=t0 data again */
        output_moveS += doutput_move;
    }
    else if(out_options.eq_config_in == true) /* use previously equilibrated configuration */
    {
        t0 = 0.0;
        lattice_vectors(a1, a2, a3, t0, params);
        xyz_to_Cart_transform(xyz_to_C, a1, a2, a3);

        read_equilib_config(xyz, xyzn, xyz0, xyz_to_C, params, pot_params, file_directory);
    }
    else/* construct initial configuration from scratch */
    {
        build_system(xyz, xyzn, xyz0, params, pot_params);
        t0 = 0.0;
    }

    /* fill contents[] */
    map(Ntot,movers,xyzn,params.L,params.sub_box_size,Lbs,params.n_box_max,contents,
        anchormain,spansmain);

    /* reinitialise movers() */
    for(i=0; i<Ntot; i++)
    {movers[i] = 0;}
    calc_R_CoM(R_short, xyz, params.Nmono, 1);


    /**---------------------------*/
    /** Evolve the system         */
    /**---------------------------*/

    /* start clock */
    run_time = (double)clock()/CLOCKS_PER_SEC;
    time_hms run_time_hms = {0, 0, 0.0};
    time_hms *prun_time_hms=&run_time_hms;

    if(out_options.system_params == true)
    {write_parameters(params, pot_params, seed, file_directory, file_suffix,
                      run_time_hms, 1.0, true);}

    /* loop over the moves */
    while(accptnum<Nmoves)
    {
        /* Update time */
        time = t0 + params.dt*accptnum;

        /* Move particle positions in the shearing reference frame if resetting the coordinate transform */
        if(time - (n_resets+1)*params.reset_time > 0)
        {remap_at_coord_reset(xyz, xyzn, Ntot, time - n_resets*params.reset_time, params, Lbs, contents);
         n_resets++;}

        /* write to file after output_move ACCEPTED moves */
        if(accptnum == (long)(output_move + 0.00001))
        {
            lattice_vectors(a1, a2, a3, time-(n_resets*params.reset_time), params);
            xyz_to_Cart_transform(xyz_to_C, a1, a2, a3);

            /* Polymer centres of mass */
            if(out_options.polymer_CoM==true)
            {R_CoM_file(xyz, params.Npoly, params.Nmono, time, xyz_to_C, file_directory, file_suffix);}

            /* Polymer diffusion */
            if(out_options.CoM_displacement==true)
            {
                calc_R_CoM(R_CoM, xyz, params.Nmono, 1);
                calc_dR_CoM_squared(R_CoM_squared, prev_R_CoM, R_CoM,
                                     time, nint(output_move/doutput_move), max_power);
            }

            /* Off-diagonal part of diffusion */
            if(out_options.off_diag==true)
            {off_diag_file(dR_acc, dR_acc_red, time, file_directory, file_suffix);}

            output_move += doutput_move;

        }

        /* Write to file after output_move ATTEMPTED moves */
        if(movenum == (long)(output_moveS + 0.00001))
        {

            /* Polymer diffusion */
            if(out_options.CoM_displacement==true)
            {
                calc_dR_CoM_squared(R_short_squared, prev_R_short, R_short,
                            t0 + params.dt*movenum, nint(output_moveS/doutput_move), max_power);
            }

            output_moveS += doutput_move;
        }

        /* Write to file after each Zimm time. This is particularly useful for static quantities of polymers */
        if(time > Zimm_count * params.Zimm_time)
        {
            lattice_vectors(a1, a2, a3, time-(n_resets*params.reset_time), params);
            xyz_to_Cart_transform(xyz_to_C, a1, a2, a3);

            /* Particle locations */
            if(out_options.unwrapped_config==true)
            {unwrapped_coords_file(xyz, Ntot, time, xyz_to_C, file_directory, file_suffix);}

            /* Particle locations in periodic coordinates */
            if(out_options.wrapped_config==true)
            {wrapped_coords_file(xyz, Ntot, params.L, time, xyz_to_C, file_directory, file_suffix);}

            /* Particle forces */
            if(out_options.forces_out==true)
            {write_forces(forces, Ntot, time, file_directory, file_suffix);}

            /* Polymer radii evolutions */
            if(out_options.polymer_radii==true)
            {get_metric(metric,a1,a2,a3);
             radii_file(xyz, params.Npoly, params.Nmono, time, metric, file_directory, file_suffix);}

            Zimm_count++;
        }

        /* Increment attempted move number. used to check if bond potentials have been calulated */
        movenum++;

        /* Make a WMCD move */
        evolve(xyz,xyzn,forces,params,pot_params,time,Lbs,contents,movers,
               moving_pcls, non_move_pcls, non_move_ptr,
               movenum,paccptnum,bondref,casesall,move_accpt,precycled_w,
               precycled_F, Ppw, pwCDF, R_short, dR_acc, dR_acc_red);


        /* Record accept/reject data */
        if(out_options.accept_reject==true)
        {add_accept_reject_stats(move_accpt, R_accpt, k_accpt);}
    }


    /* stop clock */
    run_time = ( (double)clock()/CLOCKS_PER_SEC ) - run_time;

    split_time(run_time, prun_time_hms);


    /**-----------------------------*/
    /** Print some info and finish  */
    /**-----------------------------*/

    lattice_vectors(a1, a2, a3, time-(n_resets*params.reset_time), params);
    xyz_to_Cart_transform(xyz_to_C, a1, a2, a3);
    get_metric(metric, a1, a2, a3);

    if(out_options.system_params == true)// && atoi(argv[2]) == 1)
    {write_parameters(params, pot_params, seed, file_directory, file_suffix,
                      run_time_hms, (double)Nmoves/movenum,false);}

    /* Build equilibrated configuration for later use */
    if(out_options.eq_config_out==true)
    {write_equilib_config(xyz, Ntot, xyz_to_C, file_directory, file_suffix);}

    /* write configuration to continue this run later */
    if(out_options.cont_config_out==true)
    {write_continue_config(xyz, Ntot, t0 + params.dt*accptnum, xyz_to_C, file_directory, file_suffix);}

    /* write CoM displacement data */
    if(out_options.CoM_displacement==true)
    {
        write_cumulative_data(R_CoM_squared, max_power, file_directory, file_suffix, "dR_CoM");
        write_cumulative_data(R_short_squared, max_power, file_directory, file_suffix, "dR_short");
    }

    /* Write accept-reject statistics file */
    if(out_options.accept_reject==true)
    {accpt_files(R_accpt, k_accpt, NRbins, Nkbins);}

    /* calculate the static structure factor for the final state */
    if(out_options.static_structure_factor==true)
    {structure_factor_file(xyz, params.Npoly, params.Nmono, metric);}


    printf("\nTotal runtime for %g %s: %dh %dm %.3lfs\n",
           (params.set_time_max==false)?(double)Nmoves:(double)Nmoves*params.dt,
           (params.set_time_max==false)?"moves":"units of time",
           run_time_hms.hours, run_time_hms.minutes, run_time_hms.seconds);
    printf("Global fraction of moves accepted: %lf\n\n",(double)Nmoves/movenum);

    /* list output files */
    printf("OUTPUTS:\n%s%s%s%s%s%s%s%s%s%s%s%s",
           (out_options.system_params==true)?   "System parameters\n":"",
           (out_options.eq_config_out==true)?   "Equilibrated configuration\n":"",
           (out_options.cont_config_out==true)? "Continue configuration\n":"",
           (out_options.unwrapped_config==true)?"Unwrapped coords\n":"",
           (out_options.wrapped_config==true)?  "Wrapped coords\n":"",
           (out_options.forces_out==true)?      "Forces\n":"",
           (out_options.accept_reject==true)?   "Accept-reject stats\n":"",
           (out_options.polymer_radii==true)?   "Polymer radii evolution\n":"",
           (out_options.polymer_CoM==true)?     "Polymer centre of mass\n":"",
           (out_options.off_diag==true)?        "Diffusion off-diagonal\n":"",
           (out_options.CoM_displacement==true)?"Polymer diffusion\n":"",
           (out_options.static_structure_factor==true)?"Static structure factor\n":"");

    printf("------------------------------\n");
    printf("------------------------------\n");

    free(contents);
    free(moving_pcls);
    free(non_move_pcls);
    free(non_move_ptr);

    return EXIT_SUCCESS;
}




void evolve(double xyz[][3],double xyzn[][3],double forces[][3],System_parameters params,
            Potentials pot_params,double time, int Lbs,int contents1[Lbs][Lbs][Lbs][params.n_box_max],
            int *movers, int moving_pcls[Lbs][Lbs][Lbs][params.n_box_max],
            int non_movers[Lbs][Lbs][Lbs][params.n_box_max], int * pnon_movers[Lbs][Lbs][Lbs],
            long movenum,long *paccptnum,long bondref[][params.Nmono*params.Npoly],const int casesall[13][3],
            double move_accpt[3],Recycle_wavelet *recycled_w,Recycle_Fourier *recycled_F,
            double Ppw, double pwCDF[][2], double R_short[3],
            double dR_acc[3], double dR_acc_reduced[3])
{
    /* Loop counters */
    int ix, iy, iz;

    /* Extract parameters from param_list */
    int L1=params.L, bs1=params.sub_box_size;
    int Ntot1=params.Ntot;

    /* Coordinate objects */
    double a1[3],a2[3],a3[3];
    lattice_vectors(a1,a2,a3,mod_lf(time, params.reset_time),params);
    double xyz_to_C[3][3], C_to_xyz[3][3], metric[3][3];
    Cart_to_xyz_transform(C_to_xyz,a1,a2,a3);
    xyz_to_Cart_transform(xyz_to_C,a1,a2,a3);
    get_metric(metric,a1,a2,a3);

    /* Fourier space coordinate objects */
    double b1[3],b2[3],b3[3],kextent[3];
    double kxyz_to_kC[3][3], kC_to_kxyz[3][3], kmetric[3][3];

    double image_vectors[8][3] = {{0,0,0},{-L1,0,0},{0,-L1,0},{0,0,-L1},{-L1,-L1,0},{0,-L1,-L1},{-L1,0,-L1},{-L1,-L1,-L1}};
    double image_dist[8];
    for(ix=0; ix<8; ix++)
    {image_dist[ix] = find_dr2(image_vectors[ix], metric);}


    /* Box and moving particle trackers */
    int inbox[3]={0,0,0}, imover, nmovers=0, nmovers_box, nnon_movers_box, part_num;
    int remap[Ntot1], iremap=0; /* list of all particles in spans, all of which get re-mapped */
    double part_loc[3], delta[3], length_delta;

    /* Energy change associated with the move */
    double dU = 0.0, F_fwd[3], F_rev[3], amp_weight, sMC_terms;
    /* Move acceptance probability */
    double Pacc, rand_num;

    /* Centre of mass vectors */
    double R_CoM_start[3], R_CoM_end[3], dR_CoM[3], dR_CoM_Cart[3], dR_CoM_det[3];
    calc_R_CoM(R_CoM_start, xyz, params.Nmono, 1);

    /* Choose a wavelet or a plane wave, reusing the previous type if it failed */
    double movetype;
    if(recycled_w->reuse == true && params.recycle == true)
    {movetype = Ppw + 0.1;}
    else if(recycled_F->reuse == true && params.recycle == true)
    {movetype = 0.5*Ppw;}
    else
    {movetype = POW2_32_INV_OPEN*(pcg32_random()+1);}
    /* Note: to get correct probability shift RNG and bring down to the range (0,1] */


    /* Quatities determining which sub-boxes could contain moving particles */
    int anchorxyz[3], spanxyz[3];
    int pwanchor[3]={0,0,0}, pwspans[3]={Lbs-1,Lbs-1,Lbs-1};

    /* structures containing all move parameters, initialised at arbitrary valid values */
    wavelet_parameters wavelet={1.0,{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},{1.0,0.0,0.0},1.0};
    wavelet_parameters *pwavelet = &wavelet;

    Fourier_parameters plane_wave={1.0,{1.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,1.0,0.0},0.0,1.0};
    Fourier_parameters *pplane_wave = &plane_wave;

    /* if pick a wavelet, generate a wavelet and find movers */
    if(movetype>Ppw)
    {
        /* generate new wavelet move */
        get_wavelet_parameters(pwavelet, params.lambda_min, params.lambda_max,
                               params.recycle, recycled_w, Ntot1, xyz, xyz_to_C);

        /* convert wavelet centre to flow coordinates */
        transform_coords(wavelet.centre_Cart, wavelet.centre, C_to_xyz);

        /* record move as a wavelet along with its R */
        move_accpt[0] = 0.5; /* <1, safe from numerical errors */
        move_accpt[1] = wavelet.radius;

        /*-------------------------------*/
        /* identify which particles move */
        /*-------------------------------*/

        /* find boxes in which the wavelet lies */
        get_wavelet_axes(wavelet.extent, wavelet.radius, C_to_xyz);
        anchor(wavelet.centre, wavelet.extent, L1, bs1, anchorxyz);
        spans(wavelet.centre, wavelet.extent, bs1, spanxyz);

        /* find particles inside spans region */
        for(ix=0; ix<=spanxyz[0]; ix++)
        {   inbox[0] = mod_int(anchorxyz[0] + ix,Lbs);

            for(iy=0; iy<=spanxyz[1]; iy++)
            {   inbox[1] = mod_int(anchorxyz[1] + iy,Lbs);

                for(iz=0; iz<=spanxyz[2]; iz++)
                {   inbox[2] = mod_int(anchorxyz[2] + iz,Lbs);

                    imover=0;
                    nmovers_box=0;
                    nnon_movers_box=0;
                    /* loop over contents until get 0 (i.e. no more particles in the box) */
                    while(contents1[inbox[0]][inbox[1]][inbox[2]][imover]>0)
                    {
                        /* check if the particle actually moves */
                        part_num = contents1[inbox[0]][inbox[1]][inbox[2]][imover];

                        delta[0] = xyz[part_num-1][0] - wavelet.centre[0];
                        delta[1] = xyz[part_num-1][1] - wavelet.centre[1];
                        delta[2] = xyz[part_num-1][2] - wavelet.centre[2];

                        /* shift to be the nearest image */
                        params.pshort_wrap(delta, L1, metric, image_vectors, image_dist);

                        /* only move if inside wavelet */
                        length_delta = find_dr2(delta, metric);

                        if(length_delta <= wavelet.radius*wavelet.radius)
                        {
                            /* add mover to the list */
                            movers[nmovers] = part_num;
                            moving_pcls[inbox[0]][inbox[1]][inbox[2]][nmovers_box] = part_num;

                            /* increment counters */
                            nmovers++;
                            nmovers_box++;
                        }
                        else
                        {
                            non_movers[inbox[0]][inbox[1]][inbox[2]][nnon_movers_box] = part_num;
                            nnon_movers_box++;
                        }
                        /* fill list of all local particles */
                        remap[iremap] = part_num;
                        /* increment counters */
                        iremap++;
                        imover++;
                    }

                    /* ensure box lists end with a 0 */
                    moving_pcls[inbox[0]][inbox[1]][inbox[2]][nmovers_box] = 0;
                    non_movers[inbox[0]][inbox[1]][inbox[2]][nnon_movers_box] = 0;
                    pnon_movers[inbox[0]][inbox[1]][inbox[2]] = non_movers[inbox[0]][inbox[1]][inbox[2]];

                }
            }
        }

        /* scale amplitude */
        wavelet.amplitude = params.A0/sqrt(nmovers);
    }
    /* if pick plane wave, generate it and set all particles as movers */
    else
    {
        /* if this is the first Fourier move or the last one was accepted, generate a new wavevector */
        if(recycled_F->reuse == false || params.recycle == false)
        {
            if(params.periodic == false)
            {planewave(pplane_wave, pwCDF);}
            else
            {
                reciprocal_vectors(b1,b2,b3,a1,a2,a3,params.L);
                Cart_to_xyz_transform(kC_to_kxyz,b1,b2,b3);
                xyz_to_Cart_transform(kxyz_to_kC,b1,b2,b3);
                get_metric(kmetric,b1,b2,b3);

                get_wavelet_axes(kextent,(double)params.Lcart, kC_to_kxyz);
                planewave_cart(pplane_wave, 2*M_PI/L1, params.lambda_max, params.Lcart, kextent, kmetric);
            }
        }
        else /* reuse previous wavevector if the last one was rejected */
        {
            pplane_wave->k = recycled_F->k;
            /* generate wave orientation (khat) */
            if(params.periodic==false)
            {direction(pplane_wave->khat);}
            else
            {permute_k(recycled_F->khat, pplane_wave);}
        }

        /* shift wavevector and polarisation to flow coords */
        transform_coords(plane_wave.khat, plane_wave.khat_xyz, C_to_xyz);

        /* generate phase in [0,2pi) */
        pplane_wave->phase = 2*M_PI*POW2_32_INV_OPEN*pcg32_random();


        /* record move as a plane wave along with its k */
        move_accpt[0] = 1.5; /* >1, safe from numerical errors */
        move_accpt[1] = plane_wave.k;

        /* automatically pass move if zero wave vector is chosen.
         * no dynamics need to be calculated for these moves */
        if(plane_wave.k==0.0)
        {*paccptnum = *paccptnum + 1;
         move_accpt[2] = 1.1;
         return;}

        /* set all particles as movers */
        nmovers = Ntot1;
        iremap = nmovers;
        for(ix=1; ix<=Ntot1; ix++)
        {movers[ix-1] = ix;
         remap[ix-1] = ix;}

        for(ix=0; ix<3; ix++)
        {anchorxyz[ix] = pwanchor[ix];
         spanxyz[ix] = pwspans[ix];}

        for(ix=0; ix<Lbs; ix++)
        {for(iy=0; iy<Lbs; iy++)
        {for(iz=0; iz<Lbs; iz++)
         {pnon_movers[ix][iy][iz] = non_movers[ix][iy][iz];}
        }}

    }

    /* reinitialise arrays for movers */
    for(ix=0; ix<nmovers; ix++)
    {
        forces[movers[ix]-1][0] = 0.0;
        forces[movers[ix]-1][1] = 0.0;
        forces[movers[ix]-1][2] = 0.0;
    }


    /** ----------------- **/
    /** Calc force and U0 **/
    /** ----------------- **/

    if(movetype>Ppw)
    {dU = -UF_move(xyz,params,forces,params.smart_MC,metric,anchorxyz,spanxyz,params.EV_fn,params.bond_fn,
                   pot_params,Lbs,moving_pcls,pnon_movers,nmovers,movenum,bondref[0],casesall,
                   image_vectors,image_dist);}
    else
    {dU = -UF_move(xyz,params,forces,params.smart_MC,metric,anchorxyz,spanxyz,params.EV_fn,params.bond_fn,
                   pot_params,Lbs,contents1,pnon_movers,nmovers,movenum,bondref[0],casesall,
                   image_vectors,image_dist);}


    /** ----------------- **/
    /** Get orientations  **/
    /** ----------------- **/

    if(movetype>Ppw)
    {
        if(params.smart_MC == true)
        {
            sum_forces_w(xyz,movers,params,wavelet,forces,F_fwd,xyz_to_C,metric,image_vectors,image_dist);

            /**-------------------------------------------------**/
            /** Temporary section. Used for deterministic <v.v> **/
            /**-------------------------------------------------**/

            /* Find deterministic part of move */
            amp_weight = 1.0/nmovers;
            for(ix=0; ix<3; ix++)
            {wavelet.polarisation[ix] = amp_weight * params.sMC_mob * F_fwd[ix];}

            wavelet.amplitude = sqrt(dot(wavelet.polarisation, wavelet.polarisation));

            for(ix=0; ix<3; ix++)
            {wavelet.polarisation[ix] /= wavelet.amplitude;}

            /* update particle positions xyzn */
            for(ix=0; ix<nmovers; ix++)
            {
                /* get postion vector */
                part_loc[0] = xyz[movers[ix]-1][0];
                part_loc[1] = xyz[movers[ix]-1][1];
                part_loc[2] = xyz[movers[ix]-1][2];

                /* find moved vector */
                movew(part_loc, params, wavelet, xyz_to_C, C_to_xyz, metric, image_vectors, image_dist);

                /* update position vector */
                xyzn[movers[ix]-1][0] = part_loc[0];
                xyzn[movers[ix]-1][1] = part_loc[1];
                xyzn[movers[ix]-1][2] = part_loc[2];
            }

            /* Find deterministic component of CoM displacement */
            calc_R_CoM(R_CoM_end, xyzn, params.Nmono, 1);
            for(ix=0; ix<3; ix++)
            {dR_CoM_det[ix] = R_CoM_end[ix] - R_CoM_start[ix];}


            /**-------------------------------------------------**/
            /**-------------------------------------------------**/
        }
        else /* Don't bias if not using smart MC */
        {
            F_fwd[0] = 0.0;
            F_fwd[1] = 0.0;
            F_fwd[2] = 0.0;
        }

        Gaussian_rand(wavelet.polarisation);

        /* shift to have appropriate mean and variance */
        amp_weight = 1.0/nmovers;
        for(ix=0; ix<3; ix++)
        {
            wavelet.polarisation[ix] *= sqrt(amp_weight * params.sMC_var);
            wavelet.polarisation[ix] += amp_weight * params.sMC_mob * F_fwd[ix];
        }

        /* extract move amplitude and polarisation */
        wavelet.amplitude = sqrt(dot(wavelet.polarisation, wavelet.polarisation));

        for(ix=0; ix<3; ix++)
        {wavelet.polarisation[ix] /= wavelet.amplitude;}

    }
    else
    {
        if(params.smart_MC == true)
        {
            sum_forces_F(xyz, Ntot1, forces, F_fwd, xyz_to_C, plane_wave);

            /**-------------------------------------------------**/
            /** Temporary section. Used for deterministic <v.v> **/
            /**-------------------------------------------------**/

            /* Find deterministic part of move */
            amp_weight = 1.5/((plane_wave.k)*(plane_wave.k)*1.25*Ntot1);
            for(ix=0; ix<3; ix++)
            {plane_wave.polarisation[ix] = amp_weight * params.sMC_mob * F_fwd[ix];}
            /* project the vector to be perpendicular to khat. use rand_num as an unsused double temporarily */
            rand_num = dot(plane_wave.polarisation, plane_wave.khat);
            for(ix=0; ix<3; ix++)
            {plane_wave.polarisation[ix] -= rand_num * plane_wave.khat[ix];}

            /* extract move amplitude and polarisation */
            plane_wave.amplitude = sqrt(dot(plane_wave.polarisation, plane_wave.polarisation));

            for(ix=0; ix<3; ix++)
            {plane_wave.polarisation[ix] /= plane_wave.amplitude;}

            transform_coords(plane_wave.polarisation, plane_wave.polarisation_xyz, C_to_xyz);

            for(ix=0; ix<3; ix++)
            {wavelet.polarisation[ix] /= wavelet.amplitude;}

            /* update particle positions xyzn */
            for(ix=0; ix<nmovers; ix++)
            {
                /* get postion vector */
                part_loc[0] = xyz[movers[ix]-1][0];
                part_loc[1] = xyz[movers[ix]-1][1];
                part_loc[2] = xyz[movers[ix]-1][2];

                /* find moved vector */
                movepw(part_loc, plane_wave, metric);

                /* update position vector */
                xyzn[movers[ix]-1][0] = part_loc[0];
                xyzn[movers[ix]-1][1] = part_loc[1];
                xyzn[movers[ix]-1][2] = part_loc[2];
            }

            /* Find deterministic component of CoM displacement */
            calc_R_CoM(R_CoM_end, xyzn, params.Nmono, 1);
            for(ix=0; ix<3; ix++)
            {dR_CoM_det[ix] = R_CoM_end[ix] - R_CoM_start[ix];}

            /**-------------------------------------------------**/
            /**-------------------------------------------------**/

        }
        else /* Don't bias if not using smart MC */
        {
            F_fwd[0] = 0.0;
            F_fwd[1] = 0.0;
            F_fwd[2] = 0.0;
        }

        Gaussian_rand(plane_wave.polarisation);

        /* shift to have appropriate mean and variance */
        amp_weight = 1.5/((plane_wave.k)*(plane_wave.k)*1.25*Ntot1);
        for(ix=0; ix<3; ix++)
        {
            plane_wave.polarisation[ix] *= sqrt(amp_weight * params.sMC_var);
            plane_wave.polarisation[ix] += amp_weight * params.sMC_mob * F_fwd[ix];
        }

        /* project the vector to be perpendicular to khat. use rand_num as an unsused double temporarily */
        rand_num = dot(plane_wave.polarisation, plane_wave.khat);
        for(ix=0; ix<3; ix++)
        {plane_wave.polarisation[ix] -= rand_num * plane_wave.khat[ix];}

        /* extract move amplitude and polarisation */
        plane_wave.amplitude = sqrt(dot(plane_wave.polarisation, plane_wave.polarisation));
        for(ix=0; ix<3; ix++)
        {plane_wave.polarisation[ix] /= plane_wave.amplitude;}

        transform_coords(plane_wave.polarisation, plane_wave.polarisation_xyz, C_to_xyz);
    }

    /*-----------------------------------------*/
    /* move particles according to the wavelet */
    /* + accept/reject with Metropolis test    */
    /*-----------------------------------------*/

    /* update particle positions xyzn */
    for(ix=0; ix<nmovers; ix++)
    {
        /* get postion vector */
        part_loc[0] = xyz[movers[ix]-1][0];
        part_loc[1] = xyz[movers[ix]-1][1];
        part_loc[2] = xyz[movers[ix]-1][2];

        /* find moved vector */
        if(movetype>Ppw)
        {movew(part_loc, params, wavelet, xyz_to_C, C_to_xyz, metric, image_vectors, image_dist);}
        else
        {movepw(part_loc, plane_wave, metric);}

        /* update position vector */
        xyzn[movers[ix]-1][0] = part_loc[0];
        xyzn[movers[ix]-1][1] = part_loc[1];
        xyzn[movers[ix]-1][2] = part_loc[2];

        forces[movers[ix]-1][0] = 0.0;
        forces[movers[ix]-1][1] = 0.0;
        forces[movers[ix]-1][2] = 0.0;
    }

    if(params.Langevin == false || params.smart_MC == false)
    {
        /* update which box the movers lie in after moving */
        for(ix=0; ix<=spanxyz[0]; ix++)
        {   inbox[0] = mod_int(anchorxyz[0] + ix,Lbs);

            for(iy=0; iy<=spanxyz[1]; iy++)
            {   inbox[1] = mod_int(anchorxyz[1] + iy,Lbs);

                for(iz=0; iz<=spanxyz[2]; iz++)
                {   inbox[2] = mod_int(anchorxyz[2] + iz,Lbs);

                    nmovers_box=0;
                    while(moving_pcls[inbox[0]][inbox[1]][inbox[2]][nmovers_box]>0)
                    {
                        moving_pcls[inbox[0]][inbox[1]][inbox[2]][nmovers_box] = 0;
                        nmovers_box++;
                    }
                }
            }
        }

        map(nmovers, movers, xyzn, L1, bs1, Lbs, params.n_box_max, moving_pcls, anchorxyz, spanxyz);

        /* get change in energy associated with this move */
        dU += UF_move(xyzn,params,forces,params.smart_MC,metric,anchorxyz,spanxyz,params.EV_fn,
                      params.bond_fn,pot_params,Lbs,moving_pcls,pnon_movers,nmovers,movenum,
                      bondref[1],casesall,image_vectors,image_dist);


        if(movetype>Ppw)
        {sum_forces_w(xyzn,movers,params,wavelet,forces,F_rev,xyz_to_C,metric,image_vectors,image_dist);}
        else
        {sum_forces_F(xyzn, Ntot1, forces, F_rev, xyz_to_C, plane_wave);}

        /* contribution to acceptance probability from using smart MC */
        if(params.smart_MC == true)
        {
            sMC_terms = -0.5*params.beta*wavelet.amplitude*(dot(wavelet.polarisation,F_fwd)
                                                             + dot(wavelet.polarisation,F_rev))
                        -0.25*params.beta*params.sMC_mob*amp_weight*(dot(F_rev,F_rev)-dot(F_fwd,F_fwd));
        }
        else
        {sMC_terms = 0.0;}

        /* use dU to calculate the acceptance probability */
        switch(params.use_Glauber)
        {
            case true : Pacc = Glauber_test(-params.beta*dU + sMC_terms); break;
            default : Pacc = Metropolis_test(-params.beta*dU + sMC_terms);
        }

        rand_num = POW2_32_INV_OPEN*pcg32_random();

    }
    else /* guarantee move acceptance if used as a Langevin solver */
    {Pacc = 1.0;
     rand_num = 0.5;}

    if( rand_num < Pacc)
    {
        /* count move as success: (int)1.1=1 without risk from numerical error */
        move_accpt[2] = 1.1;

        /* count towards accepted moves */
        *paccptnum = *paccptnum + 1;

        /**-------------------------------------------------**/
        /** Temporary section. Used for deterministic <v.v> **/
        /**-------------------------------------------------**/

        calc_R_CoM(R_CoM_end, xyzn, params.Nmono, 1);
        /* use true=1 and false=0 to get acceptance test dependent factor */
        if(params.use_Glauber==true)
        {length_delta = 4.0 - (2.0/Pacc);}
        else
        {length_delta = 2.0 * (1.0 - exp(0.5 * params.beta * dU));}

        for(ix=0; ix<3; ix++)
        {
            dR_CoM[ix] = R_CoM_end[ix] - R_CoM_start[ix];
            transform_coords(dR_CoM, dR_CoM_Cart, xyz_to_C);

            R_short[ix] += dR_CoM[ix];
            if(params.smart_MC == true)
            {
                dR_acc[ix] += 2.0*dR_CoM_det[ix];
            }
            else
            {
                dR_acc[ix] += dR_CoM[ix]*8.0*(Pacc - 0.5);
                dR_acc_reduced[ix] += dR_CoM[ix]*(2*Pacc - 1);
            }

        }
        /**-------------------------------------------------**/
        /**-------------------------------------------------**/

        /* ensure the move parameters are not recycled */
        if(movetype>Ppw)
        {recycled_w->reuse = false;}
        else
        {recycled_F->reuse = false;}

        /* empty spans region of contents */
        for(ix=0; ix<=spanxyz[0]; ix++)
        {   inbox[0] = mod_int(anchorxyz[0] + ix,Lbs);

            for(iy=0; iy<=spanxyz[1]; iy++)
            {   inbox[1] = mod_int(anchorxyz[1] + iy,Lbs);

                for(iz=0; iz<=spanxyz[2]; iz++)
                {
                    inbox[2] = mod_int(anchorxyz[2] + iz,Lbs);

                    imover=0;
                    while(contents1[inbox[0]][inbox[1]][inbox[2]][imover]>0)
                    {
                        contents1[inbox[0]][inbox[1]][inbox[2]][imover]=0;
                        imover++;
                    }
                }
            }
        }

        /* remap movers */
        map(iremap, remap, xyzn, L1, bs1, Lbs, params.n_box_max, contents1, anchorxyz, spanxyz);


        /* update movers' positions in xyz */
        for(ix=0; ix<nmovers; ix++)
        {
            xyz[movers[ix]-1][0] = xyzn[movers[ix]-1][0];
            xyz[movers[ix]-1][1] = xyzn[movers[ix]-1][1];
            xyz[movers[ix]-1][2] = xyzn[movers[ix]-1][2];
        }

    }
    else /* rejected move */
    {
        /* count move as fail: (int)0.1=0 without risk from numerical error */
        move_accpt[2] = 0.1;

        /**-------------------------------------------------**/
        /** Temporary section. Used for deterministic <v.v> **/
        /**-------------------------------------------------**/
        calc_R_CoM(R_CoM_end, xyzn, params.Nmono, 1);

        for(ix=0; ix<3; ix++)
        {
            dR_CoM[ix] = R_CoM_end[ix] - R_CoM_start[ix];
            transform_coords(dR_CoM, dR_CoM_Cart, xyz_to_C);

            R_short[ix] += dR_CoM_Cart[ix];
        }
        /**-------------------------------------------------**/
        /**-------------------------------------------------**/

        /* mark parameters for reusing if move recycling is on*/
        if(params.recycle == true)
        {
            if(movetype>Ppw)
            {
                recycled_w->reuse = true;
                recycled_w->radius = wavelet.radius;

                recycled_w->centre[0] = wavelet.centre[0];
                recycled_w->centre[1] = wavelet.centre[1];
                recycled_w->centre[2] = wavelet.centre[2];
            }
            else
            {
                recycled_F->reuse = true;
                recycled_F->k = plane_wave.k;

                recycled_F->khat[0] = plane_wave.khat[0];
                recycled_F->khat[1] = plane_wave.khat[1];
                recycled_F->khat[2] = plane_wave.khat[2];
            }
        }


        /* return particles to their old locations */
        for(ix=0; ix<nmovers; ix++)
        {
            xyzn[movers[ix]-1][0] = xyz[movers[ix]-1][0];
            xyzn[movers[ix]-1][1] = xyz[movers[ix]-1][1];
            xyzn[movers[ix]-1][2] = xyz[movers[ix]-1][2];
        }
    }

    /* reinitialise arrays for movers */
    for(ix=0; ix<nmovers; ix++)
    {movers[ix] = 0;}

    for(ix=0; ix<=spanxyz[0]; ix++)
    {   inbox[0] = mod_int(anchorxyz[0] + ix,Lbs);

        for(iy=0; iy<=spanxyz[1]; iy++)
        {   inbox[1] = mod_int(anchorxyz[1] + iy,Lbs);

            for(iz=0; iz<=spanxyz[2]; iz++)
            {   inbox[2] = mod_int(anchorxyz[2] + iz,Lbs);

                pnon_movers[inbox[0]][inbox[1]][inbox[2]] = contents1[inbox[0]][inbox[1]][inbox[2]];
                nmovers_box=0;
                while(moving_pcls[inbox[0]][inbox[1]][inbox[2]][nmovers_box]>0 ||
                      non_movers[inbox[0]][inbox[1]][inbox[2]][nmovers_box]>0)
                {
                    moving_pcls[inbox[0]][inbox[1]][inbox[2]][nmovers_box] = 0;
                    non_movers[inbox[0]][inbox[1]][inbox[2]][nmovers_box] = 0;
                    nmovers_box++;
                }
            }
        }
    }

}


/* generate wavelet parameters, without requiring it to move any particles */
void get_wavelet_parameters(wavelet_parameters *pwavelet,double Rmin,double Rmax,
             bool recycle, Recycle_wavelet *recycled_w, int Ntot, double xyz[][3], double transform[3][3])
{
    double probR;
    double rCart[3];

    /* generate particle array index number in [0,Ntot-1] */
    int part_num = (int)pcg32_boundedrand((uint64_t)Ntot);

    if(recycled_w->reuse == false || recycle == false) /* generate a new radius and wavelet centre */
    {
        /* radius of wavelet, with prob density ~ R^-4 */
        probR = POW2_32_INV*pcg32_random();
        pwavelet->radius = Rmin/cbrt(probR + (1-probR)*pow(Rmin/Rmax,3));

    }
    else /* reuse parameters of previously rejected wavelet */
    {pwavelet->radius = recycled_w->radius;}

    /* get direction and scaled distance of wavelet centre from chosen particle */
    pwavelet->centre[0] = 1.0;
    pwavelet->centre[1] = 1.0;
    pwavelet->centre[2] = 1.0;

    /* generate random vector within init sphere */
    while(dot( pwavelet->centre , pwavelet->centre )>1)
    {
        pwavelet->centre[0] = 2*(POW2_32_INV*pcg32_random()) - 1;
        pwavelet->centre[1] = 2*(POW2_32_INV*pcg32_random()) - 1;
        pwavelet->centre[2] = 2*(POW2_32_INV*pcg32_random()) - 1;
    }

    transform_coords(xyz[part_num],rCart,transform);

    pwavelet->centre_Cart[0] = rCart[0] + (pwavelet->radius)*(pwavelet->centre[0]);
    pwavelet->centre_Cart[1] = rCart[1] + (pwavelet->radius)*(pwavelet->centre[1]);
    pwavelet->centre_Cart[2] = rCart[2] + (pwavelet->radius)*(pwavelet->centre[2]);

    /* orientation of the wavelet */
    direction(pwavelet->polarisation);

}


void planewave(Fourier_parameters *pplane_wave, double CDF[][2])
{
    int i=0;
    double randnum;

    /* k=0 should be impossible so shift RNG and bring down to the range (0,1] */
    randnum = POW2_32_INV_OPEN*(pcg32_random()+1);
    /* find value in CDF array above randnum */
    while(randnum>CDF[i][0])
    {i++;}
    /* linearly interpolate between the two neighbouring values */
    pplane_wave->k = CDF[i-1][1] + (CDF[i][1]-CDF[i-1][1])*(randnum - CDF[i-1][0])/(CDF[i][0]-CDF[i-1][0]);

    /* generate wave orientation (khat) */
    direction(pplane_wave->khat);
}


void planewave_cart(Fourier_parameters *pplane_wave, double k0, double Rmax, int Lcart,
                    double kextent[3], double kmetric[3][3])
{
    int accpt=0;
    int lmax=(int)ceil(kextent[0]), mmax=(int)ceil(kextent[1]), nmax=(int)ceil(kextent[2]);
    double lmn[3];

    while(accpt==0)
    {
        /* generate uniformly distributed integer vectors inside cube
         * of side length 2Lcart+1 centred on 0 */
        lmn[0] = (double)pcg32_boundedrand( (uint64_t)(2*lmax + 1) ) - lmax;
        lmn[1] = (double)pcg32_boundedrand( (uint64_t)(2*mmax + 1) ) - mmax;
        lmn[2] = (double)pcg32_boundedrand( (uint64_t)(2*nmax + 1) ) - nmax;

        pplane_wave->k = sqrt(find_dr2(lmn,kmetric));

        /* accept choice of k with appropriate probability. Otherwise pick another value.
         * Note that pwProb(0.0) is used instead of the sum over the whole Cartesian region
         * since only relative prob.s matter. By using the largest single value, the acceptance
         * rate is maximised without altering the distribution.
         */
        if( pplane_wave->k < k0*Lcart &&
           (POW2_32_INV_OPEN*pcg32_random()) < pwProb((pplane_wave->k) * Rmax)/pwProb(0.0) )
        {accpt = 1;}

    }

    /* automatically pass k=0 mode
     * (this will have no effect on dynamics so the other detials are irrelevant) */
    if(pplane_wave->k == 0.0)
    {return;}

    /* extract unit vector */
    pplane_wave->khat[0] = k0*lmn[0]/(pplane_wave->k);
    pplane_wave->khat[1] = k0*lmn[1]/(pplane_wave->k);
    pplane_wave->khat[2] = k0*lmn[2]/(pplane_wave->k);

}


/* move particles according to a wavelet */
void movew(double *position, System_parameters params, wavelet_parameters wavelet, double xyz_to_C[3][3],
           double C_to_xyz[3][3], double metric[3][3], double images[8][3],double image_dist[8])
{
    double delta_xyz[3], delta[3], length_delta;
    double dpos[3], dpos_xyz[3];
    double omega, deltat[3];
    double *cross_part;

    /* find absolute vector between wrapped particle and wavelet centre */
    delta_xyz[0] = position[0] - wavelet.centre[0];
    delta_xyz[1] = position[1] - wavelet.centre[1];
    delta_xyz[2] = position[2] - wavelet.centre[2];

    params.pshort_wrap(delta_xyz,params.L,metric,images,image_dist);

    transform_coords(delta_xyz, delta, xyz_to_C);
    /* scalar distance from centre of wavelet */
    length_delta = sqrt(dot(delta,delta));

    /* rotation angle(/velocity) at particle */
    omega = (wavelet.amplitude)*(1.0 - length_delta/(wavelet.radius) );

    /* component of delta transverse to rotation axis */
    deltat[0] = delta[0] - (wavelet.polarisation[0])*dot(wavelet.polarisation ,delta);
    deltat[1] = delta[1] - (wavelet.polarisation[1])*dot(wavelet.polarisation ,delta);
    deltat[2] = delta[2] - (wavelet.polarisation[2])*dot(wavelet.polarisation ,delta);

    /* component transverse to both rotation axis and deltat */
    cross_part = cross(wavelet.polarisation ,deltat);

    /* update position */
    dpos[0] = deltat[0]*(cos(omega)-1) + cross_part[0]*sin(omega);
    dpos[1] = deltat[1]*(cos(omega)-1) + cross_part[1]*sin(omega);
    dpos[2] = deltat[2]*(cos(omega)-1) + cross_part[2]*sin(omega);

    transform_coords(dpos, dpos_xyz, C_to_xyz);

    position[0] += dpos_xyz[0];
    position[1] += dpos_xyz[1];
    position[2] += dpos_xyz[2];

}


/* move particles according to a plane wave */
void movepw(double *position, Fourier_parameters plane_wave, double metric[3][3])
{
    double delta = (plane_wave.amplitude)*cos( (plane_wave.k) * find_drdR(plane_wave.khat_xyz, position, metric)
                                               + plane_wave.phase);

    position[0] = position[0] + delta*(plane_wave.polarisation_xyz[0]);
    position[1] = position[1] + delta*(plane_wave.polarisation_xyz[1]);
    position[2] = position[2] + delta*(plane_wave.polarisation_xyz[2]);
}


/* calculate change in potential energy after a move */
double UF_move(double xyz[][3],System_parameters params,double forces[][3],bool find_forces,
               double metric[3][3],int *anchorxyz,int *spanxyz,
               double (*EV_pot)(double, Potentials, double *), double (*bond_pot)(double, Potentials, double *),
               Potentials pot_params,int Lbs,
               int movers[Lbs][Lbs][Lbs][params.n_box_max],int *non_movers[Lbs][Lbs][Lbs],
               int nmovers,long movenum,long *bondref,const int casesall[13][3],
               double images[8][3],double image_dist[8])
{
    double Utot=0.0/* total interaction energy */,F /* scalar force */,separation[3],sep_squared;
    double *pF = &F;
    int inbox[3], neibox[3]; /* box coords */
    int ix, iy, iz, idir, iin, inei; /* loop indices */
    int inpart, neipart; /* particle numbers */

    int nboxmax = params.n_box_max, Nmono2 = params.Nmono;

    /* loop over region containing the wavelet */
    for(ix=0; ix<=spanxyz[0]; ix++)
    {   inbox[0] = mod_int(anchorxyz[0] + ix,Lbs);

        for(iy=0; iy<=spanxyz[1]; iy++)
        {   inbox[1] = mod_int(anchorxyz[1] + iy,Lbs);

            for(iz=0; iz<=spanxyz[2]; iz++)
            {   inbox[2] = mod_int(anchorxyz[2] + iz,Lbs);

                /* don't waste time on this sub-box if it doesn't contain any moving particles */
                if(movers[inbox[0]][inbox[1]][inbox[2]][0]==0)
                {continue;}

/*-------------------------------------------------------------------------------*/
                iin=0;

                 /* loop over pairs of particles in the 'in' box */
                while((movers[inbox[0]][inbox[1]][inbox[2]][iin] > 0) && (iin < nboxmax)/* don't try any elements that shouldn't exist */)
                {
                    /* get particle number in box, shifted by -1 for C index to start at 0 */
                    inpart = movers[inbox[0]][inbox[1]][inbox[2]][iin] - 1;

                    /*-------------------*/
                    /* bond potentials   */
                    /*-------------------*/
                    /* Force with preceding monomer */
                    if((mod_int(inpart,Nmono2) != 0) /* skip if first monomer in the chain */
                       && (bondref[inpart-1]<movenum)/* skip if this bond counted already this move*/)
                    {
                        /* get bond vector */
                        separation[0] = xyz[inpart][0] - xyz[inpart-1][0];
                        separation[1] = xyz[inpart][1] - xyz[inpart-1][1];
                        separation[2] = xyz[inpart][2] - xyz[inpart-1][2];
                        sep_squared = find_dr2(separation,metric);

                        /* get PE */
                        Utot += (*bond_pot)(sep_squared, pot_params, pF);

                        if(find_forces==true)
                        {
                            forces[inpart][0] += F * separation[0];
                            forces[inpart][1] += F * separation[1];
                            forces[inpart][2] += F * separation[2];

                            forces[inpart-1][0] -= F * separation[0];
                            forces[inpart-1][1] -= F * separation[1];
                            forces[inpart-1][2] -= F * separation[2];
                        }

                    }

                    /* Force with next monomer */
                    if((mod_int(inpart+1,Nmono2) != 0) /* skip if last monomer in the chain */
                       && (bondref[inpart+1]<movenum))
                    {
                        /* get old bond vector */
                        separation[0] = xyz[inpart][0] - xyz[inpart+1][0];
                        separation[1] = xyz[inpart][1] - xyz[inpart+1][1];
                        separation[2] = xyz[inpart][2] - xyz[inpart+1][2];
                        sep_squared = find_dr2(separation,metric);

                        /* subtract old PE */
                        Utot += (*bond_pot)(sep_squared, pot_params, pF);

                        if(find_forces==true)
                        {
                            forces[inpart][0] += F * separation[0];
                            forces[inpart][1] += F * separation[1];
                            forces[inpart][2] += F * separation[2];

                            forces[inpart+1][0] -= F * separation[0];
                            forces[inpart+1][1] -= F * separation[1];
                            forces[inpart+1][2] -= F * separation[2];
                        }

                    }

                    /* record that bonds linked to ipart have been accounted for this move */
                    bondref[inpart] = movenum;

                    if(strcmp(pot_params.pot_type,"free")==0)
                    {iin++; continue;}

                    /*-------------------------*/
                    /* inter-bead potentials   */
                    /*-------------------------*/

                    inei = iin + 1;

                    /** interactions with movers in the box */
                    while((movers[inbox[0]][inbox[1]][inbox[2]][inei] > 0) && (inei < nboxmax))
                    {
                        /* get particle number */
                        neipart = movers[inbox[0]][inbox[1]][inbox[2]][inei] - 1;

                        /* get old separation vector */
                        separation[0] = xyz[inpart][0] - xyz[neipart][0];
                        separation[1] = xyz[inpart][1] - xyz[neipart][1];
                        separation[2] = xyz[inpart][2] - xyz[neipart][2];

                        params.pshort_wrap(separation,params.L,metric,images,image_dist);
                        sep_squared = find_dr2(separation,metric);

                        /* add to PE */
                        Utot += (*EV_pot)(sep_squared, pot_params, pF);

                        if(find_forces==true)
                        {
                            forces[inpart][0] += F * separation[0];
                            forces[inpart][1] += F * separation[1];
                            forces[inpart][2] += F * separation[2];

                            forces[neipart][0] -= F * separation[0];
                            forces[neipart][1] -= F * separation[1];
                            forces[neipart][2] -= F * separation[2];
                        }

                        /* increment counter for neighbouring box */
                        inei++;
                    }

                    inei=0;

                    /** interactions with non-movers */
                    while(( non_movers[inbox[0]][inbox[1]][inbox[2]][inei] > 0) && (inei < nboxmax))
                    {
                        /* get particle number in nei box */
                        neipart = non_movers[inbox[0]][inbox[1]][inbox[2]][inei] - 1;

                        /* get old separation vector */
                        separation[0] = xyz[inpart][0] - xyz[neipart][0];
                        separation[1] = xyz[inpart][1] - xyz[neipart][1];
                        separation[2] = xyz[inpart][2] - xyz[neipart][2];

                        params.pshort_wrap(separation,params.L,metric,images,image_dist);
                        sep_squared = find_dr2(separation,metric);

                        /* add to PE */
                        Utot += (*EV_pot)(sep_squared, pot_params, pF);

                        if(find_forces==true)
                        {
                            forces[inpart][0] += F * separation[0];
                            forces[inpart][1] += F * separation[1];
                            forces[inpart][2] += F * separation[2];

                            /* Note there is no need to get the force on the non-mover */
                        }

                        /* increment counter for neighbouring box */
                        inei++;
                    }

/*-------------------------------------------------------------------------------*/

                    /* interactions in 'forwards' directions: movers -> movers and non_movers */
                    for(idir=0; idir<13; idir++)
                    {
                        /* get coords of neighbouring box */
                        neibox[0] = mod_int(inbox[0] + casesall[idir][0],Lbs);
                        neibox[1] = mod_int(inbox[1] + casesall[idir][1],Lbs);
                        neibox[2] = mod_int(inbox[2] + casesall[idir][2],Lbs);

                        inei=0;

                        /** interactions with movers in the neighbour box */
                        while((movers[neibox[0]][neibox[1]][neibox[2]][inei] > 0) && (inei < nboxmax))
                        {
                            /* get particle number in nei box */
                            neipart = movers[neibox[0]][neibox[1]][neibox[2]][inei] - 1;

                            /* get old separation vector */
                            separation[0] = xyz[inpart][0] - xyz[neipart][0];
                            separation[1] = xyz[inpart][1] - xyz[neipart][1];
                            separation[2] = xyz[inpart][2] - xyz[neipart][2];

                            params.pshort_wrap(separation,params.L,metric,images,image_dist);
                            sep_squared = find_dr2(separation,metric);

                            /* add to PE */
                            Utot += (*EV_pot)(sep_squared, pot_params, pF);

                            if(find_forces==true)
                            {
                                forces[inpart][0] += F * separation[0];
                                forces[inpart][1] += F * separation[1];
                                forces[inpart][2] += F * separation[2];

                                forces[neipart][0] -= F * separation[0];
                                forces[neipart][1] -= F * separation[1];
                                forces[neipart][2] -= F * separation[2];
                            }

                            /* increment counter for neighbouring box */
                            inei++;
                        }

                        /* avoid double counting for moves which contain all sub-boxes */
                        if(nmovers == params.Ntot)
                        {continue;}

                        inei=0;

                        /** interactions with non-movers */
                        while((non_movers[neibox[0]][neibox[1]][neibox[2]][inei] > 0) && (inei < nboxmax))
                        {
                            /* get particle number in nei box */
                            neipart = non_movers[neibox[0]][neibox[1]][neibox[2]][inei] - 1;

                            /* get old separation vector */
                            separation[0] = xyz[inpart][0] - xyz[neipart][0];
                            separation[1] = xyz[inpart][1] - xyz[neipart][1];
                            separation[2] = xyz[inpart][2] - xyz[neipart][2];

                            params.pshort_wrap(separation,params.L,metric,images,image_dist);
                            sep_squared = find_dr2(separation,metric);

                            /* add to PE */
                            Utot += (*EV_pot)(sep_squared, pot_params, pF);

                            if(find_forces==true)
                            {
                                forces[inpart][0] += F * separation[0];
                                forces[inpart][1] += F * separation[1];
                                forces[inpart][2] += F * separation[2];

                                /* Note there is no need to get the force on the non-mover */
                            }

                            /* increment counter for neighbouring box */
                            inei++;
                        }

    /*---------------------------------------------------------------------------------*/

                        neibox[0] = mod_int(inbox[0] - casesall[idir][0],Lbs);
                        neibox[1] = mod_int(inbox[1] - casesall[idir][1],Lbs);
                        neibox[2] = mod_int(inbox[2] - casesall[idir][2],Lbs);

                        inei=0;

                        /** interactions with non-movers in negative directions */
                        while((non_movers[neibox[0]][neibox[1]][neibox[2]][inei] > 0) && (inei < nboxmax))
                        {
                            /* get particle number in nei box */
                            neipart = non_movers[neibox[0]][neibox[1]][neibox[2]][inei] - 1;

                            /* get old separation vector */
                            separation[0] = xyz[inpart][0] - xyz[neipart][0];
                            separation[1] = xyz[inpart][1] - xyz[neipart][1];
                            separation[2] = xyz[inpart][2] - xyz[neipart][2];

                            params.pshort_wrap(separation,params.L,metric,images,image_dist);
                            sep_squared = find_dr2(separation,metric);

                            /* add to PE */
                            Utot += (*EV_pot)(sep_squared, pot_params, pF);

                            if(find_forces==true)
                            {
                                forces[inpart][0] += F * separation[0];
                                forces[inpart][1] += F * separation[1];
                                forces[inpart][2] += F * separation[2];

                                /* Note there is no need to get the force on the non-mover */
                            }

                            /* increment counter for neighbouring box */
                            inei++;
                        }

                    }


                    /* increment counter for current box */
                    iin++;
                }

            /* closing braces for loop over boxes */
            }
        }
    }

    return Utot;
}


