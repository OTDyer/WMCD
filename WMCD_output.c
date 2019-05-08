/*
 * Functions for outputting information from WMCD.
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
//#include "pcg_basic.h"

//extern pcg32_random_t pcg32_global;


/* construct directory and system dependent parts of file names */
void get_directory_and_suffix(char *directory, char *suffix, System_parameters params,
                              Potentials pot_params, int batch, int run_ID, bool periodic)
{
    sprintf(directory,"./");
    // modify this as appropriate

    sprintf(suffix,
            "%s_%d_%d_%s_%g_%s_%d_%d.txt",
            (periodic==true)?"_P":"",
            params.Nmono, params.Npoly,
            pot_params.solvent, params.A0,
            (params.smart_MC == false)?(/*f1*/ (params.use_Glauber == true)?"G":"M"):
                    (/*t1*/(params.Langevin == true)?"L":
                     (/*t2*/(params.use_Glauber == true)?"sG":"sM")),
            batch, run_ID);
}


/* function that writes the list of parameters that actually affect the system */
void write_parameters(System_parameters params, Potentials pot_params, int seed[2],
                      char *directory, char *suffix, time_hms run_time_hms, double Pacc, bool start)
{
    char file_name[150];
    sprintf(file_name, "%sSystem_parameters/system%s", directory, suffix);

    FILE *system_file = fopen(file_name,"a");

    /* check file opened */
    if((system_file!=NULL))
    {
        if(start == true)
        {
            /* print input parameters */
            fprintf(system_file,
                 "L = %d\nNmono = %d\nNpoly = %d\n\
                 \nlambda_min = %g\nlambda_max = %g\nlambda_a = %g\
                 \nkT = %g\nA_0 = %g\nFourier moves are %s\
                 \n%s test is used\
                 \nmove recycling is %s\n\
                 \nPot type = %s\n\tlength scale = %g\n\tenergy scale = %g\n\
                 \nBond type = %s\n\tlength scale = %g\n\tspring const = %g\n\
                 \nshear rate = %g\n\
                 \nseed = %d\nseed sequence = %d\n",
                 params.L, params.Nmono, params.Npoly,
                 params.lambda_min, params.lambda_max, params.lambda_a,
                 params.kT, params.A0,
                 (params.use_Fourier == true)?"on":"off",
                 (params.smart_MC == false)?(/*f1*/ (params.use_Glauber == true)?"Glauber":"Metropolis"):
                    (/*t1*/(params.Langevin == true)?"Langevin":
                     (/*t2*/(params.use_Glauber == true)?"smart Glauber":"smart Metropolis")),
                 (params.recycle == true)?"on":"off",
                 pot_params.pot_type, pot_params.pot_length, pot_params.pot_energy,
                 pot_params.bond_type, pot_params.bond_length, pot_params.spring_const,
                 params.shear_rate,
                 seed[0], seed[1]);
        }


        if(start == false)
        {
            /* print run-time specific details */
            fprintf(system_file,
                 "\n-----------------------------------\n\
                 \nTotal runtime for %g %s: %dh %dm %.3lfs\
                 \nGlobal fraction of moves accepted: %.8lf\n",
                 (params.set_time_max==false)?params.movenum_max:params.time_max,
                 (params.set_time_max==false)?"moves":"units of time",
                 run_time_hms.hours, run_time_hms.minutes, run_time_hms.seconds,
                 Pacc);
        }


        /* close files */
        fclose(system_file);
    }
    else /* otherwise print an error */
    {printf("Error: system parameters file failed to open\n");}
}


/* write equilibrated configuration to file */
void write_equilib_config(double xyz[][3], int Ntot, double xyz_to_C[3][3], char *directory, char *suffix)
{
    int i;
    double Cart_coords[3];

    /* build file name from relevant system parameters */
    char file_name[150];
    sprintf(file_name, "%sEquilib_configs/config%s", directory, suffix);

    FILE *config_file = fopen(file_name,"w");

    /* check file opened */
    if((config_file!=NULL))
    {
        for(i=0; i<Ntot; i++)
        {
            transform_coords(xyz[i], Cart_coords, xyz_to_C);
            fprintf(config_file,"%lf\t%lf\t%lf\n",Cart_coords[0],Cart_coords[1],Cart_coords[2]);
        }

        /* close files */
        fclose(config_file);
    }
    else /* otherwise print an error */
    {printf("Error: file for equlibrated configuration failed to open\n");}

}

/* write configuration to file to continue from in a future run */
void write_continue_config(double xyz[][3], int Ntot, double time, double xyz_to_C[3][3], char *directory, char *suffix)
{
    int i;
    double Cart_coords[3];

    /* build file name from relevant system parameters */
    char file_name[150];
    sprintf(file_name, "%sContinue_configs/config%s", directory, suffix);

    FILE *config_file = fopen(file_name,"w");

    /* check file opened */
    if((config_file!=NULL))
    {
        /* print time */
        fprintf(config_file, "%lf\n", time);

        for(i=0; i<Ntot; i++)
        {
            transform_coords(xyz[i], Cart_coords, xyz_to_C);
            fprintf(config_file,"%lf\t%lf\t%lf\n",Cart_coords[0],Cart_coords[1],Cart_coords[2]);
        }

        /* close files */
        fclose(config_file);
    }
    else /* otherwise print an error */
    {printf("Error: file for equlibrated configuration failed to open\n");}

}


/* routine ouputting particle coordinates */
void unwrapped_coords_file(double xyz[][3], int Ntot, double time, double xyz_to_C[3][3], char *directory, char *suffix)
{
    FILE *extended_file; /* file with unwrapped coordinates */

    int i; /* loop indices */
    double r_Cart[3];

    /* declare strings for filenames */
    char ext_str[150];

    sprintf(ext_str,"%sEvolving_configs/unwrapped%s", directory, suffix);

    extended_file = fopen(ext_str,"a");

    /* check file opened */
    if((extended_file!=NULL))
    {
        /* write the time */
        //fprintf(extended_file,"%lf\n",time);
        fprintf(extended_file,"%lf\n",time);

        /* write coords to file */
        for(i=0; i<Ntot; i++)
        {
            transform_coords(xyz[i], r_Cart, xyz_to_C);
            fprintf(extended_file,"%lf\t%lf\t%lf\n",r_Cart[0],r_Cart[1],r_Cart[2]);
        }

        /* close files */
        fclose(extended_file);
    }
    else /* otherwise print an error */
    {printf("Error: unwrapped coords file failed to open\n");}

}


/* routine ouputting particle coordinates after being wrapped back into the periodic box */
void wrapped_coords_file(double xyz[][3], int Ntot, int L, double time, double xyz_to_C[3][3],
                         char *directory, char *suffix)
{
    FILE *wrapped_file; /* file with wrapped coordinates */

    int i; /* loop indices */
    double *wrp, coords[3], Cart_coords[3]; /* temporary storage of wrapped coords */

    /* declare strings for filenames */
    char wrp_str[150];

    sprintf(wrp_str,"%sEvolving_configs/wrapped%s", directory, suffix);

    wrapped_file = fopen(wrp_str,"a");

    /* check file opened */
    if((wrapped_file!=NULL))
    {
        /* write the time in the simulation */
        fprintf(wrapped_file,"%lf\n",time);

        for(i=0; i<Ntot; i++)
        {
            coords[0] = xyz[i][0];
            coords[1] = xyz[i][1];
            coords[2] = xyz[i][2];

            /* wrap the coords in the shearing frame */
            wrp = wrap(coords,L);

            transform_coords(wrp, Cart_coords, xyz_to_C);

            /* wrap again in the Cartesian frame */
            wrp = wrap(Cart_coords,L);

            /* write coords to file */
            fprintf(wrapped_file,"%lf\t%lf\t%lf\n", wrp[0],wrp[1],wrp[2]);
        }

        /* close file */
        fclose(wrapped_file);
    }
    else /* otherwise print an error */
    {printf("Error: wrapped coords file failed to open\n");}

}

/* output the forces on each particle at time */
void write_forces(double forces[][3], int Ntot, double time, char *directory, char *suffix)
{
    FILE *force_file; /* file with unwrapped coordinates */

    int i; /* loop indices */

    /* declare strings for filenames */
    char f_str[150];

    sprintf(f_str,"%sEvolving_configs/forces%s", directory, suffix);

    force_file = fopen(f_str,"a");

    /* check file opened */
    if((force_file!=NULL))
    {
        /* write the time */
        //fprintf(extended_file,"%lf\n",time);
        fprintf(force_file,"%lf\n",time);

        /* write coords to file */
        for(i=0; i<Ntot; i++)
        {
            fprintf(force_file,"%lf\t%lf\t%lf\n",forces[i][0],forces[i][1],forces[i][2]);
        }

        /* close files */
        fclose(force_file);
    }
    else /* otherwise print an error */
    {printf("Error: forces file failed to open\n");}

}

/* output mean-square end-to-end and gyration radii at chosen evolution time */
void radii_file(double xyz[][3], int Npoly, int Nmono, double time,  double metric[3][3],
                char *directory, char *suffix)
{
    int i,j,k;
    double vector[3];

    /* declare mean-square end-to-end vector over polymers */
    double Re2e2_mean=0.0;
    /* declare mean-square radius of gyration over polymers */
    double Rgy2_mean=0.0;

    FILE *rad_file;

    /* declare strings for filenames */
    char file_str[150];
    sprintf(file_str, "%sradius_files/radii%s", directory, suffix);

    rad_file = fopen(file_str,"a");

    /* check file opened properly */
    if(rad_file!=NULL)
    {
        /* print file column information if at the start of a new file */
        if(ftell(rad_file) == 0)
        {fprintf(rad_file, "time\t\tR_g^2\t\tR_e^2\t\tX_e\t\tY_e\t\tZ_e \n");}


        /** Find mean square radius of gyration */
        for(i=0; i<Npoly; i++)
        {
            /* sum over upper (or lower, they're equivalent here) half triangle of (rj-rk)^2 */
            for(j=0; j<Nmono; j++)
            {
                for(k=j+1; k<Nmono; k++)
                {
                    vector[0] = xyz[i*Nmono + j][0] - xyz[i*Nmono + k][0];
                    vector[1] = xyz[i*Nmono + j][1] - xyz[i*Nmono + k][1];
                    vector[2] = xyz[i*Nmono + j][2] - xyz[i*Nmono + k][2];

                    Rgy2_mean += find_dr2(vector, metric);
                }
            }
        }
        /* average over polymers. Nmono is squared as each of the j and k sums contribute one factor */
        Rgy2_mean /= (Npoly*Nmono*Nmono);

        /** Find mean square end-to-end separation */
        for(i=0; i<Npoly; i++)
        {
            /* calculate end-to-end vectors */
            vector[0] = xyz[(i+1)*Nmono - 1][0] - xyz[i*Nmono][0];
            vector[1] = xyz[(i+1)*Nmono - 1][1] - xyz[i*Nmono][1];
            vector[2] = xyz[(i+1)*Nmono - 1][2] - xyz[i*Nmono][2];

            /* sum their square to the mean */
            Re2e2_mean += find_dr2(vector, metric);
        }
        /* average over the polymers */
        Re2e2_mean /= Npoly;

        /* print radii to file */
        fprintf(rad_file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf \n",
                time, Rgy2_mean, Re2e2_mean, vector[0], vector[1], vector[2]);

        fclose(rad_file);
    }
    else
    {printf("Error: rad_file did not open\n");}
}

/** This function name wants changing **/
/* calculate CoM displacement squared at logarithmically increassing time lags */
void calc_dR_CoM_squared(cumulative_data R_CoM_squared[], prev_quantities prev_R_CoM[][3],
                         double R_CoM[3], double time, int output_num, int pow_max)
{
    int i,j;
    bool test = true;
    int power = 0, index;

    double dR[3];

    /** for single time lag data **/
    /* only start when previous data exists */
    if(prev_R_CoM[power][0].overlap_test != 0)
    {
        /* find the CoM displacement */
        for(i=0; i<3; i++)
        {dR[i] = R_CoM[i] - prev_R_CoM[0][2].r[i];}

        /* add temporal displacement */
        R_CoM_squared[0].sum_dt = R_CoM_squared[0].sum_dt + time - prev_R_CoM[0][2].time;

        /* add spacial displacement squared */
        R_CoM_squared[0].sum_dR = R_CoM_squared[0].sum_dR + dot(dR, dR);

        /* add to count of contributions */
        R_CoM_squared[0].num_contributions = R_CoM_squared[0].num_contributions + 1;
    }

    /** for multiple time lag data **/
    while(test==true && power<=pow_max)
    {
        if( mod_long(output_num, nint(pow(2,power))) == 0)
        {
            /* only start when previous data exists */
            if(prev_R_CoM[power][0].overlap_test != 0)
            {
                for(j=0; j<=1; j++)
                {
                    index = 2*(power + 1) - j;
                    /* avoid overlapping lags */
                    if(mod_long(prev_R_CoM[power][0].overlap_test,3-j) == 0)
                    {
                         /* find the CoM displacement */
                        for(i=0; i<3; i++)
                        {dR[i] = R_CoM[i] - prev_R_CoM[power][j].r[i];}

                        /* add temporal displacement */
                        R_CoM_squared[index].sum_dt =
                            R_CoM_squared[index].sum_dt + time - prev_R_CoM[power][j].time;

                        /* add spacial displacement squared */
                        R_CoM_squared[index].sum_dR =
                            R_CoM_squared[index].sum_dR + dot(dR, dR);

                        /* add to count of contributions */
                        R_CoM_squared[index].num_contributions =
                            R_CoM_squared[index].num_contributions + 1;
                    }

                }
            }

            /* replace oldest stored times */
            prev_R_CoM[power][0].time = prev_R_CoM[power][1].time;
            prev_R_CoM[power][1].time = prev_R_CoM[power][2].time;
            prev_R_CoM[power][2].time = time;
            /* replace oldest stored vectors */
            for(i=0; i<3; i++)
            {
                prev_R_CoM[power][0].r[i] = prev_R_CoM[power][1].r[i];
                prev_R_CoM[power][1].r[i] = prev_R_CoM[power][2].r[i];
                prev_R_CoM[power][2].r[i] = R_CoM[i];
            }

            /* increment overlap test */
            prev_R_CoM[power][0].overlap_test = prev_R_CoM[power][0].overlap_test + 1;

            power++;

        }
        else
        {test = false;}
    }
}

/* write CoM displacement squared data to file */
void write_cumulative_data(cumulative_data data[], int pow_max,
                           char *directory, char *suffix, char *prefix)
{
    FILE *write_file; /* file with wrapped coordinates */

    int i;

    /* declare strings for filenames */
    char file_str[150];
    sprintf(file_str,"%sR_CoM_files/%s%s", directory, prefix, suffix);

    write_file = fopen(file_str,"w");

    if(write_file != NULL)
    {
        for(i=0; i<=2*pow_max; i++)
        {
            fprintf(write_file,"%d\t%ld\t%lf\t%lf \n",
                    data[i].num_lags,
                    data[i].num_contributions,
                    data[i].sum_dt,
                    data[i].sum_dR);
        }

        fclose(write_file);
    }
    else
    {printf("Error: cumulative data file failed to open");}

}


/* output the set of chain centre of mass vectors at chosen evolution time */
void R_CoM_file(double xyz[][3], int Npoly, int Nmono, double time, double xyz_to_C[3][3],
                char *directory, char *suffix)
{
    FILE *Rcom_file; /* file with wrapped coordinates */

    int i;
    double Rcom[3], Cart_coords[3];

    /* declare strings for filenames */
    char file_str[150];
    sprintf(file_str, "%sR_CoM_files/R_CoM%s", directory, suffix);

    Rcom_file = fopen(file_str,"a");


    /* check file opened properly */
    if(Rcom_file!=NULL)
    {
        /* if at the start of the file wirte column information */
        if(ftell(Rcom_file)==0)
        {
            fprintf(Rcom_file,"time");

            for(i=1; i<=Npoly; i++)
            {fprintf(Rcom_file,"\t\tx%d\t\ty%d\t\tz%d",i,i,i);}

            fprintf(Rcom_file," \n");
        }

        fprintf(Rcom_file,"%lf",time);

        for(i=1; i<=Npoly; i++)
        {
            calc_R_CoM(Rcom, xyz, Nmono, i);
            transform_coords(Rcom, Cart_coords, xyz_to_C);

            fprintf(Rcom_file,"\t%lf\t%lf\t%lf",Cart_coords[0],Cart_coords[1],Cart_coords[2]);
        }

        fprintf(Rcom_file," \n");

        fclose(Rcom_file);
    }
    else
    {printf("Error: diff_file did not open\n");}
}


void off_diag_file(double dR_acc[3], double dR_acc_red[3], double time, char *directory, char *suffix)
{
    FILE *od_file; /* file with wrapped coordinates */

    /* declare strings for filenames */
    char file_str[150];
    sprintf(file_str, "%sR_CoM_files/off%s", directory, suffix);

    od_file = fopen(file_str,"a");


    /* check file opened properly */
    if(od_file!=NULL)
    {
        /* if at the start of the file wirte column information */
        if(ftell(od_file)==0)
        {fprintf(od_file,"time\t\tx\t\ty\t\tz\n");}

        fprintf(od_file,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n",time, dR_acc[0], dR_acc[1], dR_acc[2]);

        fclose(od_file);
    }
    else
    {printf("Error: off_diag_file did not open\n");}
}

/* output the difference between CoM diffusion of rejected and recycled moves  */
void write_dR_acc_rej(double dR_acc_rej, double time, char *directory, char *suffix)
{
    FILE *dR_acc_rej_file; /* file with wrapped coordinates */

    /* declare strings for filenames */
    char file_str[150];

    sprintf(file_str,"%sR_CoM_files/dR_acc_rej%s", directory, suffix);

    dR_acc_rej_file = fopen(file_str,"a");

    /* check file opened properly */
    if(dR_acc_rej_file!=NULL)
    {
        fprintf(dR_acc_rej_file,"%lf\t",time);
        fprintf(dR_acc_rej_file,"%lf \n",dR_acc_rej);

        fclose(dR_acc_rej_file);
    }
    else
    {printf("Error: dR_acc_rej_file did not open\n");}
}


/* output static structure factor */
void structure_factor_file(double xyz[][3], int Npoly, int Nmono, double metric[3][3])
{
    int i, j, m;
    /* k=wavenumber, S=structure factor (variable is k)*/
    double k=0.1, S, dS, vector[3], rij;

    FILE *struc_file;

    struc_file = fopen("structure_factor.txt","w");

    if(struc_file!=NULL)
    {
        /* loop over k values --> find S(k) at each value */
        while(k<10)
        {
            /* reinitialise S */
            S = 0.0;

            /* sum over all polymers */
            for(m=0; m<Npoly; m++)
            {
                /* reinitialise dS */
                dS=0.0;

                /* sum over pairs of different particles */
                for(i=0; i<Nmono; i++)
                {
                    for(j=i+1; j<Nmono; j++)
                    {
                        vector[0] = xyz[m*Nmono + i][0] - xyz[m*Nmono + j][0];
                        vector[1] = xyz[m*Nmono + i][1] - xyz[m*Nmono + j][1];
                        vector[2] = xyz[m*Nmono + i][2] - xyz[m*Nmono + j][2];

                        rij = sqrt(find_dr2(vector,metric));

                        dS += sin(k*rij)/(k*rij);
                    }
                }
                /* double count i!=j terms */
                dS *= 2;
                /* add contribution from i=j terms */
                S += dS + Nmono;

            }
            /* normalise average over all monomers and all chains */
            S /= (Npoly*Nmono);

            fprintf(struc_file,"{%lf,%lf}\n",k,S);

            /* increment k */
            k += 0.1;
        }


        fclose(struc_file);
    }
    else
    {printf("Error: struc_file did not open\n");}
}

/* output dynamics structure factor at chosen time
 * it follows mostly the same as the static structure factor function above
 * but needs to measure multiple times
 */
void dynamic_struc_fact_file(double xyz[][3], double xyz0[][3], int Npoly, int Nmono, double time, double metric[3][3])
{
    int i, j, m;

    double k=0.1, S, vector[3], ritj0;
    /* ritj0 = ri(time) - rj(0) */

    FILE *dyn_struc_file;

    dyn_struc_file = fopen("dynamic_structure_factor.txt","a");

    if(dyn_struc_file!=NULL)
    {
        /* loop over k values --> find S(k) at each value */
        while(k<10)
        {
            /* reinitialise S */
            S = 0.0;

            /* sum over all polymers */
            for(m=0; m<Npoly; m++)
            {
                /* sum over pairs of different particles */
                for(i=0; i<Nmono; i++)
                {
                    for(j=0; j<Nmono; j++)
                    {
                        vector[0] = xyz[m*Nmono + i][0] - xyz0[m*Nmono + j][0];
                        vector[1] = xyz[m*Nmono + i][1] - xyz0[m*Nmono + j][1];
                        vector[2] = xyz[m*Nmono + i][2] - xyz0[m*Nmono + j][2];

                        ritj0 = sqrt(find_dr2(vector,metric));

                        S += sin(k*ritj0)/(k*ritj0);
                    }
                }

            }
            /* normalise average over all monomers and all chains */
            S /= (Npoly*Nmono);

            fprintf(dyn_struc_file,"{%lf,%lf,%lf}\n",k,time,S);

            /* increment k */
            k += 0.1;
        }


        fclose(dyn_struc_file);
    }
    else
    {printf("Error: dyn_struc_file did not open\n");}
}


/* write files containing move acceptance statistics */
void accpt_files(double R_accpt[][3], double k_accpt[][3], int NR, int Nk)
{
    int i;

    /* Write file for wavelet statistics */

    FILE *R_file;
    R_file = fopen("R_accpt_stats.txt","w");

    if(R_file!=NULL)
    {
        for(i=0; i<NR; i++)
        {fprintf(R_file,"{%lf, %.0lf, %.0lf}\n",R_accpt[i][0],R_accpt[i][1],R_accpt[i][2]);}

        fclose(R_file);
    }
    else
    {printf("Error: wavelet accept/reject file did not open\n");}


    FILE *k_file;
    k_file = fopen("k_accpt_stats.txt","w");

    if(k_file!=NULL)
    {
        for(i=0; i<Nk; i++)
        {fprintf(k_file,"{%lf, %.0lf, %.0lf}\n",k_accpt[i][0],k_accpt[i][1],k_accpt[i][2]);}

        fclose(k_file);
    }
    else
    {printf("Error: Fourier accept/reject file did not open\n");}

}

