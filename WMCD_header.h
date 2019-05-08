/*
 * Header file declaring functions and structures for WMCD.
 * The function definitions can be found in WMCD_input.c, WMCD_output.c, WMCD_maths.c,
 *      WMCD_phys.c and WMCD_transforms.c.
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



#ifndef WAVESIM_HEADER_H_INCLUDED
#define WAVESIM_HEADER_H_INCLUDED

#include <stdbool.h>


/*-----------------------------*/
/* struct definitions          */
/*-----------------------------*/

typedef struct
{
    int hours;
    int minutes;
    double seconds;
}time_hms;

typedef struct
{
    double radius;
    double extent[3]; /* used to give size when transformed to flow coordinates */
    double centre[3];
    double centre_Cart[3];
    double polarisation[3];
    double amplitude;
}wavelet_parameters;

typedef struct
{
    double k;
    double khat[3];
    double khat_xyz[3];
    double polarisation[3];
    double polarisation_xyz[3];
    double phase;
    double amplitude;
}Fourier_parameters;

/* objects to store data for recycling failed moves */
typedef struct
{
    bool reuse;
    double radius;
    double centre[3];
}Recycle_wavelet;

typedef struct
{
    bool reuse;
    double k;
    double khat[3];
}Recycle_Fourier;


/* object to store length and energy scales for potentials,
 * as well as which type of potential to use */
typedef struct
{
    char solvent[10];
    char pot_type[10];
    double pot_length;
    double pot_energy;
    char bond_type[10];
    double bond_length;
    double spring_const;
}Potentials;

/* Collection of system parameters used throughout the code */
typedef struct
{
    bool periodic;
    int Lcart;
    double * (*pwrap)(double *, int);
    void (*pshort_wrap)(double *, int, double array1[][3], double array2[][3], double *);
    int sub_box_size;
    int L;
    int n_box_max;
    int Nmono;
    int Npoly;
    int Ntot;
    double lambda_min;
    double lambda_max;
    double lambda_a;
    double hydro_rad;
    double A0;
    double dt;
    double kT;
    double beta;
    double (*EV_fn)(double, Potentials, double *);
    double (*bond_fn)(double, Potentials, double *);
    double shear_rate;
    double reset_time;
    bool use_Glauber;
    bool use_Fourier;
    bool recycle;
    bool set_time_max;
    bool smart_MC;
    bool Langevin;
    double sMC_var;
    double sMC_mob;
    double time_max;
    long movenum_max;
    double Zimm_time;
}System_parameters;

/* Options for what to output and when */
typedef struct
{
    int N_out;
    bool system_params;
    bool eq_config_out;
    bool eq_config_in;
    bool cont_config_out;
    bool cont_config_in;
    bool wrapped_config;
    bool unwrapped_config;
    bool forces_out;
    bool accept_reject;
    bool polymer_radii;
    bool polymer_CoM;
    bool off_diag;
    bool CoM_displacement;
    bool static_structure_factor;
}In_out_options;


typedef struct
{
    double time;
    double r[3];
    long overlap_test;
}prev_quantities;


typedef struct
{
    int num_lags;
    double sum_dt;
    double sum_dR;
    long num_contributions;
}cumulative_data;


/**-----------------------------*/
/** input.c functions           */
/**-----------------------------*/

/* get list of input parameters */
void input(/*double *params_input, */System_parameters *params_input,
           Potentials *pot_input, In_out_options *out_input, double arg_option);

/* read in equilibrated particle configurations */
void read_equilib_config(double xyz[][3], double xyzn[][3], double xyz0[][3], double C_to_xyz[3][3],
                         System_parameters params, Potentials pot_params, char *directory);

/* read in configuration from an unfinished run */
void read_continue_config(double xyz[][3], double xyzn[][3], double xyz0[][3], int Ntot,
                    double *time, System_parameters params, char *directory, char *suffix);


/**-----------------------------*/
/** output.c functions          */
/**-----------------------------*/

/* construct directory and system dependent parts of file names */
void get_directory_and_suffix(char *directory, char *suffix, System_parameters params,
                              Potentials pot_params, int batch, int run_ID, bool periodic);

void write_parameters(System_parameters params, Potentials pot_params, int seed[2],
                      char *directory, char *suffix, time_hms run_time_hms, double Pacc, bool start);

/* write equilibrated configuration to file */
void write_equilib_config(double xyz[][3], int Ntot, double xyz_to_C[3][3], char *directory, char *suffix);

/* write configuration to file to continue from in a future run */
void write_continue_config(double xyz[][3], int Ntot, double time, double xyz_to_C[3][3], char *directory, char *suffix);

/* routines ouputting particle coordinates */
void unwrapped_coords_file(double xyz[][3], int Ntot, double time, double xyz_to_C[3][3], char *directory, char *suffix);

void wrapped_coords_file(double xyz[][3], int Ntot, int L, double time, double xyz_to_C[3][3],
                         char *directory, char *suffix);

/* output the forces on each particle at time */
void write_forces(double forces[][3], int Ntot, double time, char *directory, char *suffix);

/* output mean-square end-to-end and gyration radii at chosen evolution time */
void radii_file(double xyz[][3], int Npoly, int Nmono, double time, double metric[3][3],
                char *directory, char *suffix);


/* calculate CoM displacement squared at logarithmically increassing time lags */
void calc_dR_CoM_squared(cumulative_data R_CoM_squared[], prev_quantities prev_R_CoM[][3],
                         double R_CoM[3], double time, int output_num, int pow_max);

/* write CoM displacement squared data to file */
void write_cumulative_data(cumulative_data data[], int pow_max,
                           char *directory, char *suffix, char *prefix);

/* output the set of chain centre of mass vectors at chosen evolution time */
void R_CoM_file(double xyz[][3], int Npoly, int Nmono, double time, double xyz_to_C[3][3],
                char *directory, char *suffix);

void off_diag_file(double dR_acc[3], double dR_acc_red[3], double time, char *directory, char *suffix);

/* output the difference between CoM diffusion of rejected and recycled moves  */
void write_dR_acc_rej(double dR_acc_rej, double time, char *directory, char *suffix);

/* output static structure factor */
void structure_factor_file(double xyz[][3], int Npoly, int Nmono, double metric[3][3]);

/* output dynamic structure factor at chosen time */
void dynamic_struc_fact_file(double xyz[][3], double xyz0[][3], int Npoly, int Nmono, double time, double metric[3][3]);

/* write files containing move acceptance statistics */
void accpt_files(double R_accpt[][3], double k_accpt[][3], int NR, int Nk);




/**-----------------------------*/
/** wavemaths.c functions       */
/**-----------------------------*/

/* split time in seconds into hours, minutes and seconds */
void split_time(double t_sec, time_hms *t_hms);

/* get the sign of a number */
int signum(double x);

/* dot product between 2 vectors */
double dot(double vector1[3], double vector2[3]);

/* dot product between 3x3 matrices */
void matrix_dot(double A[3][3], double B[3][3], double C[3][3]);

/* cross product between 2 vectors */
double *cross(double vector1[3], double vector2[3]);

/* arbitrary vector transverse to a given vector */
void transverse_vector(double vector[3], double trans_vector[3]);

/* round number to nearest integer */
int nint(double number);
long long_nint(double number);

/* reliable modulus operator for floating point variables */
double mod_lf(double x, double y);

/* reliable modulus operator for integer variables */
int mod_int(int x, int y);
/* reliable modulus operator for long integer variables */
long mod_long(long x, long y);

/* min function: gives smallest of the two doubles */
double min_lf(double x, double y);

/* max function: gives smallest of the two doubles */
double max_lf(double x, double y);

/* efficinet pow. use whenever a (cumbersome) POSITIVE INTEGER power is required */
double pow_int(double x, int y);

/* find unit vector in a random direction */
void direction(double *vector);

/* Generate a random number with a Gaussian distribution with mean 0 and standard deviation 1 */
void Gaussian_rand(double *vector);

/* function numerically integrating cos(t)/t from x to 2x
 * this is equivalent to CosIntegral(x)-CosIntegral(2x)
 */
double Ci_x_2x(double x);




/**-----------------------------*/
/** wavephys.c functions        */
/**-----------------------------*/
/* block of potentials */

double EV_WCA(double dist2, Potentials params, double *F);
double EV_free(double dist2, Potentials params, double *F);
double EV_NGauss(double dist2, Potentials params, double *F);
double EV_Exp(double dist2, Potentials params, double *F);

double bond_FENE(double dist2, Potentials params, double *F);
double bond_quad(double dist2, Potentials params, double *F);


/* inter-particle potential */
double pot(double dist2, Potentials parameters);
/* bond potential */
double bondpot(double dist2, Potentials parameters);

double Metropolis_test(double V);
double Glauber_test(double V);

/* calculate the CoM position of a chosen polymer */
void calc_R_CoM(double R_CoM[3], double xyz[][3], int Nmono, int chosen_poly);

/* generate CDF for thermal distribution of bondlengths at set-up */
void thermalise(int Ns, double sCDF[Ns][2], double kT, Potentials pot_params);
/* build initial configuration */
void build_system(double xyz[][3], double xyzn[][3], double xyz0[][3],
                  System_parameters params, Potentials pot_params);

/* Generate wavelet parameters */
void get_wavelet_parameters(wavelet_parameters *pwavelet, double Rmin, double Rmax,
             bool recycle, Recycle_wavelet *recycled_w, int Ntot, double xyz[][3], double transform[3][3]);

/* permute wavevector components for recycling Fourier moves */
void permute_k(double prev_khat[3], Fourier_parameters *pplane_wave);

/* move particles according to plane wave vector field */
void movepw(double *position, Fourier_parameters plane_wave, double metric[3][3]);

/* wrapped coordinates */
/* NOTE: this is currently not the shortest disance from the origin, but from the centre of the box */
double *wrap(double position[3], int L_wrap);

/* find the shortest wrapped distance */
void short_wrap(double vector[3], int L_swrap, double metric[3][3], double basis_shifts[8][3], double basis_dist[8]);

/* find the shortest wrapped distance in a Cartesian (cubic) box */
void short_wrap_Cart(double vector[3], int L_swrap);

/* null versions of wrap and short_wrap for use when periodic BC are not used, but so the
 * code can be written using a pointer to either the wrap or dont_wrap versions
 */
double *dont_wrap(double position[3], int L_wrap);

void dont_short_wrap(double vector[3], int L_swrap, double metric[3][3],
                     double basis_shifts[8][3], double basis_dist[8]);

/* box in which a particle lies */
void box(double position[3], int L_box, int bs_box, int xyz_box[3]);

/* map the particles into? the grid. use sizeof to get listlength */
void map(int listlength,int part_list[],double xyzn[][3],int L,int bs,int Lbs,
         int nboxmax,int contents[Lbs][Lbs][Lbs][nboxmax],
         int anchor[3], int spans[3]);

/* find 'lowest' corner box moved by wavelet */
void anchor(double position[3], double extent[3], int L, int bs, int anchor_vector[3]);

/* find the span of a wavelet (dist between 'lowest' and 'highest' corners) */
void spans(double position[3], double extent[3], int bs, int span_vector[3]);


/* total of forces felt by particles in a wavelet move */
void sum_forces_w(double xyz[][3], int movers[], System_parameters params, wavelet_parameters wavelet,
                  double forces[][3], double total_force[3], double xyz_to_C[3][3], double metric[3][3],
                  double images[8][3], double image_dist[8]);

/* total of forces felt by particles in a Fourier move */
void sum_forces_F(double xyz[][3], int Ntot, double forces[][3], double total_force[3],
                  double xyz_to_C[3][3], Fourier_parameters plane_wave);

void add_accept_reject_stats(double move_accpt[3], double R_accpt[][3], double k_accpt[][3]);

/* initialise arrays involved in analysing data during the run */
void initialise_cumulative_variable(cumulative_data data[], prev_quantities prev_data[][3], int max_power1);

/* function defining the unnormalised probability for a plane wave with x=|k|*Rmax */
double pwProb(double x);

/* Sum pwProb at discrete values
 * Used for normailisation and to determine the number of plane waves needed
 */
double pwProbSum(double k0, int Lmax);

double pwProbSumk2(double k0, double Rmax_pwSum, int Lmax);

/* function finding CDF for planewaves via numerical integration */
void getCDF(int Nsteps, double CDF[][2], double kmin, double kmax, double Rmax_cdf);


/**-----------------------------*/
/** transforms.c functions      */
/**-----------------------------*/
/* basis vectors for box under solvent flows */
void lattice_vectors(double a1[3], double a2[3], double a3[3], double time, System_parameters params);

/* lattice vectors in Fourier space */
void reciprocal_vectors(double b1[3], double b2[3], double b3[3], double a1[3], double a2[3], double a3[3], int L);

/* metric of coordinates under flow */
void get_metric(double metric[3][3], double a1[3], double a2[3], double a3[3]);

/* transformation matrix from skewed to Cartesian coordinates */
void xyz_to_Cart_transform(double transform[3][3], double a1[3], double a2[3], double a3[3]);

/* transformation matrix from Cartesian to skewed coordinates */
void Cart_to_xyz_transform(double transform[3][3], double a1[3], double a2[3], double a3[3]);

void transform_coords(double r[3], double R[3], double r_to_R[3][3]);

/* use the metric to obtain the scalar length squared of a vector */
double find_dr2(double dr[3], double metric[3][3]);

/* use the metric to obtain the scalar product of two vectors */
double find_drdR(double dr[3], double dR[3], double metric[3][3]);

/* find the maximal extent of the wavelet, which is ellipsoidal in flow coordinates */
void get_wavelet_axes(double axes[3], double radius, double C_to_xyz[3][3]);

/* move particles in the flow coordinates when reset to t=0 basis vectors */
void remap_at_coord_reset(double xyz[][3], double xyzn[][3], int Ntot, double time,
                          System_parameters params, int Lbs, int contents[Lbs][Lbs][Lbs][params.n_box_max]);


#endif // WAVESIM_HEADER_H_INCLUDED
