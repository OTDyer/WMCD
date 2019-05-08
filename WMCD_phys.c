/*
 * Functions for WMCD.
 * The functions defined in this file calculate, initialise or otherwise manipulate
 *      physical properties of the system.
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

#define RAND_MAX_INV 1.0/RAND_MAX
#define POW2_32_INV 1.0/(pow(2,32)-1)

/** start of block of potentials and force functions (EV = Excluded Volume) */

/* Weeks-Chandler-Anderson potential */
double EV_WCA(double dist2, Potentials params, double *F)
{
    double U, sig_r6;
    double dist2_inv = 1.0/dist2;

    if(dist2>1.259921*params.pot_length*params.pot_length)
    {
        U = 0.0;
        *F = 0.0;
    }
    else
    {
        sig_r6 = pow_int(params.pot_length*params.pot_length*dist2_inv,3); /* = (sigma/r)^6 */
        U = params.pot_energy*(4.0*sig_r6*(sig_r6-1.0) + 1.0);

        *F = 48.0*params.pot_energy*dist2_inv*sig_r6*(sig_r6-0.5);
    }

    return U;
}


double EV_free(double dist2, Potentials params, double *F)
{
    *F = 0.0;
    return 0.0;
}

/* Narrow Gaussian potential ???? */
double EV_NGauss(double dist2, Potentials params, double *F)
{
    double U;
    double width_inv = 1.0/params.pot_length, s2 = params.pot_length*params.pot_length;
    double Gauss = params.pot_energy * exp(-0.5*dist2*width_inv*width_inv);

    if(dist2 < 3.0*s2)
    {
        U = - (dist2 - s2)*Gauss + 2.0*0.22313*params.pot_energy*s2;
        *F = - width_inv*width_inv * (dist2 - 3.0*params.pot_length*params.pot_length) * Gauss;
    }
    else
    {
        U = 0.0;
        *F = 0.0;
    }
//    /** this may want a cut-off threshold to be added */
//    U = params.pot_energy * pow_int(width_inv, 3) * exp(-0.5*dist2*width_inv*width_inv);
//    *F = width_inv * width_inv * U;

    return U;
}

/* Exponential potential */
double EV_Exp(double dist2, Potentials params, double *F)
{
    double U;
    double r = sqrt(dist2), rc=0.5, s=1.0/params.pot_length;

    if(r < rc)
    {
        U = params.pot_energy * (exp(-r*s) - exp(-rc*s));
        *F = s * U / r;
    }
    else
    {
        U = 0.0;
        *F = 0.0;
    }

    return U;

}


/* Finitely Extensible Non-linear Elastic potential */
double bond_FENE(double dist2, Potentials params, double *F)
{
    double U, arg;

    if(dist2<params.bond_length*params.bond_length)
    {
        arg = 1.0 - dist2/(params.bond_length*params.bond_length);

        U = -0.5*params.spring_const*params.bond_length*params.bond_length*log(arg);
        *F = -params.spring_const/arg;
    }
    else
    {
        U = 99999999.0;
        *F = -99999.0;
    }

    return U;
}

/* Quadratic (Hookean) potential */
double bond_quad(double dist2, Potentials params, double *F)
{
    double U;

    U = 0.5*params.spring_const*(dist2 - 2*sqrt(dist2)*params.bond_length + params.bond_length*params.bond_length);
    *F = -params.spring_const*(1.0 - params.bond_length/sqrt(dist2));

    return U;
}

/** end of potentials block */

/** pot and bondpot are basically defunt, but have not yet been unpicked from the thermalise and build_system functions **/
/* inter-particle potential */
double pot(double dist2, Potentials parameters)
{
    double U;
    double length_scale = parameters.pot_length;
    double energy_scale = parameters.pot_energy;

    if(strcmp(parameters.pot_type,"WCA")==0) // WCA potential
    {
        if(dist2>1.259921*length_scale*length_scale)
        {U = 0.0;}
        else
        {U = pow_int(length_scale*length_scale/dist2,3); /* = (sigma/r)^6 */
         U = energy_scale*(4.0*U*(U-1.0) + 1.0);}
    }
    else // free particles (theta solvent)
    {U = 0.0;}

    return U;
}
/* bond potential between adjacent monomers */
double bondpot(double dist2, Potentials parameters)
{
    double U;
    double length_scale = parameters.bond_length;
    double spring_constant = parameters.spring_const;

    if(strcmp(parameters.bond_type,"FENE")==0) // FENE potential
    {
        if(dist2<length_scale*length_scale)
        {U = -0.5*spring_constant*length_scale*length_scale*log(1.0-dist2/(length_scale*length_scale));}
        else
        {U = 99999999.0;printf("stretched!\t");}
    }
    else // Hookean (quadratic) potential
    {U = 1.5*spring_constant*(dist2 - 2*sqrt(dist2)*length_scale + length_scale*length_scale);}

    return U;
}

double Metropolis_test(double V)
{return min_lf(exp(V), 1.0);}

double Glauber_test(double V)
{return 0.5*(1-tanh(-0.5*V));}


/* calculate the CoM position of a chosen polymer */
void calc_R_CoM(double R_CoM[3], double xyz[][3], int Nmono, int chosen_poly)
{
    int i,j;
    int poly_start=(chosen_poly-1)*Nmono;
    double Nmono_inv=1.0/Nmono;

    for(i=0; i<3; i++)
    {
        R_CoM[i] = 0.0;

        for(j=0; j<Nmono; j++)
        {R_CoM[i] += xyz[poly_start+j][i];}

        R_CoM[i] *= Nmono_inv;
    }
}

/* build CDF for thermal distribution of bond lengths at set-up */
void thermalise(int Ns, double sCDF[][2], double kT, Potentials pot_params)
{
    int i;

    /*----------------------------*/
    /* Find bounds on bond length */
    /*----------------------------*/
    double smin=0.8, smax=1.0, Umax=10.0;

    /* find minimum bond length from inverted pot() */
    if(strcmp(pot_params.pot_type,"free") == 0)
    {smin = 0.5;}
    else if(strcmp(pot_params.pot_type,"WCA") == 0)
    {smin = (pot_params.pot_length)
            * pow(2.0/(1+sqrt(Umax/(pot_params.pot_energy) )) ,1.0/6);}

    /* find maximum bond length from inverted bondpot() */
    if(strcmp(pot_params.bond_type,"quadratic") == 0)
    {
        smax = (pot_params.bond_length) + sqrt(2*Umax/(3*pot_params.spring_const) );
        /* failsafe against pathological parameters */
        if(smax < smin)
        {smax = 2.0*smin;}
    }
    else if(strcmp(pot_params.bond_type,"FENE") == 0)
    {smax = 0.75*(pot_params.bond_length) * sqrt(1.0 - exp(-2*Umax/pot_params.spring_const));}


    /*----------------------------*/
    /* Build CDF within bounds    */
    /*----------------------------*/

    /* separation step, sep. at loop i, 1/kT to be more efficient */
    double ds=(smax-smin)/(Ns-1), si;
    /* 1/quatities to be more efficient */
    double kTinv=1.0/kT, CDFnorm;

    /* fill sCDF with cumulative probabilities ...[0]
     * and associated bond lengths ...[1]
     */
    sCDF[0][0] = exp(-(pot(smin, pot_params)+bondpot(smin, pot_params))*kTinv);
    sCDF[0][1] = smin;

    for(i=1; i<Ns; i++)
    {
        si = smin + i*ds;

        sCDF[i][0] = sCDF[i-1][0] + exp(-(pot(si, pot_params)+bondpot(si, pot_params))*kTinv);
        sCDF[i][1] = si;
    }

    /* Normalise probabilities */
    CDFnorm = 1.0/sCDF[Ns-1][0];
    for(i=0; i<Ns; i++)
    {sCDF[i][0] *= CDFnorm;}
}


/* build initial configuration */
void build_system(double xyz[][3], double xyzn[][3], double xyz0[][3], System_parameters params,
                  Potentials pot_params)
{
    int i,j,k,l;

    /* extract relevant parameters */
    int L=params.L;
    int Nmono=params.Nmono, Npoly=params.Npoly;
    double kT = 0.1 * params.kT; /* reduce this to make initial system close(r) to equilibrium */

    /* variables containing bond vector between monomers along the polymer */
    double bond_hat[3], bond_length;

    double init_sep[3];
    int init_test=0, init_test_count;

    /* bond vector for straight initial condistion */
/*    const double bond_theta = 0.2, bond_phi = 0.2;
    const double bond_vector[3] = {cos(bond_phi)*sin(bond_theta),sin(bond_phi)*sin(bond_theta),cos(bond_theta)};
*/
    /* CDF for Boltzmann distribution between neighbouring monomers along a polymer */
    double sepCDF[100][2];
    thermalise(100,sepCDF,kT,pot_params);

    /* loop over polymers */
    for(i=0; i<Npoly; i++)
    {
        /* generate location of first monomer */
        for(k=0; k<3; k++)
        {xyz[i*Nmono][k] = L*POW2_32_INV*pcg32_random();
         xyzn[i*Nmono][k] = xyz[i*Nmono][k];
         xyz0[i*Nmono][k] = xyz[i*Nmono][k];}

        /* loop over monomers */
        for(j=1; j<Nmono; j++)
        {
            init_test_count=0;
            init_test=0;

            /* find bond vector until sufficiently far from all other particles */
            while(init_test!=2)
            {
                /* find direction */
                direction(bond_hat);
                /* get bond length from thermal (Boltzmann) dist for chosen pot'ls */
                k=0;

                bond_length = POW2_32_INV*pcg32_random();
                while(sepCDF[k][0]<bond_length){k++;}
                bond_length = sepCDF[k][1];

                for(k=0; k<3; k++)
                {
                    xyz[i*Nmono+j][k] = xyz[i*Nmono+j-1][k] + bond_length*bond_hat[k];
                }

                /* test if this location collides with another particle */
                l=0;
                init_test=0;

                /* init_test_count is used to make sure the code doesn't get stuck if collision is unavoidable.
                 * in these cases the system has to sort it out when it evolves
                 */
                while(l<(i*Nmono+j-1) && init_test_count<15)
                {

                    if(params.periodic == true) // Periodic box --> need to wrap coords
                    {
                        /* find the wrapped vector between particles */
                        init_sep[0] = mod_lf( xyz[i*Nmono+j][0] - xyz[l][0] ,L);
                        init_sep[1] = mod_lf( xyz[i*Nmono+j][1] - xyz[l][1] ,L);
                        init_sep[2] = mod_lf( xyz[i*Nmono+j][2] - xyz[l][2] ,L);

                        /* make sure it is the shortest vector by shifting to be the nearest image */
                        short_wrap_Cart(init_sep,L);
                    }
                    else // Infinite box --> don't wrap coords
                    {
                        /* find the wrapped vector between particles */
                        init_sep[0] = xyz[i*Nmono+j][0] - xyz[l][0];
                        init_sep[1] = xyz[i*Nmono+j][1] - xyz[l][1];
                        init_sep[2] = xyz[i*Nmono+j][2] - xyz[l][2];
                    }

                    /* find and note a collision with an already existing particle */
                    if(pot(dot(init_sep,init_sep), pot_params)>3*kT)
                    {
                        init_test=1;
                        break;
                    }

                    l++;
                }

                /* if there are no collisions, allow this position */
                if(init_test==0)
                {init_test = 2;}

                /* increment counter for unphysical position choices */
                init_test_count++;
            }


            /* copy to new positions */
            xyzn[i*Nmono+j][0] = xyz[i*Nmono+j][0];
            xyzn[i*Nmono+j][1] = xyz[i*Nmono+j][1];
            xyzn[i*Nmono+j][2] = xyz[i*Nmono+j][2];
            /* copy to original positions */
            xyz0[i*Nmono+j][0] = xyz[i*Nmono+j][0];
            xyz0[i*Nmono+j][1] = xyz[i*Nmono+j][1];
            xyz0[i*Nmono+j][2] = xyz[i*Nmono+j][2];
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


/* permute wavevector components for recycling Fourier moves */
void permute_k(double prev_khat[3], Fourier_parameters *pplane_wave)
{
    int i, rand1, rand2, randsign;
    double k_old[3] = {prev_khat[0],prev_khat[1],prev_khat[2]};

    /* generate a random number in {0,1,2} */
    rand1 = (int)pcg32_boundedrand(3);

    /* generate a random number in {-1,+1} */
    randsign = 2*(int)pcg32_boundedrand(2) - 1;
    /* use the sign to generate a random number in {0,1,2}\{rand1} */
    rand2 = mod_int(rand1 + randsign, 3);

    /* use rand1 and rand2 to permute k components */
    pplane_wave->khat[rand1] = k_old[0];
    pplane_wave->khat[rand2] = k_old[1];
    pplane_wave->khat[3 - (rand1 + rand2)] = k_old[2];

    /* apply a random sign change to each component (probability 0.5 in each case) */
    for(i=0; i<3; i++)
    {
        randsign = 2*(int)pcg32_boundedrand(2) - 1;
        pplane_wave->khat[i] = randsign*pplane_wave->khat[i];
    }


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


/* wrapped coordinates */
double *wrap(double position[3], int L_wrap)
{
    static double xyz_wrapped[3];

    /* wrap in each direction */
    xyz_wrapped[0] = mod_lf(position[0],(double)L_wrap);
    xyz_wrapped[1] = mod_lf(position[1],(double)L_wrap);
    xyz_wrapped[2] = mod_lf(position[2],(double)L_wrap);

    return xyz_wrapped;
}

/* find the shortest wrapped distance */
void short_wrap(double vector[3], int L_swrap, double metric[3][3], double basis_shifts[8][3], double basis_dist[8])
{
    int i, min_shift=0;
    double *wrap000;
    double min_dist=0.0, dist;

    wrap000 = wrap(vector, L_swrap);

    /* check which image is closest from sum of basis and vector-basis cross terms */
    /* no need to do this for {0,0,0} shifts as always equals 0 */
    for(i=1; i<8; i++)
    {
        dist = 2*find_drdR(wrap000, basis_shifts[i], metric) + basis_dist[i];

        if(dist < min_dist)
        {
            min_dist = dist;
            min_shift = i;
        }
    }

    /* shift to be the nearest image */
    for(i=0; i<3; i++)
    {
        vector[i] = wrap000[i] + (double)basis_shifts[min_shift][i];
    }

}

/* find the shortest wrapped distance in a Cartesian (cubic) box */
void short_wrap_Cart(double vector[3], int L_swrap)
{
    int i;

    /* shift to be the nearest image */
    for(i=0; i<3; i++)
    {
        while(fabs(vector[i])>0.5*L_swrap)
        {
            //* shift in correct direction /
            if(vector[i]>0)
            {vector[i] -= L_swrap;}
            else
            {vector[i] += L_swrap;}
        }
    }

}


/* null versions of wrap and short_wrap for use when periodic BC are not used, but so the
 * code can be written using a pointer to either the wrap or dont_wrap versions
 */
double *dont_wrap(double position[3], int L_wrap)
{return position;}

void dont_short_wrap(double vector[3], int L_swrap, double metric[3][3],
                     double basis_shifts[8][3], double basis_dist[8])
{/* this is intentionally empty as the function should not do anything*/}

/* box in which a particle lies */
void box(double position[3], int L_box, int bs_box, int xyz_box[3])
{
    double wrp[3];
    double bs_inv = 1.0/bs_box;

    wrp[0] = mod_lf(position[0],(double)L_box);
    wrp[1] = mod_lf(position[1],(double)L_box);
    wrp[2] = mod_lf(position[2],(double)L_box);

    xyz_box[0] = floor( wrp[0]*bs_inv);
    xyz_box[1] = floor( wrp[1]*bs_inv);
    xyz_box[2] = floor( wrp[2]*bs_inv);
}


/* map the particles into? the grid. use sizeof to get listlength */
void map(int listlength,int part_list[],double xyzn[][3],int L,int box_size,int Lbs,
         int nboxmax,int contents[Lbs][Lbs][Lbs][nboxmax],
         int anchor[3], int spans[3])
{
    int i,j,k,ip,ix,iy,iz; /* loop indices and box counters */
    int box_loc[3]; /*box coords for ith particle */
    int contentscount[spans[0]+1][spans[1]+1][spans[2]+1];/* used to check boxes aren't overfull */
    double position[3];
    /* used if spans crosses the box boundary (-1 in Lbs-1 is because C starts indices at 0) */
    int b[3] = {(anchor[0] + spans[0] > Lbs-1) ? Lbs - anchor[0] : 0,
                (anchor[1] + spans[1] > Lbs-1) ? Lbs - anchor[1] : 0,
                (anchor[2] + spans[2] > Lbs-1) ? Lbs - anchor[2] : 0};


    /* initialise contentscount */
    for(i=0; i<=spans[0]; i++)
    {for(j=0; j<=spans[1]; j++)
        {for(k=0; k<=spans[2]; k++)
            {contentscount[i][j][k] = 0;}
        }
    }

    /* map each of the particles in plist[] */
    for(ip=0; ip<listlength; ip++)
    {
        /* get position vector of ith particle */
        position[0] = xyzn[part_list[ip]-1][0];
        position[1] = xyzn[part_list[ip]-1][1];
        position[2] = xyzn[part_list[ip]-1][2];

        /* get box coords of ith particle */
        box(position,L,box_size,box_loc);
        ix = box_loc[0];
        iy = box_loc[1];
        iz = box_loc[2];

        if(ix >= anchor[0])
        {i = ix - anchor[0];}
        else
        {i = ix + b[0];}

        if(iy >= anchor[1])
        {j = iy - anchor[1];}
        else
        {j = iy + b[1];}

        if(iz >= anchor[2])
        {k = iz - anchor[2];}
        else
        {k = iz + b[2];}

        /* add particle to contents map */
        contents[ix][iy][iz][contentscount[i][j][k]] = part_list[ip];

        /* count particles in this box */
        contentscount[i][j][k]++;

        /* run an error if a box is full */
        if(contentscount[i][j][k] == nboxmax)
        {
            printf("Error: box (%d,%d,%d) is full\n",ix,iy,iz);
            abort();
        }
    }
}


/* find 'lowest' corner box moved by wavelet */
void anchor(double position[3], double extent[3], int L, int bs, int anchor_vector[3])
{
    double anchor_pos[3];
    /* find anchor location as double */
    anchor_pos[0] = position[0] - extent[0];
    anchor_pos[1] = position[1] - extent[1];
    anchor_pos[2] = position[2] - extent[2];

    box(anchor_pos, L, bs, anchor_vector);
}


/* find the span of a wavelet (dist between 'lowest' and 'highest' corners) */
void spans(double position[3], double extent[3], int bs, int span_vector[3])
{
    double bs_inv = 1.0/bs;

    span_vector[0] = floor((position[0]+extent[0])*bs_inv) - floor((position[0]-extent[0])*bs_inv);
    span_vector[1] = floor((position[1]+extent[1])*bs_inv) - floor((position[1]-extent[1])*bs_inv);
    span_vector[2] = floor((position[2]+extent[2])*bs_inv) - floor((position[2]-extent[2])*bs_inv);
}

/* total of forces felt by particles in a wavelet move, weighted by the
 * (relative) displacement ampliditude at their location
 */
void sum_forces_w(double xyz[][3], int movers[], System_parameters params, wavelet_parameters wavelet,
                  double forces[][3], double total_force[3], double xyz_to_C[3][3], double metric[3][3],
                  double images[8][3], double image_dist[8])
{
    int i=0, ipart;
    double Cart_coords[3], r[3], r_scalar, weight, *torque, radius_inv=1.0/wavelet.radius;

    total_force[0] = 0.0;
    total_force[1] = 0.0;
    total_force[2] = 0.0;

    while(movers[i] > 0 && i < params.Ntot)
    {
        ipart = movers[i] - 1;

        transform_coords(xyz[ipart], Cart_coords, xyz_to_C);

        r[0] = Cart_coords[0] - wavelet.centre_Cart[0];
        r[1] = Cart_coords[1] - wavelet.centre_Cart[1];
        r[2] = Cart_coords[2] - wavelet.centre_Cart[2];

        params.pshort_wrap(r,params.L,metric,images,image_dist);
        r_scalar = sqrt(dot(r,r));

        /* scale by relative displacement at the particle */
        weight = (1.0 - r_scalar*radius_inv);

        torque = cross(r, forces[ipart]);

        total_force[0] += torque[0] * weight;
        total_force[1] += torque[1] * weight;
        total_force[2] += torque[2] * weight;

        i++;
    }

}

/* total of forces felt by particles in a Fourier move, weighted by the
 * (relative) displacement ampliditude at their location
 */
void sum_forces_F(double xyz[][3], int Ntot, double forces[][3], double total_force[3],
                  double xyz_to_C[3][3], Fourier_parameters plane_wave)
{
    int i;
    double Cart_coords[3], weight, Fkhat;

    total_force[0] = 0.0;
    total_force[1] = 0.0;
    total_force[2] = 0.0;

    for(i=0; i<Ntot; i++)
    {
        transform_coords(xyz[i], Cart_coords, xyz_to_C);

        weight = cos(dot(plane_wave.khat, Cart_coords)*plane_wave.k + plane_wave.phase);

        /* sum the component of forces perpendicular to khat */
        total_force[0] = total_force[0] + weight * forces[i][0];
        total_force[1] = total_force[1] + weight * forces[i][1];
        total_force[2] = total_force[2] + weight * forces[i][2];
    }

    /* remove the component parallel to k */
    Fkhat = dot(total_force, plane_wave.khat);

    total_force[0] = total_force[0] - Fkhat*plane_wave.khat[0];
    total_force[1] = total_force[1] - Fkhat*plane_wave.khat[1];
    total_force[2] = total_force[2] - Fkhat*plane_wave.khat[2];
}



void add_accept_reject_stats(double move_accpt[3], double R_accpt[][3], double k_accpt[][3])
{
    int i;

    if(move_accpt[0]<1) /* add wavelet statistics */
    {
        i=0;
        while(R_accpt[i][0]<move_accpt[1]){i++;}

        R_accpt[i][1]++; /* count the move in the bin */
        R_accpt[i][2] += (int)move_accpt[2]; /* add 1/0 for success/fail */
    }
    else /* add Fourier statistics */
    {
        i=0;
        while(k_accpt[i][0]<move_accpt[1]){i++;}

        k_accpt[i][1]++; /* count the move in the bin */
        k_accpt[i][2] += (int)move_accpt[2]; /* add 1/0 for success/fail */
    }
}


/* initialise arrays involved in analysing data during the run */
void initialise_cumulative_variable(cumulative_data data[], prev_quantities prev_data[][3], int max_power1)
{
    int i,j,k;

    /* set all elements of all structs to 0, except the number of lags */
    j=3;
    for(i=0; i<=max_power1; i++)
    {
        if(i==0)
        {
            for(k=0; k<3; k++)
            {
                data[k].num_lags = k+1;
                data[k].sum_dt = 0.0;
                data[k].sum_dR = 0.0;
                data[k].num_contributions = 0;
            }
        }
        else
        {
            data[j].num_lags = 2*pow_int(2,i);
            data[j+1].num_lags = 3*pow_int(2,i);
            data[j].sum_dt = 0.0;
            data[j+1].sum_dt = 0.0;
            data[j].sum_dR = 0.0;
            data[j+1].sum_dR = 0.0;
            data[j].num_contributions = 0;
            data[j+1].num_contributions = 0;
            j+=2;
        }
    }

    for(i=0; i<=max_power1; i++)
    {for(j=0; j<3; j++)
     {
         prev_data[i][j].time = 0.0;
         prev_data[i][j].r[0] = 0.0;
         prev_data[i][j].r[1] = 0.0;
         prev_data[i][j].r[2] = 0.0;
         prev_data[i][j].overlap_test = 0;
     }
    }
}


/* function defining the unnormalised probability for a plane wave with x=|k|*Rmax
 * it is wavelet dependent, so will need changing along with wavelet choices
 */
double pwProb(double x)
{
    double Probvalue=0.0, x2=x*x;

    /* Use Taylor series accurate until asymptotic behviour */
    if(x<=6.15)
    {
        Probvalue += 4.112883994E-3;//4.11288399466719E-3
        if(x!=0.0) /* save time on the commonly used x=0.0 */
        {
            Probvalue -= 3.086419753E-5 * pow_int(x2,2);
            Probvalue += 2.204585537E-6 * pow_int(x2,3);
            Probvalue -= 8.103164105E-8 * pow_int(x2,4);
            Probvalue += 1.945846302E-9 * pow_int(x2,5);
            Probvalue -= 3.369850151E-11 * pow_int(x2,6);
            Probvalue += 4.453909307E-13 * pow_int(x2,7);
            Probvalue -= 4.665325271E-15 * pow_int(x2,8);
            Probvalue += 3.980495643E-17 * pow_int(x2,9);
            Probvalue -= 2.825217277E-19 * pow_int(x2,10);
            Probvalue += 1.696385229E-21 * pow_int(x2,11);
            Probvalue -= 8.736647028E-24 * pow_int(x2,12);
            Probvalue += 3.904219287E-26 * pow_int(x2,13);
        }
    }
    else /* Asymptotic form, ignores small oscillations */
    {
        Probvalue = 0.1239081775267918/(x2*x2);
    }

    /* Mulitipy by 16\pi^2 so fill in missing (4\pi)^2 in the series */
    Probvalue *= 16*M_PI*M_PI;

    return Probvalue;
}

/* Sum pwProb at discrete values
 * Used for normailisation and to determine the number of plane waves needed
 */
double pwProbSum(double k0, int Lmax)
{
    int i, j, L_layer=0;
    double k, Sumn, dSum;

    Sumn = pwProb(0.0);

    /* add contribution from cubic layers out to Lmax
     * if Lmax large (approximating infinity), just sum until contributions are insignificant
     */
    while(/*(Sumn-Sum > 0.01*pwProb(0.0)) &&*/ (L_layer<Lmax))
    {
        L_layer++; /* increment to next layer */

        /* add contribution from the 6 faces, not including any edges
         * for this the face is split into 4 tiles with equal contribution
         */
        dSum = 0.0;
        for(i=0; i<=L_layer-1; i++)
        {for(j=1; j<=L_layer-1; j++)
         {
             k = sqrt(i*i + j*j + L_layer*L_layer)*k0;
             dSum += pwProb(k);
         }
        }
        dSum *= 4; /* multiply up to all 4 tiles */
        dSum += pwProb(L_layer*k0); /* add in missing i=j=0 point */

        Sumn += 6*dSum;

        /* add contribution from the 12 edges, not including corners
         * each edge is split into 2 halves, with equal contribution
         */
        dSum = 0.0;
        for(i=1; i<=L_layer-1; i++)
        {
            k = sqrt(i*i + 2*L_layer*L_layer)*k0;
            dSum += pwProb(k);
        }
        dSum *= 2;
        dSum += pwProb(sqrt(2.0)*L_layer*k0); /* add in missing i=0 point */

        Sumn += 12*dSum;

        /* add contribution from the 8 corners */
        Sumn += 8*pwProb(sqrt(3.0)*L_layer*k0);

    }

    Sumn -= pwProb(0.0); // Take out the k=0 mode

    return Sumn;
}


/* Sum pwProb*k^{-2} at discrete values
 * Used to determine the number of plane waves needed
 */
double pwProbSumk2(double k0, double Rmax_pwSum, int Lmax)
{
    int i, j, L_layer=0;
    double k, Sumn=0.0, dSum;

    /* add contribution from cubic layers out to Lmax
     * if Lmax large (approximating infinity), just sum until contributions are insignificant
     */
    while(L_layer<Lmax)
    {
        L_layer++; /* increment to next layer */

        /* add contribution from the 6 faces, not including any edges
         * for this the face is split into 4 tiles with equal contribution
         */
        dSum = 0.0;
        for(i=0; i<=L_layer-1; i++)
        {for(j=1; j<=L_layer-1; j++)
         {
             k = sqrt(i*i + j*j + L_layer*L_layer)*k0;
             dSum += pwProb(k*Rmax_pwSum)/(k*k);
         }
        }
        dSum *= 4; /* multiply up to all 4 tiles */

        k = L_layer*k0;
        dSum += pwProb(k*Rmax_pwSum)/(k*k); /* add in missing i=j=0 point */

        Sumn += 6*dSum;

        /* add contribution from the 12 edges, not including corners
         * each edge is split into 2 halves, with equal contribution
         */
        dSum = 0.0;
        for(i=1; i<=L_layer-1; i++)
        {
            k = sqrt(i*i + 2*L_layer*L_layer)*k0;
            dSum += pwProb(k*Rmax_pwSum)/(k*k);
        }
        dSum *= 2;

        k = sqrt(2.0)*L_layer*k0;
        dSum += pwProb(k*Rmax_pwSum)/(k*k); /* add in missing i=0 point */

        Sumn += 12*dSum;

        /* add contribution from the 8 corners */
        k = sqrt(3.0)*L_layer*k0;
        Sumn += 8*pwProb(k*Rmax_pwSum)/(k*k);
    }

    return Sumn;
}


/* function finding CDF for planewaves via numerical integration.
 * this integrates only crudely (with trapezium rule) as results are later discretised anyway
 */
void getCDF(int Nsteps, double CDF[][2], double kmin, double kmax, double Rmax_cdf)
{
    int i;
    double widthk=(kmax-kmin)/Nsteps, width=widthk*Rmax_cdf, kRmin=kmin*Rmax_cdf;
    double kRi, kRi_;
    /* width of each slice. is multiplied by Rmax
     * due to a more convenient choice of variable for the integration.
     * widthk is still used to give the corresponding k value
     */

    CDF[0][0] = 0.0;
    CDF[0][1] = kmin;

    for(i=1; i<=Nsteps; i++)
    {
        kRi = kRmin + width*i;
        kRi_ = kRmin + width*(i-1);

        /* weight PDFs by k^2 to account for the area of sphere */
        CDF[i][0] = CDF[i-1][0] + 0.5*width*( kRi_*kRi_*pwProb(kRi_) + kRi*kRi*pwProb(kRi) );

        CDF[i][1] = kmin + widthk*i;

    }

    /* normalise elements */
    double norm = 1.0/CDF[Nsteps][0];

    for(i=1; i<=Nsteps; i++)
    {CDF[i][0] *= norm;}
}





