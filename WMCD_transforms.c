/*
 * Functions for coordinate transforms in WMCD.
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

#include "WMCD_header.h"
#include "pcg_basic.h"

#define RAND_MAX_INV 1.0/RAND_MAX
#define POW2_32_INV 1.0/(pow(2,32)-1)
// 'max' random numbers to give open set [0,1) upon multiplication of RNG
#define POW2_32_INV_OPEN 1.0/pow(2,32)



/* basis vectors for box under solvent flows
 * these are entered explicitly and used to obtain the metric and coordinate transforms
 */
void lattice_vectors(double a1[3], double a2[3], double a3[3], double time, System_parameters params)
{
    a1[0] = 1.0;
    a1[1] = 0.0;
    a1[2] = 0.0;

    a2[0] = -0.5*signum(params.shear_rate) + params.shear_rate*time;
    a2[1] = 1.0;
    a2[2] = 0.0;

    a3[0] = 0.0;
    a3[1] = 0.0;
    a3[2] = 1.0;
}


/* lattice vectors in Fourier space */
void reciprocal_vectors(double b1[3], double b2[3], double b3[3], double a1[3], double a2[3], double a3[3], int L)
{
    double norm = 2*M_PI/(L * dot(a1,cross(a2,a3)) );
    int i;

    for(i=0; i<3; i++)
    {
        b1[i] = norm*cross(a2,a3)[i];
        b2[i] = norm*cross(a3,a1)[i];
        b3[i] = norm*cross(a1,a2)[i];
    }
}

/* metric of coordinates under flow */
void get_metric(double metric[3][3], double a1[3], double a2[3], double a3[3])
{
    metric[0][0] = dot(a1,a1);
    metric[1][1] = dot(a2,a2);
    metric[2][2] = dot(a3,a3);

    metric[0][1] = dot(a1,a2);
    metric[1][2] = dot(a2,a3);
    metric[0][2] = dot(a3,a1);
}

/* transformation matrix from skewed to Cartesian coordinates */
void xyz_to_Cart_transform(double xyz_to_C[3][3], double a1[3], double a2[3], double a3[3])
{
    int i;

    for(i=0; i<3; i++)
    {
        xyz_to_C[i][0] = a1[i];
        xyz_to_C[i][1] = a2[i];
        xyz_to_C[i][2] = a3[i];
    }
}

/* transformation matrix from Cartesian to skewed coordinates */
void Cart_to_xyz_transform(double C_to_xyz[3][3], double a1[3], double a2[3], double a3[3])
{
    int i,j,k;

    double inv_det = 1.0/dot(a1,cross(a2,a3));

    for(i=0; i<3; i++)
    {
        j = mod_int(i+1,3);
        k = mod_int(i+2,3);

        C_to_xyz[0][i] = (a2[j]*a3[k] - a2[k]*a3[j])*inv_det;
        C_to_xyz[1][i] = (a3[j]*a1[k] - a3[k]*a1[j])*inv_det;
        C_to_xyz[2][i] = (a1[j]*a2[k] - a1[k]*a2[j])*inv_det;
    }
}


void transform_coords(double r[3], double R[3], double r_to_R[3][3])
{
    R[0] = r_to_R[0][0]*r[0] + r_to_R[0][1]*r[1] + r_to_R[0][2]*r[2];
    R[1] = r_to_R[1][0]*r[0] + r_to_R[1][1]*r[1] + r_to_R[1][2]*r[2];
    R[2] = r_to_R[2][0]*r[0] + r_to_R[2][1]*r[1] + r_to_R[2][2]*r[2];
}

/* use the metric to obtain the scalar length squared of a vector */
double find_dr2(double dr[3], double metric[3][3])
{
    double dr2;

    dr2 = metric[0][0]*dr[0]*dr[0] + metric[1][1]*dr[1]*dr[1] + metric[2][2]*dr[2]*dr[2];
    dr2 += 2*(metric[0][1]*dr[0]*dr[1] + metric[1][2]*dr[1]*dr[2] + metric[0][2]*dr[2]*dr[0]);

    return dr2;
}

/* use the metric to obtain the scalar product of two vectors */
double find_drdR(double dr[3], double dR[3], double metric[3][3])
{
    double drdR;

    drdR = metric[0][0]*dr[0]*dR[0] + metric[1][1]*dr[1]*dR[1] + metric[2][2]*dr[2]*dR[2];
    drdR += metric[0][1]*dr[0]*dR[1] + metric[1][2]*dr[1]*dR[2] + metric[0][2]*dr[2]*dR[0];
    drdR += metric[0][1]*dr[1]*dR[0] + metric[1][2]*dr[2]*dR[1] + metric[0][2]*dr[0]*dR[2];

    return drdR;
}

/* find the maximal extent of the wavelet, which is ellipsoidal in flow coordinates */
void get_wavelet_axes(double axes[3], double radius, double C_to_xyz[3][3])
{
    int i;
    /* trig functions of spherical angles with maximum coordinate value */
    double ctheta, cphi, sphi, stheta, arg;

    for(i=0; i<3; i++)
    {
        if(C_to_xyz[i][0]==0.0)
        {cphi = 0.0;
         sphi = 1.0;}
        else
        {
            arg = C_to_xyz[i][1]/C_to_xyz[i][0];
            cphi = pow(1.0 + arg*arg, -0.5);
            sphi = arg*cphi;
        }

        if(C_to_xyz[i][2]==0)
        {ctheta = 0.0;
         stheta = 1.0;}
        else
        {
            arg = (C_to_xyz[i][0]*C_to_xyz[i][0] + C_to_xyz[i][1]*C_to_xyz[i][1])/(C_to_xyz[i][2]*C_to_xyz[i][2]);
            ctheta = pow(1.0 + arg, -0.5);
            stheta = arg*ctheta;
        }

        axes[i] = C_to_xyz[i][0]*cphi*stheta + C_to_xyz[i][1]*sphi*stheta + C_to_xyz[i][2]*ctheta;
        axes[i] = fabs(axes[i])*radius;
    }
}

/* move particles in the flow coordinates when reset to t=0 basis vectors */
void remap_at_coord_reset(double xyz[][3], double xyzn[][3], int Ntot, double time,
                          System_parameters params, int Lbs, int contents[Lbs][Lbs][Lbs][params.n_box_max])
{
    int i,j,k,ipart;

    int particle_list[Ntot], inbox[3];
    for(i=1; i<=Ntot; i++)
    {particle_list[i-1] = i;}

    int anchorxyz[3] = {0,0,0}, spansxyz[3] = {Lbs-1,Lbs-1,Lbs-1};

    /* t=0 basis vectors and transform matrix */
    double a10[3], a20[3], a30[3], C_to_xyz0[3][3];
    lattice_vectors(a10, a20, a30, 0.0, params);
    Cart_to_xyz_transform(C_to_xyz0, a10, a20, a30);

    /* basis vectors and transform matrix before reset */
    double a1t[3], a2t[3], a3t[3], xyz_to_Ct[3][3];
    lattice_vectors(a1t, a2t, a3t, params.reset_time, params);
    xyz_to_Cart_transform(xyz_to_Ct, a1t, a2t, a3t);

    /* transform matrix to reset coord system */
    double xyz_to_xyz0[3][3];
    matrix_dot(C_to_xyz0, xyz_to_Ct, xyz_to_xyz0);

    for(i=0; i<Ntot; i++)
    {
        transform_coords(xyz[i], xyzn[i], xyz_to_xyz0);
        xyz[i][0] = xyzn[i][0];
        xyz[i][1] = xyzn[i][1];
        xyz[i][2] = xyzn[i][2];
    }

    /* empty spans region of contents */
    for(i=0; i<=Lbs-1; i++)
    {   for(j=0; j<=Lbs-1; j++)
        {   for(k=0; k<=Lbs-1; k++)
            {
                inbox[0] = mod_int(i,Lbs);
                inbox[1] = mod_int(j,Lbs);
                inbox[2] = mod_int(k,Lbs);

                ipart = 0;
                while(contents[inbox[0]][inbox[1]][inbox[2]][ipart]>0)
                {
                    contents[inbox[0]][inbox[1]][inbox[2]][ipart] = 0;
                    ipart++;
                }
            }
        }
    }

    /* remap movers */
    map(Ntot, particle_list, xyzn, params.L, params.sub_box_size, Lbs, params.n_box_max, contents,
        anchorxyz, spansxyz);
}


