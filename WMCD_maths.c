/*
 * Functions for standard mathematical operations in WMCD.
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

/* split time in seconds into hours, minutes and seconds */
void split_time(double t_sec, time_hms *t_hms)
{
    double t_reduced;

    t_hms->hours = (int)(t_sec/3600);

    t_reduced = t_sec - 3600*(t_hms->hours);

    t_hms->minutes = (int)(t_reduced/60);
    t_hms->seconds = t_reduced - 60*(t_hms->minutes);
}

/* get the sign of a number */
int signum(double x)
{return ((x > 0) ? 1 : ((x < 0) ? -1 : 0));}


/* dot product between 2 vectors */
double dot(double vector1[3], double vector2[3])
{
    double xx,yy,zz;
    /* contributions to dot product */
    xx = vector1[0]*vector2[0];
    yy = vector1[1]*vector2[1];
    zz = vector1[2]*vector2[2];

    return (xx+yy+zz);
}

/* dot product between 3x3 matrices */
void matrix_dot(double A[3][3], double B[3][3], double C[3][3])
{
    int i, j;

    for(i=0; i<3; i++)
    {for(j=0; j<3; j++)
     {C[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j];}
    }
}


/* cross product between 2 vectors */
double *cross(double vector1[3], double vector2[3])
{
    static double vector3[3];

    vector3[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
    vector3[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
    vector3[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];

    return vector3;
}

/* arbitrary vector transverse to a given vector */
void transverse_vector(double vector[3], double trans_vector[3])
{
    /* use khat to get displacement vector direction */
    /* first calc trig for vector (with spherical polar angles (theta,phi)) */
    double costheta = vector[2];
    double sintheta = sqrt(1-costheta*costheta);
    double cosphi, sinphi, randnum;
    /* avoid problems with phi being undefined */
    if( (costheta == 1.0) || (costheta == -1.0) )
    {cosphi = 1.0;
     sinphi = 0.0;}
    else
    {cosphi = (vector[0])/sintheta;
     sinphi = (vector[1])/sintheta;}

    /* next generate an angle in [0,2pi) in the plane perpendicular to khat */
    randnum = 2*M_PI*POW2_32_INV_OPEN*pcg32_random();

    /* finally calculate ehat */
    trans_vector[0] = -sinphi*cos(randnum) - costheta*cosphi*sin(randnum);
    trans_vector[1] =  cosphi*cos(randnum) - costheta*sinphi*sin(randnum);
    trans_vector[2] =  sintheta*sin(randnum);

}

/* nearest integer functions */
int nint(double number)
{return (int)floor(number+0.5);}

long long_nint(double number)
{return (long)floor(number+0.5);}

/* reliable modulus operator for floating point variables */
double mod_lf(double x, double y)
{return (x - y*floor(x/y));}


/* reliable modulus operator for integer variables */
int mod_int(int x, int y)
{return (int)(x - y*floor((float)x/y));}

/* reliable modulus operator for long integer variables */
long mod_long(long x, long y)
{return (long)(x - y*floor((float)x/y));}

/* min function: gives smallest of the two doubles */
double min_lf(double x, double y)
{
    if(x<y)
    {return x;}
    else
    {return y;}
}

/* max function: gives smallest of the two doubles */
double max_lf(double x, double y)
{
    if(x>y)
    {return x;}
    else
    {return y;}
}


/* efficient pow. use whenever a (cumbersome) POSITIVE INTEGER power is required */
double pow_int(double x, int y)
{
    int i;
    double xpowy=x;

    for(i=1; i<y; i++)
    {xpowy *= x;}

    return xpowy;
}


/* find unit vector in a random direction */
void direction(double *vector)
{
    /* variable for normalising to unit vectors */
    double norm;

    /* initialise outside unit sphere */
    vector[0] = 1.0;
    vector[1] = 1.0;
    vector[2] = 1.0;

    /* generate random vector within init sphere to ensure all directions are equally likely */
    while(dot(vector,vector)>1)
    {
        vector[0] = 2*(POW2_32_INV*pcg32_random()) - 1;
        vector[1] = 2*(POW2_32_INV*pcg32_random()) - 1;
        vector[2] = 2*(POW2_32_INV*pcg32_random()) - 1;
    }

    /* rescale to be a unit vector */
    norm = 1/sqrt(dot(vector,vector));
    vector[0] *= norm;
    vector[1] *= norm;
    vector[2] *= norm;
}

/* Generate a random number with a Gaussian distribution with mean 0 and standard deviation 1
 * using the Box-Muller method */
void Gaussian_rand(double *vector)
{
    double U, V;

    U = POW2_32_INV_OPEN*(1.0 + pcg32_random());
    V = POW2_32_INV_OPEN*(1.0 + pcg32_random());

    vector[0] = sqrt(-2.0 * log(U)) * cos(2*M_PI*V);
    vector[1] = sqrt(-2.0 * log(U)) * sin(2*M_PI*V);

    U = POW2_32_INV_OPEN*(1.0 + pcg32_random());
    V = POW2_32_INV_OPEN*(1.0 + pcg32_random());

    vector[2] = sqrt(-2.0 * log(U)) * cos(2*M_PI*V);
}


/* function numerically integrating cos(t)/t from x to 2x
 * this is equivalent to CosIntegral(x)-CosIntegral(2x)
 */
double Ci_x_2x(double x)
{
    int Nstrips=1000, i;
    double width=(x/Nstrips), Integral=0.0;

    /* add value at 'left' of first strip */
    Integral += cos(x)/(2*x);
    /* add twice the values at strip boundaries inside 'bulk' */
    for(i=1; i<Nstrips; i++)
    {Integral += cos(x+i*width)/(x+i*width);}
    /* add value at 'right' of last strip */
    Integral += cos(2*x)/(4*x);

    Integral *= -width;

    return Integral;
}

