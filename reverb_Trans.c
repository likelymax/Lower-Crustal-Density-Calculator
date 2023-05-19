// ===================================
// Lower crustal density Calculator
// ===================================
// This is a tool for calculating lower crustal density by comparing theoretical and 
// actual P-to-Ps transmission coefficients. The theoretical coefficient is derived from
// Aki & Richards (2002) Eqns 5.39 and 5.40, while the actual coefficient is sourced 
// from receiver functions.

// Usage:
// To use this code, execute the following command:
// ./reverb_Trans -V$vel -C$ratio -D$diff -R$limit/$mdensity/$bound/$inv -S$file -Ooutput.txt

// Parameters
// The script accepts the following parameters:
// -V : The input velocity model.
// -C : The actual P-to-Ps transmission coefficient.
// -D : The threshold that can be used to decide the output.
// -R : Lower band of lower crustal density / Upper mantle density / The upper bound of lower crustal density / Interval.
// -S : Input receiver function.
// -O : Output result.

// Output
// The output (output.txt) is the estimated lower crustal density.

//Dependencies
// This tool is designed to run on Unix-like systems and requires a compatible environment.

//Contact
// Please reach out if you have any questions or encounter issues using the tool. Contributions are also welcome.
// You can reach me via email at yitanwang@ufl.edu


#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "sac.h"
#include "sacio.h"
#include "grd3d.h"
#include "Complex.h"
#include "RSVRTmatrix.h"
#include "Tran.h"
#define MAX 150

// double gaussrand(float sigma, float mu)
// {
//     static double U, V;
//     static int phase = 0;
//     double z;
//     if(phase == 0)
//     {
//         U = rand() / (RAND_MAX + 1.0);
//         V = rand() / (RAND_MAX + 1.0);
//         z = sqrt(-2.0 * log(U))* sin(2.0 * PI * V);
//     }
//     else
//     {
//     z = sqrt(-2.0 * log(U)) * cos(2.0 * PI * V);
//     }

//     phase = 1 - phase;
//     z = mu + (z*sigma);
//     return z;
// }

int main(int argc, char **argv){
    int i, error = 0, dir, tot;
    char sac[128], output[128],velocity[128], longitude[128], latitude[128];
    float az, lon, lat, diff_range,rho2, inv, rho1, limit, shallow, tmp_T = 0, deep, lon_min, lon_max, lat_min, lat_max, ratio;
    float mean_vp, mean_vs, std_vp, std_vs;
    char line[128] = {0};
    int type, round, certainty;
    SACHEAD	hd;
    FILE *fo, *fa, *fv, *fn, *fr, *fs, *ft;
    for (i=1; !error && i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {
                case 'V':
                    strcpy(velocity,&argv[i][2]);
                    fv = fopen(velocity,"r");       
                break;
                case 'C':
                    sscanf(&argv[i][2],"%f",&ratio);
                break;
                case 'D':
                    sscanf(&argv[i][2],"%f",&diff_range);            
                break;
                case 'R':
                    sscanf(&argv[i][2],"%f/%f/%f/%f",&rho1, &rho2, &limit, &inv);            
                break;
                case 'S':
                    strcpy(sac,&argv[i][2]);
                    fs = fopen(sac,"r");                 
                break;
                case 'O':
                        strcpy(output,&argv[i][2]);
                        fo = fopen(output,"w");
                break;
            }
        }
    }
    float *r;
    float vel_p[MAX], vel_vpvs[MAX];
    float (*output_T)[7];
    float p, dp, vp, vpvs, low_mantle;
    int incid_p = 0, incid_s = 0, trans_p = 0, trans_s = 0;
    low_mantle = 3500;
    if (ratio <= 0.0){
        printf("wrong actual trancoeff\n");
        return -1;
    } 
    if (( r=read_sac(sac,&hd)) == NULL) { 
        printf("can't open %s\n", sac);
        return -1;
    }
    // read the information of the stations and events
    p = hd.user0;
    lon = hd.stlo;
    lat = hd.stla;
    i = 0;
    // read velocity model
    while (fscanf(fv, "%f %f %f\n", &dp, &vp, &vpvs) == 3){
        vel_p[i] = vp;
        vel_vpvs[i] = vpvs;
        i++;
    }
    float ve_sinc = vel_p[i-2]/vel_vpvs[i-2];
    float ve_stra = vel_p[i-3]/vel_vpvs[i-3];
    // if (p > 1/vel_p[i-2] && p < 1/ve_sinc) incid_p = 1;
    // if (p > 1/ve_sinc){
    //     incid_p = 1;
    //     incid_s = 1;
    // }
    // if (p > 1/vel_p[i-3] && p < 1/ve_stra) trans_p = 1;
    // if (p > 1/ve_stra){
    //     trans_p = 1;
    //     trans_s = 1;
    // }
    if (p > 1/vel_p[i-2] && p < 1/ve_sinc){
        printf("wrong incident P\n");
        return -1;
    }
    if (p > 1/ve_sinc){
        printf("wrong incident angles\n");
        return -1;
    }
    if (p > 1/vel_p[i-3] && p < 1/ve_stra){
        printf("wrong incident P\n");
        return -1;
    }
    if (p > 1/ve_stra){
        printf("wrong incident angles\n");
        return -1;
    }
    output_T = Trancoeff(vel_p, vel_vpvs, i, ratio, diff_range, rho1, rho2, limit, inv, r, p);
    float dif_tmp, den1, T1, den2, T2, mantle = rho2, max_dense, output_mantle, mantle_high, diff1, diff2;
    T1 = (*output_T)[0];
    den1 = (*output_T)[1];
    // diff1 = (*output_T)[2];
    // T2 = (*output_T)[3]; 
    // den2 = (*output_T)[4];
    // diff2 = (*output_T)[5];
    // max_dense = (*output_T)[6];
    // mantle_high = 3700;
    // mantle=rho2;
    // float dif_tmp1;
    fprintf(fo,"%f\n", den1);
 
    // fprintf(fo,"%f %f %f %f %f %f %d %d %d %d %f %f\n", den1, T1, den2, T2, output_mantle, p, incid_p, incid_s, trans_p, trans_s, diff1, diff2);
    return 0;
}