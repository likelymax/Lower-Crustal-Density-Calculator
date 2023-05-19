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
#define max 150
#define RDO 1638400

// ========================================
// this code is used to compare the theoretical P-to-Ps trancoeff with the actual trancoeff

int MAXIMUM(float r[RDO]){
    int i, out;
    float tmp_r = 0;
    for (i = 0; i< RDO; i++){
        if (r[i] > tmp_r){
            out = i;
            tmp_r = r[i];
        }
    }
    return out;
}
float *Trancoeff(float vel_p[max], float vel_vpvs[max], int end,  float ratio, float diff_range, float rho1, float rho2, float limit, float inv, float r[RDO], float p)
{
    float vel_s[max], mi[3], mt[3];
    static float output[7];
    static float density[RDO], trans[RDO];
    float T,  tmp_T = 0, diff, tmp_diff = diff_range;
    int i = 0, out, judge = 0;
    FILE *ft;
    output[0] = output[1] = output[2] = output[3]= output[4] = output[5] = output[6] = 0;
    for (i = 0; i< end; i++){
        vel_s[i] = vel_p[i]/vel_vpvs[i];
    }

    // input the lcvp lcvs umvp umvs into the mi matrix
    mi[0] = vel_p[i-2];
    mi[1] = vel_s[i-2];
    mt[0] = vel_p[i-3];
    mt[1] = vel_s[i-3];
    mi[2] = rho2;
    i = 0;

    // run the loop of rho1 according to the input range. here I used the range from 2600 to 3500 kg/m^3.
    // Once the difference between the theoritcal and actual trancoeff is smaller than the threshold (input)
    // the loop will be stopped and output the rho1
    while (rho1 < limit){
        mt[2] = rho1;
        RSVRTmatrix(p, mi, mt, &T);
        fprintf(ft, "%f %f\n", rho1, T);
        if ( T < 0 ){
            rho1 += inv;
            continue;
        }else{
            if (tmp_T == 0){
                diff = fabs(T - ratio);
                if (diff <= diff_range){
                    output[0] = T;
                    output[1] = rho1;
                    tmp_diff = diff;
                    output[2] = tmp_diff;
                }
            }else if (T != tmp_T) {
                if (slope(T, tmp_T, inv) < 0 ){
                    diff = fabs(T - ratio);
                    if (diff <= diff_range && diff <= tmp_diff){
                        output[0] = T;
                        output[1] = rho1;
                        tmp_diff = diff;
                        output[2] = tmp_diff;
                    }
                }else if (slope(T, tmp_T, inv) >= 0 ){
                    diff = fabs(T - ratio);
                    if (diff <= diff_range && diff <= tmp_diff){
                        output[3] = T;
                        output[4] = rho1;
                        tmp_diff = diff;
                        output[5] = tmp_diff;
                    }              
                }
            }
            tmp_T = T;  
            density[i] = rho1;
            trans[i] = T;
            rho1 += inv;
            i++;
        }
    }
    out = MAXIMUM(trans);
    output[6] = out;
    fclose(ft);
    return output;
    done:
        output[0] = -1000;
        output[1] = -1000;
        output[2] = -1000;
        output[3] = -1000;
        output[4] = -1000;
        output[5] = -1000;
        output[6] = -1000;
        return output;
}
