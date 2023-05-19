#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "Complex.h"
#include "sac.h"
#include "sacio.h"
#include "grd3d.h"
#define MAX 200

// =======================================
// This code is used to calculate the theoretical P-to-Ps transmission coefficient, 
// based on the established formulas in Aki & Richards (2002) Eqn 5.39 and 5.40

void RSVRTmatrix(float p, float m2[3], float m1[3], float *T){
    // p = ray parameter 
    // m2 = matrix of incident wave [vp vs rho]
    // m1 = matrix of transmitted wave [vp vs rho]

    // vertical slowness
    float f_ai, f_at, f_bi, f_bt;
    complex fai_cplx, fat_cplx, fbi_cplx, fbt_cplx;
    float alpha1, beta1, alpha2, beta2, rho1, rho2;
    alpha1 = m1[0];
    alpha2 = m2[0];
    beta1 = m1[1];
    beta2 = m2[1];
    rho1 = m1[2];
    rho2 = m2[2];

    if (p > 1/alpha2){
        fai_cplx = cmplx(0, -1*sqrt(p * p - 1/pow(alpha2,2)));
    }else{
        fai_cplx = cmplx(sqrt(1/pow(alpha2,2) - p * p), 0);
    }

    if (p > 1/alpha1){
        fat_cplx = cmplx(0, -1*sqrt(p * p - 1/pow(alpha1,2)));
    }else{
        fat_cplx = cmplx(sqrt(1/pow(alpha1,2) - p * p), 0);
    }    
    
    if (p > 1/beta2){
        fbi_cplx = cmplx(0, -1*sqrt(p * p - 1/pow(beta2,2)));
    }else{
        fbi_cplx = cmplx(sqrt(1/pow(beta2,2) - p * p), 0);
    }

    if (p > 1/beta1){
        fbt_cplx = cmplx(0, -1*sqrt(p * p - 1/pow(beta1,2)));
    }else{
        fbt_cplx = cmplx(sqrt(1/pow(beta1,2) - p * p), 0);
    }

    float a, b, c, d;
    a = rho2 * (1 - 2 * pow(beta2,2) * p * p) - rho1 * (1 - 2*pow(beta1,2)*p * p);
    b = rho2 * (1 - 2*pow(beta2,2)*p * p) + 2*rho1*pow(beta1,2) * p * p;
    c = rho1 * (1 - 2 * pow(beta1,2)*p * p) + 2 * rho2 * pow(beta2,2) * p * p;
    d = 2*(rho2 * pow(beta2,2) - rho1 * pow(beta1,2));

    complex E_cplx, F_cplx, G_cplx, H_cplx, D_cplx, D_conj;
    E_cplx = cplus(dmltp(b, fat_cplx), dmltp(c, fai_cplx));
    F_cplx = cplus(dmltp(b, fbt_cplx), dmltp(c, fbi_cplx));
    G_cplx = cplus(cmplx(a,0), cngtv(dmltp(d, cmltp(fat_cplx, fbi_cplx))));
    H_cplx = cplus(cmplx(a, 0), cngtv(dmltp(d, cmltp(fai_cplx, fbt_cplx))));
    D_cplx = cplus(cmltp(E_cplx, F_cplx), dmltp(p*p ,cmltp(G_cplx, H_cplx)));
    D_conj = conjg(D_cplx);
    float D;
    D = cmltp(D_cplx, D_conj).x;
    complex Tp;

    Tp = dmltp(-2 * rho2 * p * alpha2/beta1, dmltp(1/D, cmltp(cmltp(fai_cplx, G_cplx), D_conj)));
    *T = Tp.x;
}