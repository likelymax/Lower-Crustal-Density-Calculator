#define MAX 150
#define RDO 1638400
#ifndef _TRAN_
    #define _TRAN_
    float *Trancoeff(float vel_p[MAX], float vel_vpvs[MAX], int end, float ratio, float diff_range, float rho1, float rho2, float limit, float inv, float r[RDO], float p);
#endif 