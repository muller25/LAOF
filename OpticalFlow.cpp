#include "OpticalFlow.h"

double OpticalFlow::psi(double s, double e)
{
    return sqrt(s + e);
}

double OpticalFlow::psi_d(double s, double e)
{
    return 0.5 / sqrt(s + e);
}

double OpticalFlow::phi(double s, double e)
{
    return sqrt(s + e);
}

double OpticalFlow::phi_d(double s, double e)
{
    return 0.5 / sqrt(s + e);
}



