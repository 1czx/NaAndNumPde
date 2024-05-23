#include "ThreeBody.h"
VectorXd ThreeBody( VectorXd u, double t ){
    VectorXd v = u;
    double mu = 0.012277471;
    double u1 = u(0);
    double u2 = u(1);
    double u3 = u(2);
    v(0) = u(3);
    v(1) = u(4);
    v(2) = u(5);
    v(3) = 2*u(4) + u1 - mu*(u1+mu-1)/(pow(u2*u2+u3*u3+pow(u1+mu-1,2),3.0/2)) - (1-mu)*(u1+mu)/(pow(u2*u2+u3*u3+pow(u1+mu,2),3.0/2));
    v(4) = -2*u(3) + u2 - mu*u2/(pow(u2*u2+u3*u3+pow(u1+mu-1,2),3.0/2)) - (1-mu)*u2/(pow(u2*u2+u3*u3+pow(u1+mu,2),3.0/2));
    v(5) = -mu*u3/(pow(u2*u2+u3*u3+pow(u1+mu-1,2),3.0/2)) - (1-mu)*u3/(pow(u2*u2+u3*u3+pow(u1+mu,2),3.0/2));
    return v;
}
