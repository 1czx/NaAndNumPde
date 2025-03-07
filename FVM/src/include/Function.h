#pragma once
#include "Matrix.h"

class TimeFunction2D{

    public:

    virtual double val( const double & x, const double & y, const double & t ) = 0;

    double operator()( const double & x, const double & y, const double & t ){
        return val(x,y,t);
    }

    double Int2D( const double & a, const double & b, const double & c, const double & d, const double & t ){
        double h = (d-c)/4;
        return 2.0*h*(7*(IntX(a,b,c,t)+IntX(a,b,d,t))+32*(IntX(a,b,c+h,t)+IntX(a,b,c+3*h,t))+12*IntX(a,b,c+2*h,t))/45;
    }

    double IntX( const double & a, const double & b, const double & y, const double & t ){
        double h = (b-a)/4;
        return 2.0*h*(7*(val(a,y,t)+val(b,y,t))+32*(val(a+h,y,t)+val(a+3*h,y,t))+12*val(a+2*h,y,t))/45;
    }

    double IntY( const double & x, const double & c, const double & d, const double & t ){
        double h = (d-c)/4;
        return 2.0*h*(7*(val(x,c,t)+val(x,d,t))+32*(val(x,c+h,t)+val(x,c+3*h,t))+12*val(x,c+2*h,t))/45;
    }

};