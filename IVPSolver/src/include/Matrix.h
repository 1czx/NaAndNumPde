#pragma once

#include<lapacke.h>
#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<algorithm>

using std::vector;
using std::string;
using std::cout;
using std::cin;
using std::ofstream;
using std::ifstream;
using std::function;
using std::endl;
using std::max;
using std::min;
using std::abs;

class VectorXd{
    public:
    
    VectorXd():n{0} {}
    VectorXd( const int & N ):n{N}, a(N,0) {}
    VectorXd( const vector<double> & x ):a{x}, n{int(x.size())} {}

    double * data(){ return a.data(); }


    double operator() ( const int & i ) const {
        return a[i];
    }

    double & operator() ( const int & i ){
        return a[i];
    }

    int size() const {
        return n;
    }

    void resize( const int & nn ){
        a.resize(nn);
        n = nn;
    }

    double l1Norm(){
        double sum = 0;
        for( int i = 0; i < n; i++ ) sum += abs(a[i]);
        return sum;
    }

    double l2Norm(){
        double sum = 0;
        for( int i = 0; i < n; i++ ) sum += pow(a[i],2);
        return sqrt(sum);
    }
    double lInfNorm(){
        double Max = abs(a[0]);
        for( int i = 1; i < n; i++ ) Max = max(Max,abs(a[i]));
        return Max;
    }

    VectorXd Abs(){
        vector<double> temp(n,0);
        for( int i = 0 ; i < n; i++ ) temp[i] = abs(a[i]);
        return temp; 
    }

    VectorXd operator+( VectorXd b ){
        if( b.size() != this->n ){
            cout << "Error: The sizes of two operated vectors are different. \n";
            return VectorXd();
        }
        for( int i = 0; i < n; i++ ) b(i) += (*this)(i);
        return b;
    }

    VectorXd operator-() const {
        VectorXd temp = *this;
        for( int i = 0; i < n; i++ ) temp(i) = -temp(i);
        return temp;
    }
    
    VectorXd operator*( const double & a ){
        VectorXd temp = *this;
        for( int i = 0; i < n; i++ ) temp(i) *= a;
        return temp;
    }

    VectorXd operator-( VectorXd b ){
        return *this+(-b);
    }

    private:
    
    vector<double> a;
    int n;
    
    friend class MatrixXd;
};

class MatrixXd{
    public:

    MatrixXd( int M, int N ):m{M}, n{N} {
        a = vector<double>(m*n,0);
    }

    double * data(){ return a.data(); }
    
    double operator() ( const int & i, const int & j ) const {
        return a[i*n+j];
    }

    double & operator() ( const int & i, const int & j ){
        return a[i*n+j];
    }

    VectorXd solve( VectorXd b ){
        if( n != b.size() || m != n ) return VectorXd();
        double * mat = a.data();
        double * vec = b.data();
        int ipiv[n];
        int nrhs = 1, info;
        info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,n,1,mat,n,ipiv,vec,1);
        return b;
    }

    private:

    vector<double> a;
    int m;//size of row of martix
    int n;//size of column of martix

};


VectorXd operator*( const double & a, VectorXd b ){
    int n = b.size();
    for( int i = 0; i < n; i++ ) b(i) *= a;
    return b;
}

