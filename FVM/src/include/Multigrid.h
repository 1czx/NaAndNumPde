#pragma once 

#include "Matrix.h"

class Multigrid{

    public:

    VectorXd solve( VectorXd f, const int & maxIter = 100, const double & eps = 1e-9 ){ 
        VectorXd v(f.size());
        for( int i = 0; i < maxIter; i++ ){
            VectorXd v1 = FMG(f,5,5,0);
            v = v + FMG(f,5,5,0);
            f = f - M[0]*v1;
            if( f.lInfNorm() < eps ) return v;
        }
        return v;
    }
    
    void setMatrix( const vector<SparseMatrix> & SPM ){
        M = SPM;
        maxLevel = M.size();
        // N = int(sqrt(M[0].colSize()));
    }

    protected:

    void relaxation( VectorXd & v, VectorXd f, int nu, const double & w , const int & level ){
        SparseMatrix & Ml = M[level];
        int n = Ml.rowSize();
        while( nu > 0 ){
            VectorXd v1(f);
            for( int i = 0; i < n; i++ ){
                const vector<SparseMatrix::element> & row = Ml.getRow(i);
                int aii = 0;
                for( const auto & e : row ){
                    if( e.idx != i ) v1(i) -= e.val*v(e.idx);
                    else aii = e.val;
                }
                v1(i) /= aii;
            }
            v = w*v1+(1-w)*v;
            nu--;
        }
    }
    
    virtual VectorXd restriction( const VectorXd & vh ) = 0;

    virtual VectorXd prolongation( const VectorXd & v2h ) = 0;
    
    void VC( VectorXd & v, VectorXd f, const int & nu1, const int & nu2, const int & level ){
        relaxation(v,f,nu1,2.0/3,level);
        if( level == maxLevel-1 ) relaxation(v,f,nu2,2.0/3,level);
        else{
            int n = N/(pow(2,level));
            VectorXd f2h = restriction(f-M[level]*v);
            VectorXd v2h;
            v2h = vector<double>((n/2)*(n/2),0);
            VC(v2h,f2h,nu1,nu2,level+1);
            v = v + prolongation(v2h);
            relaxation(v,f,nu2,2.0/3,level);
        }
    }

    VectorXd FMG( VectorXd f, const int & nu1, const int & nu2, const int & level ){
        VectorXd v;
            if( level == maxLevel-1 ){
                int n = f.size();
                v = vector<double>(n,0);
                VC(v,f,nu1,nu2,level);
            }else{
                VectorXd f2h = restriction(f);
                VectorXd v2h = FMG(f2h,nu1,nu2,level+1);
                v = prolongation(v2h);
                VC(v,f,nu1,nu2,level);
            }
            return v;
    }

    vector<SparseMatrix> M;
    int N;
    int maxLevel;

};