#pragma once 

#include "Matrix.h"
#include "Multigrid.h"
#include "Function.h"

enum BoundaryConditionType{ Dirichlet, Neumann };

struct BoundaryCondition{
    BoundaryConditionType BCT;
    TimeFunction2D * f;

    BoundaryCondition(): BCT{Dirichlet}, f{nullptr} {}
};

class MultigridForFVM : public Multigrid{
    
    protected:
    VectorXd restriction( const VectorXd & vh ){

    }

    VectorXd prolongation( const VectorXd & v2h ){
        
    }
    
};

class ADESolver{

    public:

    void solve( ){

    }

    void setGrid( const int & N ){
        n = N;
        Vphi.resize(n*n);
        FphiX.resize(n*(n+1));
        FphiY.resize(n*(n+1));
        for( int i = 0; i < 4; i++ ) boundary[i].resize(2*n);
    }
    void setGrid( int && N ){
        n = std::move(N);
        Vphi.resize(n*n);
        FphiX.resize(n*(n+1));
        FphiY.resize(n*(n+1));
        for( int i = 0; i < 4; i++ ) boundary[i].resize(2*n);
    }

    void setDiffusivity( const int & Nu ){ nu = Nu; }
    void setDiffusivity( int && Nu ){ nu = std::move(Nu); }

    void setCourantNum( const int & CrN ){ Cr = CrN; }
    void setCourantNum( int && CrN ){ Cr = std::move(CrN); }

    ADESolver(){}
    ADESolver( const int & N, const double & Nu, const double & CrN ):n(N), nu(Nu), Cr(CrN) {}

    private:

    void setInitialVal( TimeFunction2D * phi ){

    }

    MultigridForFVM LinerSysSolver;
    int n;
    double nu;
    double Cr;
    VectorXd Vphi;
    VectorXd FphiX;
    VectorXd FphiY;
    double corner[4];
    VectorXd boundary[4];
    BoundaryCondition BC[4];


};