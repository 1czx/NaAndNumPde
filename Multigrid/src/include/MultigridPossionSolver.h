#pragma once

#include"Matrix.h"
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<algorithm>
#include<cmath>
#include<functional>
#include"json/json.h"
#include<chrono>
#include<time.h>
#include<iomanip>

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

enum BCT{Dirichlet,Nuemann};
enum restrictionType{FullWeight,Injection};
enum prolongationType{linear,quadratic};
enum cycleType{Vcycle,FullMG};

template<int dim>
struct userInput{

};

template<>
struct userInput<2>{
    vector<BCT> boundaryConditonType;//in the order of bottom, right, top, left, circle
    restrictionType rt;
    prolongationType pt;
    cycleType ct;
    int maxIter;
    double eps;
    string initialGauss;

    function<double(const double &,const double &)> f;//RHS
    function<double(const double &,const double &)> u;//Dirichlet boundary condition
    function<double(const double &,const double &)> gx;//Nuemann condition for vertical boundary 
    function<double(const double &,const double &)> gy;//Nuemann condition for horizontal boundary
    string outputfile;

    void readFile( const string & filename ){
        Json::Reader reader;
	    Json::Value root;
        ifstream in(filename, std::ios::binary);
        if( !in.is_open() ){
	    	cout << "Error opening file\n";
	    	return;
	    }
        if( reader.parse(in,root) ){
            if( root["restriction"].asString() == "FullWeight" ) rt = FullWeight;
            else rt = Injection;
            if( root["prolongation"].asString() == "linear" ) pt = linear;
            else pt = quadratic;
            if( root["cycles"].asString() == "VC" ) ct = Vcycle;
            else ct = FullMG;
            maxIter = root["maxIter"].asInt();
            eps = root["RA"].asDouble();
            initialGauss = root["initialGauss"].asString();
            outputfile = root["opfile"].asString();
            Json::Value BC = root["BoundaryCondition"];
            boundaryConditonType.clear();
            for( int i = 0; i < BC.size(); i++ ){
                if( BC[i].asString() == "Dirichlet" ) boundaryConditonType.push_back(Dirichlet);
                else boundaryConditonType.push_back(Nuemann);
            }
        }else{
            cout << "Parse Error\n";
        }
        in.close();
    }
    userInput() {}
    userInput(
        function<double(const double &,const double &)> _f,
        function<double(const double &,const double &)> _u,
        function<double(const double &,const double &)> _gx,
        function<double(const double &,const double &)> _gy
    ): f{_f}, u{_u}, gx{_gx}, gy{_gy} {}
};

template<>
struct userInput<1>{
    vector<BCT> boundaryConditonType;//in the order of bottom, right, top, left, circle
    restrictionType rt;
    prolongationType pt;
    cycleType ct;
    int maxIter;
    double eps;
    string initialGauss;

    function<double(const double &)> f;//RHS
    function<double(const double &)> u;//Dirichlet boundary condition
    function<double(const double &)> g;//Nuemann boundary condition
    string outputfile;

    void readFile( const string & filename ){
        Json::Reader reader;
	    Json::Value root;
        ifstream in(filename, std::ios::binary);
        if( !in.is_open() ){
	    	cout << "Error opening file\n";
	    	return;
	    }
        if( reader.parse(in,root) ){
            if( root["restriction"].asString() == "FullWeight" ) rt = FullWeight;
            else rt = Injection;
            if( root["prolongation"].asString() == "linear" ) pt = linear;
            else pt = quadratic;
            if( root["cycles"].asString() == "VC" ) ct = Vcycle;
            else ct = FullMG;
            maxIter = root["maxIter"].asInt();
            eps = root["RA"].asDouble();
            initialGauss = root["initialGauss"].asString();
            outputfile = root["opfile"].asString();
            Json::Value BC = root["BoundaryCondition"];
            boundaryConditonType.clear();
            for( int i = 0; i < BC.size(); i++ ){
                if( BC[i].asString() == "Dirichlet" ) boundaryConditonType.push_back(Dirichlet);
                else boundaryConditonType.push_back(Nuemann);
            }
        }else{
            cout << "Parse Error\n";
        }
        in.close();
    }
    
    userInput() {}
    userInput(
        function<double(const double &)> _f,
        function<double(const double &)> _u,
        function<double(const double &)> _g
    ): f{_f}, u{_u}, g{_g} {}
};

template<int dim>
struct Result{
        VectorXd U;//numerical solution
        VectorXd U_cap;//ture solution
        VectorXd E;//abs error |U-U_cap|
        int iter;//the number of iteration
        double time;//
        int n;
        double L1NormErr(){
            double sum = 0;
            sum = E.l1Norm();
            sum = sum/pow(n,dim);
            return sum;
        }
        double L2NormErr(){
            double sum = 0;
            sum = E.l2Norm();
            sum = sum/pow(n,dim/2);
            return sum;
        }
        double LInfNormErr(){
            return E.lInfNorm();
        }
        void outputU( const string & filename ){
            ofstream fout(filename);
            double h = 1.0/n;
            if( dim == 1 ){
                for( int i = 0; i < n+1; i++ ){
                    fout << i*h << " " << U(i) << endl;
                }
            }else{
                for( int j = 0; j < n+1; j++ ){
                    for( int i = 0; i < n+1; i++ ){
                        fout << U(j*(n+1)+i) << " ";
                    }
                    fout << endl;
                }
            }
            fout.close();
        }
        void outputE( const string & filename ){
            ofstream fout(filename);
            double h = 1.0/n;
            if( dim == 1 ){
                for( int i = 0; i < n+1; i++ ){
                    fout << i*h << " " << E(i) << endl;
                }
            }else{
                for( int j = 0; j < n+1; j++ ){
                    for( int i = 0; i < n+1; i++ ){
                        fout << E(j*(n+1)+i) << " ";
                    }
                    fout << endl;
                }
            }
            fout.close();
        }
};

template<int dim>
class MultigridPossionSolver{
    public:
        MultigridPossionSolver():N{0}, maxLevel{1} {}
        void setUserInput( const userInput<dim> & uI ){ 
            input = uI;
            // inputValid = 1;
        }

        void setGridSpacing( const int & n ){ N = n; }

        Result<dim> getResult(){ return result; }

        void solve( ofstream & fout ){
            int s;
            if( dim == 1 ) s = N+1;
            else s = (N+1)*(N+1);
            VectorXd v = vector<double>(s,0);
            if( input.initialGauss != "ZERO" ){
                ifstream fin(input.initialGauss);
                for( int i = 0; i < s; i++ ) fin >> v(i);
                fin.close();
            }  
            maxLevel = floor(log2(N));
            VectorXd f = getRhsF();
            auto start = std::chrono::high_resolution_clock::now();
            for( int i = 0; i < input.maxIter; i++ ){
                VectorXd v1 = vector<double>(s,0);
                if( input.ct == FullMG ) v1 = FMG(f,5,5,1);
                else v1 = VC(v1,f,5,5,1);
                // relaxation(v,f,1,2.0/3,1);
                v = v + v1;
                f = f - MMV(v1,1);
                fout << f.lInfNorm() << " ";
                if( f.lInfNorm() < input.eps ){
                    result.iter = i+1;
                    break;
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> tm = end - start;
            result.time = tm.count();
            result.n = N;
            result.U = v;
            result.U_cap = v;
            getTureSolution();
            result.E = result.U_cap - result.U;
            result.E = result.E.Abs();
        }

        void solve( ){
            int s;
            if( dim == 1 ) s = N+1;
            else s = (N+1)*(N+1);
            VectorXd v = vector<double>(s,0);
            if( input.initialGauss != "ZERO" ){
                ifstream fin(input.initialGauss);
                for( int i = 0; i < s; i++ ) fin >> v(i);
                fin.close();
            }  
            maxLevel = floor(log2(N));
            VectorXd f = getRhsF();
            auto start = std::chrono::high_resolution_clock::now();
            for( int i = 0; i < input.maxIter; i++ ){
                VectorXd v1 = vector<double>(s,0);
                if( input.ct == FullMG ) v1 = FMG(f,5,5,1);
                else v1 = VC(v1,f,5,5,1);
                // relaxation(v,f,1,2.0/3,1);
                v = v + v1;
                f = f - MMV(v1,1);
                if( f.lInfNorm() < input.eps ){
                    result.iter = i+1;
                    cout << std::setw(12) << std::left << f.lInfNorm();
                    break;
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> tm = end - start;
            result.time = tm.count();
            result.n = N;
            result.U = v;
            result.U_cap = v;
            getTureSolution();
            result.E = result.U_cap - result.U;
            result.E = result.E.Abs();
        }

    private:
        void getTureSolution()requires(dim == 1){
            double h = 1.0/N;
            for( int i = 0; i <= N; i++ ) result.U_cap(i) = input.u(i*h);
        }

        void getTureSolution()requires(dim == 2){
            double h = 1.0/N;
            for( int j = 0; j <= N; j++ ){
                for(int i = 0; i <= N; i++ ) result.U_cap(j*(N+1)+i) = input.u(i*h,j*h);
            }
        }

        VectorXd getRhsF()requires(dim == 1){
            VectorXd f;
            double h = 1.0/N;
            f=vector<double>(N+1,0);
            for( int i = 1; i < N; i++ ) f(i) = input.f(i*h);
            if( input.boundaryConditonType[0] == Nuemann ) f(0) = input.g(0)-h*input.f(0)/2;
            else f(0) = input.u(0);
            if( input.boundaryConditonType[1] == Nuemann ) f(N) = input.g(1)+h*input.f(1)/2;
            else f(N) = input.u(1);
            return f;
        }

        VectorXd getRhsF()requires(dim == 2){
            VectorXd f;
            double h = 1.0/N;
            double h2 = h*h;
            f=vector<double>((N+1)*(N+1),0);
            f(0) = input.u(0,0);
            f(N) = input.u(1,0);
            f(N*(N+1)) = input.u(0,1);
            f(N*(N+2))= input.u(1,1);
            for( int j = 1; j < N; j++ ){
                for( int i = 1; i < N; i++ ){
                    f(j*(N+1)+i) = input.f(i*h,j*h);
                }
            }
            int NuemannCount = 0;
            //bottom boundary
            if( input.boundaryConditonType[0] == Dirichlet ){
                for( int i = 1; i < N; i++ ){
                    f(i) = input.u(i*h,0);
                }
            }else{
                NuemannCount++;
                for( int i = 1; i < N; i++ ){
                    f(i) = h*input.f(i*h,0)-2*input.gy(i*h,0);
                }
                // if(input.boundaryConditonType[1] == Nuemann ){
                //     f(N) = h*input.f(1,0)-2*(input.gy(1,0)-input.gx(1,0));
                // }
            }
            //right boundary
            if( input.boundaryConditonType[1] == Dirichlet ){
                for( int j = 1; j < N; j++ ){
                    f(j*(N+1)+N) = input.u(1,j*h);
                }
            }else{
                NuemannCount++;
                for( int j = 1; j < N; j++ ){
                    f(j*(N+1)+N) = h*input.f(1,j*h)+2*input.gx(1,j*h);
                }
            //     if(input.boundaryConditonType[2] == Nuemann ){
            //         f(N*(N+2)) = h*input.f(1,1)+2*(input.gy(1,1)+input.gx(1,1));
            //    }
            }
            //top boundary
            if( input.boundaryConditonType[2] == Dirichlet ){
                for( int i = 1; i < N; i++ ){
                    f(N*(N+1)+i) = input.u(i*h,1);
                }
            }else{
                NuemannCount++;
                for( int i = 1; i < N; i++ ){
                    f(N*(N+1)+i) = h*input.f(i*h,1)+2*input.gy(i*h,1);
                }
                // if(input.boundaryConditonType[3] == Nuemann ){
                //     f(N*(N+1)) = h*input.f(0,1)+2*(input.gy(0,1)-input.gx(0,1));
                // }
            }
            //left boundary
            if( input.boundaryConditonType[3] == Dirichlet ){
                for( int j = 1; j < N; j++ ){
                    f(j*(N+1)) = input.u(0,j*h);
                }
            }else{
                for( int j = 1; j < N; j++ ){
                    f(j*(N+1)) = h*input.f(0,j*h)-2*input.gx(0,j*h);
                }
                // if(input.boundaryConditonType[0] == Nuemann ){
                //     if( NuemannCount == 3 ){
                //         f(0) = input.u(0,0);
                //         return f;
                //     }
                //     f(0) = h*input.f(0,0)-2*(input.gy(0,0)+input.gx(0,0));
                // }
            }
            return f;
        }

        //matrix multiply vector specialized for multigrid to slove possion
        VectorXd MMV( const VectorXd & v, const int & level ){
            int n = N/(pow(2,level-1));
            double h = 1.0/n;
            double h2 = h*h;
            VectorXd v1 = v;
            switch( dim ){
                case 1:{
                    for( int i = 1; i < n; i++ ) v1(i) = (2*v(i)-v(i-1)-v(i+1))/h2;
                    if( input.boundaryConditonType[0] == Nuemann ) v1(0) = (v(1)-v(0))/h;
                    else v1(0) = v(0);
                    if( input.boundaryConditonType[1] == Nuemann ) v1(n) = (v(n)-v(n-1))/h;
                    else v1(n) = v(n);
                    break;
                }
                case 2:{
                    for( int j = 1; j < n; j++ ){
                        for( int i = 1; i < n; i++ ){
                            v1(j*(n+1)+i) = (4*v(j*(n+1)+i)-(v(j*(n+1)+i-1)+v(j*(n+1)+i+1)+v((j+1)*(n+1)+i)+v((j-1)*(n+1)+i)))/h2;
                        }
                    }

                    v1(0) = v(0);
                    v1(n) = v(n);
                    v1(n*(n+1)) = v(n*(n+1));
                    v1(n*(n+2))= v(n*(n+2));

                    if( input.boundaryConditonType[0] == Dirichlet ) for( int i = 1; i < n; i++ ) v1(i) = v(i);
                    else{ 
                        for( int i = 1; i < n; i++ ) v1(i) = (4*v(i)-(v(i-1)+v(i+1)+2*v(n+1+i)))/h;
                        // if(input.boundaryConditonType[1] == Nuemann ) v1(n) =(4*v(n)-(2*v(n-1)+2*v(2*n+1)))/h;
                    }

                    if( input.boundaryConditonType[1] == Dirichlet ) for( int j = 1; j < n; j++ ) v1(j*(n+1)+n) = v(j*(n+1)+n);
                    else{ 
                        for( int j = 1; j < n; j++ ) v1(j*(n+1)+n) =(4*v(j*(n+1)+n)-(v((j+1)*(n+1)+n)+v((j-1)*(n+1)+n)+2*v(j*(n+1)+n-1)))/h;
                        // if(input.boundaryConditonType[2] == Nuemann ) v1(n*(n+2)) = (4*v(n*(n+2))-(2*v(n*(n+2)-1)+2*v(n*(n+2)-n-1)))/h;
                    }

                    if( input.boundaryConditonType[2] == Dirichlet ) for( int i = 1; i < n; i++ ) v1(n*(n+1)+i) = v(n*(n+1)+i);
                    else{ 
                        for( int i = 1; i < n; i++ ) v1(n*(n+1)+i) = (4*v(n*(n+1)+i)-(v(n*(n+1)+i-1)+v(n*(n+1)+i+1)+2*v((n-1)*(n+1)+i)))/h;
                        // if(input.boundaryConditonType[3] == Nuemann ) v1(n*(n+1)) = (4*v(n*(n+1))-(2*v(n*(n+1)+1)+2*v(n*(n+1)-n-1)))/h;
                    }

                    if( input.boundaryConditonType[3] == Dirichlet ) for( int j = 1; j < n; j++ ) v1(j*(n+1)) = v(j*(n+1));
                    else{ 
                        for( int j = 1; j < n; j++ ) v1(j*(n+1)) = (4*v(j*(n+1))-(v((j+1)*(n+1))+v((j-1)*(n+1))+2*v(j*(n+1)+1)))/h;
                        // if(input.boundaryConditonType[0] == Nuemann ){
                        //     int NuemannCount = 0;
                        //     for( int i = 0; i < 4; i++ ) if( input.boundaryConditonType[i] == Nuemann ) NuemannCount++;
                        //     if( NuemannCount == 4 ) v1(0) = v(0);
                        //     else v1(0) = (4*v(0)-(2*v(1)+2*v(n+1)))/h;
                        // }
                    }
                    break;
                }
                default:
                    break;
            }
            return v1;
        }

        void relaxation( VectorXd & v, VectorXd f, int nu, const double & w , const int & level ){
            int n = N/(pow(2,level-1));
            double h = 1.0/n;
            double h2 = pow(h,2);
            VectorXd v1 = v;
            switch( dim ){
                case 1:{
                    while(nu > 0){
                        for( int i = 1; i < n; i++ ) v1(i) = 0.5*(v(i-1)+v(i+1)+h2*f(i));
                        if( input.boundaryConditonType[0] == Nuemann ) v1(0) = v(1)-h*f(0);
                        else v1(0) = f(0);
                        if( input.boundaryConditonType[1] == Nuemann ) v1(n) = v(n-1)+h*f(n);
                        else v1(n) = f(n);
                        v = w*v1+(1-w)*v;
                        nu--;
                    }
                    break;
                }
                case 2:{
                    while(nu > 0){
                        for( int j = 1; j < n; j++ ){
                            for( int i = 1; i < n; i++ ){
                                v1(j*(n+1)+i) = (v(j*(n+1)+i-1)+v(j*(n+1)+i+1)+v((j+1)*(n+1)+i)+v((j-1)*(n+1)+i)+h2*f(j*(n+1)+i))/4;
                            }
                        }

                        v1(0) = f(0);
                        v1(n) = f(n);
                        v1(n*(n+1)) = f(n*(n+1));
                        v1(n*(n+2))= f(n*(n+2));

                        if( input.boundaryConditonType[0] == Dirichlet ) for( int i = 1; i < n; i++ ) v1(i) = f(i);
                        else{ 
                            for( int i = 1; i < n; i++ ) v1(i) = (v(i-1)+v(i+1)+2*v(n+1+i)+h*f(i))/4;
                            // if(input.boundaryConditonType[1] == Nuemann ) v1(n) = (2*v(n-1)+2*v(2*n+1)+h*f(n))/4;
                        }

                        if( input.boundaryConditonType[1] == Dirichlet ) for( int j = 1; j < n; j++ ) v1(j*(n+1)+n) = f(j*(n+1)+n);
                        else{ 
                            for( int j = 1; j < n; j++ ) v1(j*(n+1)+n) = (v((j+1)*(n+1)+n)+v((j-1)*(n+1)+n)+2*v(j*(n+1)+n-1)+h*f(j*(n+1)+n))/4;
                            // if(input.boundaryConditonType[2] == Nuemann ) v1(n*(n+2)) = (2*v(n*(n+2)-1)+2*v(n*(n+2)-n-1)+h*f(n*(n+2)))/4;
                        }

                        if( input.boundaryConditonType[2] == Dirichlet ) for( int i = 1; i < n; i++ ) v1(n*(n+1)+i) = f(n*(n+1)+i);
                        else{ 
                            for( int i = 1; i < n; i++ ) v1(n*(n+1)+i) = (v(n*(n+1)+i-1)+v(n*(n+1)+i+1)+2*v((n-1)*(n+1)+i)+h*f(n*(n+1)+i))/4;
                            // if(input.boundaryConditonType[3] == Nuemann ) v1(n*(n+1)) = (2*v(n*(n+1)+1)+2*v(n*(n+1)-n-1)+h*f(n*(n+1)))/4;
                        }

                        if( input.boundaryConditonType[3] == Dirichlet ) for( int j = 1; j < n; j++ ) v1(j*(n+1)) = f(j*(n+1));
                        else{ 
                            for( int j = 1; j < n; j++ ) v1(j*(n+1)) = (v((j+1)*(n+1))+v((j-1)*(n+1))+2*v(j*(n+1)+1)+h*f(j*(n+1)))/4;
                            // if(input.boundaryConditonType[0] == Nuemann ){
                            //     int NuemannCount = 0;
                            //     for( int i = 0; i < 4; i++ ) if( input.boundaryConditonType[i] == Nuemann ) NuemannCount++;
                            //     if( NuemannCount == 4 ) v1(0) = f(0);
                            //     else v1(0) = (2*v(1)+2*v(n+1)+h*f(0))/4;
                            // }
                        }

                        v = w*v1+(1-w)*v;
                        nu--;
                    }
                    break;
                }
                default:
                    break;
            }
        }

        //G-S interation
        // void relaxation( VectorXd & v, VectorXd f, int nu, const double & w , const int & level ){
        //     int n = N/(pow(2,level-1));
        //     double h = 1.0/n;
        //     double h2 = pow(h,2);
        //     VectorXd v1 = v;
        //     switch( dim ){
        //         case 1:{
        //             while(nu > 0){
        //                 for( int i = 1; i < n; i++ ) v1(i) = 0.5*(v(i-1)+v(i+1)+h2*f(i));
        //                 if( input.boundaryConditonType[0] == Nuemann ) v1(0) = v(1)-h*f(0);
        //                 else v1(0) = f(0);
        //                 if( input.boundaryConditonType[1] == Nuemann ) v1(n) = v(n-1)+h*f(n);
        //                 else v1(n) = f(n);
        //                 v = w*v1+ (1-w)*v;
        //                 nu--;
        //             }
        //             break;
        //         }
        //         case 2:{
        //             while(nu > 0){
        //                 for( int j = 1; j < n; j++ ){
        //                     for( int i = 1; i < n; i++ ){
        //                         v(j*(n+1)+i) = (v(j*(n+1)+i-1)+v(j*(n+1)+i+1)+v((j+1)*(n+1)+i)+v((j-1)*(n+1)+i)+h2*f(j*(n+1)+i))/4;
        //                     }
        //                 }
        //                 if( input.boundaryConditonType[0] == Dirichlet ) for( int i = 1; i < n; i++ ) v(i) = f(i);
        //                 else{ 
        //                     for( int i = 1; i < n; i++ ) v(i) = (v(i-1)+v(i+1)+2*v(n+1+i)+h*f(i))/4;
        //                     if(input.boundaryConditonType[1] == Nuemann ) v(n) = (2*v(n-1)+2*v(2*n+1)+h*f(n))/4;
        //                 }
        //                 if( input.boundaryConditonType[1] == Dirichlet ) for( int j = 1; j < n; j++ ) v(j*(n+1)+n) = f(j*(n+1)+n);
        //                 else{ 
        //                     for( int j = 1; j < n; j++ ) v(j*(n+1)+n) = (v((j+1)*(n+1)+n)+v((j-1)*(n+1)+n)+2*v(j*(n+1)+n-1)+h*f(j*(n+1)+n))/4;
        //                     if(input.boundaryConditonType[2] == Nuemann ) v(n*(n+2)) = (2*v(n*(n+2)-1)+2*v(n*(n+2)-n-1)+h*f(n*(n+2)))/4;
        //                 }
        //                 if( input.boundaryConditonType[2] == Dirichlet ) for( int i = 1; i < n; i++ ) v(n*(n+1)+i) = f(n*(n+1)+i);
        //                 else{ 
        //                     for( int i = 1; i < n; i++ ) v(n*(n+1)+i) = (v(n*(n+1)+i-1)+v(n*(n+1)+i+1)+2*v((n-1)*(n+1)+i)+h*f(n*(n+1)+i))/4;
        //                     if(input.boundaryConditonType[3] == Nuemann ) v(n*(n+1)) = (2*v(n*(n+1)+1)+2*v(n*(n+1)-n-1)+h*f(n*(n+1)))/4;
        //                 }
        //                 if( input.boundaryConditonType[3] == Dirichlet ) for( int j = 1; j < n; j++ ) v(j*(n+1)) = f(j*(n+1));
        //                 else{ 
        //                     for( int j = 1; j < n; j++ ) v(j*(n+1)+n) = (v((j+1)*(n+1))+v((j-1)*(n+1))+2*v(j*(n+1)+1)+h*f(j*(n+1)))/4;
        //                     if(input.boundaryConditonType[1] == Nuemann ){
        //                         int NuemannCount = 0;
        //                         for( int i = 0; i < 4; i++ ) if( input.boundaryConditonType[i] == Nuemann ) NuemannCount++;
        //                         if( NuemannCount == 4 ) v(0) = f(0);
        //                         else v(0) = (2*v(1)+2*v(n+1)+h*f(0))/4;
        //                     }
        //                 }
        //                 v = (1-w)*v1+ w*v;
        //                 v1 = v;
        //                 nu--;
        //             }
        //             break;
        //         }
        //         default:
        //             break;
        //     }
        // }

        VectorXd restriction( const VectorXd & vh ){
            int nh;
            if( dim == 1 ) nh = vh.size()-1;
            else nh = sqrt(vh.size())-1;
            int n2h = nh/2;
            VectorXd v2h;
            switch(dim){
                case 1:{
                    v2h.resize(n2h+1);
                    if( input.rt == FullWeight ){
                        for( int i = 1; i < n2h; i++ ) v2h(i) =  (vh(2*i-1)+2*vh(2*i)+vh(2*i+1))/4;
                        v2h(0)=vh(0);
                        v2h(n2h)=vh(nh);
                    }else{
                        for( int i = 0; i <= n2h; i++ ) v2h(i) = vh(2*i);
                    }
                    break;
                }
                case 2:{
                    v2h.resize((n2h+1)*(n2h+1));
                    if( input.rt == FullWeight ){
                        for( int j = 1; j < n2h; j++ ){
                            for( int i = 1; i < n2h; i++ ){
                                v2h(i+j*(n2h+1)) = vh(2*i+2*j*(nh+1));
                                for( int ii = -1; ii <= 1 ; ii+=2 ){
                                    v2h(i+j*(n2h+1)) += (vh(2*i+ii+2*j*(nh+1))+vh(2*i+(2*j+ii)*(nh+1)))/2;
                                    for( int jj = -1; jj <= 1; jj+=2 ){
                                        v2h(i+j*(n2h+1)) += vh(2*i+ii+(2*j+jj)*(nh+1))/4;
                                    }
                                }
                                v2h(i+j*(n2h+1)) /= 4;
                            }
                        }
                        for( int i = 1; i < n2h; i++ ){
                            v2h(i) = (2*vh(2*i)+vh(2*i-1)+vh(2*i+1))/4;
                            v2h(i*(n2h+1)+n2h) = (2*vh((2*i)*(nh+1)+nh)+vh((2*i-1)*(nh+1)+nh)+vh((2*i+1)*(nh+1)+nh))/4;
                            v2h(n2h*(n2h+1)+i) = (2*vh(nh*(nh+1)+2*i)+vh(nh*(nh+1)+2*i-1)+vh(nh*(nh+1)+2*i+1))/4;
                            v2h(i*(n2h+1)) = (2*vh((2*i)*(nh+1))+vh((2*i-1)*(nh+1))+vh((2*i+1)*(nh+1)))/4;
                        }
                        v2h(0) = vh(0);
                        v2h(n2h) = vh(nh);
                        v2h(n2h*(n2h+1)) = vh(nh*(nh+1));
                        v2h(n2h*(n2h+1)+n2h) = vh(nh*(nh+1)+nh);
                    }else{
                        for( int j = 0; j <= n2h; j++ ){
                            for( int i = 0; i <= n2h; i++ ) v2h(i+j*(n2h+1)) = vh(2*i+2*j*(nh+1));
                        }
                    }
                    break;
                }
                default:
                    break;
            }
            return v2h;
        }

        VectorXd prolongation( const VectorXd & v2h ){
            int n2h;
            if( dim == 1 ) n2h = v2h.size()-1;
            else n2h = sqrt(v2h.size())-1;
            int nh = 2*n2h;
            VectorXd vh;
            switch(dim){
                case 1:{
                    vh.resize(nh+1);
                    if( input.pt == linear ){
                        for( int i = 0; i < n2h; i++ ){
                            vh(2*i) = v2h(i);
                            vh(2*i+1) = (v2h(i)+v2h(i+1))/2;
                        }
                        vh(nh) = v2h(n2h);
                    }else{
                        for( int i = 2; i <= n2h; i+=2 ){
                            vh(2*i) = v2h(i);
                            vh(2*(i-1)) = v2h(i-1);
                            vh(2*i-3) = (3*v2h(i-2)+6*v2h(i-1)-v2h(i))/8;
                            vh(2*i-1) = (3*v2h(i)+6*v2h(i-1)-v2h(i-2))/8;
                        }
                        vh(0) = v2h(0);
                        // for( int i = 1; i < n2h-1; i++ ){
                        //         vh(2*i) = v2h(i);
                        //         vh(2*i+1) = (-v2h(i-1)+9*v2h(i)+9*v2h(i+1)-v2h(i+2))/16; 
                        // }
                        // vh(0) = v2h(0);
                        // vh(1) = (3*v2h(0)+6*v2h(1)-v2h(2))/8;
                        // vh(nh) = v2h(n2h);
                        // vh(nh-1) = (3*v2h(n2h)+6*v2h(n2h-1)-v2h(n2h-2))/8;
                    }
                    break;
                }
                case 2:{
                    vh.resize((nh+1)*(nh+1));
                        for( int j = 0; j < n2h; j++ ){
                            for( int i = 0 ; i < n2h; i++ ){
                                vh(2*j*(nh+1)+2*i) = v2h(j*(n2h+1)+i);
                                vh(2*j*(nh+1)+2*i+1) = (v2h(j*(n2h+1)+i)+v2h(j*(n2h+1)+i+1))/2;
                                vh((2*j+1)*(nh+1)+2*i) = (v2h(j*(n2h+1)+i)+v2h((j+1)*(n2h+1)+i))/2;
                                vh((2*j+1)*(nh+1)+2*i+1) = (v2h(j*(n2h+1)+i)+v2h(j*(n2h+1)+i+1)+v2h((j+1)*(n2h+1)+i)+v2h((j+1)*(n2h+1)+i+1))/4;
                            }
                        }
                        for( int k = 0; k < n2h; k++ ){
                            vh(nh*(nh+1)+2*k+1) = (v2h(n2h*(n2h+1)+k)+v2h(n2h*(n2h+1)+k+1))/2;
                            vh((2*k+1)*(nh+1)+nh) = (v2h(k*(n2h+1)+n2h)+v2h((k+1)*(n2h+1)+n2h))/2;
                        }
                        vh(nh*(nh+1)+nh) = v2h(n2h*(n2h+1)+n2h);
                    break;
                }
                default:
                    break;
            }   
            return vh;
        }

        VectorXd VC( VectorXd v, VectorXd f, const int & nu1, const int & nu2, const int & level ){
            relaxation(v,f,nu1,2.0/3,level);
            if( level == maxLevel ) relaxation(v,f,nu2,2.0/3,level);
            else{
                int n = N/(pow(2,level-1));
                VectorXd f2h = restriction(f-MMV(v,level));
                VectorXd v2h;
                if( dim == 1 ) v2h = vector<double>(n/2+1,0);
                else v2h = vector<double>((n/2+1)*(n/2+1),0);
                v2h = VC(v2h,f2h,nu1,nu2,level+1);
                v = v + prolongation(v2h);
                relaxation(v,f,nu2,2.0/3,level);
            }
            return v;
        }

        VectorXd FMG( VectorXd f, const int & nu1, const int & nu2, const int level ){
            VectorXd v;
            if( level == maxLevel ){
                int n = f.size();
                v = vector<double>(n,0);
                v = VC(v,f,nu1,nu2,level);
            }else{
                VectorXd f2h = restriction(f);
                VectorXd v2h = FMG(f2h,nu1,nu2,level+1);
                v = prolongation(v2h);
                v = VC(v,f,nu1,nu2,level);
            }
            return v;
        }

        userInput<dim> input;
        int N;//Grid Spacing
        Result<dim> result;
        // bool inputValid;
        int maxLevel;
};