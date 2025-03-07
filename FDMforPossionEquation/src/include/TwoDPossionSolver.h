#pragma once

#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<algorithm>
#include<functional>
#include"Eigen/Dense"
#include"json/json.h"

using namespace Eigen;
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

struct userInput{
    vector<BCT> boundaryConditonType;//in the order of bottom, right, top, left, circle
    bool domainType; //0: (0,1)^2; 1: (0,1)^2\D
    double centerX;
    double centerY;
    double radius;
    function<double(const double &,const double &)> f;//RHS
    function<double(const double &,const double &)> u;//Dirichlet boundary condition
    function<double(const double &,const double &)> gx;//Nuemann condition for vertical boundary 
    function<double(const double &,const double &)> gy;//Nuemann condition for horizontal boundary
    function<double(const double &,const double &)> g;//Nuemann condition for circle boundary
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
            domainType = root["DomainType"].asBool();
            outputfile = root["opfile"].asString();
            Json::Value BC = root["BoundaryCondition"];
            boundaryConditonType.clear();
            for( int i = 0; i < BC.size(); i++ ){
                if( BC[i].asString() == "Dirichlet" ) boundaryConditonType.push_back(Dirichlet);
                else boundaryConditonType.push_back(Nuemann);
            }
            if( domainType == 1 ){
                centerX = root["centerX"].asDouble();
                centerY = root["centerY"].asDouble();
                radius = root["radius"].asDouble();
            }
        }else{
            cout << "Parse Error\n";
        }
        in.close();
    }

    userInput():centerX{0.5},centerY{0.5},radius{0} {}
    userInput(
        const function<double(const double &,const double &)> & _f,//RHS
        const function<double(const double &,const double &)> & _u,//Dirichlet boundary condition
        const function<double(const double &,const double &)> & _gx,//Nuemann condition for vertical boundary 
        const function<double(const double &,const double &)> & _gy,//Nuemann condition for horizontal boundary
        const function<double(const double &,const double &)> & _g //Nuemann condition for circle boundary
    ):f{_f},u{_u},gx{_gx},gy(_gy),g(_g),centerX{0.5},centerY{0.5},radius{0} {}
};

struct Result{
        VectorXd U;//numerical solution
        VectorXd U_cap;//ture solution
        VectorXd E;//abs error |U-U_cap|
        int n;
        double L1NormErr(){
            double sum = 0;
            sum = E.lpNorm<1>();
            sum = sum/(n*n);
            return sum;
        }
        double L2NormErr(){
            double sum = 0;
            sum = E.lpNorm<2>();
            sum = sum/n;
            return sum;
        }
        double LInfNormErr(){
            return E.lpNorm<Infinity>();
        }
        void outputU( const string & filename ){
            ofstream fout(filename);
            for( int j = 0; j < n+1; j++ ){
                for( int i = 0; i < n+1; i++ ){
                    fout << U(j*(n+1)+i) << " ";
                }
                fout << endl;
            }
        }
        void outputE( const string & filename ){
            ofstream fout(filename);
            for( int j = 0; j < n+1; j++ ){
                for( int i = 0; i < n+1; i++ ){
                    fout << E(j*(n+1)+i) << " ";
                }
                fout << endl;
            }
        }
};

class TwoDPossionSolver{
    public:
        TwoDPossionSolver():inputValid{0}, n{0}{ }
        void setUserInput( const userInput & uI ){ 
            input = uI;
            // inputValid = 1;
        }
        
        
        void setGridSpacing( const int & N ){ n = N; }

        void solve(){
            if( !isInputValid() ){
                cout << "Error: The user-specified parameter is invalid!" << endl;
                return;
            }
            double h = 1.0/n;
            double h2= pow(h,2);
            VectorXd f = VectorXd::Zero((n+1)*(n+1));
            MatrixXd A = MatrixXd::Zero((n+1)*(n+1),(n+1)*(n+1));
            //left bottom
            A(0,0) = 1;
            f(0) = input.u(0,0);
            //right bottom
            A(n,n) = 1;
            f(n) = input.u(1,0);
            //left top
            A(n*(n+1),n*(n+1)) = 1;
            f(n*(n+1)) = input.u(0,1);
            //right top
            A(n*(n+2),n*(n+2)) = 1;
            f(n*(n+2))= input.u(1,1);
            switch (input.domainType)
            {
            case 0:{
                getRegularFDMatrix(A,f);      
                break;
            }
            case 1:{ 
                getIrregularFDMatrix(A,f);
                break;
            }
            default:
                return;
            }
            result.U = A.partialPivLu().solve(f);
            result.U_cap = VectorXd::Zero(((n+1)*(n+1)));
            for( int j = 0; j < n+1; j++){
                for( int i = 0; i < n+1; i++ ){
                    if( (input.domainType == 0 && isEDP( i,j )) || (!isInsideDisk( i,j ) && isEDP(i,j)) ) result.U_cap(j*(n+1)+i) = input.u(i*h,j*h);
                    else result.U_cap(j*(n+1)+i) = result.U(j*(n+1)+i);
                }
            }
            result.E = result.U_cap - result.U;
            result.E = result.E.array().abs();
            result.n = n;
        }

        Result getResult(){
            return result;
        }

    private:

        bool isInputValid(){
            if( n == 0 ) return 0;

            if( input.f == nullptr ) return 0;

            if( input.domainType == 1 && input.boundaryConditonType.size() != 5 ) return 0;

            if( input.domainType == 0 && input.boundaryConditonType.size() != 4 ) return 0;
            
            bool typeD = 0, typeN = 0;
            for( auto & t: input.boundaryConditonType ){
                if( t == Dirichlet ) typeD = 1;
                else typeN = 1;
            }

            if( typeD = 1 && (input.u == nullptr ) ) return 0; 
            if( typeN = 1 && ( (input.g == nullptr ) || (input.gx == nullptr ) || (input.gy == nullptr ) ) ) return 0; 

            if( input.centerX < 1 && input.centerX > 0 && input.centerY < 1 && input.centerY > 0 ){
                int count = 0;
                for( int i = 0; i <= 1; i++ ){
                    for( int j =0; j <= 1; j++ ){
                        if( isInsideDisk(i*n,j*n) ) continue;
                        if( abs(i-input.centerX) < input.radius &&abs (j-input.centerY) < input.radius ) count++;
                    }
                }
                if( count >= 2 ) return 0;
            }

            if( input.domainType == 0 ) return 1;

            bool flag = 1;
            double h = 1.0/n;
            int NCI = floor(input.centerX/h+0.5);//nearest grid coordinate i to center
            int NCJ = floor(input.centerY/h+0.5);//nearest grid coordinate j to center
            NCI = min(max(1,NCI),n-1);
            NCJ = min(max(1,NCI),n-1);

            if( isInsideDisk(NCI,NCJ) ){
                int count = 1;
                for( int i = 1; i <= 4; i++ ){
                    for( int t = -1; t < 1; t+=2 ){
                        if( isInsideDisk(NCI+t*i,NCJ) && isEDP(NCI+t*i,NCJ) ) count++ ;
                        if( isInsideDisk(NCI,NCJ+t*i) && isEDP(NCI,NCJ+t*i) ) count++ ;
                        if( count >= 4 ) return 1;
                    }
                    for( int i = -1; i <= 1; i+=2 ){
                        for( int j = -1; j <= 1; j+=2 ){
                            if( isInsideDisk(NCI+i,NCJ+j) && isEDP(NCI+i,NCJ+j) ) count++ ;
                            if( count >= 4 ) return 1;
                        }
                    }
                }
                if( count < 4 ) return 0;
            }else return 0;
            
            return flag;
        }
        
        void getRegularFDMatrix( MatrixXd & A, VectorXd & f ){
            double h = 1.0/n;
            double h2= pow(h,2);
            for( int j = 1; j < n; j++ ){
                for( int i = 1; i < n; i++ ){
                    A(j*(n+1)+i,j*(n+1)+i) = 4;
                    A(j*(n+1)+i,(j-1)*(n+1)+i) = -1;
                    A(j*(n+1)+i,(j+1)*(n+1)+i) = -1;
                    A(j*(n+1)+i,j*(n+1)+(i-1)) = -1;
                    A(j*(n+1)+i,j*(n+1)+(i+1)) = -1;
                    f(j*(n+1)+i) = h2*input.f(i*h,j*h);
                }
            }
            int NuemannCount = 0;
            //bottom boundary
            if( input.boundaryConditonType[0] == Dirichlet ){
                for( int i = 1; i < n; i++ ){
                    A(i,i) = 1;
                    f(i) = input.u(i*h,0);
                }
            }else{
                NuemannCount++;
                for( int i = 1; i < n; i++ ){
                    A(i,i) = 4;
                    A(i,i-1) = -1;
                    A(i,i+1) = -1;
                    A(i,n+1+i) = -2;
                    f(i) = h2*input.f(i*h,0)-2*h*input.gy(i*h,0);
                }
                if(input.boundaryConditonType[1] == Nuemann ){
                    A(n,n) = 4;
                    A(n,n-1) = -2;
                    A(n,n+n+1) = -2;
                    f(n) = h2*input.f(1,0)-2*h*(input.gy(1,0)-input.gx(1,0));
                    // A(n,n) = -1;
                    // A(n,2*n+1) = 1;
                    // f(n) = h*input.gy(1,0);
                }
            }
            //right boundary
            if( input.boundaryConditonType[1] == Dirichlet ){
                for( int j = 1; j < n; j++ ){
                    A(j*(n+1)+n,j*(n+1)+n) = 1;
                    f(j*(n+1)+n) = input.u(1,j*h);
                }
            }else{
                NuemannCount++;
                for( int j = 1; j < n; j++ ){
                    A(j*(n+1)+n,j*(n+1)+n) = 4;
                    A(j*(n+1)+n,(j-1)*(n+1)+n) = -1;
                    A(j*(n+1)+n,(j+1)*(n+1)+n) = -1;
                    A(j*(n+1)+n,j*(n+1)+n-1) = -2;
                    f(j*(n+1)+n) = h2*input.f(1,j*h)+2*h*input.gx(1,j*h);
                }
                if(input.boundaryConditonType[2] == Nuemann ){
                    A(n*(n+2),n*(n+2)) = 4;
                    A(n*(n+2),n*(n+2)-1) = -2;
                    A(n*(n+2),n*(n+2)-n-1) = -2;
                    f(n*(n+2)) = h2*input.f(1,1)+2*h*(input.gy(1,1)+input.gx(1,1));
                    // A(n*(n+2),n*(n+2)) = 1;
                    // A(n*(n+2),n*(n+2)-1) = -1;
                    // f(n*(2+n)) = h*input.gx(1,1);
               }
            }
            //top boundary
            if( input.boundaryConditonType[0] == Dirichlet ){
                for( int i = 1; i < n; i++ ){
                    A(n*(n+1)+i,n*(n+1)+i) = 1;
                    f(n*(n+1)+i) = input.u(i*h,1);
                }
            }else{
                NuemannCount++;
                for( int i = 1; i < n; i++ ){
                    A(n*(n+1)+i,n*(n+1)+i) = 4;
                    A(n*(n+1)+i,n*(n+1)+i+1) = -1;
                    A(n*(n+1)+i,n*(n+1)+i-1) = -1;
                    A(n*(n+1)+i,n*(n+1)+i-n-1) = -2;
                    f(n*(n+1)+i) = h2*input.f(i*h,1)+2*h*input.gy(i*h,1);
                }
                if(input.boundaryConditonType[3] == Nuemann ){
                    A(n*(n+1),n*(n+1)) = 4;
                    A(n*(n+1),n*(n+1)+1) = -2;
                    A(n*(n+1),n*(n+1)-n-1) = -2;
                    f(n*(n+1)) = h2*input.f(0,1)+2*h*(input.gy(0,1)-input.gx(0,1));
                    // A(n*(n+1),n*(n+1)) = 1;
                    // A(n*(n+1),(n-1)*(n+1)) = -1;
                    // f(n*(n+1)) = h*input.gy(0,1);
                }
            }
            //left boundary
            if( input.boundaryConditonType[3] == Dirichlet ){
                for( int j = 1; j < n; j++ ){
                    A(j*(n+1),j*(n+1)) = 1;
                    f(j*(n+1)) = input.u(0,j*h);
                }
            }else{
                for( int j = 1; j < n; j++ ){
                    A(j*(n+1),j*(n+1)) = 4;
                    A(j*(n+1),(j-1)*(n+1)) = -1;
                    A(j*(n+1),(j+1)*(n+1)) = -1;
                    A(j*(n+1),j*(n+1)+1) = -2;
                    f(j*(n+1)) = h2*input.f(0,j*h)-2*h*input.gx(0,j*h);
                }
                if(input.boundaryConditonType[0] == Nuemann ){
                    if( NuemannCount == 3 ){
                        A(0,0) = 1;
                        f(0) = input.u(0,0);
                        return;
                    }
                    A(0,0) = 4;
                    A(0,1) = -2;
                    A(0,n+1) = -2;
                    f(0) = h2*input.f(0,0)-2*h*(input.gy(0,0)+input.gx(0,0));
                    // A(0,0) = -1;
                    // A(0,1) = 1;
                    // f(0) = h*input.gx(0,0);
                }
            }
        }

        void getIrregularFDMatrix( MatrixXd & A, VectorXd & f ){
            double h = 1.0/n;
            double h2= pow(h,2);
            for( int j = 1; j < n; j++ ){
                for( int i = 1; i < n; i++ ){
                    if( !isInsideDisk(i,j) ){
                        A(j*(n+1)+i,j*(n+1)+i) = 4;
                        A(j*(n+1)+i,(j-1)*(n+1)+i) = -1;
                        A(j*(n+1)+i,(j+1)*(n+1)+i) = -1;
                        A(j*(n+1)+i,j*(n+1)+(i-1)) = -1;
                        A(j*(n+1)+i,j*(n+1)+(i+1)) = -1;
                        f(j*(n+1)+i) = h2*input.f(i*h,j*h);
                        double alpha{1}, theta{1};
            
                        for( int t = -1; t <= 1; t+=2 ){
                            
                            if( isInsideDisk( i+t, j ) ){
                                if( input.boundaryConditonType[4] == Dirichlet ){
                                    double dx = sqrt(pow(input.radius,2)-pow(j*h-input.centerY,2));
                                    alpha = (abs(i*h-input.centerX)-dx)/h;
                                    A(j*(n+1)+i,j*(n+1)+i-t) = -2/(1+alpha);
                                    f(j*(n+1)+i) += 2*input.u((i+t*alpha)*h,j*h)/(alpha*(alpha+1));
                                }else setNCondID(A,f,i+t,j);
                            }

                            if( isInsideDisk( i, j+t ) ){
                                if( input.boundaryConditonType[4] == Dirichlet ){
                                    double dy = sqrt(pow(input.radius,2)-pow(i*h-input.centerX,2));
                                    theta = (abs((j+t)*h-input.centerY)-dy)/h;
                                    A(j*(n+1)+i,(j-t)*(n+1)+i) = -2/(1+theta);
                                    f(j*(n+1)+i) += 2*input.u(i*h,(j+t*theta)*h)/(theta*(theta+1));
                                }else setNCondID(A,f,i,j+t);
                            }

                        }
                        A(j*(n+1)+i,j*(n+1)+i) = 2/alpha + 2/theta;
                    
                    }else{
                        if( A(j*(n+1)+i,j*(n+1)+i) == 0 ) A(j*(n+1)+i,j*(n+1)+i) = 1;
                    }
                }
            }
            int NuemannCount = 0;
            for( int i = 0 ;i < 4; i++ ) if( input.boundaryConditonType[i] == Nuemann ) NuemannCount++;
            //bottom boundary
            if( input.boundaryConditonType[0] == Dirichlet ){
                for( int i = 1; i < n; i++ ){
                    A(i,i) = 1;
                    f(i) = input.u(i*h,0);
                }
            }else{
                for( int i = 1; i < n; i++ ){
                    if( !isInsideDisk( i, 0 ) ){
                        if( isInsideDisk( i-1, 0 ) || isInsideDisk( i+1, 0 ) ){
                            A(i,i) = -1;
                            A(i,i+n+1) = 1;
                            f(i) = h*input.gy(i*h,0);
                        }else{
                            A(i,i) = 4;
                            A(i,i-1) = -1;
                            A(i,i+1) = -1;
                            A(i,n+1+i) = -2;
                            f(i) = h2*input.f(i*h,0)-2*h*input.gy(i*h,0);
                            if( isInsideDisk( i, 1 ) ) setNCondID(A,f,i,1);
                        }
                    }else A(i,i) = 1;
                }
                if(input.boundaryConditonType[1] == Nuemann ){
                    if( isInsideDisk( n, 0 ) ){
                        A(n,n) = 1;
                        f(n) = 0;
                    }else if( NuemannCount != 4 ){
                        if( isInsideDisk( n-1, 0 ) ){
                            A(n,n) = -1;
                            A(n,2*n+1) = 1;
                            f(n) = h*input.gy(1,0);
                        }else if( isInsideDisk( n, 1 ) ){
                            A(n,n) = 1;
                            A(n,n-1) = -1;
                            f(n) = h*input.gx(1,0);
                        }else {
                            A(n,n) = 4;
                            A(n,n-1) = -2;
                            A(n,n+n+1) = -2;
                            f(n) = h2*input.f(1,0)-2*h*(input.gy(1,0)-input.gx(1,0));
                        }
                    }else NuemannCount = 0;
                }
            }
            //right boundary
            if( input.boundaryConditonType[1] == Dirichlet ){
                for( int j = 1; j < n; j++ ){
                    A(j*(n+1)+n,j*(n+1)+n) = 1;
                    f(j*(n+1)+n) = input.u(1,j*h);
                }
            }else{
                for( int j = 1; j < n; j++ ){
                    if( !isInsideDisk( n, j ) ){
                        if( isInsideDisk( n, j-1 ) || isInsideDisk( n, j+1 ) ){
                            A(j*(n+1)+n,j*(n+1)+n) = 1;
                            A(j*(n+1)+n,j*(n+1)+n-1) = -1;
                            f(j*(n+1)+n) = h*input.gx(1,j*h);
                        }
                        else{
                            A(j*(n+1)+n,j*(n+1)+n) = 4;
                            A(j*(n+1)+n,(j-1)*(n+1)+n) = -1;
                            A(j*(n+1)+n,(j+1)*(n+1)+n) = -1;
                            A(j*(n+1)+n,j*(n+1)+n-1) = -2;
                            f(j*(n+1)+n) = h2*input.f(1,j*h)+2*h*input.gx(1,j*h);
                            if( isInsideDisk( n-1, j ) ) setNCondID(A,f,n-1,j);
                        }
                    }else A(j*(n+1)+n,j*(n+1)+n) = 1;
                }
                if(input.boundaryConditonType[2] == Nuemann ){
                    if( isInsideDisk( n, n ) ){
                        A(n*(n+2),n*(n+2)) = 1;
                        f(n*(n+2)) = 0;
                    }else if( NuemannCount != 4 ){
                        if( isInsideDisk( n-1, n ) ){
                            A(n*(n+2),n*(n+2)) = 1;
                            A(n*(n+2),(n-1)*(n+1)+n) = -1;
                            f(n*(n+2)) = h*input.gy(1,1);
                        }else if( isInsideDisk( n, n-1 ) ){
                            A(n*(n+2),n*(n+2)) = 1;
                            A(n*(n+2),n*(n+2)-1) = -1;
                            f(n*(n+2)) = h*input.gx(1,1);
                        }else{
                            A(n*(n+2),n*(n+2)) = 4;
                            A(n*(n+2),n*(n+2)-1) = -2;
                            A(n*(n+2),n*(n+2)-n-1) = -2;
                            f(n*(n+2)) = h2*input.f(1,1)+2*h*(input.gy(1,1)+input.gx(1,1));
                        }
                    }else NuemannCount = 0;
                }
            }
            //top boundary
            if( input.boundaryConditonType[0] == Dirichlet ){
                for( int i = 1; i < n; i++ ){
                    A(n*(n+1)+i,n*(n+1)+i) = 1;
                    f(n*(n+1)+i) = input.u(i*h,1);
                }
            }else{
                for( int i = 1; i < n; i++ ){
                    if( !isInsideDisk( i, n ) ){
                        if( isInsideDisk( i-1, n ) || isInsideDisk( i+1, n ) ){
                            A(n*(n+1)+i,n*(n+1)+i) = 1;
                            A(n*(n+1)+i,(n-1)*(n+1)+i) = -1;
                            f(n*(n+1)+i) = h*input.gy(i*h,1);
                        }else{
                            A(n*(n+1)+i,n*(n+1)+i) = 4;
                            A(n*(n+1)+i,n*(n+1)+i+1) = -1;
                            A(n*(n+1)+i,n*(n+1)+i-1) = -1;
                            A(n*(n+1)+i,n*(n+1)+i-n-1) = -2;
                            f(n*(n+1)+i) = h2*input.f(i*h,1)+2*h*input.gy(i*h,1);
                            if( isInsideDisk( i, n-1 ) ) setNCondID(A,f,i,n-1);
                        }
                    }else A(n*(n+1)+i,n*(n+1)+i) = 1;
                }
                if(input.boundaryConditonType[3] == Nuemann ){
                    if( isInsideDisk( 0, n ) ){
                        A(n*(n+1),n*(n+1)) = 1;
                        f(n*(n+1)) = 0;
                    }else if( NuemannCount != 4 ){
                        if( isInsideDisk( 0, n-1 ) ){
                            A(n*(n+1),n*(n+1)) = -1;
                            A(n*(n+1),n*(n+1)+1) = 1;
                            f(n*(n+1)) = h*input.gx(0,1);
                        }else if( isInsideDisk( 1, n ) ){
                            A(n*(n+1),n*(n+1)) = 1;
                            A(n*(n+1),(n-1)*(n+1)) = -1;
                            f(n*(n+1)) = h*input.gy(0,1);
                        }else{
                            A(n*(n+1),n*(n+1)) = 4;
                            A(n*(n+1),n*(n+1)+1) = -2;
                            A(n*(n+1),n*(n+1)-n-1) = -2;
                            f(n*(n+1)) = h2*input.f(0,1)+2*h*(input.gy(0,1)-input.gx(0,1));
                        }
                    }else NuemannCount = 0;
                }
            }
            //left boundary
            if( input.boundaryConditonType[3] == Dirichlet ){
                for( int j = 1; j < n; j++ ){
                    A(j*(n+1),j*(n+1)) = 1;
                    f(j*(n+1)) = input.u(0,j*h);
                }
            }else{
                for( int j = 1; j < n; j++ ){
                    if( !isInsideDisk( 0, j ) ){
                        if( isInsideDisk( 0, j-1 ) || isInsideDisk( 0, j+1 ) ){
                            A(j*(n+1),j*(n+1)) = -1;
                            A(j*(n+1),j*(n+1)+1) = 1;
                            f(j*(n+1)) = h*input.gx(0,j*h);
                        }else{
                            A(j*(n+1),j*(n+1)) = 4;
                            A(j*(n+1),(j-1)*(n+1)) = -1;
                            A(j*(n+1),(j+1)*(n+1)) = -1;
                            A(j*(n+1),j*(n+1)+1) = -2;
                            f(j*(n+1)) = h2*input.f(0,j*h)-2*h*input.gx(0,j*h);
                            if( isInsideDisk( 1, j ) ) setNCondID(A,f,1,j);
                        }
                    }else A(j*(n+1),j*(n+1)) = 1;
                }
                if(input.boundaryConditonType[0] == Nuemann ){
                    if( isInsideDisk( 0, 0 ) ){
                        A(0,0) = 1;
                        f(0) = 0;
                    }else if( NuemannCount != 4 ){
                        if( isInsideDisk( 0, 1 ) ){
                            A(0,0) = -1;
                            A(0,1) = 1;
                            f(0) = h*input.gx(0,0);
                        }else if( isInsideDisk( 1, 0 ) ){
                            A(0,0) = -1;
                            A( 0, n+1 ) = 1;
                            f(0) = h*input.gy(0,0);
                        }else{
                            A(0,0) = 4;
                            A(0,1) = -2;
                            A(0,n+1) = -2;
                            f(0) = h2*input.f(0,0)-2*h*(input.gy(0,0)+input.gx(0,0));
                        }
                    }else NuemannCount = 0;
                }
            }
        }

        int sgn( const double & x ){
            if( x > 0 ) return 1;
            else if( x < 0 ) return -1;
            else return 0;
        }

        bool isInsideDisk( const int & i, const int & j ){
            double h = 1.0/n;
            if( pow(i*h-input.centerX,2)+pow(j*h-input.centerY,2) <= pow(input.radius,2) ) return 1;
            else return 0;
        }
        //if the point is the eqution-discreated point
        bool isEDP( const int & i, const int & j ){
            if( i >= 1 && i <= n-1 && j >= 1 && j <= n-1 ) return 1;
            else return 0;
        }
        //set Nuemann boundary conditon on irregular boundary
        void setNCondID( MatrixXd & A, VectorXd & f, const int & i, const int & j ){
            double h = 1.0/n;
            double dx = i*h-input.centerX;
            double dy = j*h-input.centerY;
            double d = sqrt(dx*dx+dy*dy);
            int signx = sgn(dx);
            int signy = sgn(dy);
            if( abs(dy) > abs(dx) ){
                A(j*(n+1)+i,j*(n+1)+i) = abs(dy)/d;
                A(j*(n+1)+i,(j+signy)*(n+1)+i+signx) = -abs(dx)/d;
                A(j*(n+1)+i,(j+signy)*(n+1)+i) = -(abs(dy)-abs(dx))/d;
                f(j*(n+1)+i) = h*input.g(input.centerX+dx*input.radius/d, input.centerY+dy*input.radius/d);
            }else if( abs(dy) < abs(dx) ) {
                A(j*(n+1)+i,j*(n+1)+i) = abs(dx)/d;
                A(j*(n+1)+i,(j+signy)*(n+1)+i+signx) = -abs(dy)/d;
                A(j*(n+1)+i,j*(n+1)+i+signx) = -(abs(dx)-abs(dy))/d;
                f(j*(n+1)+i) = h*input.g(input.centerX+dx*input.radius/d, input.centerY+dy*input.radius/d);
            }else{
                A(j*(n+1)+i,j*(n+1)+i) = 1;
                A(j*(n+1)+i,(j+signy)*(n+1)+i+signx) = -1;
                f(j*(n+1)+i) = sqrt(2)*h*input.g(input.centerX+dx*input.radius/d, input.centerY+dy*input.radius/d);
            }
        }

        userInput input;
        int n;//Grid Spacing
        Result result;
        bool inputValid;

};