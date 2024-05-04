#include"MultigridPossionSolver.h"
#include<chrono>
using namespace std;

string outputDir{"../result/"};
string inputDir{"../input/"};

//test1
double testU1( const double & x,const double & y ){
    return exp(y+sin(x));
}

double testF1( const double & x,const double & y ){
    return -exp(y+sin(x))*(1-sin(x)+cos(x)*cos(x));
}

double testGx1( const double & x,const double & y ){
    return cos(x)*exp(y+sin(x));
}

double testGy1( const double & x,const double & y ){
    return exp(y+sin(x));
}

//test2
double testU2( const double & x,const double & y ){
    return sin(M_PI*x)*sin(M_PI*y);
}

double testF2( const double & x,const double & y ){
    return 2.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
}

double testGx2( const double & x,const double & y ){
    return M_PI*cos(M_PI*x)*sin(M_PI*y);
}

double testGy2( const double & x,const double & y ){
    return M_PI*sin(M_PI*x)*cos(M_PI*y);
}

//test3
double testU3( const double & x,const double & y ){
    return pow(x,3)+pow(y,3);
}

double testF3( const double & x,const double & y ){
    return -6*(x+y);
}

double testGx3( const double & x,const double & y ){
    return 3*x*x;
}

double testGy3( const double & x,const double & y ){
    return 3*y*y;
}


void test(
        const function<double(const double &,const double &)> & f,//RHS
        const function<double(const double &,const double &)> & u,//Dirichlet boundary condition
        const function<double(const double &,const double &)> & gx,//Nuemann condition for vertical boundary 
        const function<double(const double &,const double &)> & gy,//Nuemann condition for horizontal boundary
        const string & testid
    ){
        cout << "\nStart of " + testid << endl;
        userInput<2> input(f,u,gx,gy);
        MultigridPossionSolver<2> slover;
        Result<2> result;
        for( int i = 1; i <= 4; i++ ){
            cout << endl;
            string filename{inputDir+"input"+to_string(i)+".json"};
            input.readFile(filename);
            slover.setUserInput(input);
            int n = 32;
            ofstream fout(outputDir+testid+input.outputfile+"ErrorNorm");
            ofstream fout2(outputDir+testid+input.outputfile+"Time");
            ofstream fout3(outputDir+testid+input.outputfile+"Res");
            for( int k = 0; k < 4; k++ ){
                string outputFile{outputDir+testid+"_n="+to_string(n)+"_"+input.outputfile};
                slover.setGridSpacing(n);
                slover.solve(fout3);
                fout3 << endl;
                result = slover.getResult();
                result.outputU(outputFile+"_U");
                result.outputE(outputFile+"_E");
                fout << result.LInfNormErr() << endl;
                fout2 << result.time << endl;
                cout << testid+"_n="+to_string(n)+"_"+input.outputfile + " : cost time is " << result.time << "ms; number of iteration is " << result.iter << endl; 
                n=n*2;
            }
            fout.close();
            fout2.close();
            fout3.close();
        }
        cout << "\nEnd of " + testid << endl;
    }

void testLUtime( userInput<2> input ){
        int n = 64;
        double h = 1.0/n;
        double h2= pow(h,2);
        VectorXd f = vector<double>((n+1)*(n+1),0);
        MatrixXd A((n+1)*(n+1),(n+1)*(n+1));
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
        auto start = std::chrono::high_resolution_clock::now();
        A.solve(f);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> tm = end - start;
        cout << "\nTest cost of LU factorization(2Dtest3 n = 64): " << tm.count() << "ms" << endl;
    }

int main(){
    test(testF1,testU1,testGx1,testGy1,"2Dtest1");
    test(testF2,testU2,testGx2,testGy2,"2Dtest2");
    test(testF3,testU3,testGx3,testGy3,"2Dtest3");
    userInput<2> input(testF3,testU3,testGx3,testGy3);
    input.readFile(inputDir+"input1.json");
    testLUtime(input);
}