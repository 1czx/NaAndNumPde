#include "TwoDPossionSolver.h"
#include <omp.h>

using namespace std;

string outputDir{"../result/"};
string inputDir{"../input/"};

//test1
double testU( const double & x,const double & y ){
    return sin(M_PI*x)*sin(M_PI*y);
}

double testF( const double & x,const double & y ){
    return 2.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
}

double testGx( const double & x,const double & y ){
    return M_PI*cos(M_PI*x)*sin(M_PI*y);
}

double testGy( const double & x,const double & y ){
    return M_PI*sin(M_PI*x)*cos(M_PI*y);
}

double testG( const double & x,const double & y ){
    return -((x-0.1)*testGx(x,y)+(y-0.1)*testGy(x,y))/0.3;
}

void test(
        const function<double(const double &,const double &)> & f,//RHS
        const function<double(const double &,const double &)> & u,//Dirichlet boundary condition
        const function<double(const double &,const double &)> & gx,//Nuemann condition for vertical boundary 
        const function<double(const double &,const double &)> & gy,//Nuemann condition for horizontal boundary
        const function<double(const double &,const double &)> & g, //Nuemann condition for circle boundary
        const string & testid
    ){
        userInput input(f,u,gx,gy,g);
        TwoDPossionSolver slover;
        Result result;
        for( int i = 1; i <= 9; i++ ){
            string filename{inputDir+"input"+to_string(i)+".json"};
            input.readFile(filename);
            slover.setUserInput(input);
            int n = 8;
            ofstream fout(outputDir+testid+input.outputfile+"ErrorNorm");
            for( int k = 0; k < 4; k++ ){
                string outputFile{outputDir+testid+"n="+to_string(n)+input.outputfile};
                slover.setGridSpacing(n);
                slover.solve();
                result = slover.getResult();
                result.outputU(outputFile+"U");
                result.outputE(outputFile+"E");
                n=n*2;
                fout << result.L1NormErr() << " " << result.L2NormErr() << " " << result.LInfNormErr() << endl;
            }
            fout.close();
        }
    }

int main(){
    test(testF,testU,testGx,testGy,testG,"test2");
}