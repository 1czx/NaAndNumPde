#include"MultigridPossionSolver.h"
#include<iomanip>
using namespace std;

string inputDir{"../input/"};
//test1D
double testU1D( const double & x ){
    return exp(x+sin(x));
}

double testF1D( const double & x ){
    return (sin(x)-pow(1+cos(x),2))*exp(x+sin(x));
}

double testG1D( const double & x ){
    return (1+cos(x))*exp(x+sin(x));
}


//test2D
double testU2D( const double & x,const double & y ){
    return exp(y+sin(x));
}

double testF2D( const double & x,const double & y ){
    return -exp(y+sin(x))*(1-sin(x)+cos(x)*cos(x));
}

double testGx2D( const double & x,const double & y ){
    return cos(x)*exp(y+sin(x));
}

double testGy2D( const double & x,const double & y ){
    return exp(y+sin(x));
}

int main(){
    cout << "\nStart of test accuracy\n" << endl;
    cout << setw(8) << left << "eps" << setw(12) << left << "1DResidual" << setw(12) << left << "2DResidual" << setw(8) << left << "1Diter" << setw(8) << left << "2Diter" << endl;
    userInput<1> ip1d(testF1D,testU1D,testG1D);
    userInput<2> ip2d(testF2D,testU2D,testGx2D,testGy2D);    
    ip1d.readFile(inputDir+"input1.json");
    ip2d.readFile(inputDir+"input1.json");
    MultigridPossionSolver<1> s1;
    s1.setUserInput(ip1d);
    MultigridPossionSolver<2> s2;
    s2.setUserInput(ip2d);
    double eps = 1e-9;
    int n = 256;
    s1.setGridSpacing(n);
    s2.setGridSpacing(n);
    while( eps > (2.2)*1e-16 ){
        ip1d.eps = eps;
        ip2d.eps = eps;
        cout << setw(8) << left << eps; 
        s1.setUserInput(ip1d);
        s2.setUserInput(ip2d);
        s1.solve();
        s2.solve();
        cout << setw(8) << left << s1.getResult().iter << setw(8) << left << s2.getResult().iter << endl;
        eps /= 10; 
    }
    eps = 2.2e-16;
    ip1d.eps = eps;
    ip2d.eps = eps;
    s1.setUserInput(ip1d);
    s2.setUserInput(ip2d);
    cout << setw(8) << left << eps; 
    s1.solve();
    s2.solve();
    cout << setw(8) << left << s1.getResult().iter << setw(8) << left << s2.getResult().iter << endl;
    cout << "\nEnd of test accuracy\n"; 
}