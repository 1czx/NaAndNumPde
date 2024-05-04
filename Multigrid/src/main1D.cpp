#include"MultigridPossionSolver.h"

using namespace std;

string outputDir{"../result/"};
string inputDir{"../input/"};

//test1
double testU1( const double & x ){
    return exp(x+sin(x));
}

double testF1( const double & x ){
    return (sin(x)-pow(1+cos(x),2))*exp(x+sin(x));
}

double testG1( const double & x ){
    return (1+cos(x))*exp(x+sin(x));
}

//test2
double testU2( const double & x ){
    return sin(M_PI*x);
}

double testF2( const double & x ){
    return M_PI*M_PI*sin(M_PI*x);
}

double testG2( const double & x ){
    return M_PI*cos(M_PI*x);
}

//test3
double testU3( const double & x ){
    return pow(x,3);
}

double testF3( const double & x ){
    return -6*x;
}

double testG3( const double & x ){
    return 3*x*x;
}



void test(
        const function<double(const double &)> & f,//RHS
        const function<double(const double &)> & u,//Dirichlet boundary condition
        const function<double(const double &)> & g,//Nuemann condition 
        const string & testid
    ){  
        cout << "\nStart of " + testid << endl;
        userInput<1> input(f,u,g);
        MultigridPossionSolver<1> slover;
        Result<1> result;
        for( int i = 1; i <= 8; i++ ){
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

int main(){
    test(testF1,testU1,testG1,"1Dtest1");
    test(testF2,testU2,testG2,"1Dtest2");
    test(testF3,testU3,testG3,"1Dtest3");
}