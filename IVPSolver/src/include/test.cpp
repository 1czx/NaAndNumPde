#include"TimeIntegrator.h"
#include"TimeIntegratorFactory.h"

VectorXd ThreeBody( VectorXd u, double t ){
    VectorXd v = u;
    double mu = 0.012277471;
    double u1 = u(0);
    double u2 = u(1);
    double u3 = u(2);
    v(0) = u(3);
    v(1) = u(4);
    v(2) = u(5);
    v(3) = 2*u(4) + u1 - mu*(u1+mu-1)/(pow(u2*u2+u3*u3+pow(u1+mu-1,2),3.0/2)) - (1-mu)*(u1+mu)/(pow(u2*u2+u3*u3+pow(u1+mu,2),3.0/2));
    v(4) = -2*u(3) + u2 - mu*u2/(pow(u2*u2+u3*u3+pow(u1+mu-1,2),3.0/2)) - (1-mu)*u2/(pow(u2*u2+u3*u3+pow(u1+mu,2),3.0/2));
    v(5) = -mu*u3/(pow(u2*u2+u3*u3+pow(u1+mu-1,2),3.0/2)) - (1-mu)*u3/(pow(u2*u2+u3*u3+pow(u1+mu,2),3.0/2));
    return v;
}

int main(){
    TimeIntegratorFactory & Fac = TimeIntegratorFactory::getInstance();
    // VectorXd u0(vector<double>{0.994,0,0,0,-2.0015851063790825224,0});
    // double T1 = 17.06521656015796;
    VectorXd u0(vector<double>{0.879779227778,0,0,0,-0.379677780949,0});
    double T1 = 19.140540691377;
    TimeIntegrator * solver;
    // vector<string> method{"ABMs","AMMs","BDFs","classicalRK","ESDIRK","GaussLegendreRKMs","Fehlberg45","DormandPrince54"};
    // vector<string> method{"GaussLegendreRKMs"};
    vector<string> method{"AMMs","BDFs","ABMs"};
    for(auto & str:method){
        solver = Fac.createTimeIntegrator(str);
        // GaussLegendreRKMs solver;
        // classicalRK solver;
        // BDFs solver;
        // ABMs solver;
        // AMMs solver;
        // Fehlberg45 solver;
        // DormandPrince54 solver;
        int k = 3000;
        Result result;
        // for( int i = 0; i < 4; i++ ){
        // solver.setSteps(k);
        // solver.solve(ThreeBody,u0,T1,5);
        // result = solver.getResult();
        // cout << (u0-result.U.back()).lInfNorm() << endl;
        // k *= 2;
        // }
        // result.output("test.txt");
        solver->setSteps(k);
        solver->solve(ThreeBody,u0,T1,2);
        result = solver->getResult();
        cout << (u0-result.U.back()).lInfNorm() << endl;
        // cout << max(abs(u0(0)-result.U.back()(0)),abs(u0(1)-result.U.back()(1)));
        result.output("test.txt");
    }
}