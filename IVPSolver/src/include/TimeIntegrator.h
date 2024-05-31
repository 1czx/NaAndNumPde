#pragma once
#include "Matrix.h"
#include <functional>

using RHS = function<VectorXd(VectorXd,double)>;

inline double Ea = 1e-5;
inline double Er = 1e-7;

enum MethodType{LLM,OneStep};
enum StepType{fixed,adaptive};

const vector<vector<double>> ABMcoef{{1}, {3.0/2, -0.5}, {23.0/12, -16.0/12, 5.0/12}, {55.0/24, -59.0/24, 37.0/24, -9.0/24}};
const vector<vector<double>> AMMcoef{{0.5, 0.5}, {5.0/12, 8.0/12, -1.0/12}, {9.0/24, 19.0/24, -5.0/24, 1.0/24}, {251.0/720, 646.0/720, -264.0/720, 106.0/720, -19.0/720}};
const vector<vector<double>> BDFcoef{{1, -1}, {2.0/3, -4.0/3, 1.0/3}, {6.0/11, -18.0/11, 9.0/11, -2.0/11}, {12.0/25, -48.0/25, 36.0/25, -16.0/25, 3.0/25}};
const vector<double> CRKc{0,0.5,0.5,1};
const vector<double> CRKb{1.0/6,1.0/3,1.0/3,1.0/6};

class Result{
    public:
    
    void clear(){
        U.clear();
    }
    
    void output( const string & filename ){
        ofstream fout(filename);
        int n = U[0].size();
        for( auto & u: U ){
            for( int i = 0; i < n; i++ ) fout << u(i) << " ";
            fout << endl;
        }
        fout.close();
    }

    int step;
    double time;
    vector<VectorXd> U;
};

class TimeIntegrator{

    public:

    TimeIntegrator():steps{0},currentStep{0},currentTime{0},MT{OneStep},ST{fixed} {}
    TimeIntegrator( const MethodType & mt ):steps{0},currentStep{0},currentTime{0},MT{mt},ST{fixed}{}
    TimeIntegrator( const StepType & st ):steps{0},currentStep{0},currentTime{0},MT{OneStep},ST{st}{}
    void setSteps( const int & s ){ steps = s; } 
    MethodType getMethodType() { return MT; }
    StepType getStepType() { return ST; }

    Result getResult(){ return result; }
    //s is step for LLM, or stage for one-step method, while s = 0 representing that step or stage for this method is not optional
    void solve( RHS f, VectorXd u0, const double & T, int s = 0 ){
        result.clear();
        currentStep = 0;
        currentTime = 0;
        result.U.push_back(u0);
        double k = T/steps;
        if( MT == LLM ){
            for( ; currentStep < s-1; currentStep++ ){
                VectorXd u(result.U.back());
                VectorXd y = f(u,currentStep*k);
                for( int i = 0; i < 4; i++ ){
                    y = f(result.U.back()+k*CRKc[i]*y,currentStep*k+CRKc[i]*k);
                    u = u + k*CRKb[i]*y;
                }
                result.U.push_back(u);
            }
        }
        if( ST == fixed ){ 
            for( ; (currentStep < steps) ; currentStep++ ){ 
                result.U.push_back( oneStep( f, k, s ) );
            }
        }else for( ; (currentTime < T) ; currentStep++ ) result.U.push_back( oneStep( f, k, s ) );
        // cout << currentStep << " " << currentTime <<endl;
        result.step = currentStep;
        result.time = currentTime;
    }

    protected:

    virtual VectorXd oneStep( RHS & f, double & k, int s ) = 0;

    Result result;
    MethodType MT;
    StepType ST;
    int steps;
    int currentStep;
    double currentTime;
};

class ABMs:public TimeIntegrator{

    public:
    ABMs():TimeIntegrator(LLM) {}

    private:
    
    VectorXd oneStep( RHS & f, double & k, int s ){
        if( s < 1 || s > 4 ){
            result.clear();
            currentStep = steps;
            std::cerr<< "Error: " << s << " step ABMs is not implemented!" << endl;
        }
        VectorXd u(result.U.back());
        for( int i = 0; i < s; i++ ) u = u + k*ABMcoef[s-1][i]*f(result.U[currentStep-i],(currentStep-i)*k);
        return u;
    }

};

class AMMs:public TimeIntegrator{

    public:
    AMMs():TimeIntegrator(LLM) {}

    private:
    VectorXd oneStep( RHS & f, double & k, int s ){
        if( s < 1 || s > 4 ){
            result.clear();
            currentStep = steps;
            std::cerr<< "Error: " << s << " step AMMs is not implemented!" << endl;
        }
        VectorXd u(result.U.back());
        for( int i = 0; i < s; i++ ) {u = u + k*ABMcoef[s-1][i]*f(result.U[currentStep-i],(currentStep-i)*k);}
        u = result.U.back() + k*AMMcoef[s-1][0]*f(u,(currentStep+1)*k);
        for( int i = 0; i < s; i++ ) {u = u + k*AMMcoef[s-1][i+1]*f(result.U[currentStep-i],(currentStep-i)*k);}
        return u;
    }

};

class BDFs:public TimeIntegrator{
    public:
    BDFs():TimeIntegrator(LLM) {}
    
    private:
    VectorXd oneStep( RHS & f, double & k, int s ){
        if( s < 1 || s > 4 ){
            result.clear();
            currentStep = steps;
            std::cerr<< "Error: " << s << " step BDFs is not implemented!" << endl;
        }
        VectorXd u(result.U.back());
        VectorXd temp(u);
        double esp = 2.2e-16;
        int iter = 0;
        while(iter < 1000){
            u = k*BDFcoef[s-1][0]*f(temp,(currentStep+1)*k);
            for( int i = 1; i <= s; i++ ) u = u - BDFcoef[s-1][i]*result.U[currentStep-i+1];
            if( (u - temp).lInfNorm() <= esp ) break;
            temp = u;
            iter++;
        }
        return u;
    }
};

class classicalRK:public TimeIntegrator{

    public:

    private:
    VectorXd oneStep( RHS & f, double & k, int s ){
        VectorXd u(result.U.back());
        VectorXd y = f(u,currentStep*k);
        for( int i = 0; i < 4; i++ ){
            y = f(result.U.back()+k*CRKc[i]*y,currentStep*k+CRKc[i]*k);
            u = u + k*CRKb[i]*y;
        }
        return u;
    }

};

class ESDIRK:public TimeIntegrator{

    public:

    ESDIRK():coef{ {1.0/2, 1.0/4, 1.0/4},
                   {83.0/250, 8611.0/62500, -1743.0/31250, 1.0/4},
                   {31.0/50, 5012029.0/34652500, -654441.0/2922500, 174375.0/388108, 1.0/4},
                   {17.0/20, 15267082809.0/155376265600, -71443401.0/120774400, 730878875.0/902184768, 2285395.0/8070912, 1.0/4},
                   {1, 82889.0/524892, 0, 15625.0/83664, 69875.0/102672, -2260.0/8211, 1.0/4} } {}

    private:

    VectorXd oneStep( RHS & f, double & k, int s ){
        VectorXd u(result.U.back());
        vector<VectorXd> y(6,f(u,currentStep*k));
        vector<VectorXd> temp(6,f(u,currentStep*k));
        VectorXd e(6);
        double eps = 2.2e-16;
        int iter = 0;
        while( iter < 1000 ){
            for( int i = 2; i <= 6; i++ ){
                VectorXd yt(u.size());
                for( int j = 1; j <= i; j++ ) yt = yt + coef[i-2][j]*y[j-1];
                y[i-1] = f(u+k*yt,currentStep*k+coef[i-2][0]*k);
                e(i-1) = (y[i-1]-temp[i-1]).lInfNorm();
            }
            if( e.lInfNorm() <= eps ) break;
            temp = y;
            iter++;
        }
        for( int i = 0; i < 6; i++ ) u = u + k*coef[4][i+1]*y[i];
        return u;
    }

    vector<vector<double>> coef;

};

class GaussLegendreRKMs:public TimeIntegrator{

    public:

    GaussLegendreRKMs():coef{{{1}, {1.0/2, 1.0/2}}} {
        double w1 = 1.0/4, w2 = sqrt(3)/6.0;
        vector<vector<double>> temp{
            {2*w1, 2*w1},
            {2*w1-w2, w1, w1-w2},
            {2*w1+w2, w1+w2, w1}
        };
        coef.push_back(temp);
        double w3, w4;
        w1 = 5.0/36; w2 = 2.0/9; w3 = sqrt(15)/30.0; w4 = sqrt(15)/24.0;
        temp = vector<vector<double>>{
            {2*w1, 2*w2, 2*w1},
            {2*w1+w2-3*w3, w1, w2-2*w3, w1-w3},
            {1.0/2, w1+w4, w2, w1-w4},
            {2*w1+w2+3*w3, w1+w3, w2+2*w3, w1}
        };
        coef.push_back(temp);
        double w5, v1, v2, v3 ,v4 ,v5;
        w1 = 1.0/8 - sqrt(30)/144.0; v1 = 1.0/8 + sqrt(30)/144.0;
        w2 = 1.0/2*sqrt((15+2*sqrt(30))/35.0);v2 = 1.0/2*sqrt((15-2*sqrt(30))/35.0);
        w3 = w2*(1.0/6+sqrt(30)/24.0); v3 = v2*(1.0/6-sqrt(30)/24.0);
        w4 = w2*(1.0/21+5*sqrt(30)/168.0); v4 = v2*(1.0/21-5*sqrt(30)/168.0);
        w5 = w2 - 2*w3; v5 = v2 - 2*v3; 
        temp = vector<vector<double>>{
            {2*w1, 2*v1, 2*v1, 2*w1},
            {1.0/2-w2, w1, v1-w3+v4, v1-w3-v4, w1-w5},
            {1.0/2-v2, w1-v3+w4, v1, v1-v5, w1-v3-w4},
            {1.0/2+v2, w1+v3+w4, v1+v5, v1, w1+v3-w4},
            {1.0/2+w2, w1+w5, v1+w3+v4, v1+w3-v4, w1}
        };
        coef.push_back(temp);

        w1 = sqrt(5.0/36+sqrt(70)/126);w2 = sqrt(5.0/36-sqrt(70)/126); 
        vector<double> c{0.5-w1,0.5-w2,0.5,0.5+w2,0.5+w1};
        VectorXd a(25);
        MatrixXd A(25,25);
        for( int m = 0; m < 5; m++ ){
            for( int i = 0; i < 5; i++ ){
                for( int j = 0; j < 5; j++ ){
                    A(i*5+m,i*5+j) = pow(c[j],m);
                }
                a(i*5+m) = pow(c[i],m+1)/(m+1);
            }
        }
        a = A.solve(a);
        temp.clear();
        vector<double> b;
        for( int i = 0; i < 5; i++ ) b.push_back(2*a(i*5+i));
        temp.push_back(b);
        for( int i = 0; i < 5; i++ ){
            vector<double> t{c[i]};
            for( int j = 0; j < 5; j++ ) t.push_back(a(i*5+j));
            temp.push_back(t);
        }
        coef.push_back(temp);
    }

    private:
    VectorXd oneStep( RHS & f, double & k, int s ){
        if( s < 1 || s > 5 ){
            result.clear();
            currentStep = steps;
            std::cerr<< "Error: " << s << " stage Gauss-Legendre RKMs is not implemented!" << endl;
        }
        vector<vector<double>> table = coef[s-1];
        VectorXd u(result.U.back());
        vector<VectorXd> y(s,f(u,currentStep*k));
        vector<VectorXd> temp(s,f(u,currentStep*k));
        VectorXd e(s);
        double eps = 2.2e-16;
        int iter = 0;
        while(iter < 1000){
            for( int i = 1; i <= s; i++ ){
                VectorXd yt(u.size());
                for( int j = 1; j <= s; j++ ) yt = yt + table[i][j]*y[j-1];
                y[i-1] = f(u+k*yt,currentStep*k+table[i][0]*k);
                e(i-1) = (y[i-1]-temp[i-1]).lInfNorm();
            }
            if( e.lInfNorm() <= eps ) break;
            temp = y;
            iter++;
        }
        for( int i = 0; i < s; i++ ) u = u + k*table[0][i]*y[i];
        return u; 
    }

    vector<vector<vector<double>>> coef;

};

class EmbeddedRK:public TimeIntegrator{
    
    public:
    EmbeddedRK():TimeIntegrator(adaptive) {}

    protected:
    VectorXd oneStep( RHS & f, double & k, int s ){
        VectorXd u1(result.U.back());
        VectorXd u2(result.U.back());
        s = coef.size()-1;
        vector<VectorXd> y;
        VectorXd e = VectorXd(vector<double>(u1.size(),Ea)) + Er*u1.Abs();
        int q = min(coef[0][0],coef[1][0]);
        double tk = k;
        while(1){
            y = vector<VectorXd>(s,vector<double>(u1.size(),0));
            y[0] = f(result.U.back(),currentTime);
            u1 = result.U.back() + coef[0][1]*k*y[0];
            u2 = result.U.back() + coef[1][1]*k*y[0];
            for(int i = 1; i < s; i++ ){
                VectorXd yt(u1.size());
                for(int j = 1; j <= i; j++ ) y[i] = y[i] + coef[i+1][j]*y[j-1];
                y[i] = f(result.U.back()+k*y[i],currentTime+coef[i+1][0]*k);
                u1 = u1 + coef[0][i+1]*k*y[i];
                u2 = u2 + coef[1][i+1]*k*y[i];
            }
            VectorXd & E = coef[0][0] > coef[1][0]? u2 : u1;
            E = coef[0][0] > coef[1][0]? (E - u1): (E - u2);
            // u2 = u1 - u2;
            for( int i = 0; i < u1.size(); i++ ) E(i) = E(i)/e(i);
            double Ei = E.lInfNorm();
            k = k*min(max(0.2, 0.9*pow(1.0/Ei, 1.0/(q+1))), 5.0);
            // if( k > maxK ){
            //     k = maxK;
            //     break;
            // }
            if(Ei <= 1) break; 
            tk = k;
        }
        currentTime += tk;
        return coef[0][0] > coef[1][0]? u1 : u2;
        // return u1;
    }

    vector<vector<double>> coef;
};



class Fehlberg45:public EmbeddedRK{

    public:

    Fehlberg45() {
        coef = vector<vector<double>>{ 
                {4, 25.0/216, 0, 1408.0/2565, 2197.0/4104, -1.0/5, 0},
                {5, 16.0/135, 0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55},
                {1.0/4, 1.0/4},
                {3.0/8, 3.0/32, 9.0/32},
                {12.0/13, 1932.0/2197, -7200.0/2197, 7296.0/2197},
                {1, 439.0/216, -8, 3680.0/513, -845.0/4104},
                {1.0/2, -8.0/27, 2, -3544.0/2565, 1859.0/4104, -11.0/40}    
               };
    }
};

class DormandPrince54:public EmbeddedRK{

    public:

    DormandPrince54() {
        coef = vector<vector<double>>{
            {5, 35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0},
            {4, 5179.0/57600, 0, 7571.0/16695, 393.0/640, -92097.0/339200, 187.0/2100, 1.0/40},
            {1.0/5, 1.0/5},
            {3.0/10, 3.0/40, 9.0/40},
            {4.0/5, 44.0/45, -56.0/15, 32.0/9},
            {8.0/9, 19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729},
            {1, 9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656},
            {1, 35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84}
        };
    }
};
