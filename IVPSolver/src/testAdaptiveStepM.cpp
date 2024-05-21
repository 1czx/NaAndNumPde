#include "TimeIntegrator.h"
#include "TimeIntegratorFactory.h"
#include "json/json.h"
#include "catch.hpp"


TEST_CASE("Adaptive Step Method", "[TimeIntegrator]"){

    auto ThreeBody = [](VectorXd u, double t) -> VectorXd{
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
    };

    cout << "Test for Adaptive Step Methods" << endl;
    string inputDir = "../input/";
    string outPutDir = "../result/";
    TimeIntegratorFactory & Fac = TimeIntegratorFactory::getInstance();
    TimeIntegrator * solver;
    Json::Reader reader;
	Json::Value root;
    ifstream in(inputDir+"Test.json", std::ios::binary);
    REQUIRE( in.is_open() == 1 );
    REQUIRE( reader.parse(in,root) == 1 );
    Json::Value & Methods = root["Adaptive"];
    int n = Methods.size();
    
    SECTION( "(10.190)" ){
        VectorXd u0(vector<double>{0.994,0,0,0,-2.0015851063790825224,0});
        double T1 = 17.06521656015796;
        for( int i = 0 ; i < n; i++ ){
            Json::Value & method = Methods[i];
            string ID = method["MethodName"].asString();
            SECTION( ID ){
                solver = Fac.createTimeIntegrator(ID);
                REQUIRE( solver != nullptr );
            }
        }
    }

    SECTION( "(10.191)" ){   
        VectorXd u0(vector<double>{0.879779227778,0,0,0,-0.379677780949,0});
        double T1 = 19.140540691377;
        for( int i = 0 ; i < n; i++ ){
            Json::Value & method = Methods[i];
            string ID = method["MethodName"].asString();
            SECTION( ID ){
                solver = Fac.createTimeIntegrator(ID);
                REQUIRE( solver != nullptr );
            }
        }
    }
}
