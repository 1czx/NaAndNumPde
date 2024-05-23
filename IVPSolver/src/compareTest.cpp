#include "TimeIntegrator.h"
#include "TimeIntegratorFactory.h"
#include "ThreeBody.h"
#include "json/json.h"
#include "catch.hpp"
#include <iomanip>

TEST_CASE("Compare Different Method", "[TimeIntegrator]"){
    string inputDir = "../input/";
    string outPutDir = "../result/";
    TimeIntegratorFactory & Fac = TimeIntegratorFactory::getInstance();
    TimeIntegrator * solver;
    Json::Reader reader;
	Json::Value root;
    ifstream in(inputDir+"CompareMethod.json", std::ios::binary);
    REQUIRE( in.is_open() == 1 );
    REQUIRE( reader.parse(in,root) == 1 );
    VectorXd u0(vector<double>{0.994,0,0,0,-2.0015851063790825224,0});
    double T1 = 17.06521656015796;

    SECTION( "Compare Performance" ){
        cout << "Compare Performance for Different Method with different Steps" << endl;
        Json::Value & methods = root["CompareStep"];
        int n = methods.size();
        int s = 0;
        for( int i = 0 ; i < n; i++ ){
            Json::Value & method = methods[i];
            string ID = method["MethodName"].asString();
            solver = Fac.createTimeIntegrator(ID);
            if( solver->getStepType() == fixed ){
                s = method["s"].asInt();
            }else{
                Ea = method["Ea"].asDouble();
                Er = method["Er"].asDouble();
            }
            solver->setSteps(method["step"].asInt());
            solver->solve(ThreeBody,u0,T1,s);
            Result result = solver->getResult();
            if( solver->getStepType() == fixed && s != 0 ) ID = ID + "(s=" + std::to_string(s) + ")";
            result.output(outPutDir+ID+"_Step="+std::to_string(result.step)+"_P");
            cout << std::setw(15) << std::left << ID << ": step = " << std::setw(5) << std::left << result.step << ", Max-norm Error = " << (result.U.back()-u0).lInfNorm() << "\n"; 
        }

    }

    SECTION( "Compare Time Performance" ){
        cout << "Compare CPU Time Performance for Different Method Achieving Max-norm Error of 1e-3" << endl;
        Json::Value & methods = root["CompareTime"];
        int n = methods.size();
        int s;
        std::chrono::duration<double, std::milli> tm;
        for( int i = 0 ; i < n; i++ ){
            Json::Value & method = methods[i];
            string ID = method["MethodName"].asString();
            solver = Fac.createTimeIntegrator(ID);
            if( solver->getStepType() == fixed ){
                s = method["s"].asInt();
            }else{
                Ea = method["Ea"].asDouble();
                Er = method["Er"].asDouble();
            }
            solver->setSteps(method["step"].asInt());
            auto start = std::chrono::high_resolution_clock::now();
            solver->solve(ThreeBody,u0,T1,s);
            auto end = std::chrono::high_resolution_clock::now();
            tm = end - start;
            Result result = solver->getResult();
            if( solver->getStepType() == fixed && s != 0 ) ID = ID + "(s=" + std::to_string(s) + ")";
            cout << std::setw(22) << std::left << ID << ": step = " << result.step << ", CPU Time =" << tm.count() << "ms, Max-norm Error = " << (result.U.back()-u0).lInfNorm() << "\n"; 
        }
    }

    cout << "\n";
}
