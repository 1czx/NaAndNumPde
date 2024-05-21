#include "TimeIntegrator.h"
#include "TimeIntegratorFactory.h"
#include "ThreeBody.h"
#include "json/json.h"
#include "catch.hpp"

TEST_CASE("Adaptive Step Method", "[TimeIntegrator]")
{
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
    vector<double> R{10,1,1e-1};
    SECTION( "(10.190)" ){
        VectorXd u0(vector<double>{0.994,0,0,0,-2.0015851063790825224,0});
        double T1 = 17.06521656015796;
        for( int i = 0 ; i < n; i++ ){
            Json::Value & method = Methods[i];
            string ID = method["MethodName"].asString();
            SECTION( ID ){
                solver = Fac.createTimeIntegrator(ID);
                REQUIRE( solver != nullptr );
                cout << "\n(10.190): " << ID << endl;
                Result result;
                bool flag = 0;
                Ea = 1e-4;
                std::chrono::duration<double, std::milli> tm;
                while( Ea >= 1e-14 ){
                    for( auto & r: R ){
                        Er = Ea*r;
                        solver->setSteps(100);
                        auto start = std::chrono::high_resolution_clock::now();
                        solver->solve(ThreeBody,u0,T1);
                        auto end = std::chrono::high_resolution_clock::now();
                        tm = end - start;
                        result = solver->getResult();
                        if( result.step > 180 && (result.U.back()-u0).lInfNorm() < 1 ){
                            break;
                            flag = 1;
                        }
                    }
                    if( flag == 1 ) break;
                    Ea /= 10;
                }
                cout << "Eabs: " << Ea << ", Erel: " << Er << ", StopTime: "<< result.time << ", step: " << result.step << ", max step length: " << T1/result.step << ", CPU Time: " << tm.count() << "ms" << endl;
                result.output(outPutDir+"(10.190)"+ID+"_Step="+std::to_string(result.step)+"_P");                
                delete solver
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
                cout << "\n(10.190): " << ID << endl;
                Result result;
                bool flag = 0;
                Ea = 1e-4;
                std::chrono::duration<double, std::milli> tm;
                while( Ea >= 1e-14 ){
                    for( auto & r: R ){
                        Er = Ea*r;
                        solver->setSteps(100);
                        auto start = std::chrono::high_resolution_clock::now();
                        solver->solve(ThreeBody,u0,T1);
                        auto end = std::chrono::high_resolution_clock::now();
                        tm = end - start;
                        result = solver->getResult();
                        if( result.step > 280 && (result.U.back()-u0).lInfNorm() < 0.01 ){
                            break;
                            flag = 1;
                        }
                    }
                    if( flag == 1 ) break;
                    Ea /= 10;
                }
                cout << "Eabs: " << Ea << ", Erel: " << Er << ", StopTime: "<< result.time << ", step: " << result.step << ", max step length: " << T1/result.step << ", CPU Time: " << tm.count() << "ms" << endl;
                result.output(outPutDir+"(10.191)"+ID+"_Step="+std::to_string(result.step)+"_P");
                delete solver;           
            }
        }
    }
}
