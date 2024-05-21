#include "TimeIntegrator.h"
#include "TimeIntegratorFactory.h"
#include "json/json.h"
#include "catch.hpp"


TEST_CASE("Fixed Step Method", "[TimeIntegrator]")
{
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
    cout << "Test for Fixed Step Methods" << endl;
    string inputDir = "../input/";
    string outPutDir = "../result/";
    TimeIntegratorFactory & Fac = TimeIntegratorFactory::getInstance();
    TimeIntegrator * solver;
    Json::Reader reader;
	Json::Value root;
    ifstream in(inputDir+"Test.json", std::ios::binary);
    REQUIRE( in.is_open() == 1 );
    REQUIRE( reader.parse(in,root) == 1 );
    Json::Value & Methods = root["Fixed"];
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
                Json::Value & pList = method["Accuracy"];
                int s = 0;
                for( int j = 0; j < pList.size(); j++ ){
                    int p = pList[j].asInt();
                    if( p >= 3 ){
                        SECTION( std::to_string(pList[j].asInt()) ){
                            if( pList.size() > 1 ) s = j+1;
                            cout << "\n(10.190): " << ID << " p = " << p << endl;
                            int step;
                            if( p == 3 ) step = 6000;
                            else if( p == 4 ) step = 3000;
                            else step = 1500;
                            Result result;
                            while( step < 200000 ){
                                solver->setSteps(step);
                                solver->solve(ThreeBody,u0,T1,s);
                                result = solver->getResult();
                                if( (u0-result.U.back()).lInfNorm() < 1 ) break;
                                step *= 2;
                            }
                            REQUIRE( step < 200000 );
                            cout << "Num of Step: "<< step <<  ", max step length: " << T1/step;
                            string filename( outPutDir+"(10.190)"+ID+"_P="+std::to_string(p));
                            result.output(filename+"_Step="+std::to_string(step)+"_P");
                            ofstream fout(filename+"_A");
                            for( int k = 0; k < 5; k++ ){
                                solver->setSteps(step);
                                auto start = std::chrono::high_resolution_clock::now();
                                solver->solve(ThreeBody,u0,T1,s);
                                auto end = std::chrono::high_resolution_clock::now();
                                std::chrono::duration<double, std::milli> tm = end - start;
                                if( k == 0 ) cout << ", CPU Time: " << tm.count() << "ms" << endl;
                                result = solver->getResult();
                                fout << step << " " << (u0-result.U.back()).lInfNorm() << " " << tm.count() << endl ;
                                step*=2;
                            }
                        }
                    }
                }
                delete solver;
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
                Json::Value & pList = method["Accuracy"];
                int s = 0;
                for( int j = 0; j < pList.size(); j++ ){
                    int p = pList[j].asInt();
                    if( p >= 2 ){
                        SECTION( std::to_string(pList[j].asInt()) ){
                            if( pList.size() > 1 ) s = j+1;
                            cout << "\n(10.191): " << ID << " p = " << p << endl;
                            int step = 125;
                            Result result;
                            while( step < 20000 ){
                                solver->setSteps(step);
                                solver->solve(ThreeBody,u0,T1,s);
                                result = solver->getResult();
                                if( (u0-result.U.back()).lInfNorm() < 0.01 ) break;
                                step *= 2;
                            }
                            REQUIRE( step < 20000 );
                            string filename( outPutDir+"(10.191)"+ID+"_P="+std::to_string(p));
                            result.output(filename+"_Step="+std::to_string(step)+"_P");
                            VectorXd Ubar = result.U.back();
                            if( step < 500 ) solver->setSteps(500);
                            auto start = std::chrono::high_resolution_clock::now();
                            solver->solve(ThreeBody,u0,T1,s);
                            auto end = std::chrono::high_resolution_clock::now();
                            result = solver->getResult();
                            std::chrono::duration<double, std::milli> tm = end - start;
                            cout << "Num of Step: "<< max(step,500) <<  ", max step length: " << T1/(max(step,500)) << ", CPU Time: " << tm.count() << "ms" << endl;
                            result.output(filename+"_Step="+std::to_string(step)+"_P");
                            ofstream fout(filename+"_A");
                            for( int k = 0; k < 4; k++ ){
                                step*=2;
                                solver->setSteps(step);
                                start = std::chrono::high_resolution_clock::now();
                                solver->solve(ThreeBody,u0,T1,s);
                                end = std::chrono::high_resolution_clock::now();
                                tm = end - start;
                                result = solver->getResult();
                                VectorXd UR = result.U.back() + (1.0/(pow(2,p)-1))*(result.U.back() - Ubar);
                                fout << step << " " << (UR-Ubar).lInfNorm() << " " << tm.count() << endl ;
                                Ubar = result.U.back();
                            }
                        }
                    }
                }
                delete solver;
            }
        }
    }
    cout << endl;
}


