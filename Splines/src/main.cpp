#include<iostream>
#include<iomanip>
#include"splines.h"
#include<fstream>
#include"json/json.h"
using namespace std;

class f1:public Function<1>{
    public:
    virtual double val( const double & x ){
        return 1.0/(1+25*x*x);
    }
};

class f2:public Function<1>{
    public:
    virtual double val( const double & x ){
        return 1.0/(1+x*x);
    }
};

class sigmoid:public Function<1>{
    private:
        double a;
    public:
    sigmoid( const double & _a ):a{_a} {} 
    virtual double val( const double & x ){
        return 1.0/(1+exp(-a*x));
    }
};

class Heart:public Function<2>{
    public:
    virtual Vector<double,2> val( const double & x ){
        Vector<double,2> p;
        p(0) = sqrt(3)*sin(x);
        p(1) = 2.0/3.0*(sqrt(abs(p(0)))+sqrt(3)*cos(x));
        return p;
    }
};

int intervalNumberForPlot = 1000;
string dir{"../result/"};

void ProblemA( const double & left, const double & right, const vector<int> & knotsNumber ){
    cout << "Problem A" << endl;
    f1 f;
    ofstream fout,fout1,fout2;
    vector<vector<double> > errors(5);
    ppSpline<3> ppS;
    double length = right - left;
    BCType bctype;
    double x;
    int i = 0;
    fout1.open(dir+"problemAErrors");
    fout2.open(dir+"problemAErrors(*N^4)");
    for( auto & N: knotsNumber ){
        fout1 << N << " ";
        fout2 << N << " ";
    }
    fout1 << endl;
    fout2 << endl;
    for( bctype = complete; bctype <= periodic; bctype = BCType(bctype+1) ){
        string type;
        switch(bctype){
            case 1:
                type = "complete";
                break;
            case 2:
                type = "nature";
                break;
            case 3:
                type = "second";
                break;
            case 4:
                type = "notAKnot";
                break;
            case 5:
                type = "periodic";
                break;
            default:
                break;
        }
        fout.open(dir+"problemA("+type+")");
        for( int i = 0; i <= intervalNumberForPlot; i++ ){
                x = left + double(i)*length/intervalNumberForPlot;
                fout << x << " ";
        }
        fout << endl;
        for( int i = 0; i <= intervalNumberForPlot; i++ ){
                x = left + double(i)*length/intervalNumberForPlot;
                fout << f(x) << " ";
        }
        fout << endl;
        for( auto & N: knotsNumber ){
            vector<double> knots;
            for( int i = 0; i < N; i++ ) knots.push_back(left+i*length/(N-1));
            ppS.setKnots(knots);
            ppS.fitCurve(f, bctype);
            for( int i = 0; i <= intervalNumberForPlot; i++ ){
                x = left + double(i)*length/intervalNumberForPlot;
                fout << ppS(x) << " ";
            }
            fout << endl;
            double maxError = 0;
            for( int i = 1; i < N; i++ ){
                x = (knots[i-1]+knots[i])/2;
                double error = abs(f(x)-ppS(x));
                maxError = error > maxError ? error : maxError;
            }
            errors[i].push_back(maxError);
        }
        for( int j = 0; j < knotsNumber.size(); j++ ){
            fout1 << errors[i][j] << " " ;
            fout2 << errors[i][j]*pow(knotsNumber[j],4) << " ";
        }
        fout1 << endl;
        fout2 << endl;
        i++;
        fout.close();
    }
    cout << "Max-norm of the interpolation error vector at mid-points of the subintervals:" << endl;
    cout << setw(8) << std::left << " " << setw(12) << std::left << "complete" << setw(12) << std::left << "nature" << setw(12) << std::left << "second" << setw(12) << std::left << "notAKnot" << setw(12) << std::left << "periodic" << endl;
    for( int i = 0; i < knotsNumber.size(); i++ ){
        cout << "N = " << setw(4) << std::left << knotsNumber[i];
        for ( int j = 0; j < 5; j++ ) cout << setw(12) << std::left << errors[j][i];
        cout << endl;
    }
    cout << "Normalized(error*N^4) max-norm of the interpolation error vector at mid-points of the subintervals:" << endl;
    cout << setw(8) << std::left << " " << setw(12) << std::left << "complete" << setw(12) << std::left << "nature" << setw(12) << std::left << "second" << setw(12) << std::left << "notAKnot" << setw(12) << std::left << "periodic" << endl;
    for( int i = 0; i < knotsNumber.size(); i++ ){
        cout << "N = " << setw(4) << std::left << knotsNumber[i];
        for ( int j = 0; j < 5; j++ ) cout << setw(12) << std::left << errors[j][i]*pow(knotsNumber[i],4);
        cout << endl;
    }

    cout << "The data for plot has been stored in the dir " << dir << endl;

    cout << endl;
}

void ProblemBCD( const double & left, const double & right, const vector<double> & testSites ){
    cout << "ProblemBCD" << endl;
    double length = right - left;
    f2 f;
    CardinalBSpline<2> CBS2;
    CardinalBSpline<3> CBS3;
    CBS2.setInterval(left,right);
    CBS3.setInterval(left,right);
    CBS2.fitCurve(f);
    CBS3.fitCurve(f);
    double x;
    ofstream fout(dir+"problemBCD");
    for( int i = 0; i <= intervalNumberForPlot; i++ ){
        x = left + double(i)*length/intervalNumberForPlot;
        fout << x << " " << f(x) << " " << CBS2(x) << " " << CBS3(x) << endl;
    }
    fout.close();
    cout << "The data for plot has been stored in the dir " << dir << endl;
    cout << "Absolute error of quadratic and cubic cardinal B-splines:" << endl;
    cout << setw(12) << std::left << "x" << setw(16) << std::left << "CBS2(x)" << " " << setw(16) << std::left << "CBS3(x)" << endl;  
    for( auto & x: testSites ){
        cout << "x = " << setw(4) << std::left << x << "    " << setw(16) << std::left << abs(f(x)-CBS2(x)) << " " << setw(16) << std::left << abs(f(x)-CBS3(x)) << endl;
    }
    fout.close();
    cout << endl;
}

void plotHeart(){
    Heart h;
    double t;
    ofstream fout(dir+"Heart");
    Vector<double,2> p;
    for( int i = 0; i <= intervalNumberForPlot; i++ ){
        t = 2*double(i)*M_PIf64/intervalNumberForPlot;
        p = h(t);
        fout << p(0) << " " << p(1) << endl; 
    }
    fout.close();
}

void ProblemE( const double & a, const vector<int> & intervalNumber ){
    cout << "ProblemE";
    string str{""};
    if( a!= 0 ){
        cout << ", alpha = " << a;
        str = "alpha="+to_string(a);
    }
    cout << endl;
    Heart h;
    sigmoid f(a);
    double length = M_PIf64;
    vector<double> pointsPara;
    vector<double> pointsPara1;
    vector<double> pointsPara2;
    BSpline<3,2> BS32;
    ofstream fout;
    for( auto & n: intervalNumber ){
        pointsPara.clear();
        pointsPara1.clear();
        pointsPara2.clear();
        double t;
        pointsPara.push_back(0);
        pointsPara1.push_back(0);
        for( int i = 1; i <= n/2-1 ; i++ ){
            if( a == 0 ) t = 2*length*double(i)/n;
            else t = length*f(2*length*double(i)/n-M_PI_2f64);
            pointsPara.push_back(t);
            pointsPara1.push_back(t);
        }
        pointsPara1.push_back(length);
        pointsPara2.push_back(length);
        pointsPara.push_back(length);
        for( int i = n/2+1; i <= n-1; i++ ){
            if( a == 0 ) t = 2*length*double(i)/n;
            else t = length+length*f(length*double(2*i-n)/n-M_PI_2f64);
            pointsPara.push_back(t);
            pointsPara2.push_back(t);
        }
        pointsPara.push_back(2*length);
        pointsPara2.push_back(2*length);
        double x, y;
        Vector<double,2> p;
        fout.open(dir+"problemE(n="+to_string(n)+",pieses=1)"+str);
        double right = BS32.fitCurve( h, pointsPara ,nature );
        for( int i = 0; i <= intervalNumberForPlot; i++ ){
            t = double(i)*right/intervalNumberForPlot;
            p = BS32(t);
            x = p(0);
            y = p(1);
            if( x != 0 && y != 0 ) fout << x << " " << y << endl;
        }
        fout.close();
        fout.open(dir+"problemE(n="+to_string(n)+",pieses=2)"+str);
        right = BS32.fitCurve(h,pointsPara1,nature);
        for( int i = 0; i <= intervalNumberForPlot; i++ ){
            t = double(i)*right/intervalNumberForPlot;
            p = BS32(t);
            x = p(0);
            y = p(1);
            if( x != 0 && y != 0 ) fout << x << " " << y << endl;
        }
        right = BS32.fitCurve(h,pointsPara2,nature);
        for( int i = 0; i <= intervalNumberForPlot; i++ ){
            t = double(i)*right/intervalNumberForPlot;
            p = BS32(t);
            x = p(0);
            y = p(1);
            if( x != 0 && y != 0 ) fout << x << " " << y << endl;
        }
        fout.close();
    }
    cout << "The data for plot has been stored in the dir " << dir << endl;
    cout << endl;
}

int main(){
    Json::Reader reader;
	Json::Value root;
    ifstream in("parameter.json", ios::binary);
    if( !in.is_open() ){
		cout << "Error opening file\n";
		return 0;
	}
    if( reader.parse(in,root)){
        dir = root["outputDir"].asString();
        intervalNumberForPlot = root["intervalNumberForPlot"].asInt();
        const Json::Value A = root["ProblemA"], BCD = root["ProblemBCD"], E = root["ProblemE"];
        if( A["isTest"].asBool() ){
            vector<int> knotsNumber;
            Json::Value temp = A["knotsNumber"];
            double left = A["left"].asDouble(), right = A["right"].asDouble();
            for( int i = 0; i < temp.size(); i++ ) knotsNumber.push_back(temp[i].asInt());
            ProblemA(left,right,knotsNumber);
        }
        if( BCD["isTest"].asBool() ){
            vector<double> testSites;
            Json::Value temp = BCD["testSites"];
            double left = BCD["left"].asDouble(), right = BCD["right"].asDouble();
            for( int i = 0; i < temp.size(); i++ ) testSites.push_back(temp[i].asDouble());
            ProblemBCD(left,right,testSites);
        }
        if( E["isTest"].asBool() ){
            vector<int> intervalNumbers;
            Json::Value temp = E["intervalNumbers"];
            for( int i = 0; i < temp.size(); i++ ) intervalNumbers.push_back(temp[i].asInt());
            double a = E["alpha"].asDouble();
            plotHeart();
            ProblemE(0,intervalNumbers);
            ProblemE(a,intervalNumbers);
        }
    }else{
        cout << "Parse Error\n";
    }

    in.close();
}