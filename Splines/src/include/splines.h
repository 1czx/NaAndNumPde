/**
 * @file splines.h
 * @author czx 3210103924
 * @brief implement arbitrary dimension liner and cubic piecewise polynomial splines and arbitrary order B-form splines and one dimension cardinal B splines
 * @version 1.0
 * @date 2024-01-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#pragma once

#include<cmath>
#include"function.h"
#include<vector>
#include<algorithm>
#include"Eigen/Dense"
#include<iostream>

using namespace std;
using namespace Eigen;
/**
 * @brief type of B-form splines
 * 
 */
enum BSplineType{ myDefault1, cardinal };
/**
 * @brief boundary condition for cubic splines
 * 
 */
enum BCType{ myDefault2, complete, nature, second, notAKnot, periodic };

// bool isComputeConditionNumber = 0;

/**
 * @brief compute 2-norm of points
 * 
 * @tparam Dim 
 * @param p a point in the Dim dimension space
 * @return double: 2-norm of points
 */
template<int Dim>
double l2Nrom( const Vector<double,Dim> & p ){
    double norm = 0;
    for( int i = 0; i < Dim ;i++ ) norm += (p(i)*p(i));
    norm = sqrt(norm);
    return norm;
}

/**
 * @brief arbitrary order BSplines for curve in arbitrary dimension
 * 
 * @tparam Order 
 * @tparam Dim 
 * @tparam t type of B-form splines
 */
template<int Order, int Dim = 1, BSplineType t = myDefault1>
class BSpline:public Function<Dim>{
    private:
        /**
         * @brief splines for component of curve
         * 
         */
        array< BSpline<Order,1,t>, Dim > vec;
    public:
        /**
         * @brief default construct a new BSpline object
         * 
         */
        BSpline(){}
        /**
         * @brief fitting a curve by points
         * 
         * @param points a series of points on the curve you want to fit
         * @param bctype boundary condition type
         * @param boundaryCondition boundary condition
         * @return double: the endpoints of cumulative chordal lengths
         */
        double fitCurve( vector<Vector<double,Dim> > & points, const BCType & bctype, const vector<Vector<double,Dim> > & boundaryCondition = vector<Vector<double,Dim> >{} ){
            vector<double> CCL{0};//cumulative chordal lengths
            int n = points.size();
            for(int i = 1; i < n; i++ ){
                Vector<double,Dim> temp = points[i]-points[i-1];
                CCL.push_back( CCL[i-1]+l2Nrom(temp) );
            }
            vec[0].setKnots(CCL);
            for( int i = 1; i < Dim; i++ ) vec[i] = vec[0];
            vector<double> component;
            vector<double> BCComponent;
            for( int i = 0; i < Dim; i++ ){
                for( auto & p: points ) component.push_back(p(i));
                for( auto & bc: boundaryCondition ) BCComponent.push_back(bc(i));
                vec[i].fitCurve( component, bctype, BCComponent );
                component.clear();
                BCComponent.clear();
            }
            return CCL.back();
        }
        /**
         * @brief fitting a curve by function
         * 
         * @param f function you want to fit
         * @param pointsPara knots of parameter of function
         * @param bctype boundary condition type
         * @return double cumulative chordal lengths
         */
        double fitCurve( Function<Dim> & f, const vector<double> & pointsPara, const BCType & bctype ){
            if( Order != 3 && bctype != 0 ){
                cout << "Error: The bound condition is only for cubic-Splines." << endl;
                return 0;
            }
            vector<Vector<double,Dim> > points;
            for( auto & x: pointsPara ) points.push_back(f(x));
            vector<Vector<double,Dim> > boundarycondition;
            switch(bctype){
                case 0:{
                    break;
                }
                case 1:{
                    boundarycondition.push_back(f.diffVal(pointsPara.front()));
                    boundarycondition.push_back(f.diffVal(pointsPara.back()));
                    break;
                }
                case 2:{
                    break;
                }
                case 3:{
                    boundarycondition.push_back(f.diffVal(pointsPara.front(),2));
                    boundarycondition.push_back(f.diffVal(pointsPara.back(),2));
                }
                case 4:{
                    break;
                }
                case 5:{
                    break;
                }
                default:{
                    cout << "Error: This boundary conditions of cubic splines haven't impletemented!" << endl;
                    return 0;
                }
            }
            return fitCurve(points,bctype,boundarycondition);
        }

        virtual Eigen::Vector<double,Dim> val( const double & x ){
            Eigen::Vector<double,Dim> v;
            for( int i = 0; i < Dim; i++ ) v(i) = (vec[i])(x);
            return v;
        }

};

/**
 * @brief specialization for one dimension B-form splines
 * 
 * @tparam Order 
 * @tparam t type of B-form splines
 */
template<int Order, BSplineType t>
class BSpline<Order,1,t>:public Function<1>{
    
    protected:
    // int Dim;
    // int Order;
    // BSplineType t;

    /**
     * @brief the Basis function of B-form splines
     * 
     */
    class BBasis:public Function<1>{
        
        private:
            /**
             * @brief knots of Basis function
             * 
             */
            vector<double> knots;
            /**
             * @brief order of Basis function
             * 
             */
            int order;
        public:
            /**
             * @brief express of function as polynomial in each interval
             * 
             */
            vector<polynomial> polys;
            /**
             * @brief Construct a B splines Basis dependent on knots, initially it is zero order
             * 
             * @param _knots 
             */
            BBasis( const vector<double> & _knots ):polys(int(_knots.size()-1),0), knots{_knots}, order(int(_knots.size()-2)) {
                sort(knots.begin(),knots.end());
                polys[0] = 1;
            }
            // explicit BBasis( const int & order ):polys(order+1,polynomial(order)), knots(order), Order{order} {}

            virtual double val( const double & x ){
                if( x < knots.front() || x > knots.back() ) return 0;
                else if( x == knots.back() ) return (polys.back())(x);
                else{
                    vector<double>::iterator iter = upper_bound(knots.begin(),knots.end(), x);
                    int index = iter - knots.begin() - 1;
                    return (polys[index])(x);
                }
            }

            double diffVal( const double & x, const int & n = 1 ){
                if( x < knots.front() || x > knots.back() ) return 0;
                else if( x == knots.back() ) return (polys.back())(x);
                else{
                    vector<double>::iterator iter = upper_bound(knots.begin(),knots.end(), x);
                    int index = iter - knots.begin() - 1;
                    return polys[index].diffVal(x,n);
                }
            }

        friend class Bspline;
    };
    /**
     * @brief a series of basis of spline on the interpolation knots
     * 
     */
    vector<BBasis> base;
    /**
     * @brief interpolation knots
     * 
     */
    vector<double> knots;
    /**
     * @brief expression of splines as the result of addition with basis multiply the coefficient which computed in the fitCurve
     * 
     */
    vector<polynomial> express;
    /**
     * @brief is the spline set the knots and compute the basis, if it's not, users can't fit the curve
     * 
     */
    bool isInitialized;
    /**
     * @brief is the spline fitting some curve, if it's not, users can't get the value at any points of splines
     * 
     */
    bool isFitted;
    /**
     * @brief extent the knots inputted with extra knots for basis computing
     * 
     * @return vector<double> 
     */
    vector<double> extentKnots(){
        int n = knots.size(); 
        vector<double> extendedKnots;
        for( int i = -Order; i < 0; i++ ) extendedKnots.push_back(knots.front()+i);
        for( int i = 0; i < n; i++ ) extendedKnots.push_back(knots[i]);
        for( int i = 1; i <= Order; i++ ) extendedKnots.push_back(knots.back()+i);
        return extendedKnots;
    }
    /**
     * @brief commpute basis on the extentedKnots, stored in the attribute,base , and set isInitialized as 1 
     * 
     * @param extentedKnots
     */
    void computeBasis( const vector<double> & extentedKnots ){
        int n = extentedKnots.size();
        for( int i = 0; i < n-1; i++ ){
            if( (i+Order+2) <= n  ) base.push_back( BBasis(vector<double>( (extentedKnots.begin() + i), (extentedKnots.begin() + i + Order + 2) ) ) );
            else base.push_back( BBasis(vector<double>( (extentedKnots.begin() + i), extentedKnots.end() ) ) );
        }
        n = base.size();
        for( int k = 1; k <= Order; k++ ){
            for( int i = 0; i < n-k; i++ ){
                    double temp1 = 1.0/(extentedKnots[i+k]-extentedKnots[i]);
                    double temp2 = 1.0/(extentedKnots[i+k+1]-extentedKnots[i+1]);
                    base[i].polys[0] = polynomial({-extentedKnots[i]*temp1,temp1})*base[i].polys[0];
                for ( int j = 1; j <= k; j++ ){
                    base[i].polys[j] = polynomial({-extentedKnots[i]*temp1,temp1})*base[i].polys[j]+polynomial({extentedKnots[i+k+1]*temp2,-temp2})*base[i+1].polys[j-1];
                }
            }
            base.pop_back();
        }

        isInitialized = 1;
    }
    /**
     * @brief creat a matrix with condition that the value of splines is equal to the value of fitted funciton at the interpolation knots
     * 
     * @return MatrixXd 
     */
    MatrixXd baseMat(){
        int N = knots.size();
        MatrixXd A = MatrixXd::Zero(N,N+Order-1);
        for( int i = 0; i < N; i++ )
            for( int j = i; j < i+Order; j++ ) A(i,j) = base[j](knots[i]);
        return A;
    }
    /**
     * @brief creat a matrix for compute coefficient of basis for cubic spline with boundary conditon of n order derivative
     * 
     * @param n order of derivative
     * @return MatrixXd 
     */
    MatrixXd MatOfOrderThree(const int & n){
        int N = knots.size();
        MatrixXd A = MatrixXd::Zero(N+Order-1,N+Order-1);
        for( int j = 0; j < Order; j++ ) A(0, j) = base[j].diffVal(knots[0],n);
        A.block(1,0,N,N+Order-1) = baseMat();
        for( int j = N-1; j < N+Order-1; j++ ) A(N+Order-2,j) = base[j].diffVal(knots.back(),n);
        return A;
    }
    
    public:
    /**
     * @brief default Construct a new BSpline object
     * 
     */
    BSpline():isInitialized{0}, isFitted{0} {}
    /**
     * @brief Construct a BSpline which have setted interpolation knots and computed basis
     * 
     * @param _knots interpolation knots
     */
    BSpline( const vector<double> & _knots ):knots{_knots}, isInitialized{0}, isFitted{0} {
        sort(knots.begin(),knots.end());
        computeBasis(extentKnots());
    }
    /**
     * @brief Set interpolation knots for splines and compute basis
     * 
     * @param _knots interpolation knots
     */
    void setKnots( const vector<double> & _knots ){
        knots = _knots;
        sort(knots.begin(),knots.end());
        base.clear();
        computeBasis(extentKnots());
        if( isFitted ){
            express.clear();
            isFitted = 0;
        }
    }
    /**
     * @brief compute the coefficient of basis and compute the express of spline to interpolate a series of points
     * 
     * @param f the value of curve at the interpolation knots
     * @param bctype boundary condtion type
     * @param bondaryCondition inputted extra condition
     */
    void fitCurve( vector<double> & f, const BCType & bctype, const vector<double> & bondaryCondition = vector<double>{} ){
        if( !isInitialized ){
            cout << "Error: This splines don't set knots!" << endl;
            return;
        }
        if( Order != 3 && bctype != 0 ){
            cout << "Error: The bound condition is only for 3-order-Splines." << endl;
            return;
        }
        int N = knots.size();
        MatrixXd A = MatrixXd::Zero(N+Order-1,N+Order-1);
        VectorXd rhs = VectorXd::Zero(N+Order-1);
        switch(bctype){
            case 0:{
                if( t != cardinal ){
                    // 可以优化 再说，这下不得不优化了, 条件数也太大了
                    // int index = 0;
                    // VectorXd rhs = VectorXd::Zero(N+Order-1);
                    // for( int i = 0; i < Order-1; i++ ){
                    //     for( int j = 0; j < Order; j++ ) A(i,j) = base[j].diffVal(knots[0],i+1);
                    //     rhs(index++) = f.diffVal(knots[0],i+1);
                    // }
                    // A.block(Order-1,0,N,N+Order-1) = baseMat();
                    // for( int i = 0; i < N; i++ ) rhs(index++) = f(knots[i]);
                    
                    //没时间就算了
                }else{
                    switch(Order){
                        case 2:{
                            for( int j = 0; j < Order; j++ ) A(0,j) = base[j](knots[0]);
                            rhs(0) = f[0];
                            for( int i = 1; i < N+Order-2; i++ ){
                                for( int j = i-1; j < i+Order; j++ ) A(i,j) = base[j](knots[i]-0.5);
                            }
                            for( int j = N-1; j < N+Order-1; j++ ) A(N+Order-2,j) = base[j](knots.back());
                            for( int i = 0; i < f.size(); i++ ) rhs(i) = f[i];
                            break;
                        }
                        case 3:{
                            A = MatOfOrderThree(1);
                            rhs(0) = bondaryCondition[0];
                            for( int i = 0; i < N; i++ ) rhs(i+1) = f[i];
                            rhs(N+1) = bondaryCondition[1];
                            break;
                        }
                        default:{
                            cout << "Error: The cardinal B splines, whose order isn't 2 or 3, haven't impletemented!" << endl;             
                            return;
                        }
                    }
                }
                break;
            }
            case 1:{
                if( bondaryCondition.size() != Order-1 ){
                    cout << "Error: The boundary conditions are insufficient!" << endl;
                    return;
                }
                A = MatOfOrderThree(1);
                rhs(0) = bondaryCondition[0];
                for( int i = 0; i < N; i++ ) rhs(i+1) = f[i];
                rhs(N+1) = bondaryCondition[1];
                break;
            }
            case 2:{
                if( bondaryCondition.size() != 0 ){
                    cout << "Error: The number of boundary conditions is wrong!" << endl;
                    return;
                }
                A = MatOfOrderThree(2);
                for( int i = 0; i < N; i++ ) rhs(i+1) = f[i];
                break;
            }
            case 3:{
                if( bondaryCondition.size() != Order-1 ){
                    cout << "Error: The number of boundary conditions is wrong!" << endl;
                    return;
                }
                A = MatOfOrderThree(2);
                rhs(0) = bondaryCondition[0];
                for( int i = 0; i < N; i++ ) rhs(i+1) = f[i];
                rhs(N+1) = bondaryCondition[1];
                break;
            }
            default:{
                cout << "Error: This boundary conditions of cubic splines haven't impletemented!" << endl;
                return;
            }
        }
        // JacobiSVD<MatrixXd> svd;
        // svd.compute(A);
        // cout << svd.singularValues().maxCoeff()/svd.singularValues().minCoeff() << endl;
        // VectorXd rhs = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(condition.data(),condition.size());
        // cout << A << endl << endl << rhs;
        VectorXd coef = A.colPivHouseholderQr().solve(rhs);
        
        express.clear();
        polynomial temp;
        for( int i = 0; i < N-1; i++ ){
            temp = 0;
            for( int j = Order; j >= 0; j-- ) temp = temp + coef(i+Order-j)*base[i+Order-j].polys[j];
            express.push_back(temp);
        }
        
        isFitted = 1;
    }
    /**
     * @brief compute the coefficient of basis and compute the express of spline to interpolate a function
     * 
     * @param f function you want to interpolate
     * @param bctype boundary condition
     */
    void fitCurve( Function<1> & f, const BCType & bctype ){
        if( !isInitialized ){
            cout << "Error: This splines don't set knots!" << endl;
            return;
        }
        if( Order != 3 && bctype != 0 ){
            cout << "Error: The bound condition is only for 3-order-Splines." << endl;
            return;
        }
        int N = knots.size();
        vector<double> fval;
        vector<double> boundaryCondition;
        switch(bctype){
            case 0:{
                if( t != cardinal ){

                }else{
                    switch(Order){
                        case 2:{
                            fval.push_back(f(knots[0]));
                            for( int i = 1; i < N; i++ ) fval.push_back(f((knots[i-1]+knots[i])/2));
                            fval.push_back(f(knots[N-1]));
                            fitCurve(fval,bctype);
                            return;
                        }
                        case 3:{
                            boundaryCondition.push_back(f.diffVal(knots[0]));
                            boundaryCondition.push_back(f.diffVal(knots[N-1]));
                            break;
                        }
                        default:{
                            cout << "Error: The cardinal B splines, whose order isn't 2 or 3, haven't impletemented!" << endl;             
                            return;
                        }
                    }
                }
                break;
            }
            case 1:{
                boundaryCondition.push_back(f.diffVal(knots[0]));
                boundaryCondition.push_back(f.diffVal(knots[N-1]));
                break;
            }
            case 2:{
                break;
            }
            case 3:{
                boundaryCondition.push_back(f.diffVal(knots[0],2));
                boundaryCondition.push_back(f.diffVal(knots[N-1],2));
                break;
            }
            default:{
                cout << "Error: This boundary conditions of cubic splines haven't impletemented!" << endl;
                return;
            }
        }
        for( int i = 0; i < N; i++ ) fval.push_back(f(knots[i]));
        fitCurve(fval,bctype,boundaryCondition);
    }

    virtual double val( const double & x ){
        if( !isFitted ){
            cout << "Error: This splines is not fitting any function!" << endl;
            return -1;
        }
        if( x < knots.front() || x > knots.back() ) return 0;
        else if( x == knots.back() ) return (express.back())(x);
        else{
            vector<double>::iterator iter = upper_bound(knots.begin(),knots.end(), x);
            int index = iter - knots.begin() - 1;
            return (express[index])(x);
        }
    }

    double diffVal( const double & x, const int & n = 1 ){
        if( !isFitted ){
            cout << "Error: This splines is not fitting any function!" << endl;
            return -1;
        }
        if( x < knots.front() || x > knots.back() ) return 0;
        else if( x == knots.back() ) return (express.back())(x);
        else{
            vector<double>::iterator iter = upper_bound(knots.begin(),knots.end(), x);
            int index = iter - knots.begin() - 1;
            return express[index].diffVal(x,n);
        }
    }

};
/**
 * @brief one dimension cardinal B-form splines
 * 
 * @tparam Order 
 */
template<int Order>
class CardinalBSpline:public BSpline<Order,1,cardinal>{
    public:
    /**
     * @brief default construct a new Cardinal B Spline object
     * 
     */
    CardinalBSpline(){}
    /**
     * @brief Construct a new Cardinal B Spline object with interpolation interval [left,right]
     * 
     * @param left 
     * @param right 
     */
    CardinalBSpline( const int & left, const int & right ){
        for( int i = left; i <= right; i++ ) this->knots.push_back(i);
        this->computeBasis(this->extentKnots());
    }
    /**
     * @brief Set the Interval [left,right]
     * 
     * @param left 
     * @param right 
     */
    void setInterval( const int & left, const int & right ){
        this->knots.clear();
        this->base.clear();
        for( int i = left; i <= right; i++ ) this->knots.push_back(i);
        this->computeBasis(this->extentKnots());
        if( this->isFitted ){
            this->express.clear();
            this->isFitted = 0;
        }
    }
    /**
     * @brief use cardinal B Splines to interpolate a function 
     * 
     * @param f function which you want to interpolate 
     */
    void fitCurve( Function<1> & f ){
        this->BSpline<Order,1,cardinal>::fitCurve( f, myDefault2 );
    }
};
/**
 * @brief arbitrary order BSplines for curve in arbitrary dimension
 * 
 * @tparam Order 
 * @tparam Dim 
 */
template<int Order, int Dim = 1 >
class ppSpline:public Function<Dim>{
    private:
        /**
         * @brief piecewise polynomial Splines for component of curve
         * 
         */
        array< ppSpline<Order>, Dim > vec;
    public:
        ppSpline(){}
        /**
         * @brief fitting a curve by points
         * 
         * @param points a series of points on the curve you want to fit
         * @param bctype boundary condition type
         * @param boundaryCondition boundary condition
         * @return double: the endpoints of cumulative chordal lengths
         */
        double fitCurve( vector<Vector<double,Dim> > & points, const BCType & bctype, const vector<Vector<double,Dim> > & boundaryCondition = vector<Vector<double,Dim> >{} ){
            vector<double> CCL{0};//cumulative chordal lengths
            int n = points.size();
            for(int i = 1; i < n; i++ ){
                Vector<double,Dim> temp = points[i]-points[i-1];
                CCL.push_back( CCL[i-1]+l2Nrom(temp) );
            }
            for( int i = 0; i < Dim; i++ ) vec[i].setKnots(CCL);
            vector<double> component;
            vector<double> BCComponent;
            for( int i = 0; i < Dim; i++ ){
                for( auto & p: points ) component.push_back(p(i));
                for( auto & bc: boundaryCondition ) BCComponent.push_back(bc(i));
                vec[i].fitCurve( component, bctype, BCComponent );
                component.clear();
                BCComponent.clear();
            }
            return CCL.back();
        }
        /**
         * @brief fitting a curve by function
         * 
         * @param f function you want to fit
         * @param pointsPara knots of parameter of function
         * @param bctype boundary condition type
         * @return double: cumulative chordal lengths
         */
        double fitCurve( Function<Dim> & f, const vector<double> & pointsPara, const BCType & bctype ){
            if( Order != 3 && bctype != 0 ){
                cout << "Error: The bound condition is only for cubic-Splines." << endl;
                return 0;
            }
            vector<Vector<double,Dim> > points;
            for( auto & x: pointsPara ) points.push_back(f(x));
            vector<Vector<double,Dim> > boundarycondition;
            if( Order == 1 ) return fitCurve( points, bctype );
            if( Order != 3 ){
                cout << "Error: The splines of " << Order << " haven't impletemented!" << endl;
                return 0;
            } 
            switch(bctype){
                case 0:{
                    break;
                }
                case 1:{
                    boundarycondition.push_back(f.diffVal(pointsPara.front()));
                    boundarycondition.push_back(f.diffVal(pointsPara.back()));
                    break;
                }
                case 2:{
                    break;
                }
                case 3:{
                    boundarycondition.push_back(f.diffVal(pointsPara.front(),2));
                    boundarycondition.push_back(f.diffVal(pointsPara.back(),2));
                }
                case 4:{

                    break;
                }
                case 5:{
                    break;
                }
                default:{
                    cout << "Error: This boundary conditions of cubic splines haven't impletemented!" << endl;
                    return 0;
                }
            }
            return fitCurve(points,bctype,boundarycondition);
        }

        virtual Eigen::Vector<double,Dim> val( const double & x ){
            Eigen::Vector<double,Dim> v;
            for( int i = 0; i < Dim; i++ ) v(i) = (vec[i])(x);
            return v;
        }
};
/**
 * @brief specialization for one dimension piecewise polynomial splines
 * 
 * @tparam Order 
 */
template<int Order>
class ppSpline<Order,1>:public Function<1>{
    private:
    // int Order;
    /**
     * @brief interpolation knots
     * 
     */
    vector<double> knots;
    /**
     * @brief expression of splines as a piecewise polynomials
     * 
     */
    vector<polynomial> express;
    /**
     * @brief is the spline set the knots, if it's not, users can't fit the curve
     * 
     */
    bool isInitialized;
    /**
     * @brief is the spline fitting some curve, if it's not, users can't get the value at any points of splines
     * 
     */
    bool isFitted;

    public:
    /**
     * @brief default construct a new pp Spline object
     * 
     */
    ppSpline():isInitialized{0}, isFitted{0} {}
    /**
     * @brief Construct a PP Spline which have setted interpolation knots
     * 
     * @param _knots interpolation knots
     */
    ppSpline( vector<double> _knots ): knots{_knots}, isInitialized{1}, isFitted{0} {
        sort(knots.begin(),knots.end());
    }
    /**
     * @brief Set interpolation knots for splines
     * 
     * @param _knots interpolation knots
     */
    void setKnots( const vector<double> & _knots ){
        knots = _knots;
        sort(knots.begin(),knots.end());
        isInitialized = 1;
        if( isFitted ){
            express.clear();
            isFitted = 0;
        }
    }
    /**
     * @brief compute polynomial on each interpolation subinterval to interpolate a series of points
     * 
     * @param fval the value of curve at the interpolation knots
     * @param bctype boundary condtion type
     * @param bondaryCondition inputted extra condition
     */
    void fitCurve( const vector<double> & fval, const BCType & bctype, const vector<double> & boundaryCondition = vector<double>{} ){
        if( !isInitialized ){
            cout << "Error: This splines don't set knots!" << endl;
            return;
        }
        int n = knots.size();
        vector<double> h;
        vector<double> deltaFval;
        for( int i = 1; i < n; i++ ){
            h.push_back(knots[i]-knots[i-1]);
            deltaFval.push_back(fval[i]-fval[i-1]);
        } 
        switch(Order){
            case 1:{
                vector<double> coef;
                for( int i = 0; i < n-1; i++ ){
                    coef[1] = deltaFval[i]/h[i];
                    coef[0] = fval[i]-knots[i]*coef[1];
                    express.push_back(coef);
                }
                break;
            }
            case 3:{
                MatrixXd A = MatrixXd::Identity(n,n);
                for( int i = 0 ; i < n-2; i++ ){
                    A(i+1,i) = h[i];
                    A(i+1,i+1) = 2*(h[i]+h[i+1]);
                    A(i+1,i+2) = h[i+1];
                }
                VectorXd rhs = VectorXd::Zero(n);
                for( int i = 1; i < n-1; i++ ) rhs(i) = 6*(deltaFval[i]/h[i]-deltaFval[i-1]/h[i-1]);
                switch(bctype){
                    case 1:{
                        A(0,0) = 2*h[0];
                        A(0,1) = h[0];
                        rhs(0) = 6*(deltaFval[0]/h[0]-boundaryCondition[0]);
                        A(n-1,n-2) = h[n-2];
                        A(n-1,n-1) = 2*h[n-2];
                        rhs(n-1) = 6*(boundaryCondition[1]-deltaFval[n-2]/h[n-2]);
                        break;
                    }
                    case 2:{
                        
                        break;
                    }
                    case 3:{
                        rhs(0) = boundaryCondition[0];
                        rhs(n-1) = boundaryCondition[1];
                        break;
                    }
                    case 4:{
                        A(0,0) = -h[1];
                        A(0,1) = h[1]+h[0];
                        A(0,2) = -h[0];
                        A(n-1,n-3) = -h[n-2];
                        A(n-1,n-2) = h[n-3]+h[n-2];
                        A(n-1,n-1) = -h[n-3];
                        break;
                    }
                    case 5:{
                        A(0,0) = 2*(h[0]+h[n-2]);
                        A(0,1) = h[0];
                        A(0,n-2) = h[n-2];
                        A(n-1,0) = 1;
                        A(n-1,n-1) = -1;
                        rhs(0) = 6*(deltaFval[0]/h[0]-deltaFval[n-2]/h[n-2]);
                        break;
                    }
                    default:{
                        cout << "Error: This boundary conditions of cubic splines haven't impletemented!" << endl;
                        return;
                    }
                }

                VectorXd M = A.colPivHouseholderQr().solve(rhs);

                vector<double> coef(4);
                for( int i = 1; i < n; i++ ){
                    coef[3] = (M(i)-M(i-1))/(6*h[i-1]);
                    coef[2] = M(i)/2-3*knots[i]*coef[3];
                    coef[1] = deltaFval[i-1]/h[i-1] - (knots[i]+knots[i-1])*coef[2] - (knots[i]*knots[i]+knots[i]*knots[i-1]+knots[i-1]*knots[i-1])*coef[3];
                    coef[0] = fval[i];
                    double temp = 0;
                    for( int j = 3; j > 0; j-- ){
                        temp += coef[j];
                        temp *= knots[i];
                    }
                    coef[0] -= temp;
                    express.push_back(coef);
                }

                break;
            }
            default:{
                cout << "Error: The splines of " << Order << " haven't impletemented!" << endl;
                return;
            }
        }
        isFitted = 1;

    }
    /**
     * @brief compute polynomial on each interpolation subinterval to interpolate a function
     * 
     * @param f function that you want to interpolate
     * @param bctype boundary condtion type
     */
    void fitCurve( Function<1> & f, const BCType & bctype ){
        if( !isInitialized ){
            cout << "Error: This splines don't set knots!" << endl;
            return;
        }
        int n = knots.size();
        vector<double> fval;
        vector<double> boundaryCondition;
        for( int i = 0; i < n; i++ ) fval.push_back(f(knots[i])); 
        switch(Order){
            case 1:{
                fitCurve(fval,bctype);
            }
            case 3:{
                switch(bctype){
                    case 1:{
                        boundaryCondition.push_back(f.diffVal(knots[0]));
                        boundaryCondition.push_back(f.diffVal(knots[n-1]));
                    }
                    case 2:{
                        break;
                    }
                    case 3:{
                        boundaryCondition.push_back(f.diffVal(knots[0],2));
                        boundaryCondition.push_back(f.diffVal(knots[n-1],2));
                        break;
                    }
                    case 4:{
                        break;
                    }
                    case 5:{
                        break;
                    }
                    default:{
                        cout << "Error: This boundary conditions of cubic splines haven't impletemented!" << endl;
                        return;
                    }
                }
                fitCurve(fval,bctype,boundaryCondition);
            }
        }
    }

    virtual double val( const double & x ){
        if( !isFitted ){
            cout << "Error: This splines is not fitting any function!" << endl;
            return -1;
        }
        if( x < knots.front() || x > knots.back() ) return 0;
        else if( x == knots.back() ) return (express.back())(x);
        else{
            vector<double>::iterator iter = upper_bound(knots.begin(),knots.end(), x);
            int index = iter - knots.begin() - 1;
            return (express[index])(x);
        }
    }

    double diffVal( const double & x, const int & n = 1 ){
        if( !isFitted ){
            cout << "Error: This splines is not fitting any function!" << endl;
            return -1;
        }
        if( x < knots.front() || x > knots.back() ) return 0;
        else if( x == knots.back() ) return (express.back())(x);
        else{
            vector<double>::iterator iter = upper_bound(knots.begin(),knots.end(), x);
            int index = iter - knots.begin() - 1;
            return express[index].diffVal(x,n);
        }
    }
};