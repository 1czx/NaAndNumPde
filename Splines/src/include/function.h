/**
 * @file function.h
 * @author czx 3210103924
 * @brief implement a function class and a polynomial class
 * @version 1.0
 * @date 2024-01-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#pragma once

#include<vector>
#include"Eigen/Dense"

/**
 * @brief A function(math) abstract class
 * 
 * @tparam Dim is dimension of function($R->R^{Dim}$)
 */
template<int Dim>
class Function{
    public:
    /**
     * @brief pure virtual function to return the value of function at x
     * 
     * @param x independent variable of function 
     * @return Eigen::Vector<double,Dim> a point in the Dim dimension space
     */
    virtual Eigen::Vector<double,Dim> val( const double & x ) = 0;
    /**
     * @brief override the operator () to use this class in a expression just like function in the math
     * 
     * @param x independent variable of function
     * @return Eigen::Vector<double,Dim> a point in the Dim dimension space
     */
    Eigen::Vector<double,Dim> operator()( const double & x ){
        return val(x);
    }

    /**
     * @brief impletement a numerical derivative for any funcition, user can choose override to make it more precise
     * 
     * @param x independent variable of function
     * @param n order of derivative
     * @return Eigen::Vector<double,Dim> a point in the Dim dimension space
     */
    Eigen::Vector<double,Dim> diffVal( const double & x, const int & n = 1 ){
        double eps = 10e-5;
        
        if( n == 1 ) return (val(x+eps) - val(x-eps))/(2*eps);
        else{
            std::vector<Eigen::Vector<double,Dim>> temp(n+1);
            for( int i = 0; i <= n; i++ ) temp[i] = val(x+(2*i-n)*eps);
            for(int k = n; k > 0; k-- ){
                for( int i = 0; i < n; i++ ) temp[i] = (temp[i+1]-temp[i])/(2*eps);
            }
            return temp[0];
        }    
    }
    
};

/**
 * @brief specialization for one dimension function, definition is as same as high dimension
 * 
 * @tparam 
 */
template<>
class Function<1>{

    private:
    
    public:
    /**
     * @brief pure virtual function to return the value of function at x
     * 
     * @param x independent variable of function 
     * @return double a real number
     */
    virtual double val( const double & x ) = 0;
    /**
     * @brief impletement a numerical derivative for any funcition, user can choose override to make it more precise
     * 
     * @param x independent variable of function
     * @param n order of derivative
     * @return double a real number
     */
    double diffVal( const double & x, const int & n = 1 ){
        double eps = 10e-4;
        
        if( n == 1 ) return (val(x+eps) - val(x-eps))/(2*eps);
        else{
            std::vector<double> temp(n+1);
            for( int i = 0; i <= n; i++ ) temp[i] = val(x+(2*i-n)*eps);
            for(int k = n; k > 0; k-- ){
                for( int i = 0; i < n; i++ ) temp[i] = (temp[i+1]-temp[i])/(2*eps);
            }
            return temp[0];
        }    
    }
    /**
     * @brief override the operator () to use this class in a expression just like function in the math
     * 
     * @param x independent variable of function
     * @return double a real number
     */
    double operator()( const double & x ){
        return val(x);
    }
};

/**
 * @brief polynomial inherited from Function<1>
 * 
 */
class polynomial:public Function<1>{
    private:
    /**
     * @brief coefficient of polynomial, and coeff[i] is coefficient of x^i
     * 
     */
    std::vector<double> coeff;
    /**
     * @brief order of polynomial
     * 
     */
    int Order;

    public:

    /**
     * @brief Construct a polynomial as a constant function
     * 
     * @param a polynomial = a
     */
    polynomial( const double & a = 0 ): Order{0}, coeff(1) {
        coeff[0] = a;
    }
    /**
     * @brief Construct a polynomial with coefficient
     * 
     * @param a a vector stored coefficient of polynomial
     */
    polynomial( const std::vector<double> a ): Order{int(a.size())-1}, coeff{a} {}
    /**
     * @brief Construct a n order polynomial
     * 
     * @param order number of order
     */
    explicit polynomial( const int & order ): Order{order}, coeff(order+1) {}

    /**
     * @brief compute the value of polynomial at x
     * 
     * @param x independent variable of polynomial
     * @return double 
     */
    virtual double val( const double & x ){
        double v{coeff[Order]};
        for( int i = Order-1; i >= 0; i-- ){
            v *= x;
            v += coeff[i];
        }
        return v;        
    }
    /**
     * @brief compute the derivative of polynomial at x
     * 
     * @param x independent variable of polynomial
     * @param n order of derivative
     * @return double derivative of polynomial at x
     */
    double diffVal( const double & x, const int & n = 1 ){
        if( n > Order ) return 0;
        else{
            int temp = 1;
            for( int i = Order-n+1; i <= Order; i++ ) temp*=i;
            double df{temp*coeff[Order]};
            for(int i = Order-1 ; i >= n; i-- ){
                df *= x;
                temp /= (i+1);
                temp *= (i-n+1);
                df += (temp*coeff[i]);
            }
            return df;
        }
    }

    /**
     * @brief override the operator + to compute addition of two polynomials
     * 
     * @param a one of polynomial to operator
     * @param b the other polynomial to operator
     * @return polynomial result of addition of two polynomials
     */
    friend polynomial operator+( const polynomial & a, const polynomial & b );
    /**
     * @brief override the operator - to compute subtraction of two polynomials
     * 
     * @param a minuend
     * @param b subtrahend
     * @return polynomial result of subtraction of two polynomials
     */
    friend polynomial operator-( const polynomial & a, const polynomial & b );
    /**
     * @brief overrider the operator - to get a polynomial whose coefficients are opposite of this
     * 
     * @return polynomial whose coefficients are opposite of this
     */
    polynomial operator-(){
        polynomial temp = *this;
        for (int i = 0; i <= Order; i++ ) temp.coeff[i] = -temp.coeff[i];
        return temp;
    }
    /**
     * @brief overrider the operator - to get a polynomial whose coefficients are opposite of this
     * 
     * @return polynomial whose coefficients are opposite of this
     */
    polynomial operator-() const{
        polynomial temp = *this;
        for (int i = 0; i <= Order; i++ ) temp.coeff[i] = -temp.coeff[i];
        return temp;
    }
    
    /**
     * @brief override the operator * to compute multiplication of two polynomials
     * 
     * @param a one of polynomial to operator
     * @param b the other polynomial to operator
     * @return polynomial result of multiplication of two polynomials
     */
    friend polynomial operator*( const polynomial & a, const polynomial & b );
    /**
     * @brief override the operator / to compute that a polynomial divide a real number 
     * 
     * @param a dividend, a polynomial
     * @param b divisor, a real number
     * @return polynomial result of that a polynomial divide a real number
     */
    friend polynomial operator/( const polynomial & a, const double & b );


};

polynomial operator+( const polynomial & a, const polynomial & b ){
    polynomial temp;
    if( b.Order > a.Order ){
        temp = b;
        for(int i = 0; i <= a.Order; i++ ) temp.coeff[i] += a.coeff[i];
    }else{
        temp = a;
        for(int i = 0; i <= b.Order; i++ ) temp.coeff[i] += b.coeff[i];
    }
    return temp;
}

polynomial operator-( const polynomial & a, const polynomial & b ){
    return a + (-b);
}

polynomial operator*( const polynomial & a, const polynomial & b ){
    polynomial temp{a.Order+b.Order};
    for( int j = 0; j <= a.Order; j++ ){
        for( int i = 0; i <= b.Order; i++ ) temp.coeff[i+j] += b.coeff[i]*a.coeff[j];    
    }
    return temp;
}

polynomial operator/( const polynomial & a, const double & b ){
    return a*(1.0/b);
}
