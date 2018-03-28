#ifndef MULLER_CPP_
#define MULLER_CPP_

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <boost/math/tools/polynomial.hpp>
#include "DeflationPolynomial.cpp"
using namespace boost;
using namespace boost::math;
using namespace boost::math::tools;
using namespace std;


template <class W>
class Muller {


public:
    Muller() = default;


    /**
     * Function responsible for calculating the coefficients of the quadratic equation
     * @param Poly Polynomial of the quadratic equation
     * @return Returns the coefficients of the quadratic equation
     */
    template <typename T>
    std::complex<T> *Constants(polynomial<T> Poly, std::complex<T> x0,  std::complex<T> x1,  std::complex<T> x2){
        std::complex<T> *constants = new  std::complex<T>[3];
        std::complex<T> alpha2 = x1 - x2;
        std::complex<T> alpha1 = x0 - x1;

        std::complex<T> pi2 = (EvaluateComplex(&Poly,x1)-EvaluateComplex(&Poly,x2))/alpha2;
        std::complex<T> pi1 = (EvaluateComplex(&Poly,x0)-EvaluateComplex(&Poly,x1))/alpha1;
        constants[0] = (pi1 -pi2)/(alpha1-alpha2);
        constants[1] = constants[0]*alpha1 + pi1;
        constants[2] = EvaluateComplex(&Poly,x0);
        return  constants; }

    /**
     * Function responsible for evaluating the complex number in the polynomial
     * @param Poly Polynomial to which the complex number is evaluated
     * @param x Complex number
     * @return Returns the result of the evaluation
     */
    template <typename T>
    std::complex<T> EvaluateComplex(polynomial<T> *Poly, std::complex<T> x)
    {
        polynomial<T> tempoly = *Poly;
        std::vector<T> vecttem = tempoly.data();
        int n = vecttem.size();
        std::complex<T> value = 0;
        int it = 0;

        while (it < n)
        {
            value = value + vecttem[it] * (ElevateComplex(x, it));
            it = it + 1;
        }

        return value;
    }


    /**
     * Function responsible for raising the grade of the entire complex number
     * @param x Complex number
     * @param n Grade
     * @return Returns the complex number elevated
     */
    template <typename T>
    std::complex<T> ElevateComplex(std::complex<T> x, int n)
    {
        int i = 0;
        std::complex<T> val = 1;
        while (i < n)
        {
            val = x * val;
            i = i + 1;
        }

        return val;
    }


    /**
     * Function in charge of evaluating the value of X in the function with the general formula, and recalculates what
     * the next value of X will be so on until an acceptable error is reached
     * @param Polyn Polynomial to which is evaluated
     * @return Returns a NaN if it did not find any root in the number of iterations
     */

    template <typename T>
    std::complex<T> getRootA(polynomial<T> Poly,std::complex<T> PoX0, std::complex<T> PoX1, std::complex<T> PoX2){
        T error=T();
        std::complex<T> *Roots= new std::complex<T>[3];;
        std::complex<T> *constants = new std::complex<T>[3];
        Roots[0] = PoX0;
        Roots[1] = PoX1;
        Roots[2] = PoX2;
        T max = std::numeric_limits<T>::digits;
        std::complex<T> r1= Roots[0];
        T tolerance = std::sqrt(std::numeric_limits<T>::epsilon());
        int i(0);
        while (i<max){
            constants = Constants(Poly,Roots[0],Roots[1],Roots[2]);
            std::complex<T>  rV(r1);
            std::complex<T>  discriminant = constants[1]*constants[1] - std::complex<T>(4)*constants[0]*constants[2];
            r1 = rV - ((std::complex<T>(2)*constants[2])/
                       (constants[1] + (constants[1].real()/std::abs(constants[1].real()))*sqrt(discriminant)));
            Roots[2] = Roots[1];
            Roots[1] = Roots[0];
            Roots[0] = r1;
            if (std::abs(r1) > std::numeric_limits<T>::epsilon()) {
                error = std::abs((abs(r1)-abs(rV))/abs(r1))*T(100);
            }
            else {
                error = T(0);}

            if (error < tolerance) {
                return r1;}
            i++;
        }

        return std::numeric_limits<T>::quiet_NaN();

    }

    /**
     * Function that is responsible for storing the obtained roots, then deflates the equation and
     * returns to obtain the next, so on until the equation reaches zero grade
     * @param Polynomial to which the roots are obtained
     * @return The roots
     */

    template <typename T>
    std::complex<T> *getRoots(polynomial<T> *Poly,std::complex<T> X0, std::complex<T> X1, std::complex<T> X2){
        polynomial<T> res = {0};
        polynomial<T> *ptrR = &res;
        polynomial<T> polynomial = *Poly;
        T tolerance = std::sqrt(std::numeric_limits<T>::epsilon());
        std::complex<T> *Roots = new std::complex<T>[Poly->size()-1];
        DeflationPolynomial<T> deflation;
        int i(0);
        while(i<(Poly->size()-1)){
            Roots[i]= getRootA(polynomial,X0,X1,X2);
            if(abs(Roots[i].imag())>tolerance){
                polynomial = deflation.deflate2(&polynomial,Roots[i],ptrR);
                Roots[i+1].real(Roots[i].real());
                Roots[i+1].imag(-Roots[i].imag());
                i++;
            }
            else {
                polynomial = deflation.deflate(&polynomial, Roots[i].real(), ptrR);
            }
            i++;
        }
        return Roots;
    }


    /**
     * Function in charge of polishing the roots, the function is re-evaluated but
     * with the roots obtained to obtain a more accurate result
     * @param Poly Polynomial that is polished
     * @return Returns the polished roots
     */
    template <typename T>
    std::complex<T> *getPolished(polynomial<T> *Poly,std::complex<T> X0, std::complex<T> X1, std::complex<T> X2,int n){
        std::complex<T> *datos = new std::complex<T>[Poly->size()-1];
        datos = getRoots<T>(Poly, X0, X1, X2);
        int i = 0;
        while(i<n){
            datos = getRoots<T>(Poly,X0, datos[1],X2);
            i++;
        }
        return  datos;
    }


    /**
     * Function responsible for calling the other functions that evaluate the polynomial,
     * with the difference of indicating whether the roots are polished or not
     * @param Polyn Polynomial to which the roots are calculated
     * @param Polish Boolean value to determine whether or not to polish the roots
     * @return Returns the polished roots or not, as the case may be
     */
    template <typename T>
    std::complex<T> *RootsPol(polynomial<T> *Poly,std::complex<T> X0, std::complex<T> X1, std::complex<T> X2,bool Polish,int n)
    {
        if(Polish){return getPolished(Poly,X0,X1,X2,n);}
        else{return getRoots(Poly,X0,X1,X2);}
    }
};

#endif

