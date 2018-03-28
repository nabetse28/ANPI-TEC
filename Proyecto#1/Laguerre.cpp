#ifndef LAGUERRE_CPP_
#define LAGUERRE_CPP_

#include <iostream>
#include <complex>
#include <vector>
#include <boost/math/tools/polynomial.hpp>
#include "DeflationPolynomial.cpp"

using namespace boost::math::tools;
template <class T>
class Laguerre
{


    T Tolerance = std::sqrt(std::numeric_limits<T>::epsilon());

public:

    /**
     * Main function responsible for calculating all roots using the Leguerre method
     * @param Poly Polynomial to which the roots are calculated
     * @return Returns the roots of the polynomial
     */
    std::complex<T>* getRoots(polynomial<T> *Poly, std::complex<T> x)
    {
        polynomial<T> polyn = *Poly;
        polynomial<T> *ptrpol = &polyn;;
        int stop = polyn.size();
        DeflationPolynomial<T> df;
        std::complex<T> xm = x;
        polynomial<T> res = {0};
        polynomial<T> *ptr = &res;
        std::complex<T> *roots = new std::complex<T>[stop];
        polynomial<T> DesDeflate = polyn;
        polynomial<T> *ptrDeflate = &DesDeflate;
        std::complex<T> root;
        // the deflation of the equation is done
        for(int i = 0; i<stop-1;i++){
            root = getRootA(ptrDeflate, xm);
            if (abs(root.imag()) > Tolerance)
            {
                roots[i] = root;
                std::complex<T> otherRoot(root.real(), -root.imag());
                roots[i + 1] = otherRoot;
                i++;
                DesDeflate = df.deflate2(ptrDeflate, root, ptr);
            }
            else
            {
                T rootreal = root.real();
                DesDeflate = df.deflate(ptrDeflate, rootreal, ptr);
                roots[i] = rootreal;
            }
        }
        return roots;
    }

    /**
     * the x obtained in the equation is evaluated and the new x is obtained
     * @param Poly Polynomial to which the x is evaluated
     * @param x to be evaluated
     * @return Returns a new x
     */
    std::complex<T> getRootA(polynomial<T> *Poly, std::complex<T> x)
    {
        polynomial<T> polyn = *Poly;
        std::complex<T> fx = EvaluateComplex(Poly, x);
        polynomial<T> TempPolyd1 = getDer(Poly);
        polynomial<T> *prtd1 = &TempPolyd1;
        std::complex<T> fd1x = EvaluateComplex(prtd1, x);
        polynomial<T> TempPolyd2 = getDer(prtd1);
        polynomial<T> *prtd2 = &TempPolyd2;
        std::complex<T> fd2x = EvaluateComplex(prtd2, x);
        int n = polyn.size();
        std::complex<T> nc = n;
        std::complex<T> xOld = x;
        std::complex<T> x1 = x;
        std::complex<T> fOld = fx;
        std::complex<T> fi = fx;

        T error(0);
        std::complex<T> variable = 0;
        for(int i=0; i<80;i++){
            fi = fOld;
            xOld = x1;
            fd1x = EvaluateComplex(prtd1, xOld);
            fd2x = EvaluateComplex(prtd2, xOld);
            variable = (getVariable(nc, fOld, fd1x, fd2x));
            x1 = (xOld) - ((nc * fOld) / ((fd1x) + operatorS(xOld, prtd1, sqrt(variable))));
            fOld = EvaluateComplex(Poly, x1);

            if (std::abs(x1) > std::numeric_limits<T>::epsilon()) {
                error = std::abs((abs(x1)-abs(fOld))/abs(x1))*T(100);
            }
            else {
                error = T(0);}

            if (error < Tolerance) {
                return x1;}

        }
        return x1;
    }



    std::complex<T> getVariable(std::complex<T> n, std::complex<T> f, std::complex<T> d1p, std::complex<T> d2p)
    {
        std::complex<T> one = 1;
        std::complex<T> result = (n - one) * (((n - one) * (d1p * d1p)) - (n * f * d2p));
        return result;
    }

    std::complex<T> operatorS(std::complex<T> xk, polynomial<T> *p, std::complex<T> r)
    {
        std::complex<T> sig = EvaluateComplex(p, xk);
        if (sig.real() >= 0)
        {

            return r;
        }
        else
        {

            return -r;
        }
    }


    polynomial<T> getDer(polynomial<T> *Polyn)
    {
        polynomial<T> TempPoly = *Polyn;
        std::vector<T> VecTemp = TempPoly.data();
        int n = VecTemp.size();
        boost::array<T, 10> d3a = {};
        int it = 1;

        while (it < n)
        {
            d3a[it - 1] = VecTemp[it] * it;
            it = it + 1;
        }
        polynomial<T> derivate(d3a.begin(), d3a.end());
        return derivate;
    }


    /**
     * Function responsible for evaluating the complex number in the polynomial
     * @param Poly Polynomial to which the x is evaluated
     * @param complex number to be evaluated
     * @return Returns the result of the evaluation
     */
    std::complex<T> EvaluateComplex(polynomial<T> *Poly, std::complex<T> x)
    {
        polynomial<T> polyn = *Poly;
        std::vector<T> vecttem = polyn.data();
        int n = polyn.size();
        std::complex<T> value = 0;
        for(int i = 0; i<n;i++){
            value = value + vecttem[i] * (ElevateComplex(x, i));
        }

        return value;
    }

    /**
    * Function responsible for raising the grade of the entire complex number
    * @param x Complex number
    * @param n Grade
    * @return Returns the complex number elevated
    */
    std::complex<T> ElevateComplex(std::complex<T> x, int n)
    {
        std::complex<T> num = 1;
        for(int i = 0; i<n;i++){
            num = x * num;
        }
        return num;
    }

    /**
  * Function responsible for calling the other functions that evaluate the polynomial,
  * with the difference of indicating whether the roots are polished or not
  * @param Poly Polynomial to which the roots are calculated
  * @param Polish Boolean value to determine whether or not to polish the roots
  * @return Returns the polished roots or not, as the case may be
  */
    std::complex<T> *evalRoot(polynomial<T> Poly, std::complex<T> x0,bool Polish, int n){
        polynomial<T> *polyn = &Poly;
        std::complex<T> *data = new std::complex<T>[polyn->size()-1];
        data = getRoots(polyn, x0);
        if (Polish) {
            for (int i = 0; i < n; i++) {
                data = getRoots(polyn, data[1]);
            }
        }
        return data;
    }

};

#endif