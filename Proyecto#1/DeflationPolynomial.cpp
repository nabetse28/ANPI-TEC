#ifndef DEFLACIONPOLINOMIAL_CPP_
#define DEFLACIONPOLINOMIAL_CPP_

#include <complex>
#include <vector>
#include <boost/math/tools/polynomial.hpp>

using namespace boost::math::tools;
template <class T>
class DeflationPolynomial
{

public:
    /**
     * Function used to deflate a polynomial of real roots
     * @param polynomial<T> *poly Polynomial pointer to perform deflation
     * @param T *root Pointer of the root of the polynomial that will be used to deflate
     * @param polynomial<T> *remainder Pointer polynomial where the remainder of the operation is stored
     * @return Returns the quotient of division in polynomial format
    **/
    polynomial<T> deflate(polynomial<T> *poly, T root, polynomial<T> *remainder)
    {
        T temroot = root; // the value of the root is assigned to a local variablez
        polynomial<T> temdiv{{-temroot, 1}}; // polynomial created to divide by the input
        polynomial<T> temquotient = (*poly) / temdiv; // polynomial created by dividing the input polynomial by the factor with root, to obtain the quotient
        polynomial<T> temremainder = (*poly) % temdiv; // polynomial created by applying the input polynomial module with the factor with root, to obtain the remainder
        *remainder = temremainder;
        return temquotient;
    }

    /**
     * Function used to deflate a polynomial of complex roots
     * @param polynomial<T> *poly Polynomial pointer to perform deflation
     * @param complex<T> *root Pointer to the root of the polynomial that will be used to deflate
     * @param polynomial<T> *remainder Pointer polynomial where the remainder of the operation is stored
     * @return Returns the quotient of division in polynomial format
    **/
    polynomial<T> deflate2(polynomial<T> *poly, std::complex<T> root, polynomial<T> *remainder)
    {
        std::complex<T> temroot1 = root;   // the value of the complex root is assigned to a local variable
        std::complex<T> conjugate(temroot1.real(), -temroot1.imag()); // the conjugate of the root is created
        std::complex<T> tempA = temroot1*conjugate;   // the value of the complex root multiplied by its conjugate
        std::complex<T> tempB = -((temroot1)+(conjugate)); // the value of the complex root plus its conjugate
        polynomial<T> tempoly = *poly; // the value of the polynomial is assigned to a local variable
        polynomial<T> temB = {tempA.real(), tempB.real(), 1}; // polynomial created by multiplying (x-root) (x + rootc), where rootc is the root conjugate
        polynomial<T> quotient = tempoly/temB; // polynomial created by dividing the input polynomial by the factor with root, to obtain the quotient
        polynomial<T> temremainder = tempoly%temB; // polynomial created by applying the input polynomial module with the root factor, to obtain remainder
        *remainder = temremainder;
        return quotient;
    }
};

#endif