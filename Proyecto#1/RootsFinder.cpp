
#include <boost/math/tools/polynomial.hpp>
#include "Muller.cpp"
#include "Laguerre.cpp"
using namespace boost;
using namespace boost::math;
using namespace boost::math::tools;
using namespace std;



class rootsFinder {


public:
    rootsFinder() {}


    template<typename T>
    std::complex<T> *getMuller(bool isPoolish, polynomial<T> poly, std::complex<T> x0, std::complex<T> x1, std::complex<T> x2,int iter) {
        Muller<T> b = Muller<T>();
        std::complex<T> *arr = new std::complex<T>[poly.size()-1];
        arr = b.getRoots(&poly, x0, x1, x2);
        if (isPoolish) {
            for (int i = 0; i < iter; i++) {
                arr = b.getRoots(&poly, x0, arr[1], x2);
            }
        }
        Print(arr, poly.size()-1);
        return arr;

    }

    template<typename T>
    std::complex<T> *getLaguerre(bool isPoolish, polynomial<T> poly, std::complex<T> x0, int iter) {
        Laguerre<T> b = Laguerre<T>();
        polynomial<T> *tem = &poly;
        std::complex<T> *arr = new std::complex<T>[4];
        arr = b.getRoots(tem, x0);
        if (isPoolish) {
            for (int i = 0; i < iter; i++) {
                arr = b.getRoots(tem, arr[1]);
            }
        }
        Print(arr, poly.size()-1);
        return arr;



    }
    template<typename T>
    void Print(std::complex<T> *arr, int n)
    {
        for (int i = 0; i < n; ++i)
        {
            std::cout << "Root " << i << " : " << arr[i] << "\n";
        }
    }

};

