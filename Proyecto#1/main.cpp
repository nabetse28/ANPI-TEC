/**
 * Technological Institute of Costa Rica
 * Computer Engineering Area
 * Teacher: Dr. Pablo Alvarado Moya
 * CE-3102 Numerical Analysis for Engineering
 * @Author: Ronny Quesada Arias and Esteban Herrera Vargas
 * @Date: 14.02.2018
 * @Description: Project 1
 */

#include <boost/math/tools/polynomial.hpp>
#include "UserInterface.cpp"

namespace bmt = boost::math::tools;
using namespace bmt;
using namespace std;


int main()
{


    /*typedef std::complex<double> cplx;
    bmt::polynomial<cplx> p{{cplx(1,0),cplx(0,1),cplx(0,5)}};
    std::cout << p << std::endl;*/

    UserInterface gui;
    gui.RunTest();
    

    return 0;
}