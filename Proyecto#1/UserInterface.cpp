
#include <boost/math/tools/polynomial.hpp>
#include "RootsFinder.cpp"

using namespace std;

/**
 * Class responsible for instantiating all methods with the given examples and the case of the teacher
 */


class UserInterface {

private:
    rootsFinder finder;
public:
    UserInterface() {}

    /**
     *Function responsible for instantiating a polynomial with at least three real roots and
     * another with real and complex roots, to exercise the Laguerre method and the Muller method
     */
    void RunTest(){
        std::complex<double> evaluate = 0;
        std::complex<float> evaluate1 = 0;
        std::complex<double> x0 = 0;
        std::complex<float> x0F = 0;
        std::complex<double> x1 = 1;
        std::complex<float> x1F = 2;
        std::complex<double> x2 = 3;
        std::complex<float> x2F = 3;
        boost::math::tools::polynomial<double> poly1{{-18, 9, 7, 1, 1}};//(x^4 + x^3 + 7x^2 + 9x -18 = 0)
        polynomial<double> poly2{{4, 0, -5, 0,1}};//(x^4 - 5x^2 + 4 = 0)
        polynomial<double> poly3{{0, 169, 10, 1}};//(x^3 + 10x^2 + 169x = 0)

        polynomial<double> *polyptrd1 = &poly1;
        polynomial<double> *polyptrd2 = &poly2;
        polynomial<double> *polyptrd3 = &poly3;

        polynomial<float> poly4{{-18, 9, 7, 1, 1}};//(x^4 + x^3 + 7x^2 + 9x -18 = 0)
        polynomial<float> poly5{{4, 0, -5, 0,1}};//(x^4 - 5x^2 + 4 = 0)
        polynomial<float> poly6{{0, 169, 10, 1}};//(x^3 + 10x^2 + 169x = 0)

        polynomial<float> *polyptrd4 = &poly4;
        polynomial<float> *polyptrd5 = &poly5;
        polynomial<float> *polyptrd6 = &poly6;
        rootsFinder finder;
        cout<<"Laguerre’s Method"<<"\n\n";
        std::cout << "Function 1: f(x) = x^4 + x^3 + 7x^2 + 9x -18 with Laguerre’s method, type double: " << endl;
        cout<<"\n";
        cout<<"Without polished"<<"\n";
        finder.getLaguerre<double>(false, poly1, evaluate, 1);
        cout<<"The roots are polished"<<"\n";
        finder.getLaguerre<double>(true, poly1, evaluate, 6);
        cout<<"\n";
        std::cout << "Function 1: f(x) = x^4 + x^3 + 7x^2 + 9x -18 with Laguerre’s method, type float: " << endl;
        cout<<"\n";
        cout<<"Without polish"<<"\n";
        finder.getLaguerre<float>(false, poly4, evaluate1, 1);
        cout<<"The roots are polished"<<"\n";
        finder.getLaguerre<float>(true, poly4, evaluate1, 6);
        cout<<"\n";
        std::cout << "Function 2: f(x) = x^4 - 5x^2 + 4  with Laguerre’s method, type double: " << endl;
        cout<<"\n";
        cout<<"Without polished"<<"\n";
        finder.getLaguerre<double>(false, poly2, evaluate, 1);
        cout<<"The roots are polished"<<"\n";
        finder.getLaguerre<double>(true, poly2, evaluate, 6);
        cout<<"\n";
        std::cout << "Function 2: f(x) = x^4 - 5x^2 + 4  with Laguerre’s method, type float:" << endl;
        cout<<"\n";
        cout<<"Without polished"<<"\n";
        finder.getLaguerre<float>(false, poly5, evaluate1, 1);
        cout<<"The roots are polished"<<"\n";
        finder.getLaguerre<float>(true, poly5, evaluate1, 6);
        cout<<"\n";
        std::cout << "Function 3: f(x) = x^3 + 10x^2 + 169x  with Laguerre’s method, type double:" << endl;
        cout<<"\n";
        cout<<"Without polished"<<"\n";
        finder.getLaguerre<double>(false, poly3, evaluate, 1);
        cout<<"The roots are polished"<<"\n";
        finder.getLaguerre<double>(true, poly3, evaluate, 6);
        cout<<"\n";
        std::cout << "Function 3: f(x) = x^3 + 10x^2 + 169x  with Laguerre’s method, type float: " << endl;
        cout<<"\n";
        cout<<"Without polished"<<"\n";
        finder.getLaguerre<float>(false, poly6, evaluate1, 1);
        cout<<"The roots are polished"<<"\n";
        finder.getLaguerre<float>(true, poly6, evaluate1, 6);
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        cout<<"\n\n";
        cout<<"Muller’s Method"<<"\n\n";
        std::cout << "Function 1: f(x) = x^4 + x^3 + 7x^2 + 9x -18 with Muller’s method, type double: " << endl;
        cout<<"\n";
        cout<<"Without polished"<<"\n";
        finder.getMuller<double>(false, poly1, x0,x1,x2,1);
        cout<<"The roots are polished "<<"\n";
        finder.getMuller<double>(true, poly1, x0,x1,x2, 6);
        cout<<"\n";
        std::cout << "Function 1: f(x) = x^4 + x^3 + 7x^2 + 9x -18 with Muller’s method, type float:" << endl;
        cout<<"\n";
        cout<<"Without polished"<<"\n";
        finder.getMuller<float>(false, poly4, x0F,x1F,x2F, 1);
        cout<<"The roots are polished"<<"\n";
        finder.getMuller<float>(true, poly4, x0F,x1F,x2F, 6);
        cout<<"\n";
        std::cout << "Function 2: f(x) = x^4 - 5x^2 + 4  with Muller’s method, type double: " << endl;
        cout<<"\n";
        cout<<"Without polished"<<"\n";
        finder.getMuller<double>(false, poly2, x0,x1,x2, 1);
        cout<<"The roots are polished"<<"\n";
        finder.getMuller<double>(true, poly2, x0,x1,x2, 6);
        cout<<"\n";
        std::cout << "Function 2: f(x) = x^4 - 5x^2 + 4  with Muller’s method, type float: " << endl;
        cout<<"\n";
        cout<<"Without polished"<<"\n";
        finder.getMuller<float>(false, poly5, x0F,x1F,x2F, 1);
        cout<<"The roots are polished"<<"\n";
        finder.getMuller<float>(true, poly5, x0F,x1F,x2F, 6);
        cout<<"\n";
        std::cout << "Function 3: f(x) = x^3 + 10x^2 + 169x  with Muller’s method, type double: " << endl;
        cout<<"\n";
        cout<<"Without polished"<<"\n";
        finder.getMuller<double>(false, poly3, x0,x1,x2, 1);
        cout<<"The roots are polished"<<"\n";
        finder.getMuller<double>(true, poly3, x0,x1,x2, 6);
        cout<<"\n";
        std::cout << "Function 3: f(x) = x^3 + 10x^2 + 169x  with Muller’s method, type float: " << endl;
        cout<<"\n";
        cout<<"Without polished"<<"\n";
        finder.getMuller<float>(false, poly6, x0F,x1F,x2F, 1);
        cout<<"The roots are polished"<<"\n";
        finder.getMuller<float>(true, poly6, x0F,x1F,x2F, 6);
        cout<<"\n\n";
    }

    /**
     * Function to show the main menu
     */
    void MainMenu(){
        int option(0);
        while(true){
            cout<<"\n";
            cout<<"Project 1: Roots of Polynomials"<<"\n";
            cout<<"Menu:\n"<<"1.Run test with specific polynomials\n"<<"2.Do your own test\n"<<"3.exit\n";
            cout<<"Enter an option number: "<<"\n";
            cin>>option;
            if(option==1){
                RunTest();
            }
            else if(option==2){
                InsertPrecision();
            }
            else if(option==3){ break;}
        }
    }

    /**
     * Function to insert the precision and the polynomial grade
     */
    void InsertPrecision(){
        int grade;
        int precision;
        while(true){
            cout<<"Enter polynomial grade: "<<"\n";
            cin>>grade;
            cout<<"Enter precision: "<<"\n";
            cout<<"Menu:\n"<<"1.Double\n"<<"2.Float\n"<<"3.exit\n";
            cin>>precision;
            if(precision==1){
                InsertPolynomial<double>(grade+1);
            }
            else if(precision==2){
                InsertPolynomial<float>(grade+1);
            }
            else if(precision==3){
                break;
            }
        }
    }


    /**
     * Function to insert the polynomial coefficients
     * @param grade of polynomial
     */
    template <typename T>
    void InsertPolynomial(int grade){
        const int i = grade;
        T coef;
        std::vector<T> vector;
        for(int n=0;n<grade;n++){
            cout<<"Enter grade coefficient"<<n<<"\n";
            cin>>coef;
            vector.push_back(coef);
        }
        polynomial<T> poly(vector.begin(),vector.end());
        SelectMethod(poly);
    }

    /**
     * Function to select the polynomial root search method
     * @param poly Polynomial to which the selected method is applied
     */
    template <typename T>
    void SelectMethod(polynomial<T> poly)
    {
        while(true){
            int method;
            cout<<"Select the method: "<<"\n";
            cout<<"Menu:\n"<<"1.Laguerre’s Method\n"<<"2.Muller’s Method\n"<<"3.exit\n";
            cin>>method;
            if(method==1){
                Laguerre(poly);
            }
            else if(method==2){
                Muller(poly);
            }
            else if(method==3){
                break;
            }
        }
    }

    /**
     * Function responsible for requesting what is necessary to execute the Laguerre method
     * @param poly Polynomial to which the Laguerre method is applied
     */
    template <typename T>
    void Laguerre(polynomial<T> poly){
        int value;
        int poolish;
        while(true){
            cout<<"Enter value: "<<"\n"<<"\n";
            cin>>value;
            cout<<"Menu:\n"<<"1.Poolish\n"<<"2.Not Poolish\n"<<"3.exit\n";
            cin>>poolish;
            std::complex<T> x0 = value;
            if(poolish==1){
                finder.getLaguerre<T>(true,poly,x0,5);
            }
            else if(poolish==2){
                finder.getLaguerre<T>(false,poly,x0,1);
            }
            else if(poolish==3){
                break;
            }

        }
    }

    /**
    * Function responsible for requesting what is necessary to execute the Muller method
    * @param poly Polynomial to which the Muller method is applied
    */
    template <typename T>
    void Muller(polynomial<T> poly){
        int value0;
        int value1;
        int value2;
        int poolish;
        while(true){
            cout<<"Enter x0: "<<"\n";
            cin>>value0;
            cout<<"Enter x1: "<<"\n";
            cin>>value1;
            cout<<"Enter x2: "<<"\n";
            cin>>value2;
            cout<<"Menu:\n"<<"1.Poolish\n"<<"2.Not Poolish\n"<<"3.exit\n";
            cin>>poolish;
            std::complex<T> x0 = value0;
            std::complex<T> x1 = value1;
            std::complex<T> x2 = value2;
            if(poolish==1){
                finder.getMuller<T>(true,poly,x0,x1,x2,5);
            }
            else if(poolish==2){
                finder.getMuller<T>(false,poly,x0,x1,x2,1);
            }
            else if(poolish==3){
                break;
            }

        }
    }
};

