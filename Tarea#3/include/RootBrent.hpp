/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */


#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

using namespace std;


#ifndef ANPI_ROOT_BRENT_HPP
#define ANPI_ROOT_BRENT_HPP

namespace anpi {

    /**
     * Find the roots of the function funct looking for it in the
     * interval [xl,xu], using the Brent's method.
     *
     * @param funct a std::function of the form "T funct(T x)"
     * @param xl lower interval limit
     * @param xu upper interval limit
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>

    T rootBrent(const std::function<T(T)>& funct,T xl,T xu,const T eps = numeric_limits<T>::epsilon()) {


        const  T maxi = numeric_limits<T>::digits*2; //Valor maximo de iteraciones que se esperan

        T a = xl, b = xu, c = xu, d, e;
        T fa = funct(a), fb = funct(b);
        T fc,p,q,r,s,error_r,xm;

        if(xl >= xu){
            throw anpi::Exception("El intervalo esta invertido");
        }


        if((funct(xl)*funct(xu)>T(0))){
            throw anpi::Exception("Unenclosed root");

        }


        fc = fb;

        for(int iter=0;iter<maxi;iter++){
            //Se ajustan los intervalos que se estan utilizando
            if((fb > 0.0 && fc >0.0) || (fb < 0.0 && fc < 0.0)){
                c=a;
                fc=fa;
                e=d=b-a;

            }

            if(abs(fc) < abs(fb)){
                a=b;
                b=c;
                c=a;
                fa=fb;
                fb=fc;
                fc=fa;

            }

            //Se asegura que haya convergencia
            error_r = 2.0*eps*abs(b)+0.5*eps;
            xm=0.5*(c-b);

            if(abs(xm) <= error_r || fb == 0.0){
                cout<< b<< endl;
                return b;
            }

            if( abs(e) >= error_r && abs(fa) > abs(fb)){
                // Prueba el metodo de interpolacion cuadratica
                s=fb/fa;
                if(a == c){
                    p = 2.0*xm*s;
                    q=1.0-s;
                }else{
                    q=fa/fc;
                    r=fb/fc;
                    p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                    q=(q-1.0)*(r-1.0)*(s-1.0);


                }

                //Revisa los intervalos que se estan utilizando

                if(p > 0.0){
                    q = -q;
                }

                p = abs(p);

                T min1 = 3.0*xm*q-abs(error_r*q);
                T min2=abs(e*q);

                if (2.0*p < (min1 < min2 ? min1 : min2)){
                    //Se continua con el metodo de interpolacion cuadratica
                    e=d;
                    d=p/q;

                }else{
                    //Utiliza el metodo de biseccion si no funciona l de interpolacion cuadratica
                    d=xm;
                    e=d;
                }


            }else{
                d= xm;
                e=d;
            }

            //Movimiento de los limites para acercarse a la raiz verdadera
            a=b;
            fa=fb;
            //Se prueba un nuevo valor inicial para encontrar la raiz
            if(abs(d) > error_r){
                b+=d;
            }else{
                b +=error_r*((0<xm) - (xm<0));
            }
            fb=funct(b);

        }

        throw anpi::Exception("Se excedio el numero de iteraciones");
        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }
}

#endif



