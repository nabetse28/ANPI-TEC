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

#ifndef ANPI_NEWTON_RAPHSON_HPP
#define ANPI_NEWTON_RAPHSON_HPP

namespace anpi {


    /**
     * Find the roots of the function funct looking by means of the
     * Newton-Raphson method
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param xi initial root guess
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */


    template<typename T>

    T rootNewtonRaphson(const std::function<T(T)>& funct,T xi,const T eps = numeric_limits<T>::epsilon()) {

        const T maxi = numeric_limits<T>::digits;

        T x=xi; // Se define el valor con el que x tiene que comenzar
        T xp=x-1; // El valor de xp el cual es parecido a x
        T h= eps; // El h a utilizar se le asigna el valor de epsilon


        for(int i=0;i<maxi;i++){

            T fxm = funct(xp - h); // Esta es para la derivada centrada
            T fxp = funct(xp + h);
            T dfxl = fxp - fxm;
            T dev = dfxl / (2 * h); // Se calcula la derivada centrada para tener el otro valor de x
            x = xp;
            T fx = funct(x); // Se evalua el dato en la funcion
            T dfx = dev; // Esta es la derivada obtenida
            xp = x - (fx / dfx); // Este es el valor de xp esperado


            if(abs(xp-x) <= eps){ // Condicion de finalizacion, si el error es menos que eps entonces ya se llego a la aproximacion
                cout<< xp<<endl; // y se retorna xp
                return xp;
            }


        }

        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }



}

#endif
