#ifndef ADAMS_MOULTON_H
#define ADAMS_MOULTON_H
#include <stdexcept>
#include <iostream>
#include <type_traits>
#include <tuple>
#include <cmath>

#if __cplusplus >= 202002L
#include <concepts>
template<class T>
concept arithmetic = std::integral<T> || std::floating_point<T>;
#endif
#include <iterator>
/***************************************************************************//**
* A class for solving a system of __linear__ differential equations using the 4 
* point Adams-Moulton method.
* N.B.:
* This solver is only valid for situations where f(x, y) = G(x)y(x)
*******************************************************************************/
class LinearAdamsMoulton4
{

    public:
#if __cplusplus >= 202002L
	template<std::random_access_iterator Y, 
        std::random_access_iterator G, 
        arithmetic Timestep = double>
#else
	template<class Y, class G, class Timestep = double>
#endif
    Y solve(Y y, Y y_end, G g, Timestep h = 1)
    const noexcept{
        double l = 9./24;
        while(y != y_end){
            double M = 1 - h*l*(*g);
            *y = *(y-1) + h/24*(*(y-3)*(*(g-3)) - 5*(*(y-2)*(*(g-2))) + 19*(*(y-1)*(*(g - 1))));
            *y++ /= M;
            g++;
        }
        return y;
    }

#if __cplusplus >= 202002L
	template<std::random_access_iterator Y, 
        std::random_access_iterator G, 
        arithmetic Timestep = double>
#else
	template<class Y, class G, class Timestep = double>
#endif
    std::tuple<Y, Y> solve(std::tuple<Y, Y> y, Y y_end, std::tuple<G, G, G, G> g, Timestep h = 1)
    const noexcept{
#if __cplusplus >= 201703L
        auto [y1, y2] = y;
        auto [g11, g12, g21, g22] = g;
#else
        Y y1 = std::get<0>(y), y2 = std::get<1>(y);
        G g11 = std::get<0>(g), g12 = std::get<1>(g), 
          g21 = std::get<2>(g), g22 = std::get<3>(g);
#endif
        double l = 9*h/24;
        while(y1 != y_end){
            /**********************************************************************
             * M = [1, 0, 0, 1] -h*l[g11, g12, g21, g22]
             * Delta = det|M| = det|[1 - h*l*g11, -h*l*g12, -h*l*g21, 1-h*l*g22]| =
             *       = (1 - h*l*g11)*(1 - h*l*g22) - (h*l)^2*g12*g21 =
             *       = 1 - h*l*(g11 + g22) + (h*l)^2(g11*g22) - (h*l)^2*g12*g21 =
             *       = 1 - h*l*(g11 + g22) + (h*l)^2(g11*g22 - g12*g21)
             * M^-1 = 1/Delta * [M22, -M21, -M12, M11] =
             *      = 1/Delta * [1 - h*l*g22, h*l*g12, h*l*g21, 1 - h*l*g11]
             * y = M^-1(y[-1] + h/24(y[-3]g[-3] - 5y[-2]g[-2] + 19y[-1]g[-1]))
             *********************************************************************/
            double Delta = (1 - *g11*l)*(1 - *g22*l) - (*g12*l)*(*g21*l);
            double Mi[4] = { 1/Delta * (1 - l*(*g22)),  1/Delta*     l*(*g12), 
                             1/Delta *      l*(*g21) ,  1/Delta*(1 - l*(*g11)) };
            /*******************************************************************
             * y = Mi ( y[-1] + h/24 (g[-3]y[-3] -5g[-2]y[-2] +19g[-1]y[-1]) )
             * |y1| = |Mi11    Mi12|x( |y1[-1]| + h/24(
             * |y2|   |Mi21    Mi22| ( |y2[-1]| +     (
             *
             *                         |g11[-3]    g12[-3]|x|y1[-3]|
             *                         |g21[-3]    g22[-3]| |y2[-3]|
             * 
             *                     - 5 |g11[-2]    g12[-2]|x|y1[-2]|
             *                         |g21[-2]    g22[-2]| |y2[-2]|
             *
             *                    + 19 |g11[-1]    g12[-1]|x|y1[-1]| ) )
             *                         |g21[-1]    g22[-1]| |y2[-1]| ) )
             ******************************************************************/
            *y1 = *(y1-1) + h/24*( *(g11-3)*(*(y1-3)) + *(g12-3)*(*(y2-3)) 
                              - 5*(*(g11-2)*(*(y1-2)) + *(g12-2)*(*(y2-2))) 
                              +19*(*(g11-1)*(*(y1-1)) + *(g12-1)*(*(y2-1))) );
            *y2 = *(y2-1) + h/24*( *(g21-3)*(*(y1-3)) + *(g22-3)*(*(y2-3)) 
                              - 5*(*(g21-2)*(*(y1-2)) + *(g22-2)*(*(y2-2))) 
                              +19*(*(g21-1)*(*(y1-1)) + *(g22-1)*(*(y2-1))) );

            *y1 = Mi[0]*(*y1) + Mi[1]*(*y2);
            *y2 = Mi[2]*(*y1) + Mi[3]*(*y2);

            y1++, y2++;
            g11++, g12++, g21++, g22++;
        }
        return {y1, y2};
    }
};
#endif //ADAMS_MOULTON_H
