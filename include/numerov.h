#ifndef NUMEROV_SOLVER_H
#define NUMEROV_SOLVER_H

#include <stdexcept>
#include <iostream>
#include <type_traits>

#if __cplusplus >= 202002L
#include <concepts>
template<class T>
concept arithmetic = std::integral<T> || std::floating_point<T>;
#endif
#include <iterator>
/***************************************************************************//**
* A class for solving linear second order differential equations using the 
* Numerov method.
*******************************************************************************/
class Numerov{
	public:
	/***********************************************************************//**
	* Use Numerov's algorithm to solve the differential equation.
	* \t__Y__ - Random access iterator to result container (e.g. std::vector)\n
	* \t__G__ - Random access iterator to g container (e.g. std::vector)\n
	* \t__S__ - Random access iterator to g container (e.g. std::vector)\n
	* Input:\n
	* \t__y__ - iterator to the first element in the result container to solve for\n
	* \t__y_end__ - iterator one step part the last element in the result container to solve for\n
	* \t__g__ - iterator to the element of the g container matching y\n
	* \t__s__ - iterator to the element of the s container matching y\n
	Returns:\n
	\t Iterator to one past the inversion point.\n
	NOTE : The elements of g must be multiplied with the step size squared (h^2).
	***************************************************************************/
#if __cplusplus >= 202002L
	template<std::random_access_iterator Y, 
        std::random_access_iterator G, 
        std::random_access_iterator S,
        arithmetic Timestep>
#else
	template<class Y, class G, class S, class Timestep>
#endif
	Y solve(Y y, Y y_end, G g, S s, Timestep h = 1)
	 const noexcept
	{
        return solve(y, y-1, y-2, y_end, g, g-1, g-2, s, s-1, s-2, h);
	}

	/***********************************************************************//**
	* Use Numerov's algorithm to solve the differential equation.
	* \t__Y__ - Forward iterator to result container\n
	* \t__G__ - Forward iterator to g container\n
	* \t__S__ - Forward iterator to g container\n
	* Input:\n
	* \t__y__ - iterator to the first element in the result container to solve for\n
	* \t__y_m1__ - iterator to the element one step before y\n
	* \t__y_m1__ - iterator to the element two steps before y\n
	* \t__y_end__ - iterator one step part the last element in the result container to solve for\n
	* \t__g__ - iterator to the element of the g container matching y\n
	* \t__g_m1__ - iterator to the element one step before g\n
	* \t__g_m1__ - iterator to the element two steps before g\n
	* \t__s__ - iterator to the element of the s container matching y\n
	* \t__s_m1__ - iterator to the element one step before s\n
	* \t__s_m1__ - iterator to the element two steps before s\n
	Returns:\n
	\t Iterator to one past the inversion point.\n
	NOTE : The elements of g must be multiplied with the step size squared (h^2).
	***************************************************************************/
#if __cplusplus >= 202002L
	template<std::forward_iterator Y, 
        std::forward_iterator G, 
        std::forward_iterator S,
        arithmetic Timestep>
#else
	template<class Y, class G, class S, class Timestep>
#endif
	Y solve(Y y, Y y_m1, Y y_m2, Y y_end, 
            G g, G g_m1, G g_m2, 
            S s, S s_m1, S s_m2, Timestep h = 1)
	 const noexcept
	{
		auto f = [=](typename G::value_type g){return 1 + g*h*h/12;};
		// Solve equation
		while(y != y_end){
			*y =  *y_m1++*(12 - 10*f(*g_m1++)) - *y_m2++*f(*g_m2++) +
				  1./12*(*s++ + 10*(*s_m1++) + *s_m2++);
			*y++ /= f(*g++);
		}
		return y;
	}

	/***************************************************************************
	* Derivative expressions found using the same approach as                  *
	* V.I. Tselyaev,                                                           *
	* A generalized Numerov method for linear second-order differential        *
	* equations involving a first derivative term,                             *
	* Journal of Computational and Applied Mathematics,                        *
	* Volume 170, Issue 1,                                                     *
	* 2004,                                                                    *
	* https://doi.org/10.1016/j.cam.2003.12.042.                               *
	* (https://www.sciencedirect.com/science/article/pii/S0377042704000159)    *
	***************************************************************************/
#if __cplusplus >= 202002L
	template<std::forward_iterator Y, 
        std::forward_iterator G, 
        std::forward_iterator S,
        arithmetic Timestep>
#else
	template<class Y, class G, class S, class Timestep>
#endif
	Y derivative(Y f_start, Y f_end, 
                 G g_start, S s_start, Timestep h = 1)
	 const noexcept
	{
		return f_start;
	}

};

#endif //NUMEROV_SOLVER_H
