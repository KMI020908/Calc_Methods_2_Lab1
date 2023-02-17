// Файл, содержащий функции для тестов
#ifndef TEST_FUNC_H
#define TEST_FUNC_H

#include<cmath>

// Вектор функция для решения ОДУ
template<typename Type>
std::vector<Type> sys1(Type t, const std::vector<Type> &U){
    Type k = 1.0;
    Type m = 1.0;
    return std::vector{U[1], -k / m * U[0]};
}
template<typename Type>
Type realSolution1(Type t){
    Type k = 1.0;
    Type m = 1.0;
    return cos(t);
}


template<typename Type>
std::vector<Type> sysVar1(Type t, const std::vector<Type> &U){
    return std::vector<Type>{U[1], 0.6 * U[1] - 0.6 * std::pow(U[0], 2.0) * U[1] - U[0]};
}


template<typename Type>
std::vector<Type> sysVar9(Type t, const std::vector<Type> &U){
    Type L = 1.0e7;
    Type C0 = 3.0e-7;
    Type R = 28;
    Type V0 = 0.08;
    Type f = 0.15;
    return std::vector<Type>{U[1],  V0 * sin(2.0 * M_PI * f * t) - R / L * U[1] - U[0] / (L * C0 * sqrt(1.0 + U[0] / 0.6))};
}


#endif