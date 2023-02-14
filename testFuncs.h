// Файл, содержащий функции для тестов
#ifndef TEST_FUNC_H
#define TEST_FUNC_H

#include<cmath>

// Вектор функция для решения ОДУ
template<typename Type>
std::vector<Type> sys1(Type t, std::vector<Type> &U){
    return std::vector{U[1], -U[0]};
}

template<typename Type>
std::vector<Type> sysVar1(Type t, std::vector<Type> &U){
    return std::vector<Type>{U[1], 0.6 * U[1] - 0.6 * std::pow(U[0], 2.0) * U[1] - U[0]};
}


template<typename Type>
std::vector<Type> sysVar9(Type t, std::vector<Type> &U){
    return std::vector<Type>{};
}


#endif