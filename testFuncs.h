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
std::vector<Type> funcVar9(Type t, std::vector<Type> &U){
    return std::vector<Type>{};
}


template<typename Type>
Type g(Type t, std::vector<Type> &x){
    return exp(x[0]) * sin(x[1] * t + x[0]); 
}


#endif