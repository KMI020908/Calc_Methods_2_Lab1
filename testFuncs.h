// Файл, содержащий функции для тестов
#ifndef TEST_FUNC_H
#define TEST_FUNC_H

#include<cmath>

// Вектор функция для решения ОДУ
template<typename Type>
std::vector<Type> func1(Type t, std::vector<Type> &U){
    return std::vector{U[1], -U[0]};
}




#endif