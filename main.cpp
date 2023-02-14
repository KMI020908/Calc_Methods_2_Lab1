#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"
#include<iomanip>
#include"testFuncs.h"
#include<algorithm>

// Процедура проверки решений уравнения методами бисекции и НьютонаЫ 
template<typename Type>
void checkTestEquations(Type (*f)(Type x), Type firstX, Type lastX, Type x0,
const std::string &B_EQ_FILE_PATH, const std::string &N_EQ_FILE_PATH, Type accuracy = 1e-6, Type h = 1e-4, std::size_t stopIt = 10000, std::size_t n = 51){
    
    
}

// Процедура проверки систем уравнений методом Ньютона
template<typename Type>
void checkTestSystem(Type (*f1)(Type x, Type y), Type (*f2)(Type x, Type y), std::vector<Type> (*getJacobiMatrixElems)(Type x, Type y), 
Type x0, Type y0, Type L1, Type L2, std::size_t N, const std::string &N_SYS_FILE_PATH, const std::string &A_N_SYS_FILE_PATH, 
const std::string &IT_SYS_FILE_PATH, const std::string &A_IT_SYS_FILE_PATH, Type h = 1e-4, Type accuracy = 1e-6, std::size_t stopIt = 10000){
    
}

template<typename Type>
void temp_main(){
    
}

int main(){
    //temp_main<double>();
    std::vector<double> xGrid;
    getUniformGrid(0.0, 1.0, 2, xGrid);
    std::cout << xGrid << '\n';

    std::vector<std::vector<double>> solution;
    double t0 = 0.0;
    double T = 2.0 * M_PI;
    std::vector<double> U0 = {1.0, 0.0};
    std::size_t numOfTimeInt = 6;
    forwardEulerMethod(func1, t0, T, U0, numOfTimeInt, solution);
    std::vector<double> tGrid;
    getUniformGrid(t0, T, numOfTimeInt, tGrid);
    std::cout << tGrid << '\n' << '\n';
    std::cout << solution[0] << '\n' << '\n';
    std::cout << solution[1] << '\n' << '\n';


    return 0;
}