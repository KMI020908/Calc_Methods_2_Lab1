#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"
#include<iomanip>
#include"testFuncs.h"
#include<algorithm>

template<typename Type>
void checkTestEuler(std::vector<Type>(*f)(Type t, std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
const std::string &FW_E_FILE_PATH, const std::string &BW_E_FILE_PATH, const std::string &SYM_E_FILE_PATH, 
Type h = 1e-4, Type eps = 1e-6, std::size_t iterParam = 1){
    std::size_t n = U0.size(); // Размерность задачи
    std::vector<double> tGrid; // Временная сетка
    Type tau = getUniformGrid(t0, T, numOfTimeInterv, tGrid); // Шаг сетки
    std::vector<Type> dataVec = {t0, T, tau};
    std::vector<std::vector<double>> solution; // Матрица решений

    forwardEulerMethod(f, t0, T, U0, numOfTimeInterv, solution); // Явный метод Эйлера
    writeScalarFile(n, FW_E_FILE_PATH);
    writeScalarFile(numOfTimeInterv, FW_E_FILE_PATH, true);
    writeVectorFile(dataVec, FW_E_FILE_PATH, true);
    writeVectorFile(tGrid, FW_E_FILE_PATH, true);
    for (std::size_t i = 0; i < n; i++){
        writeVectorFile(solution[i], FW_E_FILE_PATH, true);
    }

    backwardEulerMethod(f, t0, T, U0, numOfTimeInterv, solution, h, eps, iterParam); // Неявный метод Эйлера
    writeScalarFile(n, BW_E_FILE_PATH);
    writeScalarFile(numOfTimeInterv, BW_E_FILE_PATH, true);
    writeVectorFile(dataVec, BW_E_FILE_PATH, true);
    writeVectorFile(tGrid, BW_E_FILE_PATH, true);
    for (std::size_t i = 0; i < n; i++){
        writeVectorFile(solution[i], BW_E_FILE_PATH, true);
    }

    symmetricScheme(f, t0, T, U0, numOfTimeInterv, solution, h, eps, iterParam); // Симметричная схема
    writeScalarFile(n, SYM_E_FILE_PATH);
    writeScalarFile(numOfTimeInterv, SYM_E_FILE_PATH, true);
    writeVectorFile(dataVec, SYM_E_FILE_PATH, true);
    writeVectorFile(tGrid, SYM_E_FILE_PATH, true);
    for (std::size_t i = 0; i < n; i++){
        writeVectorFile(solution[i], SYM_E_FILE_PATH, true);
    }
}

template<typename Type>
void checkSpeedEst(std::vector<Type>(*f)(Type t, std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
const std::string &SPEED_FW_E_FILE_PATH, const std::string &SPEED_BW_E_FILE_PATH, const std::string &SPEED_SYM_E_FILE_PATH, 
Type h = 1e-4, Type eps = 1e-6, std::size_t iterParam = 1){
    std::vector<Type> speedEstimate;
    getSpeedEstimateDiffSystem(f, t0, T, U0, numOfTimeInterv, FW_EULER, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_FW_E_FILE_PATH);
    getSpeedEstimateDiffSystem(f, t0, T, U0, numOfTimeInterv, BW_EULER, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_BW_E_FILE_PATH);
    getSpeedEstimateDiffSystem(f, t0, T, U0, numOfTimeInterv, SYM_SCHEME, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_SYM_E_FILE_PATH);
}

template<typename Type>
void temp_main(){
    Type t0 = 0.0;
    Type T = 2.0 * 3.14;
    std::vector<Type> U0 = {1.0, 0.0};
    std::size_t numOfTimeIntervals = 6;
    Type h = 1e-4;
    Type eps = 1e-6;
    std::size_t iterParam = 1;
    checkTestEuler(sys1, t0, T, U0, numOfTimeIntervals, FW_E_FILE_PATH_1, BW_E_FILE_PATH_1, SYM_E_FILE_PATH_1, h, eps, iterParam);
    checkSpeedEst(sys1, t0, T, U0, numOfTimeIntervals, SPEED_FW_E_FILE_PATH_1, SPEED_BW_E_FILE_PATH_1, SPEED_SYM_E_FILE_PATH_1, h, eps, iterParam);

    t0 = 0.0;
    T = 200.0;
    U0[0] = 0.1;
    U0[1] = 0.1;
    numOfTimeIntervals = 2000;
    h = 1e-4;
    eps = 1e-6;
    iterParam = 1;
    checkTestEuler(sysVar1, t0, T, U0, numOfTimeIntervals, FW_E_FILE_PATH_2, BW_E_FILE_PATH_2, SYM_E_FILE_PATH_2, h, eps, iterParam);
    checkSpeedEst(sys1, t0, T, U0, numOfTimeIntervals, SPEED_FW_E_FILE_PATH_2, SPEED_BW_E_FILE_PATH_2, SPEED_SYM_E_FILE_PATH_2, h, eps, iterParam);

    
}

int main(){
    temp_main<double>();

    /*
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
    */


    return 0;
}