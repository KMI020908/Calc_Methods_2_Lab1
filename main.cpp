#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"
#include<iomanip>
#include"testFuncs.h"
#include<algorithm>

// Проверка методов Эйлера
template<typename Type>
void checkTestEuler(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
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

// Проверка методов Рунге-Кутты
template<typename Type>
void checkTestRungeKutta(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
const std::string &TW_RG_FILE_PATH, const std::string &FO_RG_FILE_PATH, bool autoStep = true, Type eps = 1e-6, Type lowEps = 1e-8){
    std::size_t n = U0.size(); // Размерность задачи
    std::vector<double> tGrid; // Временная сетка
    Type tau = getUniformGrid(t0, T, numOfTimeInterv, tGrid); // Шаг сетки
    std::vector<Type> dataVec = {t0, T, tau};
    std::vector<std::vector<double>> solution; // Матрица решений

    RungeKuttaMethod2(f, t0, T, U0, numOfTimeInterv, solution, autoStep, eps, lowEps); // Метод Рунге-Кутты 2-ого порядка
    writeScalarFile(n, TW_RG_FILE_PATH);
    writeScalarFile(numOfTimeInterv, TW_RG_FILE_PATH, true);
    writeVectorFile(dataVec, TW_RG_FILE_PATH, true);
    writeVectorFile(tGrid, TW_RG_FILE_PATH, true);
    for (std::size_t i = 0; i < n; i++){
        writeVectorFile(solution[i], TW_RG_FILE_PATH, true);
    }

    RungeKuttaMethod4(f, t0, T, U0, numOfTimeInterv, solution, autoStep, eps, lowEps); // Метод Рунге-Кутты 4-ого порядка 
    writeScalarFile(n, FO_RG_FILE_PATH);
    writeScalarFile(numOfTimeInterv, FO_RG_FILE_PATH, true);
    writeVectorFile(dataVec, FO_RG_FILE_PATH, true);
    writeVectorFile(tGrid, FO_RG_FILE_PATH, true);
    for (std::size_t i = 0; i < n; i++){
        writeVectorFile(solution[i], FO_RG_FILE_PATH, true);
    }
}

// Проверка методов Адамса
template<typename Type>
void checkTestAdams(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
const std::string &FO_AD_FILE_PATH, const std::string &PC_AD_FILE_PATH){
    std::size_t n = U0.size(); // Размерность задачи
    std::vector<double> tGrid; // Временная сетка
    Type tau = getUniformGrid(t0, T, numOfTimeInterv, tGrid); // Шаг сетки
    std::vector<Type> dataVec = {t0, T, tau};
    std::vector<std::vector<double>> solution; // Матрица решений

    AdamsMethod(f, t0, T, U0, numOfTimeInterv, solution); // Метод Aдамса 4-х шаговый
    writeScalarFile(n, FO_AD_FILE_PATH);
    writeScalarFile(numOfTimeInterv, FO_AD_FILE_PATH, true);
    writeVectorFile(dataVec, FO_AD_FILE_PATH, true);
    writeVectorFile(tGrid, FO_AD_FILE_PATH, true);
    for (std::size_t i = 0; i < n; i++){
        writeVectorFile(solution[i], FO_AD_FILE_PATH, true);
    }

    predicCorrect(f, t0, T, U0, numOfTimeInterv, solution); // Метод предсказание-коррекция 
    writeScalarFile(n, PC_AD_FILE_PATH);
    writeScalarFile(numOfTimeInterv, PC_AD_FILE_PATH, true);
    writeVectorFile(dataVec, PC_AD_FILE_PATH, true);
    writeVectorFile(tGrid, PC_AD_FILE_PATH, true);
    for (std::size_t i = 0; i < n; i++){
        writeVectorFile(solution[i], PC_AD_FILE_PATH, true);
    }
}

// Оценка скорости, когда неизвестно решение
template<typename Type>
void checkSpeedEst(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
const std::string &SPEED_FW_E_FILE_PATH, const std::string &SPEED_BW_E_FILE_PATH, const std::string &SPEED_SYM_E_FILE_PATH,
const std::string &SPEED_RG_2_FILE_PATH, const std::string &SPEED_RG_4_FILE_PATH, const std::string &SPEED_AD_4_FILE_PATH, const std::string &SPEED_PC_FILE_PATH,
Type h = 1e-4, Type eps = 1e-6, std::size_t iterParam = 1){
    std::vector<Type> speedEstimate;
    getSpeedEstimateDiffSystem(f, t0, T, U0, numOfTimeInterv, FW_EULER, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_FW_E_FILE_PATH);
    getSpeedEstimateDiffSystem(f, t0, T, U0, numOfTimeInterv, BW_EULER, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_BW_E_FILE_PATH);
    getSpeedEstimateDiffSystem(f, t0, T, U0, numOfTimeInterv, SYM_SCHEME, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_SYM_E_FILE_PATH);
    getSpeedEstimateDiffSystem(f, t0, T, U0, numOfTimeInterv, TWICE_RG, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_RG_2_FILE_PATH);
    getSpeedEstimateDiffSystem(f, t0, T, U0, numOfTimeInterv, FOURTH_RG, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_RG_4_FILE_PATH);
    getSpeedEstimateDiffSystem(f, t0, T, U0, numOfTimeInterv, FOURTH_AD, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_AD_4_FILE_PATH);
    getSpeedEstimateDiffSystem(f, t0, T, U0, numOfTimeInterv, PREDICT_CORRECT, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_PC_FILE_PATH);
}

// Оценка скорости, когда известно решение
template<typename Type>
void checkSpeedEst(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type(*realSolution)(Type t), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
const std::string &SPEED_FW_E_FILE_PATH, const std::string &SPEED_BW_E_FILE_PATH, const std::string &SPEED_SYM_E_FILE_PATH,
const std::string &SPEED_RG_2_FILE_PATH, const std::string &SPEED_RG_4_FILE_PATH, const std::string &SPEED_AD_4_FILE_PATH, const std::string &SPEED_PC_FILE_PATH,
Type h = 1e-4, Type eps = 1e-6, std::size_t iterParam = 1){
    std::vector<Type> speedEstimate;
    getSpeedEstimateDiffSystem(f, realSolution, t0, T, U0, numOfTimeInterv, FW_EULER, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate,  SPEED_FW_E_FILE_PATH);
    getSpeedEstimateDiffSystem(f, realSolution, t0, T, U0, numOfTimeInterv, BW_EULER, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_BW_E_FILE_PATH);
    getSpeedEstimateDiffSystem(f, realSolution, t0, T, U0, numOfTimeInterv, SYM_SCHEME, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_SYM_E_FILE_PATH);
    getSpeedEstimateDiffSystem(f, realSolution, t0, T, U0, numOfTimeInterv, TWICE_RG, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_RG_2_FILE_PATH);
    getSpeedEstimateDiffSystem(f, realSolution, t0, T, U0, numOfTimeInterv, FOURTH_RG, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_RG_4_FILE_PATH);
    getSpeedEstimateDiffSystem(f, realSolution, t0, T, U0, numOfTimeInterv, FOURTH_AD, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_AD_4_FILE_PATH);
    getSpeedEstimateDiffSystem(f, realSolution, t0, T, U0, numOfTimeInterv, PREDICT_CORRECT, speedEstimate, h, eps, iterParam);
    writeVectorFile(speedEstimate, SPEED_PC_FILE_PATH);
}

// Фазовые траектории
template<typename Type>
void checkPhaseTraces(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, std::size_t numOfTimeInterv, Type L, std::size_t N,
const std::string &PHASE_FW_E_FILE_PATH, const std::string &PHASE_BW_E_FILE_PATH, const std::string &PHASE_SYM_E_FILE_PATH,
const std::string &PHASE_RG_2_FILE_PATH, const std::string &PHASE_RG_4_FILE_PATH, const std::string &PHASE_AD_4_FILE_PATH, const std::string &PHASE_PC_FILE_PATH,
Type h = 1e-4, Type eps = 1e-6, std::size_t iterParam = 1){
    std::vector<std::vector<Type>> dataMatrix;
    getPhaseTraces(f, t0, T, numOfTimeInterv, FW_EULER, L, N, dataMatrix, h, eps, iterParam);
    writeMatrixFile(dataMatrix, PHASE_FW_E_FILE_PATH);
    getPhaseTraces(f, t0, T, numOfTimeInterv, BW_EULER, L, N, dataMatrix, h, eps, iterParam);
    writeMatrixFile(dataMatrix, PHASE_BW_E_FILE_PATH);
    getPhaseTraces(f, t0, T, numOfTimeInterv, SYM_SCHEME, L, N, dataMatrix, h, eps, iterParam);
    writeMatrixFile(dataMatrix, PHASE_SYM_E_FILE_PATH);
    getPhaseTraces(f, t0, T, numOfTimeInterv, TWICE_RG, L, N, dataMatrix, h, eps, iterParam);
    writeMatrixFile(dataMatrix, PHASE_RG_2_FILE_PATH);
    getPhaseTraces(f, t0, T, numOfTimeInterv, FOURTH_RG, L, N, dataMatrix, h, eps, iterParam);
    writeMatrixFile(dataMatrix, PHASE_RG_4_FILE_PATH);
    getPhaseTraces(f, t0, T, numOfTimeInterv, FOURTH_AD, L, N, dataMatrix, h, eps, iterParam);
    writeMatrixFile(dataMatrix, PHASE_AD_4_FILE_PATH);
    getPhaseTraces(f, t0, T, numOfTimeInterv, PREDICT_CORRECT, L, N, dataMatrix, h, eps, iterParam);
    writeMatrixFile(dataMatrix, PHASE_PC_FILE_PATH);
}

// Анализ методов Рунге-Кутты с автоматическим шагом
template<typename Type>
void RungeKuttaStepAnalys(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
const std::string &SA_RG_2_FILE_PATH, const std::string &SA_RG_4_FILE_PATH, Type eps = 1e-6, Type lowEps = 1e-8){
    std::size_t n = U0.size(); // Размерность задачи
    std::vector<double> tGrid; // Временная сетка
    Type tau = getUniformGrid(t0, T, numOfTimeInterv, tGrid); // Шаг сетки
    std::vector<Type> dataVec = {t0, T, tau};
    std::vector<std::vector<double>> dataMatrix; // Матрица решений

    RungeKuttaMethodStepAnalys2(f, t0, T, U0, numOfTimeInterv, dataMatrix, eps, lowEps); // Анализ Рунге-Кутты 2-ого порядка
    writeMatrixFile(dataMatrix, SA_RG_2_FILE_PATH);

    RungeKuttaMethodStepAnalys4(f, t0, T, U0, numOfTimeInterv, dataMatrix, eps, lowEps); // Анализ Рунге-Кутты 4-ого порядка
    writeMatrixFile(dataMatrix, SA_RG_4_FILE_PATH);
}

template<typename Type>
void temp_main(){
    Type t0 = 0.0;
    Type T = 100.0;
    //Type T = 6.0;
    std::vector<Type> U0 = {1.0, 0.0};
    std::size_t numOfTimeIntervals = 1000;
    //std::size_t numOfTimeIntervals = 6;
    Type h = 1e-4;
    Type eps = 1e-6;
    bool autoStep = true;
    Type lowEps = 1e-8;
    std::size_t iterParam = 15;
    checkTestEuler(sys1, t0, T, U0, numOfTimeIntervals, FW_E_FILE_PATH_1, BW_E_FILE_PATH_1, SYM_E_FILE_PATH_1, h, eps, iterParam);
    checkTestRungeKutta(sys1, t0, T, U0, numOfTimeIntervals, TW_RG_FILE_PATH_1, FO_RG_FILE_PATH_1, autoStep, eps, lowEps);
    checkTestAdams(sys1, t0, T, U0, numOfTimeIntervals, FO_AD_FILE_PATH_1, PC_AD_FILE_PATH_1);
    checkSpeedEst(sys1, realSolution1, t0, T, U0, numOfTimeIntervals, SPEED_FW_E_FILE_PATH_1, SPEED_BW_E_FILE_PATH_1, SPEED_SYM_E_FILE_PATH_1, 
    SPEED_RG_2_FILE_PATH_1, SPEED_RG_4_FILE_PATH_1, SPEED_AD_4_FILE_PATH_1, SPEED_PC_FILE_PATH_1, h, eps, iterParam);
    //RungeKuttaStepAnalys(sys1, t0, T, U0, numOfTimeIntervals, SA_RG_2_FILE_PATH_1, SA_RG_4_FILE_PATH_1, eps, lowEps);
    Type L = 20.0;
    std::size_t N = 4;
    //checkPhaseTraces(sys1, t0, T, numOfTimeIntervals, L, N, PHASE_FW_E_FILE_PATH_1, PHASE_BW_E_FILE_PATH_1, 
    //PHASE_SYM_E_FILE_PATH_1, PHASE_RG_2_FILE_PATH_1, PHASE_RG_4_FILE_PATH_1, PHASE_AD_4_FILE_PATH_1, PHASE_PC_FILE_PATH_1, h, eps, iterParam);

    t0 = 0.0;
    T = 200.0;
    U0[0] = 0.1;
    U0[1] = 0.1;
    numOfTimeIntervals = 200;
    h = 1e-3;
    eps = 1e-6;
    autoStep = true;
    lowEps = 1e-7;
    iterParam = 10;
    //checkTestEuler(sysVar1, t0, T, U0, numOfTimeIntervals, FW_E_FILE_PATH_2, BW_E_FILE_PATH_2, SYM_E_FILE_PATH_2, h, eps, iterParam);
    //checkTestRungeKutta(sysVar1, t0, T, U0, numOfTimeIntervals, TW_RG_FILE_PATH_2, FO_RG_FILE_PATH_2, autoStep, eps, lowEps);
    checkTestAdams(sysVar1, t0, T, U0, numOfTimeIntervals, FO_AD_FILE_PATH_2, PC_AD_FILE_PATH_2);
    //checkSpeedEst(sysVar1, t0, T, U0, numOfTimeIntervals, SPEED_FW_E_FILE_PATH_2, SPEED_BW_E_FILE_PATH_2, SPEED_SYM_E_FILE_PATH_2,
    //SPEED_RG_2_FILE_PATH_2, SPEED_RG_4_FILE_PATH_2, SPEED_AD_4_FILE_PATH_2, SPEED_PC_FILE_PATH_2, h, eps, iterParam);
    //RungeKuttaStepAnalys(sysVar1, t0, T, U0, numOfTimeIntervals, SA_RG_2_FILE_PATH_2, SA_RG_4_FILE_PATH_2, eps, lowEps);
    L = 20.0;
    N = 4;
    //checkPhaseTraces(sysVar1, t0, T, numOfTimeIntervals, L, N, PHASE_FW_E_FILE_PATH_2, PHASE_BW_E_FILE_PATH_2, 
    //PHASE_SYM_E_FILE_PATH_2, PHASE_RG_2_FILE_PATH_2, PHASE_RG_4_FILE_PATH_2, PHASE_AD_4_FILE_PATH_2, PHASE_PC_FILE_PATH_2, h, eps, iterParam);

    t0 = 0.0;
    T = 100.0;
    U0[0] = 0.4;
    U0[1] = 0.0;
    numOfTimeIntervals = 1000;
    h = 1e-4;
    eps = 1e-6;
    autoStep = false; 
    lowEps = 1e-8;
    iterParam = 1;
    //checkTestEuler(sysVar9, t0, T, U0, numOfTimeIntervals, FW_E_FILE_PATH_3, BW_E_FILE_PATH_3, SYM_E_FILE_PATH_3, h, eps, iterParam);
    checkTestRungeKutta(sysVar9, t0, T, U0, numOfTimeIntervals, TW_RG_FILE_PATH_3, FO_RG_FILE_PATH_3, autoStep, eps, lowEps);
    checkTestAdams(sysVar9, t0, T, U0, numOfTimeIntervals, FO_AD_FILE_PATH_3, PC_AD_FILE_PATH_3);
    //checkSpeedEst(sysVar9, t0, T, U0, numOfTimeIntervals, SPEED_FW_E_FILE_PATH_3, SPEED_BW_E_FILE_PATH_3, SPEED_SYM_E_FILE_PATH_3,
    //SPEED_RG_2_FILE_PATH_3, SPEED_RG_4_FILE_PATH_3, SPEED_AD_4_FILE_PATH_3, SPEED_PC_FILE_PATH_3, h, eps, iterParam);
    //RungeKuttaStepAnalys(sysVar9, t0, T, U0, numOfTimeIntervals, SA_RG_2_FILE_PATH_3, SA_RG_4_FILE_PATH_3, eps, lowEps);
    L = 20.0;
    N = 4;
    //checkPhaseTraces(sysVar9, t0, T, numOfTimeIntervals, L, N, PHASE_FW_E_FILE_PATH_3, PHASE_BW_E_FILE_PATH_3, 
    //PHASE_SYM_E_FILE_PATH_3, PHASE_RG_2_FILE_PATH_3, PHASE_RG_4_FILE_PATH_3, PHASE_AD_4_FILE_PATH_3, PHASE_PC_FILE_PATH_3, h, eps, iterParam);

}

int main(){
    temp_main<double>();
    return 0;
}