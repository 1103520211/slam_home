// #include <iostream>
// #include <ctime>
// // Eigen 核心部分
// #include <Eigen/Core>
// // 稠密矩阵的代数运算（逆，特征值等）
// #include <Eigen/Dense>

// using namespace std;
// using namespace Eigen;

// #define MATRIX_SIZE 100

// int main(int argc, char **argv) {

//     MatrixXd matrix_NN = MatrixXd::Random(MATRIX_SIZE, MATRIX_SIZE);
//     matrix_NN = matrix_NN * matrix_NN.transpose();  // 对称阵，保证半正定
// //    cout << "matrix_NN:\n" << matrix_NN << endl;

//     // 实对称矩阵可以正交相似对角化，保证对角化成功
//     SelfAdjointEigenSolver<MatrixXd> eigen_solver(matrix_NN);
//     c
// //    The eigenvalues are sorted in increasing order. 所有特征值均大于0，保证正定
//     cout << "Eigen values min = \n" << eigen_solver.eigenvalues()(0) << endl;

//     clock_t time_stt = clock(); // 计时
//     Matrix<double, MATRIX_SIZE, 1> v_Nd = MatrixXd::Random(MATRIX_SIZE, 1);
//     Matrix<double, MATRIX_SIZE, 1> x;

//     // 通常用矩阵分解来求，例如QR分解，速度会快很多
//     time_stt = clock();
//     x = matrix_NN.colPivHouseholderQr().solve(v_Nd);
//     cout << "time of Qr decomposition is "
//          << 1000 * (clock() - time_stt) / (double) CLOCKS_PER_SEC << "ms" << endl;
//     cout << "x = " << x.transpose() << endl;

//     // 对于正定矩阵，还可以用cholesky分解来解方程
//     time_stt = clock();
//     x = matrix_NN.ldlt().solve(v_Nd);
//     cout << "time of ldlt decomposition is "
//          << 1000 * (clock() - time_stt) / (double) CLOCKS_PER_SEC << "ms" << endl;
//     cout << "x = " << x.transpose() << endl;

//     return 0;
// }

#include <iostream>
#include <ctime>
// Eigen 核心部分
#include <eigen3/Eigen/Core>
// 稠密矩阵的代数运算（逆，特征值等）
#include <Eigen/Dense>
//xiugai
using namespace std;
using namespace Eigen;

#define Matrix_Size 100

int main(int argc, char **argv) {
    MatrixXd matrix_NN = MatrixXd::Random(Matrix_Size, Matrix_Size);
    Matrix<double,Matrix_Size,1> matrix_N = MatrixXd::Random(Matrix_Size, 1);
    matrix_NN = matrix_NN * matrix_NN.transpose();
    cout << "matrix_NN:\n" << matrix_NN << endl;
    SelfAdjointEigenSolver<MatrixXd> eigen_slover(matrix_NN);
    cout << "eigen value:/n" << eigen_slover.eigenvalues()(0) << endl;
    clock_t time1 = clock();
    Matrix<double,Matrix_Size,1> x; 
    x = matrix_NN.colPivHouseholderQr().solve(matrix_N);//qr
    cout << "x = " << x.transpose() << endl;
    
    x = matrix_NN.ldlt().solve(matrix_N);
    cout << "x = " << x.transpose() << endl;



}