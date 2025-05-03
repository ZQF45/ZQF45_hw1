#include "algebra.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>  //引入标准库函数

#define EPSILON 1e-8  //引入一个很小的阈值，用于判断浮点数是否为0，因为浮点数直接比较可能不准确


//创建矩阵
Matrix create_matrix(int rows, int cols) {   
    Matrix m;
    m.rows = rows;
    m.cols = cols;
    // 初始化为0
    memset(m.data, 0, sizeof(m.data));  //清零，memset函数是按字节操作的，这里是把m.data中的前sizeof(m.data)个字节赋为零
    return m;
}


//加法
Matrix add_matrix(Matrix a, Matrix b) {              
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return result;
}


//减法
Matrix sub_matrix(Matrix a, Matrix b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return result;
}

//乘法
Matrix mul_matrix(Matrix a, Matrix b) {
    if (a.cols != b.rows) {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, b.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < b.cols; j++) {
            for (int k = 0; k < a.cols; k++) {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return result;
}


//数乘
Matrix scale_matrix(Matrix a, double k) {
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] * k;
        }
    }
    return result;
}

//转置
Matrix transpose_matrix(Matrix a) {
    Matrix result = create_matrix(a.cols, a.rows);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[j][i] = a.data[i][j];
        }
    }
    return result;
}

//高斯消元法化下三角行列式（AI写的）
static void gaussian_elimination(Matrix *a, int *sign) {
    int n = a->rows;
    
    for (int i = 0; i < n; i++) {
        int max_row = i;
        //寻找首列中绝对值最大的行
        for (int k = i; k < n; k++) {
            if (fabs(a->data[k][i]) > fabs(a->data[max_row][i])) {
                max_row = k;
            }
        }
        //如果找到的主元不在当前行，就交换行
        if (max_row != i) {
            for (int j = 0; j < n; j++) {
                double tmp = a->data[i][j];
                a->data[i][j] = a->data[max_row][j];
                a->data[max_row][j] = tmp;
            }
            *sign *= -1;  //行交换会改变行列式符号
        }
        //如果主元很接近零，行列式就为零
        if (fabs(a->data[i][i]) < EPSILON) {
            return;
        }
        //消去下方行的第i个元素
        for (int k = i + 1; k < n; k++) {
            double factor = a->data[k][i] / a->data[i][i];
            for (int j = i; j < n; j++) {
                a->data[k][j] -= factor * a->data[i][j];
            }
        }
    }
}

//求行列式，先高斯消元法化成下三角行列式，再把主对角线元素相乘
double det_matrix(Matrix a) {
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    int n = a.rows;
    Matrix temp = a;
    int sign = 1;
    gaussian_elimination(&temp, &sign);
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        det *= temp.data[i][i];
    }
    return det * sign;
}

//求逆（AI写的），创建增广矩阵，再用高斯消元法求解
Matrix inv_matrix(Matrix a) {
    Matrix inv = create_matrix(0, 0);
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return inv;
    }
    int n = a.rows;
    double det = det_matrix(a);
    if (fabs(det) < EPSILON) {
        printf("Error: The matrix is singular.\n");
        return inv;
    }
    Matrix aug = create_matrix(n, 2 * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug.data[i][j] = a.data[i][j];
            aug.data[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (int i = 0; i < n; i++) {
        int max_row = i;
        for (int k = i; k < n; k++) {
            if (fabs(aug.data[k][i]) > fabs(aug.data[max_row][i])) {
                max_row = k;
            }
        }
        if (max_row != i) {
            for (int j = 0; j < 2 * n; j++) {
                double tmp = aug.data[i][j];
                aug.data[i][j] = aug.data[max_row][j];
                aug.data[max_row][j] = tmp;
            }
        }
        double pivot = aug.data[i][i];
        if (fabs(pivot) < EPSILON) {
            printf("Error: The matrix is singular.\n");
            return create_matrix(0, 0);
        }
        for (int j = 0; j < 2 * n; j++) {
            aug.data[i][j] /= pivot;
        }
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = aug.data[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    aug.data[k][j] -= factor * aug.data[i][j];
                }
            }
        }
    }
    inv = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv.data[i][j] = aug.data[i][j + n];
        }
    }
    return inv;
}

//求秩（AI写的），也使用了高斯消元法
int rank_matrix(Matrix a) {
    int rank = 0;
    Matrix temp = a;
    int rows = a.rows, cols = a.cols;
    for (int col = 0; col < cols && rank < rows; col++) {
        int pivot = -1;
        for (int i = rank; i < rows; i++) {
            if (fabs(temp.data[i][col]) > EPSILON) {
                pivot = i;
                break;
            }
        }
        if (pivot == -1) continue;
        if (pivot != rank) {
            for (int j = 0; j < cols; j++) {
                double tmp = temp.data[rank][j];
                temp.data[rank][j] = temp.data[pivot][j];
                temp.data[pivot][j] = tmp;
            }
        }
        for (int i = rank + 1; i < rows; i++) {
            double factor = temp.data[i][col] / temp.data[rank][col];
            for (int j = col; j < cols; j++) {
                temp.data[i][j] -= factor * temp.data[rank][j];
            }
        }
        rank++;
    }
    return rank;
}

//求矩阵的迹
double trace_matrix(Matrix a) {
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    double trace = 0;
    for (int i = 0; i < a.rows; i++) {
        trace += a.data[i][i];
    }
    return trace;
}

//打印矩阵 
void print_matrix(Matrix a) {
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}



/*#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    // ToDo
    return create_matrix(0, 0);
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    // ToDo
    return create_matrix(0, 0);
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    // ToDo
    return create_matrix(0, 0);
}

Matrix scale_matrix(Matrix a, double k)
{
    // ToDo
    return create_matrix(0, 0);
}

Matrix transpose_matrix(Matrix a)
{
    // ToDo
    return create_matrix(0, 0);
}

double det_matrix(Matrix a)
{
    // ToDo
    return 0;
}

Matrix inv_matrix(Matrix a)
{
    // ToDo
    return create_matrix(0, 0);
}

int rank_matrix(Matrix a)
{
    // ToDo
    return 0;
}

double trace_matrix(Matrix a)
{
    // ToDo
    return 0;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}
    */


