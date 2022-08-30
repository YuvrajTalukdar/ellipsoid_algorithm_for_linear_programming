#include<iostream>
#include<cblas.h>

struct blas_matrix
{
    double *matrix;
    int rows,cols;
    double get(int row,int col);//ok tested
    void set(int row,int col,double value);//ok tested
    void clear();//ok tested
    blas_matrix(double *array,int rows1,int cols1);//ok tested
    blas_matrix(int rows1,int cols1,bool initialize_with_1);//ok tested
};

blas_matrix matmul(blas_matrix x,double factor);//ok tested

blas_matrix matmul(blas_matrix x,blas_matrix y);//ok tested

blas_matrix transpose(blas_matrix x);//ok tested

blas_matrix matsub(blas_matrix x,blas_matrix y);//x=x-y//ok check

blas_matrix matadd(blas_matrix x,blas_matrix y);//x=x-y//ok check

void print_mat(blas_matrix x);//ok check