#include<iostream>
#include<cblas.h>

using namespace std;

struct blas_matrix
{
    double *matrix;
    int rows,cols;
    double get(int row,int col)//ok tested
    {   return matrix[cols*row+col];}
    void set(int row,int col,double value)//ok tested
    {   matrix[cols*row+col]=value;}
    void clear()//ok tested
    {   free(matrix);rows=0;cols=0;}
    blas_matrix(double *array,int rows1,int cols1)//ok tested
    {
        matrix=array;
        rows=rows1;
        cols=cols1;
    }
    blas_matrix(int rows1,int cols1,bool initialize_with_1)//ok tested
    {
        matrix=(double*)malloc(sizeof(double)*rows1*cols1);
        rows=rows1;
        cols=cols1;
        if(initialize_with_1)
        {
            for(int a=0;a<rows;a++)
            {
                for(int b=0;b<cols;b++)
                {
                    if(a==b)
                    {   matrix[a*cols+b]=1;}
                    else
                    {   matrix[a*cols+b]=0;}
                }
            }
        }
    }
};

blas_matrix matmul(blas_matrix x,double factor)//ok tested
{
    cblas_dscal(x.rows*x.cols,factor,x.matrix,1);//for multiplying an array by constant. array_length,constant,array,1(just put 1 here, i dont know why.)   
    return x;
}

blas_matrix matmul(blas_matrix x,blas_matrix y)//ok tested
{
    blas_matrix z(x.rows,y.cols,false);
    int m=x.rows;
    int k=x.cols;
    int n=y.cols;
    double alpha=1.0,beta=0.0;
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,alpha,x.matrix,k,y.matrix,n,beta,z.matrix,n);
    return z;
}

blas_matrix transpose(blas_matrix x)//ok tested
{
    blas_matrix y(x.cols,x.rows,false);//identity matrix
    for(int a=0;a<y.rows;a++)
    {
        for(int b=0;b<y.cols;b++)
        {   y.set(a,b,x.get(b,a));}
    }
    return y;
}

blas_matrix matsub(blas_matrix x,blas_matrix y)//x=x-y//ok check
{
    for(int a=0;a<x.rows;a++)
    {
        for(int b=0;b<x.cols;b++)
        {   x.set(a,b,x.get(a,b)-y.get(a,b));}
    }
    return x;
}

blas_matrix matadd(blas_matrix x,blas_matrix y)//x=x+y//ok check
{
    for(int a=0;a<x.rows;a++)
    {
        for(int b=0;b<x.cols;b++)
        {   x.set(a,b,x.get(a,b)+y.get(a,b));}
    }
    return x;
}

void print_mat(blas_matrix x)//ok check
{
    for(int a=0;a<x.rows;a++)
    {
        for(int b=0;b<x.cols;b++)
        {   cout<<x.get(a,b)<<",";}
        cout<<"\n";
    }
}

void test()
{
    double x[6]={1,2,3,4,5,6};
    blas_matrix x_mat(x,2,3);
    double factor=1;
    x_mat=matmul(x_mat,factor);
    cout<<"\nx_mat: \n";
    print_mat(x_mat);

    double y[6]={1,2,3,4,5,6};
    blas_matrix y_mat(y,3,2);
    cout<<"\ny_mat:\n";
    print_mat(y_mat);

    blas_matrix z_mat=matmul(x_mat,y_mat);
    cout<<"\nz_mat: \n";
    print_mat(z_mat);

    blas_matrix x_transpose=transpose(x_mat);
    cout<<"\nx_transpose:\n";
    print_mat(x_transpose);

    blas_matrix I_mat(4,4,true);
    cout<<"\nI_mat:\n";
    print_mat(I_mat);
    
    matadd(x_mat,x_mat);
    cout<<"\nx_mat after add:\n";
    print_mat(x_mat);

    matsub(x_mat,x_mat);
    cout<<"\nx_mat after sub:\n";
    print_mat(x_mat);

    I_mat.clear();
    //x_mat.clear();//no need to free this one
    //y_mat.clear();//no need to free this one
    z_mat.clear();
    x_transpose.clear();
}

int main()
{
    test();

    return 0;
}