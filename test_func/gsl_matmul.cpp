#include<iostream>
#include<gsl/gsl_blas.h>
//#include<cblas.h>

using namespace std;

gsl_matrix* matmul(double *A,int row_a,int col_a,double *B,int row_b,int col_b,int &row_c,int &col_c) 
{
    row_c=row_a;
    col_c=col_b;

    gsl_matrix_view A_mat = gsl_matrix_view_array(A,row_a,col_a);
    gsl_matrix_view B_mat = gsl_matrix_view_array(B,row_b,col_b);
    gsl_matrix *C_mat = gsl_matrix_alloc (row_c,col_c);
    
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,
    1.0,&A_mat.matrix,&B_mat.matrix,0.0,C_mat);
    
    return C_mat;
}

int main()
{
    double A[6]={1,2,3,4,5,6};
    double B[8]={1,2,3,4,5,6,7,8};
    int row_a=3,col_a=2;
    int row_b=2,col_b=4;
    int row_c,col_c;
    gsl_matrix *C=matmul(A,row_a,col_a,B,row_b,col_b,row_c,col_c);
    cout<<"\nrow_c: "<<row_c<<" col_c: "<<col_c<<endl;
    //int gh;cin>>gh;
    for(int a=0;a<row_c;a++)
    {
        for(int b=0;b<col_c;b++)
        {
            cout<<gsl_matrix_get(C,a,b)<<",";
        }
        cout<<"\n";
    }
    return  0;
}

//sudo apt-get install libblas-dev libblas64-dev libatlas-base-dev liblapack-dev libopenblas-dev libgsl-dev
//-lblas
//-lgsl