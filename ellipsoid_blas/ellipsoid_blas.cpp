#include<iostream>
#include<math.h>
#include<vector>
#include<chrono>
#include"blas_wrapper.h"

using namespace std;
using namespace std::chrono;

struct converted_data_pack
{
    blas_matrix *firing_data;
    blas_matrix *not_firing_data;
    vector<double> weight;
    double firing_barrier;//must be more than this
    double not_firing_barrier;//must be less than this
};

blas_matrix calc_M(converted_data_pack &cdp,blas_matrix M_old,int firing_data_index,int not_firing_data_index)//ok check
{
    double n=cdp.firing_data->rows+cdp.not_firing_data->rows;
    double ai[cdp.firing_data->cols];
    if(firing_data_index!=-1)
    {
        for(int a=0;a<cdp.firing_data->cols;a++)
        {   ai[a]=cdp.firing_data->get(firing_data_index,a);}
    }
    if(not_firing_data_index!=-1)
    {
        for(int a=0;a<cdp.not_firing_data->cols;a++)
        {   ai[a]=cdp.not_firing_data->get(not_firing_data_index,a);}
    }
    //cout<<"\nn: "<<n;
    double part1=pow(n,2)/(pow(n,2)-1);
    //cout<<"\npart1: "<<part1;
    blas_matrix ai_mat(ai,cdp.firing_data->cols,1);
    //cout<<"\nai_mat: \n";
    //print_mat(&ai_mat.matrix);
    //cout<<"\nM_old: \n";
    //print_mat(M_old);
    blas_matrix part2=matmul(M_old,ai_mat);
    //cout<<"\npart2: \n";
    //print_mat(part2);
    blas_matrix part3=transpose(part2);
    //cout<<"\npart3: \n";
    //print_mat(part3);
    //done
    blas_matrix ai_mat_t=transpose(ai_mat);
    blas_matrix part4=matmul(ai_mat_t,part2);
    //cout<<"\npart4: \n";
    //print_mat(part4);

    part2=matmul(part2,1.0/part4.get(0,0));
    //cout<<"\nnew_part2:\n";
    //print_mat(part2);
    blas_matrix part5=matmul(part2,part3);
    part5=matmul(part5,2.0/(1+n));
    //cout<<"\npart5:\n";
    //print_mat(part5);
    //gsl_matrix *M=gsl_matrix_alloc(M_old->size1,M_old->size2);
    M_old=matsub(M_old,part5);
    M_old=matmul(M_old,part1);
    //cout<<"\nM:\n";
    //print_mat(M_old);

    return M_old;
}

blas_matrix calc_X(converted_data_pack &cdp,blas_matrix M,blas_matrix X_old,int firing_data_index,int not_firing_data_index)//ok check
{
    //int gh;
    double n=cdp.firing_data->rows+cdp.not_firing_data->rows;
    //cout<<"\nn: "<<n;
    //cin>>gh;
    double ai[cdp.firing_data->cols];
    if(firing_data_index!=-1)
    {
        for(int a=0;a<cdp.firing_data->cols;a++)
        {   ai[a]=cdp.firing_data->get(firing_data_index,a);}
    }
    if(not_firing_data_index!=-1)
    {
        for(int a=0;a<cdp.not_firing_data->cols;a++)
        {   ai[a]=cdp.not_firing_data->get(not_firing_data_index,a);}
    }
    double part1=1.0/(n+1);
    //cout<<"\nPart 1: "<<part1;
    //cin>>gh;
    blas_matrix ai_mat(ai,cdp.firing_data->cols,1);
    //cout<<"\nai_mat:\n";
    //print_mat(&ai_mat.matrix);
    //cin>>gh;
    blas_matrix part2=matmul(M,ai_mat);
    //cout<<"\npart2: \n";
    //print_mat(part2);
    //cin>>gh;
    blas_matrix ai_mat_t=transpose(ai_mat);
    //cout<<"\nai_mat_t: \n";
    //print_mat(ai_mat_t);
    blas_matrix part3=matmul(ai_mat_t,part2);
    //cout<<"\npart3: \n";
    //print_mat(part3);
    //cin>>gh;
    //part4
    part2=matmul(part2,part1*(1.0/sqrt(part3.get(0,0))));
    //cout<<"\npart4: \n";
    //print_mat(part2);
    //cin>>gh;

    for(int a=0;a<X_old.rows;a++)
    {   X_old.set(a,0,X_old.get(a,0)-part2.get(a,0));}
    //cout<<"\nX_new:\n";
    //print_mat(X_old);
    //cin>>gh;

    return X_old;
}

void check_if_solution_found(converted_data_pack &cdp,blas_matrix X,int &firing_data_index,int &not_firing_data_index)//ok check
{
    firing_data_index=-1;
    not_firing_data_index=-1;
    for(int a=0;a<cdp.firing_data->rows;a++)
    {
        double result=0;
        for(int b=0;b<cdp.firing_data->cols;b++)
        {
            result+=cdp.firing_data->get(a,b)*X.get(b,0);
        }
        if(result>cdp.firing_barrier)//<
        {
            firing_data_index=a;
            break;
        }
    }
    if(firing_data_index==-1)
    {
        for(int a=0;a<cdp.not_firing_data->rows;a++)
        {
            double result=0;
            for(int b=0;b<cdp.not_firing_data->cols;b++)
            {
                result+=cdp.not_firing_data->get(a,b)*X.get(b,0);
            }
            if(result>cdp.not_firing_barrier)
            {
                not_firing_data_index=a;
                break;
            }
        }
    }
}

int calc_l(converted_data_pack &cdp)//ok check
{
    int l=0;
    //section2
    for(int a=0;a<cdp.firing_data->rows;a++)
    {   l+=ceil(log2(1+abs(cdp.firing_barrier)));}
    for(int a=0;a<cdp.not_firing_data->rows;a++)
    {   l+=ceil(log2(1+abs(cdp.not_firing_barrier)));}
    //section1
    for(int a=0;a<cdp.firing_data->rows;a++)
    {
        for(int b=0;b<cdp.firing_data->cols;b++)
        {   l+=ceil(log2(1+abs(cdp.firing_data->get(a,b))));}
    }
    for(int a=0;a<cdp.not_firing_data->rows;a++)
    {
        for(int b=0;b<cdp.not_firing_data->cols;b++)
        {   l+=ceil(log2(1+abs(cdp.not_firing_data->get(a,b))));}
    }
    int n=cdp.firing_data->rows+cdp.not_firing_data->rows;
    int m=cdp.firing_data->cols;
    l+=ceil(log2(n));
    l+=ceil(log2(m));
    l+=(n*m+m);

    return l;
}

void run_ellipsoid_method(converted_data_pack& cdp)
{
    //setting up initial weight
    for(int a=0;a<cdp.firing_data->cols;a++)
    {   cdp.weight.push_back(0);}
    //setting up initial matrix M, X, 2^l
    int L=calc_l(cdp);
    cout<<"\nL: "<<L;
    L=L/6;
    cout<<"\nnew L: "<<L;
    double two_power_l=pow(2,L);
    cout<<"\ntwo_power_l: "<<two_power_l;
    blas_matrix M(cdp.weight.size(),cdp.weight.size(),false);
    blas_matrix X(cdp.weight.size(),1,false);
    for(int a=0;a<cdp.weight.size();a++)
    {
        for(int b=0;b<cdp.weight.size();b++)
        {
            if(a==b)
            {   M.set(a,b,two_power_l);}
            else
            {   M.set(a,b,0);}
        }
        X.set(a,0,0);
    }
    unsigned long int max_iteration=6*pow((cdp.firing_data->rows+cdp.not_firing_data->cols+1),2)*L;
    cout<<"\nmax_iteration: "<<max_iteration;
    int firing_data_index,not_firing_data_index;
    unsigned long int k=0;
    while(k<max_iteration)//volume checking condition
    {
        check_if_solution_found(cdp,X,firing_data_index,not_firing_data_index);
        cout<<"\niteration: "<<k;
        cout<<"\nfiring_data_index: "<<firing_data_index<<" not_firing_data_index: "<<not_firing_data_index;
        cout<<"\nX:";
        print_mat(X);
        cout<<"\nM:";
        print_mat(M);
        //cout<<"\ncin>>";
        //int gh;cin>>gh;
        if(firing_data_index==-1 && not_firing_data_index==-1)//solution_found
        {   break;}
        X=calc_X(cdp,M,X,firing_data_index,not_firing_data_index);
        M=calc_M(cdp,M,firing_data_index,not_firing_data_index);
        k++;
    }
    if(k==max_iteration)
    {   cout<<"Solution not found. Ran for "<<k<<" iterations";}
    else
    {
        cout<<"\nsolution (X): \n";
        print_mat(X);
    }
    //int gh;cin>>gh;
    //testing
    //cout<<"\nl: "<<calc_l(cdp);
    //cout<<"\nM fx called: \n";
    //calc_M(cdp,M,0,-1);
    //calc_X(cdp,M,X,0,-1);
}

int main()
{
    converted_data_pack cdp;
    //double fd[2]={-1,0};
    //double nfd[2]={0,-1};
    //cdp.firing_barrier=-1;
    //cdp.firing_data=gsl_matrix_view_array(fd,1,2).matrix;
    //cdp.not_firing_barrier=-1;
    //cdp.not_firing_data=gsl_matrix_view_array(nfd,1,2).matrix;
    
    //double nfd[3]={2,3,6};
    //double fd[3]={-1,-4,-5};
    //cdp.firing_barrier=-40;
    //cdp.not_firing_barrier=60;
    //cdp.firing_data=gsl_matrix_view_array(fd,1,3).matrix;
    //cdp.not_firing_data=gsl_matrix_view_array(nfd,1,3).matrix;

    double nfd[4]={-3,1,-2,0};
    double fd[4]={-1,-1,+4,-2};
    cdp.firing_barrier=-4;
    cdp.not_firing_barrier=6;
    cdp.firing_data=new blas_matrix(fd,1,4);
    cdp.not_firing_data=new blas_matrix(nfd,1,4);
    auto start = high_resolution_clock::now();
    run_ellipsoid_method(cdp);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start); 
    cout<<"\n\nduration= "<<duration.count()<<" microseconds";
}
//g++ blas_wrapper.cpp ellipsoid_blas.cpp -std=c++17 -lblas -O2