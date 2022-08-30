#include<iostream>
#include<math.h>
#include<vector>
#include<chrono>
#include<gsl/gsl_blas.h>

using namespace std;
using namespace std::chrono;

struct converted_data_pack
{
    gsl_matrix firing_data;
    gsl_matrix not_firing_data;
    vector<float> weight;
    float firing_barrier;//must be more than this
    float not_firing_barrier;//must be less than this
};

void print_mat(gsl_matrix *mat)//ok check
{
    for(int a=0;a<mat->size1;a++)
    {
        for(int b=0;b<mat->size2;b++)
        {   cout<<gsl_matrix_get(mat,a,b)<<",";}
        cout<<"\n";
    }
}

gsl_matrix* calc_M(converted_data_pack &cdp,gsl_matrix *M_old,int firing_data_index,int not_firing_data_index)//ok check
{
    float n=cdp.firing_data.size1+cdp.not_firing_data.size1;
    double ai[cdp.firing_data.size2];
    if(firing_data_index!=-1)
    {
        for(int a=0;a<cdp.firing_data.size2;a++)
        {   ai[a]=gsl_matrix_get(&cdp.firing_data,firing_data_index,a);}
    }
    if(not_firing_data_index!=-1)
    {
        for(int a=0;a<cdp.not_firing_data.size2;a++)
        {   ai[a]=gsl_matrix_get(&cdp.not_firing_data,not_firing_data_index,a);}
    }
    //cout<<"\nn: "<<n;
    float part1=pow(n,2)/(pow(n,2)-1);
    //cout<<"\npart1: "<<part1;
    gsl_matrix *part2=gsl_matrix_alloc(cdp.firing_data.size2,1);//Mk*ai. Mk includes the identity matrix multiplied by the constant.
    gsl_matrix_view ai_mat = gsl_matrix_view_array(ai,cdp.firing_data.size2,1);
    //cout<<"\nai_mat: \n";
    //print_mat(&ai_mat.matrix);
    //cout<<"\nM_old: \n";
    //print_mat(M_old);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,M_old,&ai_mat.matrix,0.0,part2);
    //cout<<"\npart2: \n";
    //print_mat(part2);
    gsl_matrix *part3=gsl_matrix_alloc(part2->size2,part2->size1);
    gsl_matrix_transpose_memcpy(part3,part2);//part3
    //cout<<"\npart3: \n";
    //print_mat(part3);
    //done
    gsl_matrix *ai_mat_t=gsl_matrix_alloc(ai_mat.matrix.size2,ai_mat.matrix.size1);
    gsl_matrix_transpose_memcpy(ai_mat_t,&ai_mat.matrix);
    gsl_matrix *part4=gsl_matrix_alloc(1,1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,ai_mat_t,part2,0.0,part4);
    //cout<<"\npart4: \n";
    //print_mat(part4);

    gsl_matrix_scale(part2,1.0/gsl_matrix_get(part4,0,0));
    //cout<<"\nnew_part2:\n";
    //print_mat(part2);
    gsl_matrix *part5=gsl_matrix_alloc(part2->size1,part3->size2);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,part2,part3,0.0,part5);
    gsl_matrix_scale(part5,2.0/(1+n));
    //cout<<"\npart5:\n";
    //print_mat(part5);
    //gsl_matrix *M=gsl_matrix_alloc(M_old->size1,M_old->size2);
    gsl_matrix_sub(M_old,part5);
    gsl_matrix_scale(M_old,part1);
    //cout<<"\nM:\n";
    //print_mat(M_old);

    return M_old;
}

gsl_matrix* calc_X(converted_data_pack &cdp,gsl_matrix *M,gsl_matrix *X_old,int firing_data_index,int not_firing_data_index)//ok check
{
    //int gh;
    float n=cdp.firing_data.size1+cdp.not_firing_data.size1;
    //cout<<"\nn: "<<n;
    //cin>>gh;
    double ai[cdp.firing_data.size2];
    if(firing_data_index!=-1)
    {
        for(int a=0;a<cdp.firing_data.size2;a++)
        {   ai[a]=gsl_matrix_get(&cdp.firing_data,firing_data_index,a);}
    }
    if(not_firing_data_index!=-1)
    {
        for(int a=0;a<cdp.not_firing_data.size2;a++)
        {   ai[a]=gsl_matrix_get(&cdp.not_firing_data,not_firing_data_index,a);}
    }
    float part1=1.0/(n+1);
    //cout<<"\nPart 1: "<<part1;
    //cin>>gh;
    gsl_matrix *part2=gsl_matrix_alloc(cdp.firing_data.size2,1);//Mk*ai. Mk includes the identity matrix multiplied by the constant.
    gsl_matrix_view ai_mat = gsl_matrix_view_array(ai,cdp.firing_data.size2,1);
    //cout<<"\nai_mat:\n";
    //print_mat(&ai_mat.matrix);
    //cin>>gh;
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,M,&ai_mat.matrix,0.0,part2);
    //cout<<"\npart2: \n";
    //print_mat(part2);
    //cin>>gh;
    gsl_matrix *part3=gsl_matrix_alloc(1,1);
    gsl_matrix *ai_mat_t=gsl_matrix_alloc(ai_mat.matrix.size2,ai_mat.matrix.size1);
    gsl_matrix_transpose_memcpy(ai_mat_t,&ai_mat.matrix);
    //cout<<"\nai_mat_t: \n";
    //print_mat(ai_mat_t);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,ai_mat_t,part2,0.0,part3);
    //cout<<"\npart3: \n";
    //print_mat(part3);
    //cin>>gh;
    //part4
    gsl_matrix_scale(part2,part1*(1.0/sqrt(gsl_matrix_get(part3,0,0))));
    //cout<<"\npart4: \n";
    //print_mat(part2);
    //cin>>gh;

    for(int a=0;a<X_old->size1;a++)
    {   gsl_matrix_set(X_old,a,0,gsl_matrix_get(X_old,a,0)-gsl_matrix_get(part2,a,0));}
    //cout<<"\nX_new:\n";
    //print_mat(X_old);
    //cin>>gh;

    return X_old;
}

void check_if_solution_found(converted_data_pack &cdp,gsl_matrix *X,int &firing_data_index,int &not_firing_data_index)//ok check
{
    firing_data_index=-1;
    not_firing_data_index=-1;
    for(int a=0;a<cdp.firing_data.size1;a++)
    {
        float result=0;
        for(int b=0;b<cdp.firing_data.size2;b++)
        {
            result+=gsl_matrix_get(&cdp.firing_data,a,b)*gsl_matrix_get(X,b,0);
        }
        if(result>cdp.firing_barrier)//<
        {
            firing_data_index=a;
            break;
        }
    }
    if(firing_data_index==-1)
    {
        for(int a=0;a<cdp.not_firing_data.size1;a++)
        {
            float result=0;
            for(int b=0;b<cdp.not_firing_data.size2;b++)
            {
                result+=gsl_matrix_get(&cdp.not_firing_data,a,b)*gsl_matrix_get(X,b,0);
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
    for(int a=0;a<cdp.firing_data.size1;a++)
    {   l+=ceil(log2(1+abs(cdp.firing_barrier)));}
    for(int a=0;a<cdp.not_firing_data.size1;a++)
    {   l+=ceil(log2(1+abs(cdp.not_firing_barrier)));}
    //section1
    for(int a=0;a<cdp.firing_data.size1;a++)
    {
        for(int b=0;b<cdp.firing_data.size2;b++)
        {   l+=ceil(log2(1+abs(gsl_matrix_get(&cdp.firing_data,a,b))));}
    }
    for(int a=0;a<cdp.not_firing_data.size1;a++)
    {
        for(int b=0;b<cdp.not_firing_data.size2;b++)
        {   l+=ceil(log2(1+abs(gsl_matrix_get(&cdp.not_firing_data,a,b))));}
    }
    int n=cdp.firing_data.size1+cdp.not_firing_data.size1;
    int m=cdp.firing_data.size2;
    l+=ceil(log2(n));
    l+=ceil(log2(m));
    l+=(n*m+m);

    return l;
}

void run_ellipsoid_method(converted_data_pack& cdp)
{
    //setting up initial weight
    for(int a=0;a<cdp.firing_data.size2;a++)
    {   cdp.weight.push_back(0);}
    //setting up initial matrix M, X, 2^l
    int L=calc_l(cdp);
    cout<<"\nL: "<<L;
    L=L/6;
    cout<<"\nnew L: "<<L;
    float two_power_l=pow(2,L);
    cout<<"\ntwo_power_l: "<<two_power_l;
    gsl_matrix *M=gsl_matrix_alloc(cdp.weight.size(),cdp.weight.size());
    gsl_matrix *X=gsl_matrix_alloc(cdp.weight.size(),1);
    for(int a=0;a<cdp.weight.size();a++)
    {
        for(int b=0;b<cdp.weight.size();b++)
        {
            if(a==b)
            {   gsl_matrix_set(M,a,b,two_power_l);}
            else
            {   gsl_matrix_set(M,a,b,0);}
        }
        gsl_matrix_set(X,a,0,0);
    }
    unsigned long int max_iteration=6*pow((cdp.firing_data.size1+cdp.not_firing_data.size2+1),2)*L;
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
    cdp.firing_data=gsl_matrix_view_array(fd,1,4).matrix;
    cdp.not_firing_data=gsl_matrix_view_array(nfd,1,4).matrix;
    auto start = high_resolution_clock::now();
    run_ellipsoid_method(cdp);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start); 
    cout<<"\n\nduration= "<<duration.count()<<" microseconds";
}
//g++ ellipsoid_gsl.cpp -std=c++17 -lgsl -O2