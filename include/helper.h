#ifndef HELPER_H
#define HELPER_H

#include "time.h"


using namespace std;
/*

extern struct timespec t0;
struct timespec elapsed (){
	struct timespec ts;
	clock_gettime (CLOCK_REALTIME, &ts);
	if (ts.tv_nsec < t0.tv_nsec){
		ts.tv_nsec = 1000000000 + ts.tv_nsec - t0.tv_nsec;
		ts.tv_sec--;
	}
	ts.tv_sec -= t0.tv_sec;
	return (ts);
}

int timelog (const char* message){
  struct timespec ts = elapsed ();
  return (printf ("[%06ld.%09ld] %s\n", ts.tv_sec, ts.tv_nsec, message));
}
*/
void * malloc_double_align(size_t n, unsigned int a /*alignment*/, double *& output){
    void *adres=NULL;
    void *adres2=NULL;
    adres=malloc(n*sizeof(double)+a);
    size_t adr=(size_t)adres;
    size_t adr2=adr+a-(adr&(a-1u)); 	// a valid address for a alignment
    adres2=(void * ) adr2;
    output=(double *)adres2;
    return adres;                		// pointer to be used in free()
}

void print_timenl () {
	clock_t c = clock();
	double t = double(c) / CLOCKS_PER_SEC;
	cout << "Time = " << t << endl ;
}

void print_time () {
	clock_t c = clock();
	double t = double(c) / CLOCKS_PER_SEC;
	cout << "Time = " << t  << " : ";
}


//solve for Ax=b efficiently, return A^{-1}b = x, where A is fixted to be sigma_g^2XX^T/M+sigma_eI
/*
void conjugate_gradient(int n, double vg, double ve,MatrixXdr &A,  MatrixXdr &b, MatrixXdr &x , int exist_ind){
        int k=0;
        double thres=0.0001;
        int max_iter=50;
        MatrixXdr r0(n, 1);
        MatrixXdr r1(n, 1);
        MatrixXdr p(n, 1);
        MatrixXdr s(n, 1);
        for(int i=0; i<n; i++)
        {       x(i,0)=0;
        }
        double temp=1;
        double beta,alpha;
        r0=b;
        r1=b;
        MatrixXdr mask_sum = A.colwise().sum();
        while(temp>thres && k<max_iter){
                k++;
                if(k==1)
                        p = b;
                else
                {
                        MatrixXdr temp1 = r0.transpose() * r0;
                        MatrixXdr temp2 = r1.transpose() * r1;
                        beta = temp2(0,0)/ temp1(0,0);
                        p = r1+ beta*p;
                }
                s=A*p ;
                MatrixXdr temp1 = r1.transpose() * r1;
                MatrixXdr temp2 = p.transpose()*s;
                alpha = temp1(0,0)/ temp2(0,0);
                x = x+alpha * p ;
                MatrixXdr r2= r1;
                r1 = r1 - alpha* s;
                r0 = r2;
                MatrixXdr z = r1.transpose() * r1;
                temp = z(0,0);
                cout<<"Iter: "<< k <<"  " << temp <<endl;
        }
}
*/

#endif
