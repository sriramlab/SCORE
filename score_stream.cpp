/** 
 Mailman algorithm is written by Aman Agrawal 
 (Indian Institute of Technology, Delhi)
 RHE_reg, SCORE is written by Yue Ariel Wu
 (UCLA)
*/
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector> 
//#include <random>

#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>
#include "time.h"

#include "genotype.h"
#include "mailman.h"
#include "arguments.h"
#include "helper.h"
#include "storage.h"
/*
#include "boost/random.hpp"
#include "boost/accumulators/statistics/stats.hpp"
#include "boost/math/distributions/chi_squared.hpp"
#include "boost/math/distributions/normal.hpp"
#include "boost/math/special_functions.hpp"
*/
#if SSE_SUPPORT==1
	#define fastmultiply fastmultiply_sse
	#define fastmultiply_pre fastmultiply_pre_sse
#else
	#define fastmultiply fastmultiply_normal
	#define fastmultiply_pre fastmultiply_pre_normal
#endif

using namespace Eigen;
using namespace std;

// Storing in RowMajor Form
typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;
//Intermediate Variables
int blocksize;
double *partialsums;
double *sum_op;		
double *yint_e;
double *yint_m;
double **y_e;
double **y_m;

//for partition
double **y_e2; 
double **yint_e2; 
struct timespec t0;

//clock_t total_begin = clock();
MatrixXdr pheno;
MatrixXdr pheno_prime; 
MatrixXdr covariate;  
genotype g;
MatrixXdr geno_matrix; //(p,n)
int MAX_ITER;
int k,p,n;
int k_orig;

MatrixXdr c; //(p,k)
MatrixXdr x; //(k,n)
MatrixXdr v; //(p,k)
MatrixXdr means; //(p,1)
MatrixXdr stds; //(p,1)
MatrixXdr sum2;
MatrixXdr sum;  


//solve Ax= b, use for b
//
MatrixXdr yy; 
MatrixXdr yKy; 
MatrixXdr Xy;
MatrixXdr all_yKy;//use all_yKy to replace Xy, that compute on the fly   
//use for covariate
MatrixXdr WW; 

//use for missing phenotype
MatrixXdr pheno_mask2; 
vector<int> pheno_mask; 

//use for genetic correlation
vector<vector<int> > gc_index; 
options command_line_opts;

MatrixXdr FunCat; 
bool debug = false;
bool check_accuracy = false;
bool var_normalize=true;
int accelerated_em=0;
double convergence_limit;
bool memory_efficient = false;
bool missing=false;
bool fast_mode = true;
bool text_version = false;
bool use_cov=false;
bool compute_gc=false;  
bool reg = true;
bool gwas=false;
bool bpheno=false;   
bool pheno_fill=false;  // if pheno_fill, estimate tr[k2] once and use for all phenotypes 
bool noh2g=false;

int total_num_blocks=0;  
vector<string> pheno_name; 

std::istream& newline(std::istream& in)
{
    if ((in >> std::ws).peek() != std::char_traits<char>::to_int_type('\n')) {
        in.setstate(std::ios_base::failbit);
    }
    return in.ignore();
}
int read_gc(std::string filename)
{
	ifstream ifs(filename.c_str(), ios::in); 
	std::string line; 
	std::istringstream in; 

	int count=0; 
	while(std::getline(ifs, line)){
		in.clear(); 
		in.str(line); 
		string temp1, temp2; 
		in>>temp1; in>>temp2; 
		gc_index.push_back(vector<int>()); 
		gc_index[count].push_back(atoi(temp1.c_str())); 
		gc_index[count].push_back(atoi(temp2.c_str()));
		count++;

	}
	for(int i=0; i<count; i++)
		cout<<gc_index[i][0]<<endl; 
	return count; 
}
int read_cov(bool std,int Nind, std::string filename, std::string covname){
	ifstream ifs(filename.c_str(), ios::in); 
	std::string line; 
	std::istringstream in; 
	int covIndex = 0; 
	std::getline(ifs,line); 
	in.str(line); 
	string b;
	vector<vector<int> > missing; 
	int covNum=0;  
	while(in>>b)
	{
		if(b!="FID" && b !="IID"){
		missing.push_back(vector<int>()); //push an empty row  
		if(b==covname && covname!="")
			covIndex=covNum; 
		covNum++; 
		}
	}
	vector<double> cov_sum(covNum, 0); 
	if(covname=="")
	{
		covariate.resize(Nind, covNum); 
		cout<< "Read in "<<covNum << " Covariates.. "<<endl;
	}
	else 
	{
		covariate.resize(Nind, 1); 
		cout<< "Read in covariate "<<covname<<endl;  
	}

	
	int j=0; 
	while(std::getline(ifs, line)){
		in.clear(); 
		in.str(line);
		string temp;
		in>>temp; in>>temp; //FID IID 
		for(int k=0; k<covNum; k++){
			
			in>>temp;
			if(temp=="NA")
			{
				missing[k].push_back(j);
				continue; 
			} 
			double cur = atof(temp.c_str()); 
			if(cur==-9)
			{
				missing[k].push_back(j); 
				continue; 
			}
			if(covname=="")
			{
				cov_sum[k]= cov_sum[k]+ cur; 
				covariate(j,k) = cur; 
			}
			else
				if(k==covIndex)
				{
					covariate(j, 0) = cur;
					cov_sum[k] = cov_sum[k]+cur; 
				}
		}
		//if(j<10) 
		//	cout<<covariate.block(j,0,1, covNum)<<endl; 
		j++;
	}
	//compute cov mean and impute 
	for (int a=0; a<covNum ; a++)
	{
		int missing_num = missing[a].size(); 
		cov_sum[a] = cov_sum[a] / (Nind - missing_num);

		for(int b=0; b<missing_num; b++)
		{
                        int index = missing[a][b];
                        if(covname=="")
                                covariate(index, a) = cov_sum[a];
                        else if (a==covIndex)
                                covariate(index, 0) = cov_sum[a];
                } 
	}
	if(std)
	{
		MatrixXdr cov_std;
		cov_std.resize(1,covNum);  
		MatrixXdr sum = covariate.colwise().sum();
		MatrixXdr sum2 = (covariate.cwiseProduct(covariate)).colwise().sum();
		MatrixXdr temp;
//		temp.resize(Nind, 1); 
//		for(int i=0; i<Nind; i++)
//			temp(i,0)=1;  
		for(int b=0; b<covNum; b++)
		{
			cov_std(0,b) = sum2(0,b) + Nind*cov_sum[b]*cov_sum[b]- 2*cov_sum[b]*sum(0,b);
			cov_std(0,b) =sqrt((Nind- 1)/cov_std(0,b)) ;
			double scalar=cov_std(0,b); 
			for(int j=0; j<Nind; j++)
			{
				covariate(j,b) = covariate(j,b)-cov_sum[b];  
				covariate(j,b) =covariate(j,b)*scalar;
			} 
			//covariate.col(b) = covariate.col(b) -temp*cov_sum[b];
			
		}
	}	
	return covNum; 
}
int read_annotation(int Nsnp, std::string filename){
        ifstream ifs(filename.c_str(), ios::in);
        std::string line;
        std::istringstream in;
        int group_count =0;

	std::getline(ifs, line);
        in.str(line);
        string b;
        in>>b; in>>b;
        while(in>>b)
                group_count++;

        FunCat.resize(group_count, Nsnp);
        int i=0;
        while(std::getline(ifs,line)){
                in.clear();
                in.str(line);
                string temp;
                in>>temp; in>>temp ;
                for(int j=0; j<group_count; j++)
                {
                        in>>temp;
                        double cur = atof(temp.c_str());
                        FunCat(j, i)  = cur;
                }
                i++;
        }
        return group_count;
}
int read_pheno2(int Nind, std::string filename,std::string pheno_idx, bool pheno_fill){
	//if pheno_idx !="", we read in a pair of phenotypes 
	std::string phenoName1 =""; 
	std::string phenoName2 ="";
	int idx1= -1; 
	int idx2 = -1;  
	if(pheno_idx!=""){
	std::string delimiter =","; 
	size_t pos = 0; 
	pos= pheno_idx.find(delimiter); 
	if(pos != std::string::npos) 
		phenoName1 = pheno_idx.substr(0, pos); 
	cout<<"Analysing phenotype "<<phenoName1; 
	pheno_idx.erase(0, pos + delimiter.length()); 

	phenoName2 = pheno_idx; 	
	cout<<", "<< phenoName2<<endl; 
	

	}
//	pheno.resize(Nind,1); 
	ifstream ifs(filename.c_str(), ios::in); 
	
	std::string line;
	std::istringstream in;  
	int phenocount=0; 
	vector<vector<int> > missing; 
//read header
	std::getline(ifs,line); 
	in.str(line); 
	string b; 
	while(in>>b)
	{
		if(b!="FID" && b !="IID"){
			if(b == phenoName1)
				idx1 = phenocount; 
			if(b == phenoName2) 
				idx2 = phenocount; 
			phenocount++;
			missing.push_back(vector<int>());  
			pheno_name.push_back(b); 
		}
	}
	//if(pheno_idx !=0)
        //		pheno_name[0] = pheno_name[pheno_idx-1];   
	vector<double> pheno_sum(phenocount,0); 
	if(pheno_idx !=""){
		pheno.resize(Nind,2);
		pheno_mask2.resize(Nind,2);
	} 
	else{
		pheno.resize(Nind, phenocount);
		pheno_mask2.resize(Nind, phenocount);
	} 
	int i=0;  
	while(std::getline(ifs, line)){
		in.clear(); 
		in.str(line); 
		string temp;
		//fid,iid
		//todo: fid iid mapping; 
		in>>temp; in>>temp; 
		for(int j=0; j<phenocount;j++) {
			in>>temp;
			if(pheno_idx !="" && j== idx1)
			{
				if(temp=="NA")
				{	
					missing[0].push_back(i); 
					pheno_mask.push_back(0); 
					pheno_mask2(i,0)=0; 
					pheno(i,0)=0; 
				}
				else
				{
					double cur = atof(temp.c_str()); 
					pheno(i,0)=cur; 
					pheno_mask.push_back(1); 
					pheno_mask2(i,0)=1;
					pheno_sum[0]=pheno_sum[j]+cur;  
				}
			}
			if(pheno_idx != "" && j == idx2)
			{
				if(temp=="NA")
                                {
                                        missing[1].push_back(i);
                                        pheno_mask.push_back(0);
                                        pheno_mask2(i,1)=0;
                                        pheno(i,1)=0;
                                }
				else
				{
					double cur = atof(temp.c_str()); 
					pheno(i, 1) = cur; 
					pheno_mask.push_back(1); 
					pheno_mask2(i,1) =1; 
					pheno_sum[1] = pheno_sum[j] + cur; 
				}


			}
			if(pheno_idx ==""){
				if(temp=="NA")
				{
					missing[j].push_back(i); 
					pheno_mask.push_back(0); 
					pheno_mask2(i,j)=0;
					pheno(i,j)=0; 
				} 
				else{
					double cur= atof(temp.c_str()); 
					pheno(i,j)=cur; 
					pheno_sum[j] = pheno_sum[j]+cur; 
					pheno_mask2(i,j)=1; 
				}
			}
		}
		i++;
	}
	//not performing phenotype imputation
	if(! pheno_fill){
		if(pheno_idx!="")
			return 2; 
		return phenocount;
	} 
	//fill missing with mean
	if (pheno_idx!="") 
		phenocount=2; 	
	for(int a=0; a<phenocount; a++)
	{
		int missing_num= missing[a].size(); 
		double	pheno_avg = pheno_sum[a]/(Nind- missing_num); 
		//cout<<"pheno "<<a<<" avg: "<<pheno_avg<<endl; 
		for(int b=0 ; b<missing_num; b++)
		{
			int index = missing[a][b]; 
			pheno(index, a)= pheno_avg; 
			pheno_mask2(index,a)=1;
		}
	}
	return phenocount; 
}

void multiply_y_pre_fast(MatrixXdr &op, int Ncol_op ,MatrixXdr &res,bool subtract_means, int exist_ind,int phenoindex){
	for(int k_iter=0;k_iter<Ncol_op;k_iter++){
		sum_op[k_iter]=op.col(k_iter).sum();		
	}

			//cout << "Nops = " << Ncol_op << "\t" <<g.Nsegments_hori << endl;
	double vg, ve; 
	#if DEBUG==1
		if(debug){
			print_time (); 
			cout <<"Starting mailman on premultiply"<<endl;
			cout << "Nops = " << Ncol_op << "\t" <<g.Nsegments_hori << endl;
			cout << "Segment size = " << g.segment_size_hori << endl;
			cout << "Matrix size = " <<g.segment_size_hori<<"\t" <<g.Nindv << endl;
			cout << "op = " <<  op.rows () << "\t" << op.cols () << endl;
		}
	#endif


	//TODO: Memory Effecient SSE FastMultipy

	for(int seg_iter=0;seg_iter<g.Nsegments_hori-1;seg_iter++){
		mailman::fastmultiply(g.segment_size_hori,g.Nindv,Ncol_op,g.p[seg_iter],op,yint_m,partialsums,y_m);
		int p_base = seg_iter*g.segment_size_hori; 
		for(int p_iter=p_base; (p_iter<p_base+g.segment_size_hori) && (p_iter<g.Nsnp) ; p_iter++ ){
			for(int k_iter=0;k_iter<Ncol_op;k_iter++) 
				res(p_iter,k_iter) = y_m[p_iter-p_base][k_iter];
		}
	}

	int last_seg_size = (g.Nsnp%g.segment_size_hori !=0 ) ? g.Nsnp%g.segment_size_hori : g.segment_size_hori;
	mailman::fastmultiply(last_seg_size,g.Nindv,Ncol_op,g.p[g.Nsegments_hori-1],op,yint_m,partialsums,y_m);		
	int p_base = (g.Nsegments_hori-1)*g.segment_size_hori;
	for(int p_iter=p_base; (p_iter<p_base+g.segment_size_hori) && (p_iter<g.Nsnp) ; p_iter++){
		for(int k_iter=0;k_iter<Ncol_op;k_iter++) 
			res(p_iter,k_iter) = y_m[p_iter-p_base][k_iter];
	}

	#if DEBUG==1
		if(debug){
			print_time (); 
			cout <<"Ending mailman on premultiply"<<endl;
		}
	#endif


	if(!subtract_means)
		return;

	for(int p_iter=0;p_iter<p;p_iter++){
 		for(int k_iter=0;k_iter<Ncol_op;k_iter++){		 
			res(p_iter,k_iter) = res(p_iter,k_iter) - (g.get_col_mean(p_iter, phenoindex)*sum_op[k_iter]);
			if(var_normalize)
				res(p_iter,k_iter) = res(p_iter,k_iter)/(g.get_col_std(p_iter,phenoindex,exist_ind));		
 		}		
 	}	

}

void multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res,bool subtract_means, int exist_ind, int phenoindex){

	MatrixXdr op;
	op = op_orig.transpose();

	if(var_normalize && subtract_means){
		for(int p_iter=0;p_iter<p;p_iter++){
			for(int k_iter=0;k_iter<Nrows_op;k_iter++)		
				op(p_iter,k_iter) = op(p_iter,k_iter) / (g.get_col_std(p_iter,phenoindex,exist_ind));		
		}		
	}

	#if DEBUG==1
		if(debug){
			print_time (); 
			cout <<"Starting mailman on postmultiply"<<endl;
		}
	#endif
	
	int Ncol_op = Nrows_op;

	//cout << "ncol_op = " << Ncol_op << endl;

	int seg_iter;
	for(seg_iter=0;seg_iter<g.Nsegments_hori-1;seg_iter++){
		mailman::fastmultiply_pre(g.segment_size_hori,g.Nindv,Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter],op,yint_e,partialsums,y_e);
	}
	int last_seg_size = (g.Nsnp%g.segment_size_hori !=0 ) ? g.Nsnp%g.segment_size_hori : g.segment_size_hori;
	mailman::fastmultiply_pre(last_seg_size,g.Nindv,Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter],op,yint_e,partialsums,y_e);

	for(int n_iter=0; n_iter<n; n_iter++)  {
		for(int k_iter=0;k_iter<Ncol_op;k_iter++) {
			res(k_iter,n_iter) = y_e[n_iter][k_iter];
			y_e[n_iter][k_iter] = 0;
		}
	}
	
	#if DEBUG==1
		if(debug){
			print_time (); 
			cout <<"Ending mailman on postmultiply"<<endl;
		}
	#endif


	if(!subtract_means)
		return;

	double *sums_elements = new double[Ncol_op];
 	memset (sums_elements, 0, Nrows_op * sizeof(int));

 	for(int k_iter=0;k_iter<Ncol_op;k_iter++){		
 		double sum_to_calc=0.0;		
 		for(int p_iter=0;p_iter<p;p_iter++)		
 			sum_to_calc += g.get_col_mean(p_iter,phenoindex)*op(p_iter,k_iter);		
 		sums_elements[k_iter] = sum_to_calc;		
 	}		
 	for(int k_iter=0;k_iter<Ncol_op;k_iter++){		
 		for(int n_iter=0;n_iter<n;n_iter++)		
 			res(k_iter,n_iter) = res(k_iter,n_iter) - sums_elements[k_iter];		
 	}


}

void multiply_y_post_fast_partition(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res, bool subtract_means, int exist_ind, int phenoindex, int cat_num, int snp_idx)
{

   MatrixXdr op;
        op = op_orig.transpose();

        if(var_normalize && subtract_means){
                for(int p_iter=0;p_iter<p;p_iter++){
                        for(int k_iter=0;k_iter<Nrows_op;k_iter++)
                                op(p_iter,k_iter) = op(p_iter,k_iter) / (g.get_col_std(p_iter,phenoindex,exist_ind));
                }
        }

	int Ncol_op = Nrows_op;

         int seg_iter;

        for(seg_iter=0;seg_iter<g.Nsegments_hori-1;seg_iter++){
                mailman::fastmultiply_pre_normal_partition(g.segment_size_hori,g.Nindv,Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter],op, yint_e2,partialsums,y_e2, FunCat, cat_num, snp_idx);
        }
        int last_seg_size = (g.Nsnp%g.segment_size_hori !=0 ) ? g.Nsnp%g.segment_size_hori : g.segment_size_hori;
        mailman::fastmultiply_pre_normal_partition(last_seg_size,g.Nindv,Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter],op,yint_e2,partialsums,y_e2,FunCat, cat_num, snp_idx);
	//set snp_idx to 0 to start from beginning

	 for(int cat_k=0; cat_k < cat_num; cat_k++){
                for(int n_iter=0; n_iter<n; n_iter++)  {
                        for(int k_iter=0;k_iter<Ncol_op;k_iter++) {
			       int it = cat_k*n + n_iter;
                                res(k_iter,it) = y_e2[it][k_iter];
                                y_e2[it][k_iter] = 0;
                        }
                }
        }

}


void multiply_y_pre_naive_mem(MatrixXdr &op, int Ncol_op ,MatrixXdr &res, int exist_ind, int phenoindex){
	for(int p_iter=0;p_iter<p;p_iter++){
		for(int k_iter=0;k_iter<Ncol_op;k_iter++){
			double temp=0;
			for(int n_iter=0;n_iter<n;n_iter++)
				temp+= g.get_geno(p_iter,n_iter,var_normalize,phenoindex, exist_ind)*op(n_iter,k_iter);
			res(p_iter,k_iter)=temp;
		}
	}
}

void multiply_y_post_naive_mem(MatrixXdr &op, int Nrows_op ,MatrixXdr &res,int exist_ind, int phenoindex){
	for(int n_iter=0;n_iter<n;n_iter++){
		for(int k_iter=0;k_iter<Nrows_op;k_iter++){
			double temp=0;
			for(int p_iter=0;p_iter<p;p_iter++)
				temp+= op(k_iter,p_iter)*(g.get_geno(p_iter,n_iter,var_normalize,phenoindex, exist_ind));
			res(k_iter,n_iter)=temp;
		}
	}
}

void multiply_y_pre_naive(MatrixXdr &op, int Ncol_op ,MatrixXdr &res){
	res = geno_matrix * op;
}

void multiply_y_post_naive(MatrixXdr &op, int Nrows_op ,MatrixXdr &res){
	res = op * geno_matrix;
}

void multiply_y_post(MatrixXdr &op, int Nrows_op ,MatrixXdr &res,bool subtract_means,int exist_ind, int phenoindex){
    if(fast_mode)
        multiply_y_post_fast(op,Nrows_op,res,subtract_means,exist_ind,phenoindex);
    else{
		if(memory_efficient)
			multiply_y_post_naive_mem(op,Nrows_op,res,exist_ind,phenoindex);
		else
			multiply_y_post_naive(op,Nrows_op,res);
	}
}

void multiply_y_pre(MatrixXdr &op, int Ncol_op ,MatrixXdr &res,bool subtract_means,int exist_ind, int phenoindex){
    if(fast_mode)
        multiply_y_pre_fast(op,Ncol_op,res,subtract_means, exist_ind, phenoindex);
    else{
		if(memory_efficient)
			multiply_y_pre_naive_mem(op,Ncol_op,res,exist_ind, phenoindex);
		else
			multiply_y_pre_naive(op,Ncol_op,res);
	}
}
pair<double,double> weightedjack(vector<double> &t, vector<double> &m, double theta){
	if(t.size() != m.size())
	{
		cerr<<"In functions::weightedjack, Mismatch in length of t and m" <<endl; 
		exit(1); 
	}
	int g=0; 
	//double n=vector::sum(m); 
	double n=0; 
	for(int i=0; i<m.size(); i++)
		n += m[i]; 
	vector <double> res(m.size(), 0); 
	double a=0; 
	for(int i=0; i<m.size(); i++)
	{
		if(m[i]<=0)
			continue; 
		g++; 
		double h = n/m[i]; 
		a += (1-1/h)*t[i]; 
		double r = theta * h - (h-1) * t[i]; 
		res[i] =r; 
	}
	if(g==0){
		cerr<<"In function::weightedjack. Nmber of bolcks == 0"<<endl; 
		exit(1); 
	}
	double tj = theta * g; 
	tj -=a ; 
	double sj =0; 
	for(int i=0; i<m.size(); i++)
	{	
		if(m[i]<=0)
			continue; 
		double h=n/m[i]; 
		sj += pow((res[i]-tj), 2)/(h-1); 
	}
	sj /= g; 
	sj = pow(sj, 0.5); 
	
	return pair<double, double> (tj,sj); 


}
double compute_jack_knife(int j, int k, double rg, int SNP_max)
{
	MatrixXdr Xy1 = Xy.block(0, j, SNP_max, 1); 
	MatrixXdr Xy2 = Xy.block(0, k, SNP_max, 1);
	MatrixXdr y1  = pheno.block(0, j, g.Nindv, 1); 
	MatrixXdr y2 = pheno.block(0,k, g.Nindv, 1);  
	

	int block_per_chrom=5;
	vector<double> jack_knife_rg(22*block_per_chrom,1); 
	int snp_start=0;
	for(int i=0; i<22; i++)
	{
		int chrom_snp_num=g.get_chrom_snp(i); 
		int block_len= chrom_snp_num/block_per_chrom; 
		for(int t=0; t< block_per_chrom; t++)
		{
			int len = block_len; 
			if(t*block_len+block_len > chrom_snp_num)
				len = chrom_snp_num - t*block_len; 
			MatrixXdr Xy1_cur  = Xy1.block(snp_start+block_len*t, 0, len, 1); 
			MatrixXdr Xy2_cur = Xy2.block(snp_start+block_len*t, 0, len,1); 

		
			double X =yKy(j,k)*SNP_max/(SNP_max-len) - yy(j,k);  
			X = X- (Xy1_cur.transpose()*Xy2_cur/(SNP_max-len)).sum(); 
			
			double  Y = yKy(j,j)*SNP_max/(SNP_max-len)- yy(j,j);  
			Y = Y-  (Xy1_cur.transpose()* Xy1_cur/(SNP_max-len)).sum(); 
			Y = sqrt(Y); 
			double Z = yKy(k,k)*SNP_max/(SNP_max-len)- yy(k,k); 
			Z = Z- (Xy2_cur.transpose()*Xy2_cur/(SNP_max-len)).sum(); 
			Z = sqrt(Z); 
			jack_knife_rg[i*block_per_chrom+t] = X / Y/Z;
		}
		snp_start += chrom_snp_num;  
		/*if(i<10){
			cout<<jack_knife_rg(i,0)<<endl; 	
			cout<<"X: " << X <<endl <<"Y: "<<Y <<endl <<"Z: "<<Z <<endl; 
		}*/
	}
	pair<double, double> jack_output; 
	vector<double> jack_weight(22*block_per_chrom,1); 
	//weighted jacknife SE
	jack_output = weightedjack(jack_knife_rg, jack_weight, rg); 

	cout<<"Weighted jackknife SE: "<<jack_output.second<<endl; 
	//unweighted
	double result = 0; 
	for(int i=0; i<22*block_per_chrom; i++)
	{	
		double temp = jack_knife_rg[i]-rg; 
		result = result + temp*temp; 
	}
	//cout<<jack_knife_rg.block(0,0,10,1); 
	result = result *(22*block_per_chrom-1) / (22*block_per_chrom); 
	cout<<"Unweighted jackknife SE: "<<sqrt(result) <<endl; 
	return sqrt(result); 

}

void compute_se_coh(MatrixXdr &Xy1,MatrixXdr &Xy2, MatrixXdr &y1,MatrixXdr &y2,MatrixXdr &se, double h2g1,double h2g2, double h2e1,double h2e2,double tr_k2, int B , int exist_ind1, int exist_ind2, double rho_g, double rho_e, double h2g1_se, double h2g2_se, double h2g1_yKy, double h2g1_yy, double h2g2_yKy, double h2g2_yy)
{
	h2g1_se = h2g1_se*h2g1_se; 
	h2g2_se = h2g2_se* h2g2_se; 
//compute tr[y_1y_1^T(K-I)(h2g2K+h2e2I)(K-I)] 
//substude with tr[y_1y_1^T (K-I) y_2y_2^T (K-I)] 
	MatrixXdr zb=Xy1;
	MatrixXdr zb2 =Xy2; 
	//compute XXy1, XXy2
	 for(int j=0; j<p; j++){
                zb(j, 0)= zb(j,0)*stds(j,0);
		zb2(j,0) = zb2(j,0)*stds(j,0); 
	}
        MatrixXdr new_zb = zb.transpose();
	MatrixXdr new_zb2 = zb2.transpose();
        MatrixXdr new_res(1,n);
	MatrixXdr new_res2(1,n); 
        multiply_y_post_fast(new_zb, 1, new_res,false, exist_ind1,0);
	multiply_y_post_fast(new_zb2,1,new_res2, false, exist_ind2, 0); 
        MatrixXdr new_resid(1,p);
	MatrixXdr new_resid2(1,p); 
	MatrixXdr zb_scale_sum = new_zb*means;
	MatrixXdr zb_scale_sum2 = new_zb2 *means; 
        new_resid= zb_scale_sum* MatrixXdr::Constant(1,n, 1);
	new_resid2 = zb_scale_sum2*MatrixXdr::Constant(1,n,1); 
        MatrixXdr alpha1 = (new_res-new_resid).transpose() ;
	MatrixXdr alpha2 = (new_res2 -new_resid2).transpose(); 
        //compute Ky1, ky2
         for(int j=0; j< n; j++){
                alpha1(j,0)= alpha1(j,0)/p ;
		alpha2(j,0) = alpha2(j,0)/p; 
	}
	//alpha1 = (K-I)y1, alpha2 = (K-I)y2
	alpha1 = alpha1-y1;
	alpha2 = alpha2-y2;  
	MatrixXdr var_A = y2.transpose() * alpha1; 
	var_A(0,0) = var_A(0,0)*var_A(0,0); 
	//var_A(0,0) = tr[y1y1^T (K-I) y2y2^T(K-I)] 
	//compute tr[y1 y2^T(K-I)(rg*K+re*I)(K-I)] = rg(alpha2^T XX^T alpha1)/M + alpha2^Talpha1 re
	//compute X^Talpha1, X^T alpha2
	MatrixXdr res1(p,1); 
	MatrixXdr res2(p,1); 
	MatrixXdr resid1(p,1);
	MatrixXdr resid2(p,1); 
	MatrixXdr inter = means.cwiseProduct(stds); 
	multiply_y_pre_fast(alpha1, 1, res1, false, exist_ind1,0); 
	multiply_y_pre_fast(alpha2, 1, res2, false, exist_ind2,0); 
	for(int j=0; j<p; j++)
	{
		res1(j,0)= res1(j,0)*stds(j,0); 
		res2(j,0) = res2(j,0)*stds(j,0); 
	}
	resid1 = inter* alpha1.sum(); 
	resid2 = inter* alpha2.sum(); 
	MatrixXdr Xalpha1(p,1); 
	MatrixXdr Xalpha2(p,1); 
	Xalpha1 = res1 -resid1; 
	Xalpha2 = res2 -resid2; 
	//tr[y1 y2^T(K-I)(rg*K+re*I)(K-I)] = result
	double yKy1 = (Xalpha1.array()*Xalpha1.array()).sum()/p; 
	double temp1 = (alpha1.array()*alpha1.array()).sum();  
	double yKy = (Xalpha1.array()*Xalpha2.array()).sum()/p; 
	double temp = (alpha1.array()*alpha2.array()).sum(); 
	double result = yKy*rho_g + temp*rho_e + yKy1*h2g2+temp1*h2e2; 
	var_A(0,0)=0; 
	double var_X =  result; 
	double cor_XY = 2*(yKy* h2g1+ temp*h2e1) ; 
	double cor_XZ = 2*(yKy*h2g2 + temp*h2g2);
	double cor_YZ = 2*(yKy*rho_g + temp*rho_e);  

	double var_Y2 = 2*(yKy1*h2g1+ temp1*h2g1); 
	double yKy2 = (Xalpha2.array()*Xalpha2.array()).sum()/p; 
	double temp2 = (alpha2.array()*alpha2.array()).sum(); 
	double var_Z2 = 2*(yKy2*h2g2 + temp2*h2g2); 
	double E_X = tr_k2 *rho_g - exist_ind1 *rho_g; 
	double E_Y2 = tr_k2 *h2g1 - exist_ind1 *h2g1; 
	double E_Z2 = tr_k2 *h2g2 -exist_ind1*h2g2; 
	
	double E_Y = sqrt(E_Y2) - var_Y2/8/E_Y2/(sqrt(E_Y2)); 
	double E_Z = sqrt(E_Z2) - var_Z2 /8/E_Z2/(sqrt(E_Z2)); 
	double var_Y = var_Y2 / 4 / E_Y2; 
	double var_Z = var_Z2 / 4 / E_Z2; 
 
	double final_result = var_X / E_Y/E_Y /E_Z/E_Z  + var_Y*E_X*E_X / E_Y/E_Y/E_Y/E_Y/E_Z/E_Z + var_Z*E_X*E_X /E_Y /E_Y/E_Z/E_Z/E_Z/E_Z ; 

	se(0,0)=final_result;	
 
	//var_A(0,0) = var_A(0,0) + result; 
	//cout <<"var_A: "<<var_A(0,0)<<endl; 
	//double mu_A = rho_g * tr_k2 - exist_ind1* rho_g; 
	//cout<<"mu_A: "<<mu_A<<endl; 
	//double mu_B = tr_k2 - exist_ind1;
	//cout<<"mu_B: "<<mu_B<<endl;  
	//double var_B = tr_k2 / B /10;
	//cout<<"var B: " <<var_B<<endl ; 
	//double var_rg = var_A(0,0)/mu_B/mu_B + mu_A*mu_A*var_B/mu_B/mu_B/mu_B/mu_B; 
	//cout<<"var_rg: "<<var_rg<<endl; 
	//double factor = tr_k2 -2*exist_ind1 + exist_ind1*exist_ind1 / tr_k2; 
	//double E_rg = rho_g + rho_g/B/10/factor; 
	//cout<<"E_rg: " <<E_rg<<endl; 
	//double E_h2g1_2 = h2g1 + h2g1/B/10/factor;
	//cout<<"E_h2g1^2: "<<E_h2g1_2<<endl;  
	//double E_h2g2_2 = h2g2 + h2g2/B/10/factor;
	//cout<<"E_h2g2^2: "<<E_h2g2_2<<endl;  

	//double E_h2g1 = sqrt(E_h2g1_2) - h2g1_se / E_h2g1_2 / sqrt(E_h2g1_2) / 8; 
	//cout<<"E_h2g1: "<<E_h2g1<<endl; 
	//double E_h2g2 = sqrt(E_h2g2_2) - h2g2_se / E_h2g2_2 / sqrt(E_h2g2_2) / 8 ; 
	//cout<<"E_h2g2: "<<E_h2g2<<endl; 
	//double var_h2g1 = h2g1_se  / 4 / E_h2g1_2;
	//cout<<"var(h2g1): "<<var_h2g1<<endl;  
	//double var_h2g2 = h2g2_se / 4 / E_h2g2_2;
	//cout<<"var(h2g2): "<<var_h2g2<<endl;   
	//double final_result = var_rg / E_h2g1/E_h2g1/E_h2g2/E_h2g2 + E_rg* E_rg * var_h2g1 / E_h2g1/E_h2g1/E_h2g1/E_h2g1 /E_h2g2 /E_h2g2 + E_rg *E_rg * var_h2g2 /E_h2g1/E_h2g1/E_h2g2/E_h2g2/E_h2g2/E_h2g2; 
	//se(0,0)=final_result;  	
}

//previous version; computing se without covariate
void compute_se1(MatrixXdr &Xy,  MatrixXdr &y,MatrixXdr &se, double h2g, double h2e,double tr_k2, int B , int exist_ind)
{
	//compute tr[yy^T(K-I)(h2gK+h2eI)(K-I)]
	//compute X^T y 
	//imput X^y[i] p*1 vector
	cout<<"p: "<<p << "  n: "<<n<<endl;
        MatrixXdr zb =Xy;
	//compute XXy
	for(int j=0; j<p; j++)
                zb(j, 0)= zb(j,0)*stds(j,0);
        MatrixXdr new_zb = zb.transpose();
        MatrixXdr new_res(1,n);
        multiply_y_post_fast(new_zb, 1, new_res,false, exist_ind,0);
        MatrixXdr new_resid(1,p);
        MatrixXdr zb_scale_sum = new_zb*means;
        new_resid= zb_scale_sum* MatrixXdr::Constant(1,n, 1);
        MatrixXdr alpha = (new_res-new_resid).transpose() ;
	//compute Ky
	 for(int j=0; j< n; j++)
                alpha(j,0)= alpha(j,0)/p ;
	//alpha =(K-I)y
	 alpha = alpha -y;
        MatrixXdr res(p,1);
        MatrixXdr resid(p,1);
        MatrixXdr inter = means.cwiseProduct(stds);
        multiply_y_pre_fast(alpha, 1, res, false, exist_ind,0);
        for(int j=0; j<p;j++)
                res(j,0)=res(j,0)*stds(j,0);
        inter = means.cwiseProduct(stds);
        resid = inter * alpha.sum();
        MatrixXdr Xalpha(p,1);
        Xalpha = res-resid;
	//Xy =res; 
	double yKy = (Xalpha.array()*Xalpha.array()).sum() / p;
        double temp =(alpha.array()*alpha.array()).sum();
        double result = yKy*h2g+temp*h2e;
        result = 2*result + h2g*h2g*tr_k2/10/B;
	result = sqrt(result) / (tr_k2-n);
        cout<<result<<endl; 
   	MatrixXdr result1(1,1);
        result1(0,0)=result;se=result1;
}

void compute_se(MatrixXdr &Xy, MatrixXdr &y,MatrixXdr &se, double h2g, double h2e,double tr_k2, int B , int exist_ind, double tr_k, int cov_num)
{
	//mu_a and mu_B are for covarate use to compute scalar
	double mu_A = h2g * (tr_k2 - tr_k * tr_k); 
	double mu_B = tr_k2 - tr_k*tr_k / (exist_ind - cov_num); 
	//compute X^T y
	//input X^y[i] p*1 vector
	cout<<"p: "<<p << "  n: "<<n<<endl; 
	MatrixXdr zb=Xy;
	//compute XXy
	for(int j=0; j<p; j++)
		zb(j, 0)= zb(j,0)*stds(j,0); 
	MatrixXdr new_zb = zb.transpose(); 
	MatrixXdr new_res(1,n); 
	multiply_y_post_fast(new_zb, 1, new_res,false, exist_ind,0); 
	MatrixXdr new_resid(1,p); 
	MatrixXdr _scale_sum = new_zb*means; 
	new_resid= _scale_sum* MatrixXdr::Constant(1,n, 1); 
	MatrixXdr alpha = (new_res-new_resid).transpose() ;
	//compute Ky
	for(int j=0; j< n; j++)
		alpha(j,0)= alpha(j,0)/p ;
	//alpha = (K-I)y 
	//if no covariate, tr_k / (exist_ind - cov_num) = 1
	if(cov_num==0)
		alpha = alpha -((tr_k)/(exist_ind-cov_num))*y; 
	if(cov_num != 0){
		alpha = alpha - covariate * WW * (covariate.transpose()*alpha);
		alpha = alpha - ((tr_k)/(exist_ind-cov_num))*y; 
		MatrixXdr alpha_prime = covariate.transpose() * alpha; 
		alpha = alpha - covariate * WW * alpha_prime; 
	}
	MatrixXdr res(p,1); 
	MatrixXdr resid(p,1); 
	MatrixXdr inter = means.cwiseProduct(stds); 
	multiply_y_pre_fast(alpha, 1, res, false, exist_ind,0); 
	for(int j=0; j<p;j++)
		res(j,0)=res(j,0)*stds(j,0); 
	inter = means.cwiseProduct(stds); 
	resid = inter * alpha.sum(); 
	MatrixXdr Xalpha(p,1); 
	Xalpha = res-resid;
	//Xy =res;
	double yKy = (Xalpha.array()*Xalpha.array()).sum() / p; 
	double temp =(alpha.array()*alpha.array()).sum();  
	double result = yKy*h2g+temp*h2e; 
	if(cov_num==0)
		result = 2*result + h2g*h2g*tr_k2/10/B; 
	else	
		result = 2*result +tr_k2* mu_A*mu_A / mu_B/mu_B/10/B; 
	//cout<<result; 
	if(cov_num==0)
		result = sqrt(result) / (tr_k2-n);
	else 
		result = sqrt(result) / mu_B;   
	MatrixXdr result1(1,1); 
	result1(0,0)=result;se=result1;  
}
//compute_b2 for genetic correlation
void compute_b2(bool use_cov,  double exist_ind, int pheno_i, int pheno_j, MatrixXdr& pheno_prime1, MatrixXdr& pheno_prime2)
{
	int pheno_num=1; 
	MatrixXdr y1 = pheno.block(0,pheno_i, g.Nindv, 1); 
	MatrixXdr y2 = pheno.block(0,pheno_j, g.Nindv, 1); 
	MatrixXdr mask1 = pheno_mask2.block(0, pheno_i, g.Nindv, 1); 
	MatrixXdr mask2 = pheno_mask2.block(0, pheno_j, g.Nindv, 1); 
	MatrixXdr y_sum1 = y1.colwise().sum(); 
	MatrixXdr y_sum2 = y2.colwise().sum(); 
	if(!use_cov)
	{
		 MatrixXdr res1(g.Nsnp, pheno_num);
		MatrixXdr res2(g.Nsnp, pheno_num); 
                multiply_y_pre(y1,pheno_num,res1,false, exist_ind,0);
		multiply_y_pre(y2, pheno_num, res2, false, exist_ind,0); 
                for(int i=0; i<pheno_num; i++){
                        MatrixXdr cur= res1.block(0,i,g.Nsnp, 1);
                        res1.block(0,i,g.Nsnp, 1)  = cur.cwiseProduct(stds);
                	cur  = res2.block(0,i,g.Nsnp, 1); 
			res2.block(0,i,g.Nsnp, 1) = cur.cwiseProduct(stds) ;
		}
		MatrixXdr resid1(g.Nsnp, pheno_num);
		MatrixXdr resid2(g.Nsnp, pheno_num); 
                for(int i=0; i<pheno_num; i++)
                {
                        resid1.block(0,i,g.Nsnp, 1) = means.cwiseProduct(stds)*y_sum1(0,i);
			resid2.block(0,i,g.Nsnp,1 ) =means.cwiseProduct(stds)*y_sum2(0,i); 
                }
		MatrixXdr Xy1 = res1-resid1; 
		MatrixXdr Xy2 = res2-resid2; 		
		yKy = Xy1.transpose() * Xy2; 
		yKy = yKy /g.Nsnp; 
	}
	if(use_cov)
	{
		MatrixXdr y_temp1 = y1 - covariate * WW *pheno_prime1; 
		MatrixXdr y_temp2 = y2-covariate*WW* pheno_prime2; 
		y_temp1 = y_temp1.cwiseProduct(mask1); 
		y_temp2 = y_temp2.cwiseProduct(mask2); 
		y_sum1 = y_temp1.colwise().sum(); 
		y_sum2 = y_temp2.colwise().sum(); 
			
		
		MatrixXdr res1(g.Nsnp, pheno_num); 
		MatrixXdr res2(g.Nsnp, pheno_num); 
		multiply_y_pre(y_temp1, pheno_num, res1, false, exist_ind,0); 
		multiply_y_pre(y_temp2, pheno_num, res2, false, exist_ind,0); 

		for(int i=0; i<pheno_num; i++){
                        MatrixXdr cur= res1.block(0,i,g.Nsnp, 1);
                        res1.block(0,i,g.Nsnp, 1)  = cur.cwiseProduct(stds);
                        cur  = res2.block(0,i,g.Nsnp, 1);
                        res2.block(0,i,g.Nsnp, 1) = cur.cwiseProduct(stds) ;
                }
                MatrixXdr resid1(g.Nsnp, pheno_num);
                MatrixXdr resid2(g.Nsnp, pheno_num);
                for(int i=0; i<pheno_num; i++)
                {
                        resid1.block(0,i,g.Nsnp, 1) = means.cwiseProduct(stds)*y_sum1(0,i);
                        resid2.block(0,i,g.Nsnp,1 ) =means.cwiseProduct(stds)*y_sum2(0,i);
                }
		MatrixXdr Xy1 = res1-resid1;
                MatrixXdr Xy2 = res2-resid2;
                yKy = Xy1.transpose() * Xy2;
                yKy = yKy /g.Nsnp;

	}
	yy= y1.transpose()* y2; 
	if(use_cov)
		yy = yy- pheno_prime1.transpose() * WW * pheno_prime2; 


}


//compute_b1 for heritability
void compute_b1 (bool use_cov, MatrixXdr& y_sum, double exist_ind, MatrixXdr& pheno_prime_cur, bool pheno_fill, int pheno_num, int snp_idx, int cov_num){
		//	double exist_ind = exist_ind_mx(0,0);
		/*	MatrixXdr pheno_cur, mask_cur; 
			if (!pheno_fill) 
			{
			pheno_cur = pheno.block(0, pheno_i, g.Nindv, 1); 
			mask_cur = pheno_mask2.block(0,pheno_i, g.Nindv, 1); 
			}
			else 
			{
			pheno_cur = pheno; 
			mask_cur = pheno_mask2; 
			}
			*/
		for (int pheno_i=0; pheno_i<pheno_num; pheno_i+=10)
		{	
				//	MatrixXdr Xy_cur; 

				int pheno_num_cur= 10;
				if(pheno_i+10 >= pheno_num)
						pheno_num_cur = pheno_num-pheno_i;
				MatrixXdr pheno_cur = pheno.block(0,pheno_i, g.Nindv,pheno_num_cur);
				MatrixXdr stds_cur = stds.block(0, pheno_i, g.Nsnp, pheno_num_cur); 
				MatrixXdr means_cur = means.block(0, pheno_i, g.Nsnp, pheno_num_cur); 

				if(!use_cov){

						MatrixXdr res(g.Nsnp, pheno_num_cur);
						multiply_y_pre(pheno_cur,pheno_num_cur,res,false, exist_ind,0);


						res = res.cwiseProduct(stds_cur);   
						MatrixXdr resid(g.Nsnp, pheno_num_cur);
						resid = means_cur.cwiseProduct(stds_cur); 
						for(int i=0; i<pheno_num_cur; i++)
						{
								resid.block(0,i,g.Nsnp, 1) = resid.block(0, i, g.Nsnp, 1)*y_sum(0,i+pheno_i);
						}
						//		Xy_cur= res-resid; 
						Xy.block(snp_idx, pheno_i, g.Nsnp, pheno_num_cur) = res-resid; 

				}
				if(use_cov)
				{
						MatrixXdr y_temp = pheno_cur; 


						for (int i=0; i<pheno_num_cur; i++) 
						{
								MatrixXdr covariate_cur = covariate; 
								MatrixXdr cur_mask = pheno_mask2.block(0, i, g.Nindv, 1) ; 
								for(int k=0; k<cov_num; k++)
										covariate_cur.block(0,k,g.Nindv, 1) = covariate_cur.block(0,k,g.Nindv,1).cwiseProduct(cur_mask);
								MatrixXdr WW_cur = WW.block(0,(pheno_i+i)*cov_num, cov_num, cov_num) ;
								y_temp.block(0, i, g.Nindv, 1) = y_temp.block(0, i, g.Nindv, 1) - covariate_cur * WW_cur * pheno_prime_cur.block(0, pheno_i+i, cov_num, 1);  

						}
						y_temp = y_temp.cwiseProduct(pheno_mask2.block(0, pheno_i, g.Nindv, pheno_num_cur)); 
						MatrixXdr y_temp_sum=y_temp.colwise().sum();

						MatrixXdr res(g.Nsnp, pheno_num_cur);
						multiply_y_pre(y_temp,pheno_num_cur,res,false,exist_ind,0);
						res= res.cwiseProduct(stds_cur); 
						/*	for(int i=0; i<pheno_num; i++){
							MatrixXdr cur= res.block(0,i,g.Nsnp, 1);
							res.block(0,i,g.Nsnp, 1)  = cur.cwiseProduct(stds);
							}
							*/
						MatrixXdr resid(g.Nsnp, pheno_num_cur);
						resid = means_cur.cwiseProduct(stds_cur); 
						for(int i=0; i<pheno_num_cur; i++)
						{
								resid.block(0,i,g.Nsnp, 1) = resid.block(0,i, g.Nsnp, 1)*y_temp_sum(0,i);
						}
						//	Xy_cur= res-resid; 
						Xy.block(snp_idx,pheno_i,g.Nsnp, pheno_num_cur) = res-resid;
				}
		}
}

//previous version: read only one phenotype, use vector for missing mask
//out dated
void rhe_reg(MatrixXdr & zb,MatrixXdr  &Xzb_total_cur, MatrixXdr &resid_total_cur,  int pheno_i, int exist_ind, int cov_num, int cat_num, int snp_idx, int block_omit, int pheno_num)
{
	MatrixXdr cur_mask = pheno_mask2.block(0, pheno_i, g.Nindv, 1);
		
		MatrixXdr cur_means  = means.block(0, pheno_i, g.Nsnp, 1); 
		MatrixXdr cur_stds = stds.block(0, pheno_i, g.Nsnp, 1); 
		MatrixXdr zb_cur = zb; 		
                for(int b=0; b<10; b++){
                        MatrixXdr temp = zb_cur.block(0,b,g.Nindv,1);
                        zb_cur.block(0,b,g.Nindv, 1) = temp.cwiseProduct(cur_mask);
                }
		MatrixXdr  covariate_cur = covariate; 
		MatrixXdr WW_cur = WW.block(0, pheno_i*cov_num, cov_num, cov_num); 
		if(use_cov)
		{
			for(int k=0; k<cov_num; k++) 
				 covariate_cur.block(0,k,g.Nindv, 1) = covariate_cur.block(0,k,g.Nindv,1).cwiseProduct(cur_mask); 
			MatrixXdr tremp = WW_cur*covariate_cur.transpose() * zb_cur; 
			tremp =  covariate_cur * tremp; 
			zb_cur  = zb_cur-tremp; 
		}

		MatrixXdr zb_sum = zb_cur.colwise().sum(); 
/*
		if(pheno_i==0){
			cout<<"zb stats: "<<endl<<"means: "<<zb_sum/exist_ind<<endl; 
			cout<<"variance:"<<endl; 
			MatrixXdr zb_stats = zb_cur.block(0, 0, g.Nindv, 1) ; 
			zb_stats = zb_stats - zb_sum(0,0)/exist_ind *MatrixXdr::Constant(g.Nindv, 1, 1); 
			cout << zb_stats.cwiseProduct(zb_stats).sum()/exist_ind<<endl; 
		}*/
                MatrixXdr res(g.Nsnp, 10);
                multiply_y_pre(zb_cur,10,res, false, exist_ind,pheno_i);
	//	if(pheno_i==0)
	//		cout<<"FUNC rhereg res: "<<endl<<res.block(0,0,4,2)<<endl;      

		for(int j=0; j<g.Nsnp; j++)
                        for(int k=0; k<10;k++){
                             res(j,k) = res(j,k)*cur_stds(j,0);
		}
	//	if(pheno_i==0)
	//		cout<<"FUNC rhereg scaled res: "<<endl<<res.block(0,0,4,2); 
	
		MatrixXdr resid(g.Nsnp, 10);
                MatrixXdr inter = cur_means.cwiseProduct(cur_stds);
                resid = inter * zb_sum;
                MatrixXdr zb1(g.Nindv,10);
		zb1 = res-resid; 
		for(int k=0; k<10; k++){
                  for(int j=0; j<g.Nsnp;j++){
                        zb1(j,k) =zb1(j,k) *cur_stds(j,0);}}
		MatrixXdr new_zb = zb1.transpose();

                MatrixXdr new_res_mc(10, g.Nindv*cat_num);
                multiply_y_post_fast_partition(new_zb, 10, new_res_mc, false, exist_ind,pheno_i,cat_num, snp_idx);
		for(int b_iter=0; b_iter<10; b_iter++){

                MatrixXdr new_res(cat_num, g.Nindv);
		for (int fun_iter=0; fun_iter<cat_num; fun_iter++){
                        new_res.block(fun_iter, 0, 1, g.Nindv) = new_res_mc.block(b_iter, fun_iter*g.Nindv, 1, g.Nindv);
                }
		/*
 		//if using matrix the order is different 
		for (int n_iter=0; n_iter<g.Nindv; n_iter++){
                        new_res.block(0, n_iter, cat_num, 1) = new_res_mc.block(b_iter, n_iter*cat_num, 1, cat_num).transpose();
                }
		*/
                MatrixXdr new_resid(cat_num, g.Nindv);

		MatrixXdr zb_scale_sum(cat_num, 1); 
		for (int fun_iter=0; fun_iter<cat_num; fun_iter++) 
		{
			MatrixXdr temp  =new_zb.block(b_iter, 0, 1, g.Nsnp).cwiseProduct(FunCat.block(fun_iter, snp_idx, 1, g.Nsnp)); 
			zb_scale_sum(fun_iter, 0) =(temp*cur_means).sum(); 
		}
                new_resid = zb_scale_sum * MatrixXdr::Constant(1,g.Nindv, 1);
		 MatrixXdr Xzb1 = new_res- new_resid;
		if(use_cov)
                {
                        MatrixXdr temp1 = WW_cur * covariate_cur.transpose() *Xzb1.transpose();
                        MatrixXdr temp = covariate_cur * temp1;
			MatrixXdr temp_resid = Xzb1.transpose() - temp;      
                       for(int b=0; b<cat_num; b++) 
			{
				MatrixXdr one_col = temp_resid.block(0,b, g.Nindv, 1);
				temp_resid.block(0, b, g.Nindv, 1) = one_col.cwiseProduct(cur_mask);  
			}
			resid_total_cur.block(0, pheno_i*(10*cat_num)+b_iter*cat_num, g.Nindv, cat_num) += temp_resid;

                        Xzb1 = Xzb1 - temp.transpose();
		}
		else{
		resid_total_cur.block(0, pheno_i*(10*cat_num)+b_iter*cat_num, g.Nindv, cat_num) += Xzb1.transpose();
		} 
	Xzb_total_cur.block(0, pheno_i*(10*cat_num)+b_iter*cat_num, g.Nindv, cat_num) += Xzb1.transpose();
	// assume cat_num for now; TODO: fix this for partition
	//for jackknife; disabled for now
	
		
	}
}

pair<double,double> get_error_norm(MatrixXdr &c, int exist_ind, int phenoindex){
	HouseholderQR<MatrixXdr> qr(c);
	MatrixXdr Q;
	Q = qr.householderQ() * MatrixXdr::Identity(p,k);
	MatrixXdr q_t(k,p);
	q_t = Q.transpose();
	MatrixXdr b(k,n);
	multiply_y_post(q_t,k,b,true,exist_ind,0);
	JacobiSVD<MatrixXdr> b_svd(b, ComputeThinU | ComputeThinV);
	MatrixXdr u_l,d_l,v_l; 
	if(fast_mode)
        u_l = b_svd.matrixU();
    else
        u_l = Q * b_svd.matrixU();
	v_l = b_svd.matrixV();
	d_l = MatrixXdr::Zero(k,k);
	for(int kk=0;kk<k; kk++)
		d_l(kk,kk) = (b_svd.singularValues())(kk);
	
	MatrixXdr u_k,v_k,d_k;
	u_k = u_l.leftCols(k_orig);
	v_k = v_l.leftCols(k_orig);
	d_k = MatrixXdr::Zero(k_orig,k_orig);
	for(int kk =0 ; kk < k_orig ; kk++)
		d_k(kk,kk)  =(b_svd.singularValues())(kk);

	MatrixXdr b_l,b_k;
    b_l = u_l * d_l * (v_l.transpose());
    b_k = u_k * d_k * (v_k.transpose());

    if(fast_mode){
        double temp_k = b_k.cwiseProduct(b).sum();
        double temp_l = b_l.cwiseProduct(b).sum();
        double b_knorm = b_k.norm();
        double b_lnorm = b_l.norm();
        double norm_k = (b_knorm*b_knorm) - (2*temp_k);
        double norm_l = (b_lnorm*b_lnorm) - (2*temp_l);	
        return make_pair(norm_k,norm_l);
    }
    else{
        MatrixXdr e_l(p,n);
        MatrixXdr e_k(p,n);
        for(int p_iter=0;p_iter<p;p_iter++){
            for(int n_iter=0;n_iter<n;n_iter++){
                e_l(p_iter,n_iter) = g.get_geno(p_iter,n_iter,var_normalize, phenoindex, exist_ind) - b_l(p_iter,n_iter);
                e_k(p_iter,n_iter) = g.get_geno(p_iter,n_iter,var_normalize,phenoindex, exist_ind) - b_k(p_iter,n_iter);
            }
        }

        double ek_norm = e_k.norm();
        double el_norm = e_l.norm();
        return make_pair(ek_norm,el_norm);
    }
}

void construct_linear_system(MatrixXdr &cur_pheno_Xzb1, MatrixXdr &cur_pheno_Xzb2, MatrixXdr &cur_pheno_resid1 ,MatrixXdr &cur_pheno_resid2,  MatrixXdr &A, MatrixXdr &b, MatrixXdr &zb, MatrixXdr &FunCat_snpcount, int cat_num, int i,int j, int SNP_max, MatrixXdr &pheno_mask1, MatrixXdr &pheno_mask2, MatrixXdr &FunCat){
	
	MatrixXdr common_mask = pheno_mask1.cwiseProduct(pheno_mask2);
	double N_cur = common_mask.sum(); 
	//i, j is the index of phenotype
	//if(use_cov){
	A.block(cat_num, 0, 1, cat_num) = MatrixXdr::Constant(1, cat_num, 0);
	A.block(0, cat_num, cat_num, 1) = MatrixXdr::Constant(cat_num, 1, 0);
   	//}
	for(int b_iter=0; b_iter<10; b_iter++)
	{
	MatrixXdr cur_Xzb_block1 = cur_pheno_Xzb1.block(0,b_iter*cat_num,g.Nindv, cat_num); 
	MatrixXdr cur_Xzb_block2 = cur_pheno_Xzb2.block(0, b_iter*cat_num, g.Nindv, cat_num); 
	MatrixXdr cur_resid_block = cur_pheno_resid1.block(0,b_iter*cat_num, g.Nindv, cat_num);
	// ===== for tr[KiKj] -----
	// cur_Xzb_block1 is N *  cat_num
	A.block(0,0,cat_num,cat_num) += (cur_Xzb_block1.transpose() * cur_Xzb_block2) / 10 ;

//	cout<<(cur_Xzb_block1.transpose() * cur_Xzb_block2) / 10 <<endl
	// ===== for residuals (tr[Ki])  ----  
	for(int fun_i=0; fun_i<cat_num; fun_i++)
	{
		cur_resid_block.block(0, fun_i, g.Nindv,1) = cur_resid_block.block(0, fun_i, g.Nindv,1).cwiseProduct(common_mask); 
		
	} 
	 MatrixXdr W = (zb.block(0, b_iter, g.Nindv, 1).transpose())* cur_resid_block;
	 W= W/10;
	 for(int fun_i=0; fun_i<cat_num; fun_i++)
                                W(0, fun_i) = W(0, fun_i) /FunCat_snpcount(fun_i, 0);
	
			A.block(cat_num, 0, 1, cat_num)  += W ; 
			A.block(0, cat_num, cat_num, 1) += W.transpose(); 
	}//end for loop for b_iter 
	
	for (int fun_i=0; fun_i < cat_num; fun_i ++)
            for (int fun_j=0; fun_j< cat_num; fun_j ++)
               A(fun_i, fun_j) = A(fun_i, fun_j) / FunCat_snpcount(fun_i, 0) / FunCat_snpcount(fun_j, 0);

	for(int fun_iter=0; fun_iter<cat_num; fun_iter++)
                        {
                                MatrixXdr cur1 = Xy.block(0, i, SNP_max,1).cwiseProduct(FunCat.block(fun_iter, 0, 1,SNP_max).transpose());
                               MatrixXdr cur2 = Xy.block(0, j, SNP_max, 1).cwiseProduct(FunCat.block(fun_iter, 0, 1, SNP_max).transpose()); 
				 b(fun_iter, 0) = (cur1.array()*cur2.array()).sum()/FunCat_snpcount(fun_iter, 0);
                        }
}

int main(int argc, char const *argv[]){

	//clock_t io_begin = clock();
    //clock_gettime (CLOCK_REALTIME, &t0);

	pair<double,double> prev_error = make_pair(0.0,0.0);
	double prevnll=0.0;

	parse_args(argc,argv);

	
	//TODO: Memory effecient Version of Mailman

	memory_efficient = command_line_opts.memory_efficient;
	text_version = command_line_opts.text_version;
	fast_mode = command_line_opts.fast_mode;
	missing = command_line_opts.missing;
	reg = command_line_opts.reg;
	noh2g = command_line_opts.noh2g;
	//gwas=false; 
	//gwas=command_line_opts.gwas; 
	//cout<<"perform GWAS: "<<gwas<<endl; 
	bpheno = command_line_opts.bpheno; 
	cout<<"Binary phenotype: "<<bpheno<<endl;
	pheno_fill = command_line_opts.pheno_fill; 
	cout<<"Filling phenotype with mean: " <<pheno_fill <<endl;  
	cout<<"Compute heritability: "<<!noh2g<<endl; 
	
        std::string output_file = command_line_opts.OUTPUT_PATH;
         cout<<"Output write to destination "<<output_file<<endl;
	FILE *fp_output;
        fp_output = fopen(output_file.c_str(), "w") ;
        fprintf(fp_output, "%s\t%s\n","Source","Value") ;
	//get number of individuals
	std::stringstream f2;
        f2 << command_line_opts.GENOTYPE_FILE_PATH << ".fam";
        g.read_fam (f2.str());
	//cout<<"ind: "<<g.Nindv<<endl; 	
	//get phenotype
	std::string  pheno_idx = command_line_opts.pheno_idx;
	std::string filename=command_line_opts.PHENOTYPE_FILE_PATH;
        int pheno_num= read_pheno2(g.Nindv, filename, pheno_idx,pheno_fill);
	// int exist_ind =0;
       // for(int i=0; i<g.Nindv; i++)
         //       exist_ind += pheno_mask[i];
	MatrixXdr exist_ind = pheno_mask2.colwise().sum(); 
	int cov_num=0 ;
        if(filename=="")
        {
                cout<<"No Phenotype File Specified"<<endl;
                return 0 ;
        }
        cout<< "Read in "<<pheno_num << " phenotypes"<<endl;
	cout<< "There are "<<exist_ind<< " individuals with no missing phenotypes"<<endl; 
	MatrixXdr VarComp(pheno_num,2);

	//if(gwas)
	//	fast_mode=false; //for now, gwas need the genotype matrix, and compute kinship constructed with one chrom leave out 
	if(!reg)
		fast_mode=false; //force save whole genome if non randomized  
	/*
	if(text_version){
		if(fast_mode)
			g.read_txt_mailman(command_line_opts.GENOTYPE_FILE_PATH,missing);
		else
			g.read_txt_naive(command_line_opts.GENOTYPE_FILE_PATH,missing);
	}
	else{
		g.read_plink(command_line_opts.GENOTYPE_FILE_PATH,missing,fast_mode, pheno_mask2, pheno_num);
	}
*/
	std::stringstream bimfile; 
	bimfile<<command_line_opts.GENOTYPE_FILE_PATH<<".bim"; 
	std::stringstream bedfile;
	g.read_bim(bimfile.str(), pheno_num);  
	bedfile<<command_line_opts.GENOTYPE_FILE_PATH<<".bed"; 
			
//g.read_bed(bedfile.str(), missing, fast_mode, pheno_mask2, pheno_num); 
	

	//TODO: Implement these codes.
	if(missing && !fast_mode){
		cout<<"Missing version works only with mailman i.e. fast mode\n EXITING..."<<endl;
		exit(-1);
	}
	if(fast_mode && memory_efficient){
		cout<<"Memory effecient version for mailman EM not yet implemented"<<endl;
		cout<<"Ignoring Memory effecient Flag"<<endl;
	}
	if(missing && var_normalize){
		cout<<"Missing version works only without variance normalization\n EXITING..."<<endl;
		exit(-1);
	}

	
        int cat_num=1;

        string annotation_file = command_line_opts.ANNOTATION_FILE_PATH;
        if (annotation_file!="")
                cat_num= read_annotation(g.Nsnp, annotation_file);
        else
                FunCat = MatrixXdr::Constant(1, g.Nsnp, 1);
        MatrixXdr FunCat_snpcount = FunCat.rowwise().sum();
        cout<<"Total number of annotation groups: "<<cat_num<<endl;

    //MAX_ITER =  command_line_opts.max_iterations ; 
	int B = command_line_opts.batchNum;
	cout<<"num of random vectors: "<<B<<endl;  
	k_orig = command_line_opts.num_of_evec ;
	debug = command_line_opts.debugmode ;
	float tr2= command_line_opts.tr2; 
	check_accuracy = command_line_opts.getaccuracy;
	var_normalize = true; 
	accelerated_em = command_line_opts.accelerated_em;
	k = k_orig + command_line_opts.l;
	k = (int)ceil(k/10.0)*10;
	command_line_opts.l = k - k_orig;
	p = g.Nsnp;
	n = g.Nindv;
	bool toStop=false;
		toStop=true;
	//fix seed for checking
	//srand((unsigned int) time(0));
	c.resize(p,k);
	x.resize(k,n);
	v.resize(p,k);
//	means.resize(p,1);
//	stds.resize(p,1);
//	sum2.resize(p,1); 
//	sum.resize(p,1); 


/*	//Ture off non-fast mode 
 	if(!fast_mode ){
		geno_matrix.resize(p,n);
		cout<<"geno resize"<<endl; 
		g.generate_eigen_geno(geno_matrix,true,0);
		cout<<geno_matrix.rows()<<endl; 
		cout<<geno_matrix.cols()<<endl;
}
*/	
		
	//clock_t io_end = clock();

	//TODO: Initialization of c with gaussian distribution
	c = MatrixXdr::Random(p,k);


	// Initial intermediate data structures
	blocksize = k;
	int hsegsize = floor(log(g.Nindv)/log(3))-2;
	hsegsize = hsegsize>0?hsegsize:1;  	// = log_3(n)
	int hsize = pow(3,hsegsize);		 
	int vsegsize = 1; 		// = log_3(p)
	int vsize = pow(3,vsegsize);		 
	//int SNP_BLOCK_SIZE = max(50000,int(log(g.Nindv)/log(3))); 
//	if(SNP_BLOCK_SIZE > g.Nsnp)
//		SNP_BLOCK_SIZE = g.Nsnp; 
	int SNP_BLOCK_SIZE; 

	int SNP_max = g.Nsnp; 
	string temp_bedfile = bedfile.str(); 
	ifstream ifs(temp_bedfile.c_str(), ios::in|ios::binary); 
	char magic[3]; 
	ifs.read(magic, 3*sizeof(char)); 

	ofstream c_file;
        if(debug){
                c_file.open((string(command_line_opts.OUTPUT_PATH)+string("cvals_orig.txt")).c_str());
                c_file<<c<<endl;
                c_file.close();
                printf("Read Matrix\n");
        }

       #if SSE_SUPPORT==1
                if(fast_mode)
                        cout<<"Using Optimized SSE FastMultiply"<<endl;
        #endif
		//center phenotypes 
	 	MatrixXdr y_sum=pheno.colwise().sum();
		MatrixXdr exist_ind_inv(1, pheno_num); 
		for(int k=0; k<pheno_num; k++)
			exist_ind_inv(0, k) = 1/exist_ind(0,k); 
		MatrixXdr prevelance = y_sum.cwiseProduct(exist_ind_inv);
                MatrixXdr y_mean = y_sum.cwiseProduct(exist_ind_inv);
                for(int i=0; i<g.Nindv; i++)
		{
                        MatrixXdr temp =(pheno.block(i,0,1,pheno_num) - y_mean); //center phenotype 
                	pheno.block(i,0,1,pheno_num) = temp.cwiseProduct(pheno_mask2.block(i,0,1,pheno_num)); 
		}
		//normalize
		for(int pheno_i=0; pheno_i<pheno_num; pheno_i++){
		MatrixXdr pheno_cur = pheno.block(0,pheno_i, g.Nindv, 1); 
		MatrixXdr pheno_sum2 = pheno_cur.transpose()*pheno_cur; 
		double pheno_variance = pheno_sum2(0,0)/ (exist_ind(0,pheno_i)-1);
		pheno_variance = sqrt(pheno_variance);  
		pheno.block(0,pheno_i, g.Nindv, 1)= pheno_cur / pheno_variance; 


}
		y_sum=pheno.colwise().sum();
	
		//read in covariate 	
		std::string covfile=command_line_opts.COVARIATE_FILE_PATH;
                std::string covname=command_line_opts.COVARIATE_NAME;
                if(covfile!=""){
           	     use_cov=true;
                	cov_num=read_cov(true,g.Nindv, covfile, covname);
		 }
                else if(covfile=="")
                        cout<<"No Covariate File Specified"<<endl;

		//unused argument
		/*
		std::string gcfile = command_line_opts.PAIR_PATH;
		int gc_pairs=0;  
		if(gcfile!=""){
			compute_gc=true; 
			gc_pairs = read_gc(gcfile);
		}
		*/
		if(pheno_num<2)
			compute_gc=false; //can not compute genetic correlation if there is only one phenotype 
		int window_size = command_line_opts.window_size; 	
	vector<double> ve_result; 
	vector<double> vg_result; 
	vector<double> vg_se; 
	double tr_k =0 ;
        double tr_k_rsid =0;
	
	yKy.resize(pheno_num,  pheno_num);
        Xy.resize(SNP_max, pheno_num);

//	all_yKy =MatrixXdr::Constant(cat_num, (pheno_num*(pheno_num -2))/2+pheno_num, 0); 
        cout<<"genotype size: "<<g.Nindv << " X "<< g.Nsnp<<endl;
        if(use_cov)
        {
			WW.resize(cov_num, cov_num*pheno_num); 
			pheno_prime.resize(cov_num, pheno_num); 
			for(int i=0; i<pheno_num; i++)
			{
				MatrixXdr covariate_cur = covariate; 
				MatrixXdr cur_mask = pheno_mask2.block(0,i, g.Nindv, 1); 
                        	for(int k=0; k<cov_num; k++) 
					covariate_cur.block(0,k,g.Nindv,1) = covariate_cur.block(0,k,g.Nindv,1).cwiseProduct(cur_mask); 
				MatrixXdr temp = covariate_cur.transpose() * covariate_cur;
				WW.block(0,cov_num*i, cov_num, cov_num) = temp.inverse();
                       		 pheno_prime.block(0,i,cov_num,1)= covariate_cur.transpose()* pheno.block(0,i, g.Nindv, 1);
       			}
	 }	
	//get means stds, same for all phenotypes
	//
	//
	//compute total number of blocks 	
	for (int chrom_i=0; chrom_i<22; chrom_i++)
                        {
                        int chrom_snp_num = g.get_chrom_snp(chrom_i);
		        int block_num = ceil(((float)chrom_snp_num)/window_size);

			total_num_blocks += block_num; 	
			}
	cout<<"Total number of blocks: "<<total_num_blocks<<endl; 
	MatrixXdr Xzb_total= MatrixXdr::Constant(g.Nindv, cat_num*pheno_num*10,0);
	MatrixXdr resid_total= MatrixXdr::Constant(g.Nindv, cat_num*pheno_num*10, 0);


	//fix random seed for debugging
	if(debug)
		srand(1);  
	MatrixXdr zb= MatrixXdr::Random(g.Nindv, 10);
        zb = zb * sqrt(3);
	//zb = MatrixXdr::Constant(g.Nindv, 10, 1); 
	// if(use_cov)
        //{
	//	MatrixXdr tremp = WW* covariate.transpose() * zb; 
	//	tremp = covariate * tremp;  
        //        zb = zb - tremp;
        //}
	cout<<"Start reading genotype round 1/2..." <<endl; 
	int snp_idx= 0 ;
	int block_idx=0; //total 118 blocks  
	for(int chrom_i=0; chrom_i<22; chrom_i++){
	SNP_BLOCK_SIZE = g.get_chrom_snp(chrom_i); 
	cout<< "Reading chromosome "<<chrom_i<<endl; 
//        int block_len = SNP_BLOCK_SIZE /5;
	//fix block_len
	int block_len = window_size; 
	int block_num = ceil(((float)SNP_BLOCK_SIZE)/window_size); 
	//previously, block_num was fixed to 5 
	for (int block_iter=0 ; block_iter<block_num; block_iter++)
	{
		int cur_len = block_len; 
		if(block_iter * block_len + block_len > SNP_BLOCK_SIZE) 
			cur_len = SNP_BLOCK_SIZE - block_iter*block_len; 	
	g.read_bed_mailman_stream(bedfile.str(), missing, pheno_mask2, pheno_num, snp_idx, ifs, cur_len);

	
        vsegsize=g.segment_size_ver;
        vsize= pow(3, vsegsize);
        partialsums = new double [blocksize];
        sum_op = new double[blocksize];
        yint_e = new double [hsize*blocksize];
        yint_m = new double [hsize*blocksize];
        memset (yint_m, 0, hsize*blocksize * sizeof(double));
        memset (yint_e, 0, hsize*blocksize * sizeof(double));

        y_e  = new double*[g.Nindv];
        for (int i = 0 ; i < g.Nindv ; i++) {
                y_e[i] = new double[blocksize];
                memset (y_e[i], 0, blocksize * sizeof(double));
        }

        y_m = new double*[hsegsize];
        for (int i = 0 ; i < hsegsize ; i++)
                y_m[i] = new double[blocksize];
	
       // yint_e_mat  = MatrixXdr::Constant(cat_num, hsize*blocksize, 0); 

        yint_e2 = new double*[cat_num];
        for (int i=0; i<cat_num; i++){
                yint_e2[i] = new double[hsize*blocksize];
                memset(yint_e2[i], 0, hsize*blocksize* sizeof(double));
        }
	//y_e_mat = MatrixXdr::Constant(g.Nindv*cat_num, blocksize, 0); 
        y_e2 = new double*[g.Nindv*cat_num];
        for (int i=0; i< (g.Nindv*cat_num); i++) {
                y_e2[i] = new double[blocksize];
                memset(y_e2[i], 0, blocksize* sizeof(double));
        }


        means.resize(g.Nsnp, pheno_num);
        stds.resize(g.Nsnp, pheno_num);
        sum.resize(g.Nsnp, pheno_num);
        sum2.resize(g.Nsnp, pheno_num);
	for (int k=0; k<pheno_num; k++){
	for(int i=0; i<g.Nsnp; i++)
	{
				means(i,k)= g.get_col_mean(i, k); 
				stds(i, k)= 1/g.get_col_std(i, k, exist_ind(0,k)); 
				sum2(i, k) = g.get_col_sum2(i, k); 
				sum(i,k) =g.get_col_sum(i, k); 
	}

	}
	//yKy Xy only compute once 
	compute_b1( use_cov,  y_sum, exist_ind(0,0), pheno_prime, pheno_fill, pheno_num, snp_idx,cov_num);
			//comptue tr[k]
///notice this opiton is only for genome wide	



	
	if(!noh2g){
	for(int i=0; i<pheno_num; i++) 
	{
			

		rhe_reg(zb, Xzb_total,resid_total, i, exist_ind(0, i), cov_num, cat_num, snp_idx, block_idx, pheno_num);  
		if(pheno_fill)
			exit;
	}
	}
	
	block_idx++;
	snp_idx += cur_len;  
	 delete[] sum_op;
       delete[] partialsums;
       delete[] yint_e;
        delete[] yint_m;

        for (int i  = 0 ; i < hsegsize; i++)
                delete[] y_m [i];
        delete[] y_m;

        for (int i  = 0 ; i < g.Nindv; i++)
                delete[] y_e[i];
        delete[] y_e;

        for (int i=0; i<(g.Nindv *cat_num); i++)
                delete[] y_e2[i] ;
        delete[] y_e2;

        for (int i=0; i<cat_num; i++)
                delete[] yint_e2[i];
        delete[] yint_e2;	
 
	g.not_O_j.clear(); 
	g.not_O_i.clear(); 
	g.p.clear(); 
	}//end loop for jackknife blocks 
	}

	yy = pheno.transpose() * pheno;
                if(use_cov)
                {
			MatrixXdr pheno_temp = pheno;
			for(int i=0; i<pheno_num; i++) 
			{	MatrixXdr temp = covariate * WW.block(0, i*cov_num, cov_num, cov_num) * pheno_prime.block(0,i, cov_num, 1);
				MatrixXdr cur_mask = pheno_mask2.block(0, i, g.Nindv, 1);  
				temp = temp.cwiseProduct(cur_mask);  
				pheno_temp.block(0, i, g.Nindv,1) = pheno_temp.block(0, i, g.Nindv,1) - temp ; 
			}
			yy  = pheno_temp.transpose() *  pheno_temp;  
		}


	//shortcut for genetic correlation
	 if(noh2g)
        {
	//noh2g option can only be used for no partition, thus all_yKy is a row 
		yKy = Xy.transpose()*Xy; 
		yKy= yKy/ SNP_max; 
		cout<<"yKy"<< yKy<<endl;
		cout<<"yy"<<yy<<endl; 
                for(int j=0; j<pheno_num; j++)
                {        for(int k=j+1; k<pheno_num; k++)
                        {
                                double X = yKy(j,k) - yy(j,k);
                                double Y = yKy(j,j) - yy(j,j);
                                double Z = yKy(k,k) - yy(k,k);
                                double rg= X/sqrt(Y)/sqrt(Z);
                                cout<<"Coheritability factor estimation for phenotype: "<<j << " , " << k <<endl;
                                cout<<"lambda_g: "<<rg<<endl;
                                double jack_knife_se = compute_jack_knife(j,k, rg, SNP_max);
                                cout<<"Jack Knife SE: "<<jack_knife_se<<endl;
                        }
                }
                return 0;
        }
//first row of h2g_estimates save point estimates of h2g
//first for of genCov_estimates save point estimates of genetic covariance
	MatrixXdr  h2g_estimates = MatrixXdr::Constant(total_num_blocks+1, pheno_num*cat_num, 0 ); 
	MatrixXd   genCov_estimates = MatrixXdr::Constant(total_num_blocks+1, pheno_num*(pheno_num-1)*cat_num/2, 0); 
	MatrixXdr A  =MatrixXdr::Constant(cat_num+1, cat_num+1, g.Nindv);
	A.block(0,0,cat_num, cat_num) = MatrixXdr::Constant(cat_num,cat_num, 0);  
	MatrixXdr b(cat_num+1,1); 	
		
	double vg, ve;
		//cout<<"pheno:"<<endl<<pheno.block(0,0,4,2)<<endl;
		//cout<<"Xy:"<<endl<<Xy.block(0,0,4,2)<<endl; 
		for(int i=0; i<pheno_num; i++)
		{
			double N_cur = exist_ind(0,i); 
			A= MatrixXdr::Constant(cat_num+1, cat_num+1, N_cur);
			 A.block(0,0,cat_num, cat_num) = MatrixXdr::Constant(cat_num,cat_num, 0);
			MatrixXdr cur_pheno_Xzb; 
			MatrixXdr cur_pheno_resid; 
			if(pheno_fill) 
			{
				cur_pheno_Xzb = Xzb_total.block(0, 0, g.Nindv, cat_num*10);
				cur_pheno_resid = resid_total.block(0,0,g.Nindv, cat_num*10); 
			}	
			if(!pheno_fill) 
			{
				cur_pheno_Xzb = Xzb_total.block(0, i*10*cat_num, g.Nindv, cat_num*10);
				cur_pheno_resid = resid_total.block(0, i*10*cat_num, g.Nindv, cat_num*10); 
			}
			MatrixXdr cur_mask=pheno_mask2.block(0, i, g.Nindv, 1); 
			for(int fix_rhe=0; fix_rhe<cat_num*10; fix_rhe++) 
			{
				cur_pheno_Xzb.block(0, fix_rhe, g.Nindv, 1) = cur_pheno_Xzb.block(0, fix_rhe, g.Nindv, 1).cwiseProduct(cur_mask); 
			}
			
			construct_linear_system(cur_pheno_Xzb,cur_pheno_Xzb,  cur_pheno_resid,cur_pheno_resid, A, b, zb, FunCat_snpcount, cat_num, i, i,SNP_max, cur_mask,cur_mask, FunCat); 
			b(cat_num,0) = yy(i,i); 
			cout<<"A: "<<endl<<A<<endl; 
			cout<<"b: "<<endl<<b<<endl; 
			MatrixXdr herit = A.colPivHouseholderQr().solve(b);
               	 	ve = herit(cat_num,0);
                	ve_result.push_back(ve);
                	//cout<<"V(e): "<<herit(cat_num,0)<<endl;
                	//cout<<"Vp "<<herit.sum()<<endl;
			for(int fun_t=0; fun_t<cat_num; fun_t++) 	
                	{			
				string pheno_idx = "Vg/Vp("+to_string(i)+")("+to_string(fun_t)+")"; 
				fprintf(fp_output, "%s\t%f\n", pheno_idx.c_str(), herit(fun_t,0)/herit.sum());
			}
			//cout<<"V(G)("<<fun_t<<"): "<<herit(fun_t, 0)/herit.sum()<<endl;
			h2g_estimates.block(0, i*cat_num, 1,cat_num) = herit.block(0,0,cat_num, 1).transpose();
			

	
                	/*if(bpheno){
                	cout<<"Prevelance: "<<prevelance(0,0)<<endl;
                	boost::math::normal m_normal(0.0, 1.0);
                	double t = quantile(m_normal,1-prevelance(0,i));
                	double c = pdf(m_normal, t);
			c = c*c;
               	 	c= 1/c;
                	c = c* prevelance(0,i) * (1-prevelance(0,i));
                	cout<<"Liability Scale: "<<herit(0,0)*c / herit.sum()<<endl;
			}
			MatrixXdr se(1,1);
               		MatrixXdr pheno_cur = pheno.block(0,i, g.Nindv, 1);
			MatrixXdr pheno_sum2 = pheno_cur.transpose() *pheno_cur;
                	double pheno_variance = pheno_sum2(0,0) / (exist_ind(0,i)-1);
			MatrixXdr Xy_cur  = Xy.block(0,i, g.Nsnp, 1); 	
			compute_se1(Xy_cur, pheno_cur, se,vg, ve, tr2,B, exist_ind(0, i));
			vg_se.push_back(se(0,0)); 
			cout<<"phenotype variance: "<<pheno_variance<<endl;
                	cout<<"sigma_g SE: "<<se<<endl;
                	cout<<"h2g SE:"<<se/pheno_variance<<endl;
			*/
		}
		if(pheno_num>1 ){
                int phenoPairIdx_temp=0; 
		for(int j=0; j<pheno_num; j++)
                        for(int k=j+1; k<pheno_num; k++)
                        {
	
			 A= MatrixXdr::Constant(cat_num+1, cat_num+1, g.Nindv);
                        MatrixXdr cur_pheno_Xzb1;
			MatrixXdr cur_pheno_Xzb2; 
                        MatrixXdr cur_pheno_resid1;
			MatrixXdr cur_pheno_resid2; 
                        if(pheno_fill)
                        {
                                cur_pheno_Xzb1 = Xzb_total.block(0, 0, g.Nindv, cat_num*10);
                                cur_pheno_resid1 = resid_total.block(0,0,g.Nindv, cat_num*10);
                        	cur_pheno_Xzb2= cur_pheno_Xzb1; 
				cur_pheno_resid2 = cur_pheno_resid1 ;
			}
                        if(!pheno_fill)
                        {
                                cur_pheno_Xzb1 = Xzb_total.block(0, j*10*cat_num, g.Nindv, cat_num*10);
				cur_pheno_Xzb2 = Xzb_total.block(0, k*10*cat_num, g.Nindv, cat_num*10);
                                cur_pheno_resid1 = resid_total.block(0, j*10*cat_num, g.Nindv, cat_num*10);
				cur_pheno_resid2 = resid_total.block(0, k*10*cat_num, g.Nindv, cat_num*10);
                        }
			MatrixXdr cur_mask1 = pheno_mask2.block(0,j , g.Nindv, 1); 
			MatrixXdr cur_mask2 = pheno_mask2.block(0,k , g.Nindv, 1); 
                        for( int fix_rhe=0; fix_rhe<cat_num*10; fix_rhe++) 			{
			cur_pheno_Xzb1.block(0, fix_rhe, g.Nindv, 1) = cur_pheno_Xzb1.block(0,fix_rhe,g.Nindv, 1).cwiseProduct(cur_mask2); 
			}
			double N_cur = cur_mask1.cwiseProduct(cur_mask2).sum(); 
			cout<<"total comman samples: "<<N_cur<<endl; 
			A.block(0,0,cat_num+1, cat_num+1) =  MatrixXdr::Constant(cat_num+1,cat_num+1, N_cur);
			construct_linear_system(cur_pheno_Xzb1,cur_pheno_Xzb1,  cur_pheno_resid1,cur_pheno_resid2, A, b, zb, FunCat_snpcount, cat_num, j,k,SNP_max, cur_mask1, cur_mask2, FunCat);
		b(cat_num,0) = yy(j,k);
                        cout<<"A: "<<endl<<A<<endl;
                        cout<<"b: "<<endl<<b<<endl;
                        MatrixXdr herit = A.colPivHouseholderQr().solve(b);
		        for(int fun_t=0; fun_t<cat_num; fun_t++)
			{
				string pheno_idx  = "rho_g("+to_string(j)+","+to_string(k)+")("+to_string(fun_t)+")"; 
				fprintf(fp_output, "%s\t%f\n", pheno_idx.c_str(),  herit(fun_t,0));
			}
			//cout <<"rho_g("<<fun_t<<"): "<<herit(fun_t, 0)<<endl;
                        fprintf(fp_output, "%s\t%f\n", "rho_e:", herit(cat_num, 0)); 
			//cout <<"rho_e: "<<herit(cat_num,0)<<endl;
			MatrixXdr rg_estimate( 1,cat_num) ; 
			for (int partition_i=0; partition_i<cat_num; partition_i++) 
			{
				rg_estimate( 0, partition_i) = herit(partition_i, 0) / sqrt(h2g_estimates(0, k*cat_num+partition_i)) / sqrt(h2g_estimates(0, j*cat_num+partition_i)); 
				string pheno_idx  = "gamma_g("+to_string(j)+","+to_string(k)+")("+to_string(partition_i)+")";
                                fprintf(fp_output, "%s\t%f\n", pheno_idx.c_str(),  rg_estimate(0, partition_i));
			//	cout<<"gamma_g("<<partition_i<<"): "<<rg_estimate(0,partition_i)<<endl; 

			}
			genCov_estimates.block(0,0,1, cat_num) = herit.block(0,0,cat_num, 1).transpose(); 
//			cout <<"gamma_g: "<<herit.block(0,0, cat_num, 1) / sqrt(h2g_estimates(0,j))/sqrt(h2g_estimates(0, k))<<endl ; 

					

			
			}
			
        }	
	

	//read the genome second time to perform jackknife 
	cout<<"Start reading genotype round 2/2..."<<endl; 
	snp_idx=0; 
	block_idx=0; 
	string temp_bedfile2 = bedfile.str(); 
	ifstream ifs2(temp_bedfile2.c_str(), ios::in|ios::binary); 
	char magic2[3]; 
	ifs2.read(magic2, 3*sizeof(char)); 


	vector<double> jack_weight(total_num_blocks, 0); 
	for(int chrom_i=0; chrom_i<22; chrom_i++){
	cout<<"reading chromosome "<<chrom_i<<endl; 
	SNP_BLOCK_SIZE = g.get_chrom_snp(chrom_i); 
	int block_len=window_size; 
	int block_num = ceil(((float)SNP_BLOCK_SIZE)/window_size); 

	for (int block_iter=0 ; block_iter<block_num; block_iter++)
        {
                int cur_len = block_len;
                if(block_iter * block_len + block_len > SNP_BLOCK_SIZE)
                        cur_len = SNP_BLOCK_SIZE - block_iter*block_len;
	
		g.read_bed_mailman_stream(bedfile.str(),missing, pheno_mask2, pheno_num, snp_idx, ifs2, cur_len);

	vsegsize=g.segment_size_ver;
        vsize= pow(3, vsegsize);
        partialsums = new double [blocksize];
        sum_op = new double[blocksize];
        yint_e = new double [hsize*blocksize];
        yint_m = new double [hsize*blocksize];
        memset (yint_m, 0, hsize*blocksize * sizeof(double));
        memset (yint_e, 0, hsize*blocksize * sizeof(double));

        y_e  = new double*[g.Nindv];
        for (int i = 0 ; i < g.Nindv ; i++) {
                y_e[i] = new double[blocksize];
                memset (y_e[i], 0, blocksize * sizeof(double));
        }

        y_m = new double*[hsegsize];
        for (int i = 0 ; i < hsegsize ; i++)
                y_m[i] = new double[blocksize];

	//yint_e_mat  = MatrixXdr::Constant(cat_num, hsize*blocksize, 0);
        yint_e2 = new double*[cat_num];
        for (int i=0; i<cat_num; i++){
                yint_e2[i] = new double[hsize*blocksize];
                memset(yint_e2[i], 0, hsize*blocksize* sizeof(double));
        }
	//y_e_mat = MatrixXdr::Constant(g.Nindv*cat_num, blocksize, 0);
        y_e2 = new double*[g.Nindv*cat_num];
        for (int i=0; i< (g.Nindv*cat_num); i++) {
                y_e2[i] = new double[blocksize];
                memset(y_e2[i], 0, blocksize* sizeof(double));
        }
	

	means.resize(g.Nsnp, pheno_num);
        stds.resize(g.Nsnp, pheno_num);
        sum.resize(g.Nsnp, pheno_num);
        sum2.resize(g.Nsnp, pheno_num);
        for (int k=0; k<pheno_num; k++) {
	for(int i=0; i<g.Nsnp; i++)
        {
                                means(i,k)= g.get_col_mean(i, k);
                                stds(i, k)= 1/g.get_col_std(i, k, exist_ind(0,k));
                                sum2(i, k) = g.get_col_sum2(i, k);
                                sum(i,k) =g.get_col_sum(i, k);
        }
	}
	MatrixXdr Xzb_curRES=MatrixXdr::Constant(g.Nindv, cat_num*pheno_num*10, 0); 
	MatrixXdr resid_curRES = MatrixXdr::Constant(g.Nindv, cat_num*pheno_num*10, 0);  	
	 for(int i=0; i<pheno_num; i++)
        {
		rhe_reg(zb, Xzb_curRES,resid_curRES, i, exist_ind(0, i), cov_num, cat_num, snp_idx, block_idx, pheno_num);
                if(pheno_fill)
                        exit;
	}
        

	 
	MatrixXdr Xzb_curBlock = Xzb_total - Xzb_curRES; 
	MatrixXdr resid_curBlock = resid_total - resid_curRES;  
	
	
	// h2g jackknife estimates
	for(int i=0; i<pheno_num; i++)
	{
		MatrixXdr cur_pheno_Xzb; 
		MatrixXdr cur_pheno_resid;
		double N_cur = exist_ind(0, i);  
		A  =MatrixXdr::Constant(cat_num+1, cat_num+1, N_cur);
                A.block(0,0,cat_num, cat_num) = MatrixXdr::Constant(cat_num,cat_num, 0);
		cur_pheno_Xzb = Xzb_curBlock.block(0, cat_num*i*10, g.Nindv, cat_num*10); 
		cur_pheno_resid = resid_curBlock.block(0, cat_num*i*10, g.Nindv, cat_num*10); 
		MatrixXdr cur_mask=pheno_mask2.block(0,i,g.Nindv, 1); 
                
		for(int fix_rhe=0; fix_rhe<cat_num*10; fix_rhe++)
                        {
                                cur_pheno_Xzb.block(0, fix_rhe, g.Nindv, 1) = cur_pheno_Xzb.block(0, fix_rhe, g.Nindv, 1).cwiseProduct(cur_mask);
                        }
		
		MatrixXdr temp_cat_mask = FunCat;
                temp_cat_mask.block(0, snp_idx, cat_num, cur_len) = MatrixXdr::Constant(cat_num, cur_len, 0);
		MatrixXdr temp_FunCat_snpcount = temp_cat_mask.rowwise().sum();
                construct_linear_system(cur_pheno_Xzb, cur_pheno_Xzb, cur_pheno_resid, cur_pheno_resid, A, b, zb ,temp_FunCat_snpcount, cat_num, i,i , SNP_max, cur_mask, cur_mask, temp_cat_mask);
	      	b(cat_num, 0 ) = yy(i,i); 
		 MatrixXdr herit_cur = A.colPivHouseholderQr().solve(b);
		h2g_estimates.block(block_idx+1,i*cat_num, 1, cat_num)   = herit_cur.block(0,0,cat_num, 1).transpose() ;
	}
	
	
	int phenoPairIdx=0; 
	for(int j=0; j<pheno_num; j++){
	for(int k=j+1; k<pheno_num; k++)
	{
		A =MatrixXdr::Constant(cat_num+1, cat_num+1, g.Nindv); 
		MatrixXdr cur_pheno_Xzb1;
               	MatrixXdr cur_pheno_Xzb2;
                MatrixXdr cur_pheno_resid1;
                MatrixXdr cur_pheno_resid2;
		 A.block(0,0,cat_num, cat_num) = MatrixXdr::Constant(cat_num,cat_num, 0);
		cur_pheno_Xzb1 = Xzb_curBlock.block(0, j*10*cat_num, g.Nindv, 10*cat_num);
                cur_pheno_Xzb2 = Xzb_curBlock.block(0, k*10*cat_num, g.Nindv, 10*cat_num);
                cur_pheno_resid1 = resid_curBlock.block(0, j*10*cat_num, g.Nindv, 10*cat_num);
                cur_pheno_resid2 = resid_curBlock.block(0, k*10*cat_num, g.Nindv, 10*cat_num);
		MatrixXdr cur_mask1 = pheno_mask2.block(0,j , g.Nindv, 1);
                        MatrixXdr cur_mask2 = pheno_mask2.block(0,k , g.Nindv, 1);
                       for( int fix_rhe=0; fix_rhe<cat_num*10; fix_rhe++)                      {
                        cur_pheno_Xzb1.block(0, fix_rhe, g.Nindv, 1) = cur_pheno_Xzb1.block(0,fix_rhe,g.Nindv, 1).cwiseProduct(cur_mask2);
			}
                       
			double N_cur = cur_mask1.cwiseProduct(cur_mask2).sum();
                        A.block(0,0,cat_num+1, cat_num+1) =  MatrixXdr::Constant(cat_num+1,cat_num+1, N_cur); 
		MatrixXdr temp_cat_mask = FunCat;
                temp_cat_mask.block(0, snp_idx, cat_num, cur_len) = MatrixXdr::Constant(cat_num, cur_len, 0);
                MatrixXdr temp_FunCat_snpcount = temp_cat_mask.rowwise().sum();
                construct_linear_system(cur_pheno_Xzb1, cur_pheno_Xzb1, cur_pheno_resid1, cur_pheno_resid2, A, b, zb ,temp_FunCat_snpcount, cat_num, j,k , SNP_max, cur_mask1, cur_mask2, temp_cat_mask);
                b(cat_num, 0 ) = yy(j,k); 
		MatrixXdr herit_cur = A.colPivHouseholderQr().solve(b);

		
//		if(block_idx==0)
//
//		cout<<Xzb_total.block(0,0,4,4)<<endl<<A<<endl<<b<<endl<<temp_FunCat_snpcount<<endl<<Xzb_curBlock.block(0,0,4,4)<<endl ; 
	
		
		genCov_estimates.block(block_idx+1, phenoPairIdx*cat_num, 1, cat_num) = herit_cur.block(0,0,cat_num, 1).transpose(); 	
		phenoPairIdx++; 
	}}
	//genetic covariance 
	
	//clear mem
	jack_weight[block_idx] = SNP_max-cur_len; 
	block_idx++;
        snp_idx += cur_len;
         delete[] sum_op;
       delete[] partialsums;
       delete[] yint_e;
        delete[] yint_m;

        for (int i  = 0 ; i < hsegsize; i++)
                delete[] y_m [i];
        delete[] y_m;

        for (int i  = 0 ; i < g.Nindv; i++)
                delete[] y_e[i];
        delete[] y_e;

        for (int i=0; i<(g.Nindv *cat_num); i++)
                delete[] y_e2[i] ;
        delete[] y_e2;

        for (int i=0; i<cat_num; i++)
                delete[] yint_e2[i];
        delete[] yint_e2;
	 g.not_O_j.clear();
        g.not_O_i.clear();
        g.p.clear();
	}//end of iter of blocks
	}//end of iter of chromosome

	for(int i=0; i<pheno_num; i++)
	{
		
		vector<double> jackknife_estimate(total_num_blocks, 1); 
		for(int fun_t=0; fun_t<cat_num; fun_t++){
		for (int block_t=0; block_t< total_num_blocks; block_t++)
		{
			jackknife_estimate[block_t] =h2g_estimates(block_t+1, i*cat_num+fun_t); 
			 
		}
		pair<double, double> jack_output; 
		double estimate=h2g_estimates(0,i*cat_num+fun_t); 
		jack_output= weightedjack(jackknife_estimate, jack_weight, estimate); 
		string pheno_idx= "SE(Vg/Vp)("+to_string(i)+")("+to_string(fun_t)+")";
		 fprintf(fp_output, "%s\t%f\n",pheno_idx.c_str(),  jack_output.second); 
		//cout<<"SE(h2g_"<<i<<"("<<fun_t<<")): "<<jack_output.second<<endl;
		}	
	}
	int phenoPairIdx=0; 
	for (int j=0; j<pheno_num; j++) {
		for(int k=j+1; k<pheno_num; k++) 
		{
		for(int fun_t=0;  fun_t<cat_num; fun_t++){
		vector<double> jackknife_estimate(total_num_blocks, 1); 
		vector<double> jackknife_covariance(total_num_blocks, 1); 
		for(int block_t=0; block_t<total_num_blocks; block_t++)
		{
			jackknife_estimate[block_t] = genCov_estimates(block_t+1, phenoPairIdx*cat_num+fun_t)/ sqrt(h2g_estimates(block_t+1, j*cat_num+fun_t)) / sqrt(h2g_estimates(block_t+1, k*cat_num+fun_t));
			jackknife_covariance[block_t] =  genCov_estimates(block_t+1, phenoPairIdx*cat_num+fun_t); 
		}
		pair<double, double> jack_output; 
		double estimate = genCov_estimates(0, phenoPairIdx*cat_num+fun_t)/sqrt(h2g_estimates(0, j*cat_num+fun_t)) / sqrt(h2g_estimates(0, k*cat_num+fun_t)); 	
		jack_output=weightedjack(jackknife_estimate, jack_weight, estimate); 
		string pheno_idx =  "SE(gamma_g)(" + to_string(j)+","+ to_string(k)+")("+to_string(fun_t)+")";
		 fprintf(fp_output, "%s\t%f\n" , pheno_idx.c_str(), jack_output.second);
		
		pair<double, double> jack_covariance; 
		double estimate_covariance = genCov_estimates(0, phenoPairIdx*cat_num+fun_t); 
		jack_covariance=weightedjack(jackknife_covariance, jack_weight, estimate_covariance);
                pheno_idx =  "SE(rho_g)(" + to_string(j)+","+ to_string(k)+")("+to_string(fun_t)+")";
                 fprintf(fp_output, "%s\t%f\n" , pheno_idx.c_str(), jack_covariance.second);
		//cout<<"SE(rg_("<<j<<","<<k<<")("<<fun_t<<")): "<<jack_output.second<<endl; 
		}
		phenoPairIdx++; 	
		}
	}
	fclose(fp_output); 
	return 0;
}
