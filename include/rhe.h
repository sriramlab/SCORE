#ifndef rhe_H
#define rhe_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include "storage.h"
#include <assert.h>
#include <emmintrin.h>
#include <iostream>
#include <pthread.h>
struct rhe_thread_args{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>   new_zb1;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  Xy;
	int exist_ind;
	int pheno_i;
	int cat_num;
	int snp_idx;
};



#endif
