# SCORE
SCalable genetic CORrelation Estimator


We propose a scalable randomized Method-of-Moments (MoM) estimator of genetic correlations in bi-variate LMMs. The method compute genetic correlation in two senerios: the samples of phenotypes are shared and partially shared. 

### Citing SCORE

If you find that this software is useful for your research project, 
please cite our paper: 
Yue Wu, Kathryn S. Burch, Andrea Ganna, Päivi ajukanta, Bogdan Pasaniuc, Sriram Sankararaman,
Fast estimation of genetic correlation for biobank-scale data,
The American Journal of Human Genetics,
Volume 109, Issue 1,
2022
## Getting Started

### Prerequisites
The following packages are required on a linux machine to compile and use the software package. 
```
g++
cmake
make
```

### Installing
Installing SCORE is fairly simple. Follow the following steps: 
```
git clone https://github.com/sriramlab/SCORE.git
cd SCORE
mkdir build 
cd build
cmake .. 
make
```

## Documentation for SCORE

After compiling the executble SCORE is present in the build directory. 
To run SCORE, use

``./SCORE <command_line arguments> ``

Notice this is a streaming version, which use less memory but requires two reads of the genome. 

### Parameters

The values in the brackets are the command line flags for running the code without parameter file. 

```
* genotype (-g) : The path of the genotype file or plink binary file prefix.
* phenotype (-p) : The path of the phenotype file.
* output (-o): The desitnation of formatted output. (Default: output)  
* covariate (-c) : The path of the covariate file.
* batch number (-b) : Number of random vectors used in the estimator % 10. 
* phenotype number (-mpheno) : The name of the pair of phenotypes to use in the phenotype file, seperated by comma. (eg. bmi,height)
  If not specified, SCORE will compute genetic correlation estimates on the first 2 phenotypes. 
* fill in missing phenotype with mean (-fill) : Fill in missing phenotypes with mean. 
* Computing genetic correlation only (noh2g): If pass in multiple phenotypes, the program will compute genetic correlation factors (rg). If the phenotypes share samples, one may compute genetic correlation without computing heritabitliy. This option must be combined with (-fill) if there is missingess in phenotypes.  
  Otherwise will be ignored. 
```


An example parameter file is provided in the example directory. 
You can run the code using the command: 
To compute genetic correlation in the ``` build ``` directory with desired output destination:  
```
./SCORE -g ../example/all -p pheno_1.pheno.plink -mpheno 1,2 -b 10 -o [OUTPUT DESTINATION DEFINED]  
```
To compute genetic correlation for shared-sample phenotypes: 
```
./SCORE -g ../example/all -p pheno_1.pheno.plink -mpheno 1,2 -fill -b 10
```
To compute genetics correlation only for shared-sample phenotypes: 
```
./SCORE -g ../example/all -p pheno_1.pheno.plink -mpheno 1,2 -fill -noh2g -b 10 
```

### Format of output
The output for the example provided is: 
```
Source  Value
Vg/Vp(0)(0)     0.732530
Vg/Vp(1)(0)     0.154493
rho_g(0,1)(0)   0.030006
rho_e:  -1.006169
gamma_g(0,1)(0) 0.088896
SE(Vg/Vp)(0)(0) 0.065486
SE(Vg/Vp)(1)(0) 0.019234
SE(gamma_g)(0,1)(0)     0.083957
```
Vg(i)(p), Vg/Vp(i)(p), and SE(Vg/Vp(i)(p) are the estimations of genetic variance component, heritability, and standard error for trait i for annotation group p.
Currently version has only one annotation group, which is the entire genome.The extension is under construction.  
rho_g(i,j)(p), gamma_g(i,j)(p) and SE(gamma_g)(i,j)(p) are genetic covariance, genetic correlation, and the standard error of genetic correlation for trait i and j for annotation group p. 
