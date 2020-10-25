# SCORE
SCalable genetic CORrelation Estimator


We propose a scalable randomized Method-of-Moments (MoM) estimator of genetic correlations in bi-variate LMMs. The method compute genetic correlation in two senerios: the samples of phenotypes are shared and partially shared. 

### Citing SCORE

If you find that this software is useful for your research project, 
please cite our paper: 

Wu, Y., Yaschenko, A., Heydary, M. H., & Sankararaman, S. (2019). Fast estimation of genetic correlation for Biobank-scale data. bioRxiv, 525055. 
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
git clone https://github.com/ariel-wu/SCORE.git
cd SCORE
mkdir build 
cd build
cmake .. 
make
```

## Documentation for RHE_reg

After compiling the executble SCORE is present in the build directory. 
To run RHE_reg, use

*``./SCORE <command_line arguments> ``

### Parameters

The values in the brackets are the command line flags for running the code without parameter file. 

```
* genotype (-g) : The path of the genotype file or plink binary file prefix.
* phenotype (-p) : The path of the phenotype file. 
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
To compute heritability and genetic correlation: 
```
../build/RHE_reg -g 200_100 -p 200_100.pheno.plink -c 200_100.cov -b 10 
```
To compute heritability and genetic correlation for shared-sample phenotypes: 
```
../build/RHE_reg -g 200_100 -p 200_100.pheno.pheno2 -fill -b 10
```
To compute genetics correlation only for shared-sample phenotypes: 
```
../build/RHE_reg -g 200_100 -p 200_100.pheno_pheno2 -fill -noh2g -b 10 
```

