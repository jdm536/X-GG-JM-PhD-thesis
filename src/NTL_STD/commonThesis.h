/*******************************************************************************
 * Header File: commonThesis.h
 *              Contains definitions of all common functions required
 *              for lattice-based cryptography
*******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h> // open function
#include <unistd.h> // close function
#include <math.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <NTL/fileio.h>
#include <NTL/vec_double.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <NTL/xdouble.h>
#include <NTL/quad_float.h>
#include <fplll.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <bits/stdc++.h>

#define _PARENT_DIR_ "../generate_SVP_bases"
#define _PARENT_INPUT_DIR_ "../generate_SVP_bases"
#define _PARENT_OUTPUT_DIR_ "../SVP_Outputs"
#define _PARENT_OUTPUT_DIR_RANDOM_ "../Random_Outputs"
#define _PARENT_OUTPUT_DIR_PREPROCESSED_ "../SVP_Outputs_Preprocessed"
#define _PARENT_OUTPUT_DIR_PREPROCESSED_MANYRUNS_ "../SVP_Outputs_Preprocessed_ManyRuns"
#define _PARENT_OUTPUT_DIR_PREPROCESSED_3LOOP_ "../SVP_Outputs_Preprocessed_3Loop"
#define _PARENT_OUTPUT_DIR_SS_GGLLL_PREPROCESSED_ "../SVP_Outputs_SSGGLLL_Preprocessed"
#define _LLL_RED_OUTPUT_DIR_ "../generate_SVP_bases/SVP_Bases_OurLLLReduced"
#define _SS_GGLLL_RED_OUTPUT_DIR_ "../generate_SVP_bases/SVP_Bases_SS_GGLLLReduced"
#define _WILLOW_PARENT_OUTPUT_DIR_ "../Willow_SVP_Outputs"
#define _INPUT_DIR_ "SVP_Bases"
#define _INPUT_DIR_FPLLL_ "SVP_Bases_FPLLL"
#define _INPUT_DIR_RANDOM_ "Random_Bases"
#define _INPUT_DIR_LLL_RED_ "SVP_Bases_OurLLLReduced"
#define _INPUT_DIR_SSGGLLL_RED_ "SVP_Bases_SS_GGLLLReduced"
#define _INPUT_DIR_FPLLL_RED_ "SVP_Bases_LLLReduced"
#define _INPUT_FILE_PREFIX_ "SVP_Basis"
#define _INPUT_FILE_PREFIX_RANDOM_ "Random_Basis"
#define _OUTPUT_FILE_PREFIX_ "out"

#define _LLL_         "LLL"
#define _DeepLLL_     "DeepLLL"
#define _DeepLLL_5_   "DeepLLL_5"
#define _DeepLLL_10_  "DeepLLL_10"
#define _DeepLLL_15_  "DeepLLL_15"
#define _DeepLLL_20_  "DeepLLL_20"
#define _LC_DeepLLL_  "LC_DeepLLL"
#define _LC_GGLLL_    "LC_GGLLL"
#define _LC_GGLLL_5_  "LC_GGLLL_5"
#define _LC_GGLLL_10_ "LC_GGLLL_10"
#define _LC_GGLLL_15_ "LC_GGLLL_15"
#define _LC_GGLLL_20_ "LC_GGLLL_20"
#define _LC_PGGLLL_   "LC_PGGLLL"
#define _LC_PGGLLL_5G_  "LC_PGGLLL_5G"
#define _LC_PGGLLL_10G_ "LC_PGGLLL_10G"
#define _LC_PGGLLL_15G_ "LC_PGGLLL_15G"
#define _LC_PGGLLL_20G_ "LC_PGGLLL_20G"
#define _LC_PGGLLL_5L_  "LC_PGGLLL_5L"
#define _LC_PGGLLL_10L_ "LC_PGGLLL_10L"
#define _LC_PGGLLL_15L_ "LC_PGGLLL_15L"
#define _LC_PGGLLL_20L_ "LC_PGGLLL_20L"
#define _SS_LLL_      "SS_LLL"
#define _SS_GGLLL_    "SS_GGLLL"
#define _Pot_LLL_     "Pot_LLL"
#define _Pot_GGLLL_   "Pot_GGLLL"
#define _Pot_F_LC_GGLLL_   "Pot_F_LC_GGLLL"
#define _SS_F_LC_GGLLL_   "SS_F_LC_GGLLL"
#define _SS_GGLLL_STD_   "SS_GGLLL_std"
#define _Pot_GGLLL_STD_   "Pot_GGLLL_std"
#define _SS_LLL_STD_   "SS_LLL_std"
#define _Pot_LLL_STD_   "Pot_LLL_std"

using namespace std;
using namespace NTL;

#define TRUE 1
#define FALSE 0

#define FILE_ERROR -1

#define NO_BITS_INT (sizeof(unsigned int)*8)
#define NO_BITS_DOUBLE (sizeof(long double)*8)
#define MASK_MSB_INT ((unsigned int)(1<<(NO_BITS_INT-1)))
#define MASK_MSB_DOUBLE ((unsigned int)(1<<(NO_BITS_DOUBLE-1)))
#define ALL_ONE_INT ((unsigned int)(((MASK_MSB_INT-1)<<1)+1))
#define ALL_ONE_DOUBLE ((unsigned int)(((MASK_MSB_DOUBLE-1)<<1)+1))

#define CLOSEST_INTEGER(X) ( (X>0) ?  ((long double)((long double)X - (int)X) > ((long double)0.51))? (int)X + 1: (int)X : ((long double)((long double)X - (int)X) < ((long double)-0.51))? ((int)X - 1): (int)X )

#define LC (Bi, Bj, delta, mu_ij) ((Bi >= ((delta - (mu_ij*mu_ij)) * Bj)) ? TRUE : FALSE)


inline int setPrecision(int m) 
{
	int x = 0;
      	if (m <= 40) {         
		x = 1100;
        	RR::SetPrecision(x);
   
         	#ifdef NTL_ZZ_NBITS 
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (600)

      	} else if (m <= 45) {
		x = 1250;
        	RR::SetPrecision(x);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (450)

      	} else if (m <= 50) {
		x = 1325;
                RR::SetPrecision(x);
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (500)
      	} else if (m <= 55) {
		x = 1375;
                RR::SetPrecision(x);
 
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (550)
      	} else if (m <= 60) {
		x = 1550;
                RR::SetPrecision(x);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (600)
      	} else if (m <= 65) {
		x = 1650;
                RR::SetPrecision(x);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (650)
      	} else if (m <= 70) {
		x = 1800;
                RR::SetPrecision(x);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (700)
      	} else if (m <= 75) {
		x = 1900;
                RR::SetPrecision(x);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (750)
      	} else if (m <= 80) {
		x = 2100;
                RR::SetPrecision(x);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (800)
      	} else if (m <= 85) {
		x = 2300;
                RR::SetPrecision(x);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (850)
      	} else if (m <= 90) {
		x = 2500;
                RR::SetPrecision(x);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (900)
      	} else if (m <= 100) {
		x = 3000;
                RR::SetPrecision(x);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1000)
      	} else if (m <= 105) {
		x = 3200;
                RR::SetPrecision(x);
   
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1100)
      	} else if (m <= 110) {
		x = 3400;
                RR::SetPrecision(x);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1100)
      	} else if (m <= 115) {
		x = 3600;
                RR::SetPrecision(x);
   
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1150)
      	} else if (m <= 120) {
		x = 3800;
                RR::SetPrecision(x);

         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1200)
      	} else if (m <= 125) {
		x = 4000;
                RR::SetPrecision(x);
 
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1250)
      	} else if (m <= 130) {
		x = 4500;
                RR::SetPrecision(x);
 
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1300)
      	} else if (m <= 140) {
		x = 5000;
                RR::SetPrecision(x);
 
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1400)
      	} else if (m <= 150) {
		x = 5700;
                RR::SetPrecision(x);
 
         	#ifdef NTL_ZZ_NBITS
         	#undef NTL_ZZ_NBITS
         	#endif
         	#define NTL_ZZ_NBITS (1500)
      	} else {
		x = 6000;
                RR::SetPrecision(x);
         	printf ("Warning: Dimension not yet tested to find precisions. \n");
      	}
	return x;
}

// Computes n choose r
inline int NCR(int n, int r) 
{
	if (r == 0) {	
		return 1;
	}

	if (r > n / 2) {
		return NCR(n, n - r); 
	}
	
	long res = 1; 

	for (int k = 1; k <= r; ++k)
				    {
	res *= n - k + 1;
	res /= k;
							        }
	return res;
}

// Compute SVP Challenge factor
inline RR SVP_Challenge_Factor(int n, ZZ volume) 
{
   	RR n_RR;
   	conv(n_RR, n);

   	RR power, volume_RR;
   	div(power, (RR) 1.0, n_RR);
   	conv(volume_RR, volume);
  
   	double gamma;
   	RR gamma_RR, nth_root_vol, nth_root_gamma; 
   	gamma = tgamma(n/2 + 1);
   	conv(gamma_RR, gamma);
  
   	nth_root_gamma = pow(gamma_RR, power); 
  
  	RR pi_RR, root_pi_RR;
   	ComputePi(pi_RR);
   	root_pi_RR = SqrRoot(pi_RR);
   
   	pow(nth_root_vol, volume_RR, power);

   	RR RRfactor;
   	RRfactor = (nth_root_gamma / root_pi_RR) * nth_root_vol;

   	return RRfactor;
}

// Generate a random integer from the given range (min, max)
inline int generateRandomInteger (int min, int max) 
{
   int r = rand();
   r = (r % (max-min))+min;
   return r;
}

// Generate a random integer vector
inline int * generateRandomIntegerVector (int m, int min, int max) 
{
   int * pRandomIntegerVector;
   pRandomIntegerVector = (int *) calloc (m, sizeof(int));
   for (int j=0; j<m; j++) {
      pRandomIntegerVector[j] = generateRandomInteger (min, max);
   }
   return pRandomIntegerVector;
}

// Generate a random basis with n integer vectors, each of dimension m
inline int ** generateRandomIntegerBasis (int m, int n, int min, int max) 
{
   int ** ppV;
   ppV = (int **) calloc (n, sizeof(int *));
   for (int i=0; i<n; i++) {
      ppV[i] = generateRandomIntegerVector (m, min, max);
   }
   return ppV;
}

// Generate a random NTL basis (mat_ZZ) with n integer vectors, each of dimension m
inline mat_ZZ generateRandomIntegerNTLBasis (int m, int n, int min, int max) 
{
   mat_ZZ NTLBasis;
   NTLBasis.SetDims(n,m);
   for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
         int r = rand();
         r = (r % (max-min))+min;
	 NTLBasis[i][j] = r;
      }
   }
   return NTLBasis;   
}

// Generate a random NTL basis (mat_ZZ) with n integer vectors, each of dimension m
inline mat_ZZ generateNTLChallengeBasis (int n) 
{
   mat_ZZ NTLBasis;
   ZZ prime;
   long length, err; 
   length = 10*n;
   GenPrime(prime, length, err = 80); 
   NTLBasis.SetDims(n,n);
   NTLBasis[0][0] = prime;
   for (int i=1; i<n; i++) {
      ZZ r;
      RandomBnd(r, prime);
      NTLBasis[i][0] = r;
      for (int j=1; j<n; j++){
         if (j == i) {
	    NTLBasis[i][j] = 1;
	 } else {
	    NTLBasis[i][j] = 0;
	 }	 
      }
   }
   return NTLBasis;   
}

// Print the basis vectors as columns
inline int printBasisAsColumns (int ** ppBasis, int m, int n) 
{
   for (int j=0; j<m; j++) {
      for (int i=0; i<n; i++) {
         printf ("%d\t", ppBasis[i][j]);
      }
      printf ("\n");
   }
   return 0;
}

// Print the basis vectors as rows
inline int printBasisAsRows (int ** ppBasis, int m, int n) 
{
   printf ("[");
   for (int i=0; i<n; i++) {
      printf ("[");
      for (int j=0; j<m-1; j++) {
         printf ("%d, ", ppBasis[i][j]);
      }
      printf ("%d]", ppBasis[i][m-1]);
      if (i<n-1) {
         printf (",\n");
      } else {
         printf ("]\n");
      }
   }
   return 0;
}

// Print the basis double vectors as columns
inline int printBasisDoubleAsColumns (long double ** ppBasisDouble, int m, int n) 
{
   for (int j=0; j<m; j++) {
      for (int i=0; i<n; i++) {
         printf ("%Lf\t", ppBasisDouble[i][j]);
      }
      printf ("\n");
   }
   return 0;
}

// Find the inner products of vectors pV1 with pV2 of differing types
template <typename T1, typename T2>
inline long double inner_product_template (T1 * pV1, T2 * pV2, int m) 
{
   long double result = 0;
   for (int i=0; i<m; i++) {
      result += (long double)pV1[i] * (long double)pV2[i];
   }
   return result;
}

// Inner Product using int
inline int inner_product_int (int * pV1, int * pV2, int m) 
{
	int result = 0;
	for (int i=0; i<m; i++) {
		result += (pV1[i] * pV2[i]);
	}
	
	return result;
}

// Inner Product using long
inline long int inner_product_long (int * pV1, int * pV2, int m) 
{
	long int result = 0;
	for (int i=0; i<m; i++) {
		result += (long int)pV1[i] * (long int)pV2[i];
	}
	
	return result;
}

// Inner Product using long long
inline long long inner_product_longlong (int * pV1, int * pV2, int m) 
{
	long long result = 0;
	for (int i=0; i<m; i++) {
		result += (long long)pV1[i] * (long long)pV2[i];
	}
	
	return result;
}

// Inner product of vec_ZZ with vec_ZZ
inline ZZ NTL_inner_product (vec_ZZ V1, vec_ZZ V2, int m) 
{
   ZZ result;
   for (int i=0; i<m; i++) {
      result += V1[i] * V2[i];
   }
   return result;   
}

// Compute norm of a vec_ZZ
inline RR NTL_vector_norm (vec_ZZ V1, int m) 
{
   ZZ norm_squared;
   InnerProduct(norm_squared, V1, V1);
   RR RR_norm_squared, norm;
   conv(RR_norm_squared, norm_squared);
   SqrRoot(norm, RR_norm_squared);

   return norm;
}

// Inner product between vec_ZZ and double/long double
template <typename T1>
inline T1 NTL_double_inner_product (vec_ZZ V1, T1 * doubleV2, int m) 
{
   T1 result = 0;
   for (int i=0; i<m; i++) {
      T1 V1_element = (T1) 0;
      conv(V1_element, V1[i]);
      result += V1_element * doubleV2[i];
   }
   return result;
}

// Compute norm of an integer vector
inline long double vector_norm (int * pV1, int m) 
{
   long double b1NormSquared = inner_product_template (pV1, pV1, m);
   long double b1Norm;
   b1Norm = sqrtl(b1NormSquared);
   return b1Norm;
}

// Size reduce vector pV2 with pV1
inline int reduce (int * pV1, int * pV2, int m, int closest_integer_to_mu) 
{
   for (int i=0; i<m; i++) {
      pV2[i] -= closest_integer_to_mu * pV1[i];
   }
   return 0;
}

// Size reduce vec_ZZ V2 with vec_ZZ V1
inline int vec_ZZ_reduce (vec_ZZ& V1, vec_ZZ& V2, int m, int closest_integer_to_mu) 
{
   ZZ closest_ZZ_to_mu;
   conv(closest_ZZ_to_mu, closest_integer_to_mu);
   for (int i=0; i<m; i++) {
      ZZ V1_mu;
      mul(V1_mu, closest_ZZ_to_mu, V1[i]);
      sub(V2[i], V2[i], V1_mu);
   }
   return 0;
}

// Size reduce vec_ZZ V2 with vec_ZZ V1 (only NTL)
inline int NTL_reduce (vec_ZZ& V1, vec_ZZ& V2, int m, ZZ closest_ZZ_to_mu) 
{
   for (int i=0; i<m; i++) {
      ZZ V1_mu;
      mul(V1_mu, closest_ZZ_to_mu, V1[i]);
      sub(V2[i], V2[i], V1_mu);
   }
   return 0;
}


// Size reduce pV1 with pV2
template <typename T1, typename T2>
inline int reduceDoubleTemplate (T1 * pV1, T2 * pV2, int m, long double mu) 
{
   if (typeid(T1) == typeid(T2)) {
      T1 mu1 = (T1) mu;
      for (int i=0; i<m; i++) {
         pV2[i] -= mu1*pV1[i];
      }
   } else {
      pV1 = (double *) pV1;
      pV2 = (double *) pV2;
      double mu1 = (double) mu;

      for (int j=0; j<m; j++) {
         pV2[j] -= mu*pV1[j];
      }
   }
   return 0;
}

// Size reduce vector pVDouble2 with pVDouble1
template <typename T1>
inline int reduceDouble (T1 * pVDouble1, T1 * pVDouble2, int m, T1 mu) 
{
   for (int i=0; i<m; i++) {
      pVDouble2[i] -= mu * pVDouble1[i];
   }
   return 0;
}

// Copy pBasis into pBasis2 and return
inline int ** copyBasis (int ** ppBasis, int m, int n) 
{
   int ** ppBasis2;
   // Create the second 2-dimensional array
   ppBasis2 = (int **) calloc (n, sizeof(int *));
   for (int i=0; i<n; i++) {
      ppBasis2[i] = (int *) calloc (m, sizeof(int));
      // Copy the values
      for (int j=0; j<m; j++) {
         ppBasis2[i][j] = ppBasis[i][j];
      }
   }
   return ppBasis2;
}

// Copy int **ppBasis to mat_ZZ BasisNTL and return
inline mat_ZZ copyBasisToNTL (int ** ppBasis, int m, int n) 
{
   mat_ZZ BasisNTL;
   BasisNTL.SetDims(n,m);
   for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
         BasisNTL[i][j] = ppBasis[i][j];
      }
   }
   return BasisNTL;
}

// Copy long double **ppBasis to mat_ZZ BasisNTL and return
inline mat_ZZ copyBasisLDToNTL (long double ** ppBasis, int m, int n) 
{
   mat_ZZ BasisNTL;
   
   BasisNTL.SetDims(n,m);
   for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
         BasisNTL[i][j] = (int) ppBasis[i][j];
      }
   }
   return BasisNTL;
}

// Copy mat_ZZ ppBasis to mat_ZZ BasisNTL and return
inline mat_ZZ copyNTLBasisToNTL (mat_ZZ ppBasis, int m, int n) 
{
   mat_ZZ BasisNTL;
   BasisNTL.SetDims(n,m);
   for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
         BasisNTL[i][j] = ppBasis[i][j];
      }
   }
   return BasisNTL;
}

// Copy mat_ZZ ppBasis to mat_ZZ BasisNTL and return
inline mat_RR copyNTLBasisToRR (mat_ZZ ppBasis, int m, int n) 
{
   mat_RR BasisRR;

   BasisRR.SetDims(n,m);
   for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
         conv(BasisRR[i][j], ppBasis[i][j]);
      }
   }
   return BasisRR;
}

// Copy int basis to long double basis return
inline long double ** copyBasisToDouble (int ** pBasis, int m, int n) 
{
   long double ** ppBasisDouble;
   // Create the 2-dimensional array
   ppBasisDouble = (long double **) calloc (n, sizeof(long double *));
   for (int i=0; i<n; i++) {
      ppBasisDouble[i] = (long double *) calloc (m, sizeof(long double));
      // Copy the values
      for (int j=0; j<m; j++) {
         ppBasisDouble[i][j] = (long double) pBasis[i][j];
      }
   }
   return ppBasisDouble;
}

// Copy int basis to double basis return
inline double ** copyBasisToDoub (int ** pBasis, int m, int n) 
{
   double ** ppBasisDouble;
   // Create the 2-dimensional array
   ppBasisDouble = (double **) calloc (n, sizeof(double *));
   for (int i=0; i<n; i++) {
      ppBasisDouble[i] = (double *) calloc (m, sizeof(double));
      // Copy the values
      for (int j=0; j<m; j++) {
         ppBasisDouble[i][j] = (double) pBasis[i][j];
      }
   }
   return ppBasisDouble;
}

// Copy NTL vector (vec_ZZ) to long double
inline long double * copyVectorToLD (int * pVec, int m) 
{
   long double * pVecDouble;
   pVecDouble = (long double *) calloc (m, sizeof(long double));
   for (int i=0; i<m; i++) {
      pVecDouble[i] = (long double) pVec[i];
   }
   return pVecDouble;
}	

// Copy NTL vector (vec_ZZ) to double
inline double * copyNTLVectorToDouble (vec_ZZ V1, int m) 
{
   double * VecDouble;
   VecDouble = (double *) calloc (m, sizeof(double));
   for (int i=0; i<m; i++) {
      conv(VecDouble[i], V1[i]);
   }
   return VecDouble;
}	

// Copy BasisNTL (mat_ZZ) into ppBasisDouble and return
inline double ** copyNTLBasisToDouble (mat_ZZ ppBasis, int m, int n) 
{
   double ** ppBasisDouble;
   // Create the 2-dimensional array
   ppBasisDouble = (double **) calloc (n, sizeof(double *));
   for (int i=0; i<n; i++) {
      ppBasisDouble[i] = (double *) calloc (m, sizeof(double));
      // Copy the values
      for (int j=0; j<m; j++) {
	 conv(ppBasisDouble[i][j], ppBasis[i][j]);
      }
   }
   return ppBasisDouble;
}

// Copy BasisNTL (mat_ZZ) to an int ** Basis and return
inline int ** copyNTLBasisToInt (mat_ZZ pBasis, int m, int n) 
{
   int ** ppBasisInt;
   // Create the 2-dimensional array
   ppBasisInt = (int **) calloc (n, sizeof(int *));
   for (int i=0; i<n; i++) {
      ppBasisInt[i] = (int *) calloc (m, sizeof(int));
      // Copy the values
      for (int j=0; j<m; j++) {
         conv(ppBasisInt[i][j], pBasis[i][j]);
      }
   }
   return ppBasisInt;
}

// Deallocate ppBasis
inline int deleteBasis (int ** ppBasis, int n) 
{
   for (int i=0; i<n; i++) {
      free (ppBasis[i]);
   }
   free (ppBasis);
   return 0;
}

// Deallocate ppBasisDouble
template <typename T1>
inline int deleteBasisDouble (T1 ** ppBasisDouble, int n) 
{
   for (int i=0; i<n; i++) {
      free (ppBasisDouble[i]);
   }
   free (ppBasisDouble);
   return 0;
}

// Deallocate ppBasisDouble
inline int deleteBasisLongDouble (long double ** ppBasisDouble, int n) 
{
   for (int i=0; i<n; i++) {
      free (ppBasisDouble[i]);
   }
   free (ppBasisDouble);
   return 0;
}

// Check equality of bases of equal dimension 
inline int basisEqualCheck(int ** ppBasis1, int ** ppBasis2, int m, int n) 
{
   int flag = 1; 
   for (int i=0; i<n; i++){
      for (int j= 0; j<m; j++) {
         if (ppBasis1[i][j] != ppBasis2[i][j]) {
	    flag = 0;
	    break;
	 }
      }
   }
   return flag;
}

// Compute matrix of delta[i][j] 
// delta[i][j] = ||\pi_j(b_i)||^2 / ||b_j^*||^2
inline mat_RR computeDeltaMatrix (mat_ZZ ppBasis, mat_RR ppM, vec_RR pB, int n, int m) 
{
   RR C;
   mat_RR ppBasisRR, ppDelta;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);
   ppDelta.SetDims(n,m);

   for (int k=1; k<n; k++) {
      C = 0;
      InnerProduct(C, ppBasisRR[k], ppBasisRR[k]);
      for (int j=0; j<k; j++) {
         div(ppDelta[k][j], C, pB[j]);
	 C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
      }  
   }
   return ppDelta;
}

// Gram-Schmidt Orthogonalisation
inline long double ** GSO (int ** pBasis, int m, int n) 
{
   long double ** ppBasisGSO;
   ppBasisGSO = copyBasisToDouble (pBasis, m, n);

   for (int i=1; i<n; i++) {
      for (int j=i-1; j>=0; j--) {

         // Find the value of mu_ij
         long double mu_ij =
                 inner_product_template (pBasis[i], ppBasisGSO[j], m) /
                 inner_product_template (ppBasisGSO[j], ppBasisGSO[j], m);

         // Reduce vector ppBasisGSO[i] with ppBasisGSO[j]
	 // Compute b_i - mu_ij x b_j*
         reduceDouble (ppBasisGSO[j], ppBasisGSO[i], m, mu_ij);
      }
   }
   return ppBasisGSO;
}

// Double Gram-Schmidt Orthogonalisation
inline double ** GSO_double (int ** pBasis, int m, int n) 
{
   double ** ppBasisGSO;
   ppBasisGSO = copyBasisToDoub (pBasis, m, n);

   for (int i=1; i<n; i++) {
      for (int j=i-1; j>=0; j--) {

         // Find the value of mu_ij
         double mu_ij =
                 inner_product_template (pBasis[i], ppBasisGSO[j], m) /
                 inner_product_template (ppBasisGSO[j], ppBasisGSO[j], m);

         // Reduce vector ppBasisGSO[i] with ppBasisGSO[j]
         // Compute b_i - mu_ij x b_j*
         reduceDouble (ppBasisGSO[j], ppBasisGSO[i], m, mu_ij);
      }
   }
   return ppBasisGSO;
}

inline mat_RR GSO_RR (mat_ZZ ppBasis, int m, int n) 
{
   mat_RR ppBasisGSO;
   mat_RR ppBasisRR;
   ppBasisGSO.SetDims(m, n);
   ppBasisRR.SetDims(m, n);
   for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
         conv(ppBasisGSO[i][j], ppBasis[i][j]);
         conv(ppBasisRR[i][j], ppBasis[i][j]);
      }
   }
   
  for (int i=0; i<n; i++) {
     for (int j=i-1; j>=0; j--) {
        RR mu_ij, num, den;
	InnerProduct(num, ppBasisRR[i], ppBasisGSO[j]);
	InnerProduct(den, ppBasisGSO[j], ppBasisGSO[j]);

	div(mu_ij, num, den);

	vec_RR product;
	mul(product, mu_ij, ppBasisGSO[j]);
	sub(ppBasisGSO[i], ppBasisGSO[i], product);
     }
  }
  return ppBasisGSO; 
}

// Gram-Schmidt Orthogonalisation for mat_ZZ input
inline double ** GSO_NTL (mat_ZZ pBasis, int m, int n) 
{
   double ** ppBasisGSO;
   ppBasisGSO = copyNTLBasisToDouble (pBasis, m, n);

   for (int i=1; i<n; i++) {
      for (int j=i-1; j>=0; j--) {

         // Find the value of mu_ij
         double mu_ij;
         double * pBasisDouble;
	 pBasisDouble = copyNTLVectorToDouble(pBasis[i],m);
	 
	 double num;
	 double den;

	 num = inner_product_template (pBasisDouble, ppBasisGSO[j], m);
	 den = inner_product_template (ppBasisGSO[j], ppBasisGSO[j], m);

	 mu_ij = num/den;

         // Reduce vector ppBasisGSO[i] with ppBasisGSO[j]
         reduceDoubleTemplate (ppBasisGSO[j], ppBasisGSO[i], m, mu_ij);
      }
   }
   return ppBasisGSO;
}

// Update the GSO information in the way required for DeepLLL
inline int deep_LLL_GSO_update (long double ** &ppM, long double * &pB, int k , int i, int n) 
{
   	long double * pP;
   	pP = (long double *) calloc (n, sizeof(long double));
   	long double * pS;
   	pS = (long double *) calloc (n, sizeof(long double));
   	long double * pD;
   	pD = (long double *) calloc (n, sizeof(long double));

   	for (int j=0; j<n; j++) {
      		pP[j] = 0.0;
      		pS[j] = 0.0;
      		pD[j] = 0.0;
   	}


   	pP[k] = pB[k];
   	pD[k] = pB[k];

   	// Lines 2-4 of Algorithm 4 in YY18
   	for (int j=k-1; j>=i; j--) {
      		pP[j] = ppM[k][j] * pB[j];
      		pD[j] = pD[j+1] + (ppM[k][j]*pP[j]);
   	}

   	// Line 5 of Algorithm 4 in YY18 - Intialising values of S
   	for (int j=i; j<n; j++) {
     		pS[j] = 0;
   	}

      	long double T;
   	
	// Lines 6 - 18 of Algorithm 4 in YY18
   	for (int j=k; j>i; j--) {
      		T = 0.0;
      		T = ppM[k][j-1] / pD[j];
      		for (int l=n-1; l>k;l--) {
         		pS[l] = pS[l] + (ppM[l][j] * pP[j]);
         		ppM[l][j] = ppM[l][j-1] - (T * pS[l]);
      		}	

      		for (int l=k; l>=j+2;l--) {
         		pS[l] = pS[l] + (ppM[l-1][j] * pP[j]);
         		ppM[l][j] = ppM[l-1][j-1] - (T * pS[l]);
      		}

      		if (j != k) {
         		pS[j+1] = pP[j];
         		ppM[j+1][j] = ppM[j][j-1] - (T * pS[j+1]);
      		}
   	}

   	// Lines 20 - 22 of Algorithm 4 in YY18
   	T = 1.0 / pD[i];


   	for (int l=n-1; l>k; l--) {
      		ppM[l][i] = T * (pS[l] + ppM[l][i] * pP[i]);
   	}
   	for (int l=k; l>=i+2; l--) {
      		ppM[l][i] = T * (pS[l] + ppM[l-1][i] * pP[i]);
   	}

   	// Lines 23 - 29 of Algorithm 4 in YY18
   	ppM[i+1][i] = T * pP[i];

   	for (int j=0; j<i; j++) {
      		long double epsilon;
      		epsilon = ppM[k][j];
      		
		for (int l=k; l>i; l--) {
         		ppM[l][j] = ppM[l-1][j];
      		}
         	
		ppM[i][j] = epsilon;
   	}

   	// Lines 31 - 32 of Algorithm 4 in YY18
   	for (int j=k; j>i;j--) {
      		pB[j] = (pD[j] * pB[j-1]) / pD[j-1];
   	}
   	
	pB[i] = pD[i];

      	free (pD);
      	free (pP);
      	free (pS);

   	return 0;
}

// Update the GSO information in the way required for DeepLLL (RR datatypes)
inline int deep_LLL_GSO_update_RR (mat_RR &ppM, vec_RR &pB, int k , int i, int n) 
{
   // Initialise arrays pP, pS, pD as described in [YY18]
   vec_RR pP;
   pP.SetLength(n);
   vec_RR pS;
   pS.SetLength(n);
   vec_RR pD;
   pD.SetLength(n);

   for (int j=0; j<n; j++) {
      pP[j] = 0;
      pS[j] = 0;
      pD[j] = 0;
   }


   pP[k] = pB[k];
   pD[k] = pB[k];

   // Lines 2-4 of Algorithm 4 in YY18
   for (int j=k-1; j>=i; j--) {
      pP[j] = ppM[k][j] * pB[j];
      pD[j] = pD[j+1] + (ppM[k][j]*pP[j]);
   }

   // Line 5 of Algorithm 4 in YY18 - Intialising values of S
   for (int j=i; j<n; j++) {
      pS[j] = 0;
   }

   // Lines 6 - 18 of Algorithm 4 in YY18
   for (int j=k; j>i; j--) {
      RR T;
      T = 0;
      T = ppM[k][j-1] / pD[j];
      for (int l=n-1; l>k;l--) {
         pS[l] = pS[l] + (ppM[l][j] * pP[j]);
         ppM[l][j] = ppM[l][j-1] - (T * pS[l]);
      }

      for (int l=k; l>=j+2;l--) {
         pS[l] = pS[l] + (ppM[l-1][j] * pP[j]);
         ppM[l][j] = ppM[l-1][j-1] - (T * pS[l]);
      }

      if (j != k) {
         pS[j+1] = pP[j];
         ppM[j+1][j] = ppM[j][j-1] - (T * pS[j+1]);
      }
   }

   // Lines 20 - 22 of Algorithm 4 in YY18
   RR T;
   inv(T, pD[i]);


   for (int l=n-1; l>k; l--) {
      ppM[l][i] = T * (pS[l] + ppM[l][i] * pP[i]);
   }
   for (int l=k; l>=i+2; l--) {
      ppM[l][i] = T * (pS[l] + ppM[l-1][i] * pP[i]);
   }

   // Lines 23 - 29 of Algorithm 4 in YY18
   ppM[i+1][i] = T * pP[i];

   for (int j=0; j<i; j++) {
      RR epsilon;
      epsilon = ppM[k][j];
      for (int l=k; l>i; l--) {
         ppM[l][j] = ppM[l-1][j];
      }
         ppM[i][j] = epsilon;
   }

   // Lines 31 - 32 of Algorithm 4 in YY18
   for (int j=k; j>i;j--) {
      pB[j] = (pD[j] * pB[j-1]) / pD[j-1];
   }
   pB[i] = pD[i];

   return 0;
}

template <typename T1>
inline T1 find_insertion_indices_LC_GG (T1 ** ppM, T1 ** ppDelta, int ** ppBasis, T1 * pB, int n, int m, int * pK, int * pJ) 
{
   int k_min = 1;
   int i_min = 0;
   T1 C = 0.0;
   T1 minDelta =  (inner_product_template (ppBasis[1], ppBasis[1],m) / pB[0]);
   T1 delta_kj = 0.0;

   for (int k=2; k<n; k++) {
      C = (T1) inner_product_template (ppBasis[k], ppBasis[k], m);
      for (int j=0; j<k; j++) {
         // Lovasz condition definition
         delta_kj = (C / pB[j]);
         if (delta_kj < minDelta) {
            k_min = k;
            i_min = j;
            minDelta = delta_kj;
         }
         // Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
         C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
      }
   }
   *pK = k_min;
   *pJ = i_min;

   return minDelta;
}

// Find insertion indices for LC-GG
inline long double find_insertion_indices_LC_GG (long double ** ppM, int ** ppBasis, long double * pB, int n, int m, int * pK, int * pJ) 
{
   	int k_min = 1;
   	int i_min = 0;
   	long double C = 0.0;
	int inner;
	inner = inner_product_template (ppBasis[1], ppBasis[1],m);
   	long double minDelta = (inner_product_template (ppBasis[1], ppBasis[1],m) / pB[0]);
   	long double delta_kj = 0.0;

   	for (int k=2; k<n; k++) {
		C = (long double) inner_product_template(ppBasis[k], ppBasis[k], m);
      		for (int j=0; j<k; j++) {
         		// Lovasz condition definition
			delta_kj = (C / pB[j]);
         		if (delta_kj < minDelta) {
            			k_min = k;
            			i_min = j;
            			minDelta = delta_kj;
         		}
         		
			// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
         		C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
      		}
   	}
   	
	*pK = k_min;
   	*pJ = i_min;

   	return minDelta;
}

// LC-PGG recursive search (locally restricted)
inline int block_LC_PGG_recursive_search_local (long double ** ppDelta, long double ** &ppM, long double * &pB, int ** &ppBasis, int i, int k, int n, int epsilon, long double delta, int * nInsertions) 
{
	int i_min = k+1;
	int i_left = k+1;

	if (i >= k) {
		return i_min;
	}
	
	int * pTemp;

	long double minDelta = delta + 1;
	int ub = 0;
	int lb = 0;
	int i1 = k+1;
	int k1 = k+1;

	// Check all ppDelta in range [i,i+epsilon-1] and [(max(epsilon,k-epsilon), k-1]
	// Find min delta
	for (int l=i+1; l<=k; l++) {
		lb = max(i+epsilon-1, l-epsilon);
		ub = min(l, i+epsilon);
		for (int j=i; j<ub; j++) {
			if (ppDelta[l][j] < minDelta) {
                                minDelta = ppDelta[l][j];
                                i1 = j;
                                k1 = l;
                        }
		}

		for (int j=lb; j<l; j++) {
			if (ppDelta[l][j] < minDelta) {
                                minDelta = ppDelta[l][j];
                                i1 = j;
                                k1 = l;
                        }
		}
	}

	if (minDelta > delta) {
		// In this case, no insertion is to be performed
		return i_min;
	} 

	(*nInsertions)++;
	// Deep insert b_{k1} in position i1 
      	pTemp = ppBasis[k1];
	for (int j=k1-1; j>=i1; j--) {
        	ppBasis[j+1] = ppBasis[j];
      	}
      	ppBasis[i1] = pTemp;

	i_min = i1;

      	deep_LLL_GSO_update (ppM, pB, k1 , i1, n);
		
	// Index search (i,i1-1), (k1+1,k)
	i_left = block_LC_PGG_recursive_search_local (ppDelta, ppM, pB, ppBasis, i, i1-1, n, epsilon, delta, nInsertions);
	// Check here if the i_left returned by the call above is less than the minimum i obtained thus far
	// If so, we set i_min = i_left (returned by function above)
	if (i_left < i_min) {
		i_min = i_left;
	}

	block_LC_PGG_recursive_search_local (ppDelta, ppM, pB, ppBasis, k1+1, k, n, epsilon, delta, nInsertions);

	return i_min;
}

// LC-PGG recursive search (globally restricted)
inline int block_LC_PGG_recursive_search_global (long double ** ppDelta, long double ** &ppM, long double * &pB, int ** &ppBasis, int i, int k, int n, int epsilon, long double delta, int * nInsertions) 
{
	int i_min = k+1;
	int i_left = k+1;

	if (i >= k) {
		return i_min;
	}
	
	int * pTemp;

	long double minDelta = delta + 1.0;
        int lb = 0;
        int ub = 0;
        int i1 = k+1;
        int k1 = k+1;

        // Check all ppDelta in range [1,epsilon] and [(max(epsilon,k-epsilon), k-1]
        // Find min delta
        for (int l=i+1; l<=k; l++) {
		ub = min(epsilon, l);
        	lb = max(epsilon, l-epsilon);
		for (int j=i; j<ub; j++) {
			if (ppDelta[l][j] < minDelta) {
                                minDelta = ppDelta[l][j];
                                i1 = j;
                                k1 = l;
                        }
		}
		for (int j=lb; j<l; j++) {
			if (ppDelta[l][j] < minDelta) {
                                minDelta = ppDelta[l][j];
                                i1 = j;
                                k1 = l;
                        }
		}	
	}	

	if (i1 == k+1 && k1 == k+1) {
		return i_min;
	}

	if (minDelta > delta) {
		// In this case, no insertion is to be performed
		return i_min;
	} 

	(*nInsertions)++;
	// Deep insert b_{k1} in position i1 
      	pTemp = ppBasis[k1];
	for (int j=k1-1; j>=i1; j--) {
        	ppBasis[j+1] = ppBasis[j];
      	}
      	ppBasis[i1] = pTemp;

	i_min = i1;
	cout << "k = " << k1 << "   i = " << i1 << endl; 

      	deep_LLL_GSO_update (ppM, pB, k1 , i1, n);
		
	// Index search (i,i1-1), (k1+1,k)
	i_left = block_LC_PGG_recursive_search_global (ppDelta, ppM, pB, ppBasis, i, i1-1, n, epsilon, delta, nInsertions);

	// Check here if the i_left returned by the call above is less than the minimum i obtained thus far
	// If so, we set i_min = i_left (returned by function above)
	if (i_left < i_min) {
		i_min = i_left;
	}

	block_LC_PGG_recursive_search_global (ppDelta, ppM, pB, ppBasis, k1+1, k, n, epsilon, delta, nInsertions);

	return i_min;
}

// LC-PGG recursive search
inline int LC_PGG_recursive_search (long double ** ppDelta, long double ** &ppM, long double * &pB, int ** &ppBasis, int i, int k, int n, long double delta, int * nInsertions) 
{
	int i_min = k+1;
	int i_left = k+1;

	if (i >= k) {
		return i_min;
	}
	
	int * pTemp;

	long double minDelta = ppDelta[i+1][i];
	int i1 = i;
	int k1 = i+1;

	// Check all ppDelta in range [i,k]
	// Find min delta
	for (int l=i+2; l<=k; l++) {
		for (int j=i; j<l; j++) {
			if (ppDelta[l][j] < minDelta) {
				minDelta = ppDelta[l][j];
				i1 = j;
				k1 = l;
			}
		}
	}

	if (minDelta > delta) {
		// In this case, no insertion is to be performed
		
		return i_min;
	} 

	(*nInsertions)++;
	// Deep insert b_{k1} in position i1 
      	pTemp = ppBasis[k1];
	for (int j=k1-1; j>=i1; j--) {
        	ppBasis[j+1] = ppBasis[j];
      	}
      	ppBasis[i1] = pTemp;

	i_min = i1;

      	deep_LLL_GSO_update (ppM, pB, k1 , i1, n);
		
	// Index search (i,i1-1), (k1+1,k)
	i_left = LC_PGG_recursive_search (ppDelta, ppM, pB, ppBasis, i, i1-1, n, delta, nInsertions);
	// Check here if the i_left returned by the call above is less than the minimum i obtained thus far
	// If so, we set i_min = i_left (returned by function above)
	if (i_left < i_min) {
		i_min = i_left;
	}

	LC_PGG_recursive_search (ppDelta, ppM, pB, ppBasis, k1+1, k, n, delta, nInsertions);

	return i_min;
}

// LC-GG index search
inline RR RR_find_insertion_indices_LC_GG (mat_RR ppM, mat_ZZ ppBasis, vec_RR pB, int n, int m, int * pK, int * pJ) 
{
   int k_min = 1;
   int i_min = 0;
   RR C;
   C = 0;

   mat_RR ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   // Lovasz condition definition
   RR inner;
   InnerProduct(inner, ppBasisRR[1], ppBasisRR[1]);
   RR minDelta;
   div(minDelta, inner, pB[0]);
   RR delta_kj;
   delta_kj = 0;

   for (int k=2; k<n; k++) {
      C = 0;
      InnerProduct(C, ppBasisRR[k], ppBasisRR[k]);
      for (int j=0; j<k; j++) {
         // Lovasz condition definition
         div(delta_kj, C, pB[j]);
         if (delta_kj <= minDelta) {
            k_min = k;
            i_min = j;
            minDelta = delta_kj;
         }
         // Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
         C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
      }
   }
   *pK = k_min;
   *pJ = i_min;

   return minDelta;
}

// Pot-Filtered LC-GG Index search
inline RR find_insertion_indices_Pot_F_GG (mat_RR ppM, mat_ZZ ppBasis, vec_RR pB, RR delta, int n, int m, int * pK, int * pJ) 
{
   int k_min = 1;
   int i_min = 0;
   RR C;
   C = 0;

   mat_RR ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   // Lovasz condition definition
   RR inner;
   InnerProduct(inner, ppBasisRR[1], ppBasisRR[1]);
   RR minDelta;
   div(minDelta, inner, pB[0]);
   RR delta_kj;
   delta_kj = 0;
   
   vec_RR P_vec;
   P_vec.SetLength(n);

   RR P; 

   for (int k=2; k<n; k++) {
      // Compute change in Pot for each insertion of vector k and store
      P = 1.0;
      for (int j=k-1; j>=0; j--) {
         RR sum; sum = 0.0;
         for (int r=j; r<k; r++) {
            sum += ppM[k][r]*ppM[k][r]*pB[r];
         }

         RR P_mult; P_mult = 0.0;
         P_mult = (pB[k] + sum) / pB[j];
         P = P * P_mult;
	 P_vec[j] = P; 
      }

      C = 0;
      InnerProduct(C, ppBasisRR[k], ppBasisRR[k]);
      for (int j=0; j<k; j++) {
         // Lovasz condition definition
         div(delta_kj, C, pB[j]);
	 // Check for min delta, but only insert if the reduction in Pot is viable
         if (delta_kj <= minDelta) {
	    if (P_vec[j] < delta) {
               k_min = k;
               i_min = j;
               minDelta = delta_kj;
	    }
         }
         // Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
         C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
      }
   }
   *pK = k_min;
   *pJ = i_min;

   return minDelta;
}

// SS-Filtered LC-GG Index search
inline RR find_insertion_indices_SS_F_GG (mat_RR ppM, mat_ZZ ppBasis, vec_RR pB, RR delta_threshold, RR eta, int n, int m, int * pK, int * pJ) 
{
   int k_min = 1;
   int i_min = 0;
   RR C;
   C = 0;
   RR eta_p;
   eta_p = 1.0 - eta;
   mat_RR ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   // Lovasz condition definition
   RR inner;
   InnerProduct(inner, ppBasisRR[1], ppBasisRR[1]);
   RR minDelta;
   div(minDelta, inner, pB[0]);
   RR delta_kj;
   delta_kj = 0;
   
   RR S_max, S_ik;
   S_max = 0;
   S_ik = 0;

   vec_RR S_vec; 
   S_vec.SetLength(n);

   // For k = 1,...,n and j = k-1,...,1,0, we find the indices k, i where we insert b_k to reduce the value of SS the most
   for (int k=1; k<n; k++) {
      // Initialise to the k-1 case
      RR delta_k_1;
      delta_k_1 = 0;
      delta_k_1 = (pB[k] + (ppM[k][k-1] * ppM[k][k-1] * pB[k-1])) / pB[k-1];

      // Initialise D_j = D_{k-1}
      RR D_j;
      D_j = 0;
      D_j = (ppM[k][k-1] * ppM[k][k-1] * pB[k-1]) + pB[k];

      // Initialise S_ik = S_{i-1,k}
      // After initialising S and D, we can update them in the 'for' loop below for general S_jk and D_j for 0 <= j < k-1
      RR S_ik;
      S_ik = 0;
      S_ik = ppM[k][k-1] * ppM[k][k-1] * (1 - delta_k_1);

      S_vec[k-1] = S_ik;

      for (int l=k-2; l>=0; l--) {
         D_j += ppM[k][l] * ppM[k][l] * pB[l];
         S_ik += ppM[k][l] * ppM[k][l] * pB[l] * ((pB[l]/D_j) - 1);
	 
	 S_vec[l] = S_ik;

      }

      C = 0;
      InnerProduct(C, ppBasisRR[k], ppBasisRR[k]);
      for (int j=0; j<k; j++) {
         // Lovasz condition definition
         div(delta_kj, C, pB[j]);
	 // Check for min delta, but only insert if the reduction in SS is viable
         if (delta_kj <= minDelta) {
	    if (S_vec[j] > eta_p) {
               k_min = k;
               i_min = j;
               minDelta = delta_kj;
	    }
         }
         // Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
         C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
      }
   }
   *pK = k_min;
   *pJ = i_min;

   return minDelta;
}

// Finding insertion indices for block-LC-DeepLLL
inline long double find_insertion_indices_blockwise_LC_DeepLLL (long double ** ppM, int ** ppBasis, long double * pB, int n, int m, int epsilon, int k, int * pJ) 
{
   	int i_min = 0;
   	int index = 0;
   	long double C = 0.0;

   	long double * pVecDouble;
   	pVecDouble = copyVectorToLD(ppBasis[k], m);

	long double inner;
	long double minDelta = 2.0;
	long double delta_kj;

    	C = inner_product_template(pVecDouble, pVecDouble, m);

	if (k <= 2*epsilon) {
       		for (int j=0; j<k; j++) {
 	        	// Lovasz condition definition
			delta_kj = C / pB[j];
            		if (delta_kj <= minDelta) {
               			i_min = j;
               			minDelta = delta_kj;
            		}
            		
			// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            		C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
         	}
      		
	} else {
		index = k - epsilon; 
	 	// Check first epsilon indices
         	for (int j=0; j<epsilon; j++) {
	    		// Lovasz condition definition
			delta_kj = C / pB[j];
            		if (delta_kj <= minDelta) {
               			i_min = j;
               			minDelta = delta_kj;
            		}
            		
			// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            		C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 	}

	 	// Update C for the indices not being checked
         	for (int j=epsilon; j<index; j++) {
	    		C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 	}

	 	// Check indices from k-epsilon+1 to k-1
	 	for (int j=index; j<k; j++) {
            		// Lovasz condition definition
			delta_kj = C / pB[j];
            		if (delta_kj <= minDelta) {
               			i_min = j;
               			minDelta = delta_kj;
            		}
            
			// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            		C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 	}
      	}

   	*pJ = i_min;
        free (pVecDouble);

   	return minDelta;
}

// Finding insertion indices for LC-GG with restricted insertions
inline long double find_insertion_indices_blockwise_LC_GG (long double ** ppM, int ** ppBasis, long double * pB, int n, int m, int epsilon, int * pK, int * pJ) 
{
	int k_min = 1;
   	int i_min = 0;
   	int index = 0;
   	long double C = 0.0;

   	long double ** ppBasisDouble;
   	ppBasisDouble = copyBasisToDouble(ppBasis, m, n);

	long double inner;
	inner = inner_product_template(ppBasisDouble[1], ppBasisDouble[1], m);
	long double minDelta;
	minDelta = inner/pB[0];
	long double delta_kj;;

   	for (int k=2; k<n; k++) {
      		C = 0.0;
      		C = inner_product_template(ppBasisDouble[k], ppBasisDouble[k], m);
      		
		if (k <= 2*epsilon) {
        		for (int j=0; j<k; j++) {
        	    		// Lovasz condition definition
				delta_kj = C / pB[j];
            			if (delta_kj <= minDelta) {
               				k_min = k;
               				i_min = j;
               				minDelta = delta_kj;
            			}
            			
				// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            			C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
         		}
      		
		} else {
	 		index = k - epsilon; 
	 		// Check first epsilon indices
         		for (int j=0; j<epsilon; j++) {
	    			// Lovasz condition definition
				delta_kj = C / pB[j];
            			if (delta_kj <= minDelta) {
               				k_min = k;
               				i_min = j;
               				minDelta = delta_kj;
            			}
            			
				// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            			C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 		}

	 		// Update C for the indices not being checked
         		for (int j=epsilon; j<index; j++) {
	    			C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 		}

	 		// Check indices from k-epsilon+1 to k-1
	 		for (int j=index; j<k; j++) {
            			// Lovasz condition definition
				delta_kj = C / pB[j];
            			if (delta_kj <= minDelta) {
              				k_min = k;
               				i_min = j;
               				minDelta = delta_kj;
            			}
            
				// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            			C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 		}
      		}
   	}
   
	*pK = k_min;
   	*pJ = i_min;

	deleteBasisDouble(ppBasisDouble, n);

   	return minDelta;
}

// Finding insertion indices for LC-GG with restricted insertions
inline RR RR_find_insertion_indices_blockwise_LC_GG (mat_RR ppM, mat_ZZ ppBasis, vec_RR pB, int n, int m, int epsilon, int * pK, int * pJ) 
{
   int k_min = 1;
   int i_min = 0;
   int index = 0;
   RR C;
   C = 0;

   mat_RR ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   // Lovasz condition definition
   RR inner;
   InnerProduct(inner, ppBasisRR[1], ppBasisRR[1]);
   RR minDelta;
   div(minDelta, inner, pB[0]);
   RR delta_kj;
   delta_kj = 0;

   for (int k=2; k<n; k++) {
      C = 0;
      InnerProduct(C, ppBasisRR[k], ppBasisRR[k]);
      if (k <= 2*epsilon) {
         for (int j=0; j<k; j++) {
            // Lovasz condition definition
            div(delta_kj, C, pB[j]);
            if (delta_kj <= minDelta) {
               k_min = k;
               i_min = j;
               minDelta = delta_kj;
            }
            // Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
         }
      } else {
	 index = k - epsilon; 
	 // Check first epsilon indices
         for (int j=0; j<epsilon; j++) {
	    // Lovasz condition definition
            div(delta_kj, C, pB[j]);
            if (delta_kj <= minDelta) {
               k_min = k;
               i_min = j;
               minDelta = delta_kj;
            }
            // Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 }

	 // Update C for the indices not being checked
         for (int j=epsilon; j<index; j++) {
	    C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 }

	 // Check indices from k-epsilon+1 to k-1
	 for (int j=index; j<k; j++) {
            // Lovasz condition definition
            div(delta_kj, C, pB[j]);
            if (delta_kj <= minDelta) {
               k_min = k;
               i_min = j;
               minDelta = delta_kj;
            }
            // Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 }
      }
   }
   *pK = k_min;
   *pJ = i_min;

   return minDelta;
}

// SS-GG Index search
template <typename T1>
inline T1 find_insertion_indices_SS_GG (T1 ** ppM, int ** ppBasis, T1 * pB, int n, int m, int * pK, int * pJ) 
{
   int k_max = 1;
   int i_max = 0;

   T1 S_max = 0.0;
   T1 S_ik = 0.0;

   // For k = 1,...,n and j = k-1,...,1,0, we find the indices k, i where we insert b_k to reduce the value of SS the most
   for (int k=1; k<n; k++) {
      // Initialise to the k-1 case
      T1 delta_k_1 = 0.0;
      delta_k_1 = (pB[k] + (ppM[k][k-1] * ppM[k][k-1] * pB[k-1])) / pB[k-1];

      // Initialise D_j = D_{k-1}
      T1 D_j = 0.0;
      D_j = (ppM[k][k-1] * ppM[k][k-1] * pB[k-1]) + pB[k];

      // Initialise S_ik = S_{i-1,k}
      // After initialising S and D, we can update them in the 'for' loop below for general S_jk and D_j for 0 <= j < k-1
      T1 S_ik = 0.0;
      S_ik = ppM[k][k-1] * ppM[k][k-1] * (1 - delta_k_1);

      // Check if k, k-1 case improves upon current best
      if (S_ik > S_max) {
         S_max = S_ik;
         k_max = k;
         i_max = k-1;
      }

      for (int l=k-2; l>=0; l--) {
         D_j += ppM[k][l] * ppM[k][l] * pB[l];
         S_ik += ppM[k][l] * ppM[k][l] * pB[l] * ((pB[l]/D_j) - 1);

         if (S_ik > S_max) {
            k_max = k;
            i_max = l;
            S_max = S_ik;
         }
      }
   }
   *pK = k_max;
   *pJ = i_max;

   return S_max;
}

// SS-GG Index search
inline RR RR_find_insertion_indices_SS_GG (mat_RR ppM, mat_ZZ ppBasis, vec_RR pB, int n, int m, int * pK, int * pJ) 
{
   int k_max = 1;
   int i_max = 0;

   RR S_max, S_ik; 
   S_max = 0;
   S_ik = 0;

   // For k = 1,...,n and j = k-1,...,1,0, we find the indices k, i where we insert b_k to reduce the value of SS the most
   for (int k=1; k<n; k++) {
      // Initialise to the k-1 case
      RR delta_k_1;
      delta_k_1 = 0;
      delta_k_1 = (pB[k] + (ppM[k][k-1] * ppM[k][k-1] * pB[k-1])) / pB[k-1];

      // Initialise D_j = D_{k-1}
      RR D_j;
      D_j = 0;
      D_j = (ppM[k][k-1] * ppM[k][k-1] * pB[k-1]) + pB[k];

      // Initialise S_ik = S_{i-1,k}
      // After initialising S and D, we can update them in the 'for' loop below for general S_jk and D_j for 0 <= j < k-1
      RR S_ik;
      S_ik = 0;
      S_ik = ppM[k][k-1] * ppM[k][k-1] * (1 - delta_k_1);

      // Check if k, k-1 case improves upon current best
      if (S_ik > S_max) {
         S_max = S_ik;
         k_max = k;
         i_max = k-1;
      }

      for (int l=k-2; l>=0; l--) {
         D_j += ppM[k][l] * ppM[k][l] * pB[l];
         S_ik += ppM[k][l] * ppM[k][l] * pB[l] * ((pB[l]/D_j) - 1);

         if (S_ik > S_max) {
            k_max = k;
            i_max = l;
            S_max = S_ik;
         }
      }
   }
   *pK = k_max;
   *pJ = i_max;

   return S_max;
}

// Pot-GG Index search (using NTL datatypes)
inline RR RR_find_insertion_indices_dynamic_PotLLL (mat_RR ppM, mat_ZZ ppBasis, vec_RR pB, int n, int m, int * pK, int * pJ) 
{
   int k_min = 0;
   int i_min = 0;

   // Initialise P = 1.0, P_min = 1.0
   // We reinitialise P to be 1.0 every time we increment k (in standard PotLLL, for each k, we start with P = 1)
   RR P, P_min;
   P = 1;
   P_min = 1;

   // Initialise P_min to be the k=1, i=0 case
   RR P_init;
   P_init = (pB[1] + ppM[1][0]*ppM[1][0]*pB[0]) / pB[0];

   if (P_init < P_min) {
      P_min = P_init;
      k_min = 1;
      i_min = 0;
   }

   for (int k=2;k<n;k++) {
      P = 1;
      for (int j=k-1; j>=0; j--) {
         RR sum;
	 sum = 0;
         for (int r=j; r<k; r++) {
            sum += ppM[k][r]*ppM[k][r]*pB[r];
         }

      RR P_mult;
      P_mult = (pB[k] + sum) / pB[j];
      P = P * P_mult;
      if (P<P_min) {
         k_min = k;
         i_min = j;
         P_min = P;
      }
      }
   }
   *pK = k_min;
   *pJ = i_min;

   return P_min;
}

// Pot-GG Index search 
inline long double find_insertion_indices_dynamic_PotLLL (long double ** ppM, int ** ppBasis, long double * pB, int n, int m, int * pK, int * pJ) 
{
   int k_min = 0;
   int i_min = 0;

   // Initialise P = 1.0, P_min = 1.0
   // We reinitialise P to be 1.0 every time we increment k (in standard PotLLL, for each k, we start with P = 1)
   long double P = 1.0;
   long double P_min = 1.0;

   // Initialise P_min to be the k=1, i=0 case
   long double P_init = 1.0;
   P_init = (pB[1] + ppM[1][0]*ppM[1][0]*pB[0]) / pB[0];

   if (P_init < P_min) {
      P_min = P_init;
      k_min = 1;
      i_min = 0;
   }
   long double sum = 0.0;
   long double P_mult = 0.0;
   for (int k=2;k<n;k++) {
      P = 1.0;
      for (int j=k-1; j>=0; j--) {
	 sum = 0.0;
         for (int r=j; r<k; r++) {
            sum += ppM[k][r]*ppM[k][r]*pB[r];
         }

         P_mult = (pB[k] + sum) / pB[j];
         P = P * P_mult;
         if (P<P_min) {
            k_min = k;
            i_min = j;
            P_min = P;
         }
      }
   }

   *pK = k_min;
   *pJ = i_min;

   return P_min;
}

// Compute  Root Hermite Factor
inline RR NTLrootHermiteFactor (mat_ZZ ppBasis, int m, int n) 
{
   RR RR_n;
   conv(RR_n, n);
   vec_RR V1;
   V1.SetLength(m);
   mat_RR ppBasisGSO;
   ppBasisGSO = GSO_RR (ppBasis, m, n);
   
   for (int i=0; i<m; i++){
      conv(V1[i], ppBasis[0][i]);
   }
   RR volume;
   volume = 1;

   for (int i=0; i<n; i++) {
      RR gsoVectorNormSquared, gsoVectorNorm;
      InnerProduct(gsoVectorNormSquared, ppBasisGSO[i],ppBasisGSO[i]);
      SqrRoot(gsoVectorNorm, gsoVectorNormSquared);

      mul(volume, volume, gsoVectorNorm);
   }
   RR power;
   inv(power, RR_n);

   RR nthRootVolume;
   pow(nthRootVolume, volume, power);

   RR b1NormSquared;
   InnerProduct (b1NormSquared, V1, V1);

   RR b1Norm;
   SqrRoot(b1Norm, b1NormSquared);

   RR hermiteFactor, RHF;
   div(hermiteFactor, b1Norm, nthRootVolume);

   pow(RHF, hermiteFactor, power);

   return RHF;   
}

// Compute standard deviation for RR datatypes
inline RR StdDevRR(vec_RR RHFList, RR mean, int length)
{
   RR stdDev, sum, variance, term, difference;
   sum = 0;
   variance = 0;
   term = 0;
   difference = 0;

   for (int i=0; i<length; i++) {
      difference = RHFList[i] - mean;
      term = sqr(difference);
      
      sum += term;
   }
   
   variance = sum / length;
   stdDev = SqrRoot(variance);

   return stdDev;
}

// Compute the Root Hermite Factor of a Lattice Basis
inline long double rootHermiteFactor (int ** ppBasis, int m, int n) 
{
   // Compute GSO
   long double ** ppBasisGSO = GSO (ppBasis, m, n);

   // Compute volume of lattice as product of GSO lengths
   long double volume = 1.0;
   
   for (int i=0; i<n; i++) {
      long double gsoVectorNormSquared = inner_product_template (ppBasisGSO[i],ppBasisGSO[i], m);
      long double gsoVectorNorm = sqrtl(gsoVectorNormSquared);
      volume *= gsoVectorNorm;
   }

   // Compute the nth root of the volume using pow function (from math.h)
   long double power = (long double) 1.0 / (long double) n;
   long double nthRootVolume;
   nthRootVolume = pow (volume, power);

   // Compute the norm of b1 using the square root of inner product
   long long b1NormSquared = inner_product_template (ppBasis[0], ppBasis[0], m);
   long double b1Norm; 
   b1Norm = sqrtl(b1NormSquared);

   // Compute the Hermite Factor as gamma = ||b_1|| / Vol(L)^(1/n)
   long double hermiteFactor;
   hermiteFactor = (b1Norm / nthRootVolume);

   // Compute the Root Hermite Factor (RHF) as gamma^(1/n)
   long double RHF;
   RHF = pow(hermiteFactor, power);

   return RHF;
}

// Compute the Orthogonality Defect
inline long double orthogonalityDefect (int ** ppBasis, int m, int n) 
{
   // Compute GSO
   long double ** ppBasisGSO = GSO(ppBasis,m,n);

   // Compute volume of lattice as product of GSO vectors
   long double volumeGSO = 1;
   for (int i=0; i<n; i++) {
      long double gsoVectorNormSquared = inner_product_template (ppBasisGSO[i],ppBasisGSO[i], m);
      long double gsoVectorNorm = sqrt (gsoVectorNormSquared);
      volumeGSO = volumeGSO * gsoVectorNorm;
   }

   // Compute product of basis vectors
   long double volumeBasis = 1;
   for (int j=0; j<n; j++) {
      long double basisVectorNormSquared = inner_product_template (ppBasis[j],ppBasis[j], m);
      long double basisVectorNorm = sqrt (basisVectorNormSquared);
      volumeBasis = volumeBasis * basisVectorNorm;
   }

   // Orthogonality Defect = (Product of basis vectors) / Volume of Lattice
   long double orthDefect = volumeBasis / volumeGSO;
   return orthDefect;
}

// Compute the Orthogonality Defect
inline RR NTLorthogonalityDefect (mat_ZZ ppBasis, int m, int n) 
{
   mat_RR ppBasisGSO;
   ppBasisGSO = GSO_RR(ppBasis, m, n);

   RR volumeGSO;
   volumeGSO = 1;
   for (int i=0; i<n; i++) {
      RR gsoVectorNormSquared;
      InnerProduct(gsoVectorNormSquared, ppBasisGSO[i], ppBasisGSO[i]);
      mul(volumeGSO, volumeGSO, gsoVectorNormSquared); 
   }
   SqrRoot(volumeGSO, volumeGSO);

   RR volumeBasis;
   ZZ volumeBasisSquared;
   volumeBasisSquared = 1;
   for (int j=0; j<n; j++) {
      ZZ vectorNormSquared;
      InnerProduct(vectorNormSquared, ppBasis[j], ppBasis[j]);
      mul(volumeBasisSquared, volumeBasisSquared, vectorNormSquared);
   }
   conv (volumeBasis, volumeBasisSquared);
   SqrRoot(volumeBasis,volumeBasis);

   RR orthDefect;
   div(orthDefect, volumeBasis, volumeGSO);
   return orthDefect;
}

// Compute the Squared Sum of GSO Vectors
inline long double squaredSum (int ** ppBasis, int m, int n) 
{
   // Compute GSO
   long double ** ppBasisGSO = GSO(ppBasis,m,n);

   // Initialise squared sum = 0
   long double squaredSumGSO = 0.0;
   for (int i=0; i<n; i++) {
      long double gsoVectorNormSquared = inner_product_template (ppBasisGSO[i],ppBasisGSO[i], m);
      squaredSumGSO += gsoVectorNormSquared;
   }
   return squaredSumGSO;
}

// Compute the Squared Sum of GSO Vectors
inline RR NTLsquaredSum (mat_ZZ ppBasis, int m, int n) 
{
   mat_RR ppBasisGSO;
   ppBasisGSO = GSO_RR(ppBasis, m, n);
   RR squaredSum;
   squaredSum = 0;

   for (int i=0; i<n; i++) {
      RR gsoNormSquared;
      InnerProduct(gsoNormSquared, ppBasisGSO[i], ppBasisGSO[i]);
      add(squaredSum, squaredSum, gsoNormSquared);
   }

   return squaredSum;
}

// Compute the log of the potential of a basis
inline long double basisPotential (int ** ppBasis, int m, int n) 
{
   // Compute GSO
   long double ** ppBasisGSO = GSO(ppBasis,m,n);
   // Initialise log of potential = 0
   long double potential = 0;
   for (int i=0; i<n; i++) {
      // The potential term for b_i is ||b_i*||^(2(n-i)|| (since indexing from 0)
      // This is the same as (||b_i*||^2)^(n-i)
      // Log of the potential term is (n-i) * log(||b_i*||^2)
      long double gsoVectorNormSquared = inner_product_template (ppBasisGSO[i],ppBasisGSO[i], m);
      long double power = n-i;
      long double potentialTerm = power * log(gsoVectorNormSquared);
      potential += potentialTerm;
   }
   return potential;
}

// Compute the log of the potential of a basis
inline RR NTLbasisPotential (mat_ZZ ppBasis, int m, int n) 
{
   mat_RR ppBasisGSO;
   ppBasisGSO = GSO_RR(ppBasis, m, n);

   RR potential;
   potential = 0;
   for (int i=0; i<n; i++) {
      // The potential term for b_i is ||b_i*||^(2(n-i)|| (since indexing from 0)
      // This is the same as (||b_i*||^2)^(n-i)
      // Log of the potential term is (n-i) * log(||b_i*||^2)
      RR gsoVectorNormSquared;
      InnerProduct(gsoVectorNormSquared, ppBasisGSO[i], ppBasisGSO[i]);
      RR RR_power;
      int power = n-i;
      conv(RR_power, power);
      RR logNormSquared;
      log(logNormSquared, gsoVectorNormSquared);
      RR potentialTerm; 
      mul(potentialTerm, RR_power, logNormSquared);
      potential += potentialTerm;
   }
   return potential;
}

// Check if a Basis is LLL reduced
inline int LLL_check (int ** ppBasis, int m, int n, long double delta_threshold) 
{
   int isLLLReduced = TRUE;

   long double ** ppBasisGSO;

   // Find the GSO of the Basis
   ppBasisGSO = GSO (ppBasis, m, n);

   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   printBasisDoubleAsColumns (ppBasisGSO, m, n);
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   long double * pB;
   pB = (long double *) calloc (n, sizeof(long double));
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      pB[i] = inner_product_template (ppBasisGSO[i], ppBasisGSO[i], m);
      #ifdef __DEBUG__
      printf ("%Lf, ", pB[i]);
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   long double ** ppM;
   // This 2-dimensional nxn array ppDelta will store all the values of delta_ij
   long double ** ppDelta;
   ppM = (long double **) calloc (n, sizeof (long double*));
   ppDelta = (long double **) calloc (n, sizeof (long double*));
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif
   for (int i=1; i<n; i++) {
      ppM[i] = (long double *) calloc (n, sizeof (long double));
      ppDelta[i] = (long double *) calloc (n, sizeof (long double));
      for (int j=0; j<i; j++) {
         ppM[i][j] = (long double)(inner_product_template (ppBasis[i], ppBasisGSO[j], m)) / (long double)pB[j];
         if (ppM[i][j]<-0.5 || ppM[i][j]>0.5) {
            isLLLReduced = FALSE;
            printf ("\nmu[%d][%d] = %Lf", i, j, ppM[i][j]);
            return isLLLReduced;
         }
         ppDelta[i][j] = ((pB[i]/pB[j]) + (ppM[i][j] * ppM[i][j]));
         #ifdef __DEBUG__
         printf ("%Lf(%Lf), ", ppM[i][j], ppDelta[i][j]);
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
      // Check the Lovasz Condition
      if (pB[i] < ((delta_threshold - (ppM[i][i-1] * ppM[i][i-1])) * pB[i-1])) {
         isLLLReduced = FALSE;
         printf ("\nLC[%d][%d]", i, i-1);
         return isLLLReduced;
      }
   }

   deleteBasisDouble(ppM, n);
   deleteBasisDouble(ppDelta, n);
   deleteBasisDouble(ppBasisGSO, n);

   return isLLLReduced;
}

// Check if a Basis is LLL reduced
inline int NTL_LLL_check (mat_ZZ ppBasis, int m, int n, RR delta_threshold) 
{
   int isLLLReduced = TRUE;

   mat_RR ppBasisGSO;

   mat_RR ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   cout << "\n\n\nThe GSO of the basis is:\n" << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i] << endl;
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   // This 2-dimensional nxn array ppDelta will store all the values of delta_ij
   mat_RR ppDelta;
   ppM.SetDims(n,m);
   ppDelta.SetDims(n,m);
   #ifdef __DEBUG__
   cout << "The values of mu_ij are:\n" << ppM << endl;
   #endif
   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
         InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);	 
         div(ppM[i][j], inner, pB[j]);
         if (ppM[i][j]<-0.5 || ppM[i][j]>0.5) {
            isLLLReduced = FALSE;
            cout << "\nmu" << i << j << "  " << ppM[i][j] << endl;
            return isLLLReduced;
         }
	 RR frac, product;
	 div(frac, pB[i], pB[j]);
         mul(product, ppM[i][j], ppM[i][j]);
	 add(ppDelta[i][j], frac, product);
         #ifdef __DEBUG__
         cout << "\nMU_ij  " << ppM[i][j] << "\nDelta_ij  " << ppDelta[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
      // Check the Lovasz Condition
      if (pB[i] < ((delta_threshold - (ppM[i][i-1] * ppM[i][i-1])) * pB[i-1])) {
         isLLLReduced = FALSE;
         printf ("\nLC[%d][%d]", i, i-1);
         return isLLLReduced;
      }
   }
    
   ppM.kill();
   ppDelta.kill();
   ppBasisGSO.kill();
   ppBasisRR.kill();

   return isLLLReduced;
}

// The standard LLL implementation
inline int ** LLL (int ** ppBasis, int m, int n, long double delta, long long int * pSwaps, long long int * pReductions) 
{
   long double ** ppBasisGSO;

   // Find the GSO of the Basis
   ppBasisGSO = GSO (ppBasis, m, n);

   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   printBasisDoubleAsColumns (ppBasisGSO, m, n);
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   long double * pB;
   pB = (long double *) calloc (n, sizeof(long double));
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      pB[i] = inner_product_template (ppBasisGSO[i], ppBasisGSO[i], m);
      #ifdef __DEBUG__
      printf ("%Lf, ", pB[i]);
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   long double ** ppM;
   ppM = (long double **) calloc (n, sizeof (long double*));
   #ifdef __DEBUG__
   printf ("\nThe values of mu_ij are:");
   #endif
   for (int i=1; i<n; i++) {
      ppM[i] = (long double *) calloc (n, sizeof (long double));
      for (int j=0; j<i; j++) {
         ppM[i][j] = (long double)(inner_product_template (ppBasis[i], ppBasisGSO[j], m)) / (long double)pB[j];
         #ifdef __DEBUG__
         printf ("%Lf, ", ppM[i][j]);
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   // The LLL task
   *pSwaps = 0;
   *pReductions = 0;
   int k = 1;
   while (k<n) {

      // Reduce b_k with all the previous basis vectors
      for (int j=k-1; j>=0; j--) {

         int closest_integer_to_mu_kj = CLOSEST_INTEGER (ppM[k][j]);

	 if (closest_integer_to_mu_kj != 0) { 
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            // mu_ki -= CLOSEST_INTEGER (mu_kj)
            for (int i=j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_to_mu_kj * ppM[j][i];
            }

            // mu_kj = mu_kj - [mu_kj]
            ppM[k][j] -= closest_integer_to_mu_kj;
	 }
      }

      // Check Lovasz condition
      if (pB[k] >= ((delta - (ppM[k][k-1] * ppM[k][k-1])) * pB[k-1])) {

         #ifdef __DEBUG__
         printf ("\nLC succeeded for k=%d, delta=%Lf", k, ((pB[k]/pB[k-1]) + (ppM[k][k-1] * ppM[k][k-1])));
         #endif
         k++;

      } else {

         #ifdef __DEBUG__
         printf ("\nLC failed for k=%d, delta=%Lf", k, ((pB[k]/pB[k-1]) + (ppM[k][k-1] * ppM[k][k-1])));
         #endif

         // Swap vectors k, k-1
         int * pTemp;
         pTemp = ppBasis[k];
         ppBasis[k] = ppBasis[k-1];
         ppBasis[k-1] = pTemp;

         // Count the number of swaps
         (*pSwaps)++;

         // Swap the first k-2 elements in rows k, k-1 of matrix M
         for (int c=k-2; c>=0; c--) {
            // Swap rows k-1 and k: ppM[k-1][c] with ppM[k][c]
            long double temp = ppM[k-1][c];
            ppM[k-1][c] = ppM[k][c];
            ppM[k][c] = temp;
         }

         // Updating B_k and B_{k-1}
         long double mu_k_k1 = ppM[k][k-1];
         long double B_k1_updated = pB[k] + (mu_k_k1*mu_k_k1)*pB[k-1];
         ppM[k][k-1] = (mu_k_k1 * pB[k-1]) / B_k1_updated;
         pB[k] = (pB[k-1] * pB[k]) / B_k1_updated;
         pB[k-1] = B_k1_updated;

         // Update columns k-1 and k in M
         for (int r=k+1; r<n; r++) {
            long double T = ppM[r][k];
            ppM[r][k] = ppM[r][k-1] - (mu_k_k1 * T);
            ppM[r][k-1] = T + (ppM[k][k-1] * ppM[r][k]);
         }

         #ifdef __DEBUG__
         printf ("\nThe Basis vectors after the swap are:\n");
         printBasisAsRows (ppBasis, m, n);
         printf ("\nThe updated values of mu_ij are:");
         for (int i=0; i<n; i++) {
            for (int j=0; j<i; j++) {
               printf ("%Lf, ", ppM[i][j]);
            }
            printf ("\n");
         }
         #endif

         // Update k
         k = (1 > (k-1)) ? 1 : k-1;
      }

   }
   
   deleteBasisDouble(ppM, n);
   deleteBasisDouble(ppBasisGSO, n);

   return ppBasis;
}

// The standard LLL implementation using NTL mat_ZZ datatype for basis
inline mat_ZZ LLL_mat_ZZ (mat_ZZ ppBasis, int m, int n, RR delta, long long int * pSwaps, long long int * pReductions) 
{
   mat_RR ppBasisGSO;
   mat_RR ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);
   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   cout <<"\n\n\nThe GSO of the basis is:\n" << ppBasisGSO << endl;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i];
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   ppM.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\nThe values of mu_ij are:");
   #endif
   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
         div(ppM[i][j], inner, pB[j]);
         #ifdef __DEBUG__
         cout << ppM[i][j];
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   // The LLL task
   *pSwaps = 0;
   *pReductions = 0;
   int fullReductions;
   fullReductions = 0;
   int k = 1;
   while (k<n) {


      // Reduce b_k with all the previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_int_RR;
	 ZZ closest_integer_to_mu_kj;
	 round(closest_int_RR, ppM[k][j]);
	 conv(closest_integer_to_mu_kj, closest_int_RR);
         
	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
	    // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);
   
            // mu_ki -= CLOSEST_INTEGER (mu_kj)
            for (int i=j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2i
	       RR term;
	       mul (term, closest_int_RR, ppM[j][i]);
               sub(ppM[k][i], ppM[k][i], term);
            }

            // mu_kj = mu_kj - [mu_kj]
            sub(ppM[k][j], ppM[k][j], closest_int_RR);
         }
      }
      (fullReductions++) ;

      // Check Lovasz condition
      if (pB[k] >= ((delta - (ppM[k][k-1] * ppM[k][k-1])) * pB[k-1])) {

         #ifdef __DEBUG__
         cout << "\nLC succeeded for k="<< k << " delta=" << ((pB[k]/pB[k-1]) + (ppM[k][k-1] * ppM[k][k-1]));
         cout << "\nLC succeeded for k="<< k << " delta=" << ((pB[k]/pB[k-1]) + (ppM[k][k-1] * ppM[k][k-1]));
         #endif
         k++;

      } else {

         #ifdef __DEBUG__
         printf ("\nLC failed for k=" << k << " delta=%Lf" << ((pB[k]/pB[k-1]) + (ppM[k][k-1] * ppM[k][k-1]));
         #endif

         // Swap vectors k, k-1
	 vec_ZZ pTemp;
         pTemp = ppBasis[k];
         ppBasis[k] = ppBasis[k-1];
         ppBasis[k-1] = pTemp;

         // Count the number of swaps
         (*pSwaps)++;

         // Swap the first k-2 elements in rows k, k-1 of matrix M
         for (int c=k-2; c>=0; c--) {
            // Swap rows k-1 and k: ppM[k-1][c] with ppM[k][c]
            RR temp;
	    temp = ppM[k-1][c];
            ppM[k-1][c] = ppM[k][c];
            ppM[k][c] = temp;
         }

         // Updating B_k and B_{k-1}
         RR mu_k_k1;
	 mu_k_k1 = ppM[k][k-1];
         RR B_k1_updated;
	 B_k1_updated = pB[k] + (mu_k_k1*mu_k_k1)*pB[k-1];
         ppM[k][k-1] = (mu_k_k1 * pB[k-1]) / B_k1_updated;
         pB[k] = (pB[k-1] * pB[k]) / B_k1_updated;
         pB[k-1] = B_k1_updated;

         // Update columns k-1 and k in M
         for (int r=k+1; r<n; r++) {
            RR T;
	    T = ppM[r][k];
            ppM[r][k] = ppM[r][k-1] - (mu_k_k1 * T);
            ppM[r][k-1] = T + (ppM[k][k-1] * ppM[r][k]);
         }

         #ifdef __DEBUG__
         printf ("\nThe Basis vectors after the swap are:\n");
         cout << ppBasis;
         printf ("\nThe updated values of mu_ij are:");
         for (int i=0; i<n; i++) {
            for (int j=0; j<i; j++) {
               cout << ppM[i][j];
            }
            printf ("\n");
         }
         #endif
         // Update k
         k = (1 > (k-1)) ? 1 : k-1;
      }

   }

   ppM.kill();
   ppBasisGSO.kill();
   ppBasisRR.kill();

   return ppBasis;
}

// NTL implementation of LLL
inline mat_ZZ LLL_NTL (mat_ZZ BasisNTL, int m, int n, long a, long b, long verbose) 
{
   ZZ det_squared;
   mat_ZZ U;
   long deep = 0;
   double delta = 0.999;
   LLL_RR (BasisNTL, U, delta, deep, 0, verbose);
   //LLL(det_squared, BasisNTL, U, a, b, verbose);

   U.kill();
   return BasisNTL;
}

// NTL implementation of BKZ
inline mat_ZZ BKZ_NTL (mat_ZZ BasisNTL, int m, int n, RR delta, long blocksize, long verbose) 
{
   mat_ZZ U;
   double delta_double;
   conv(delta_double, delta);
   long prune;
   prune = 0;

   // If not FP, use QP, QP1, XD then RR
   BKZ_RR (BasisNTL, U, delta_double, blocksize, prune, 0,  verbose);

   U.kill();
   return BasisNTL;

}

// LC-DeepLLL
inline mat_ZZ LC_DeepLLL_NTL (mat_ZZ ppBasis, int m, int n, RR delta, long long int * pSwaps, long long int * pReductions) 
{
   	mat_RR ppBasisGSO, ppBasisRR;
   	ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);
   	// Find the GSO of the Basis
   	ppBasisGSO = GSO_RR (ppBasis, m, n);

	int * insertion_positions;
        insertion_positions = (int *) calloc (n, sizeof(int));	
	for (int i=0; i<n; i++)
	{
	    insertion_positions[i] = 0;
	}
   	#ifdef __DEBUG__
   	cout << "\n\n\nThe GSO of the basis is:\n" << ppBasisGSO;
   	#endif

   	// This 1-dimensional array pB will store all the values of ||b_i*||^2
   	vec_RR pB;
   	pB.SetLength(n);
   	#ifdef __DEBUG__
   	printf ("\nThe values of ||b_i*||^2 are:\n");
   	#endif
   	for (int i=0; i<n; i++) {
      		InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      		#ifdef __DEBUG__
      		cout << pB[i];
      		#endif
   	}

   	// This 2-dimensional nxn array ppM will store all the values of mu_ij
   	mat_RR ppM;
   	ppM.SetDims(n,m);
   	#ifdef __DEBUG__
   	printf ("\n\nThe values of mu_ij are:");
   	#endif
   	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
	 		RR inner;
	 		InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
         		ppM[i][j] = inner / pB[j];
         		#ifdef __DEBUG__
         		cout << ppM[i][j] << endl; 
         		#endif
      		}
      		#ifdef __DEBUG__
      		printf ("\n");
      		#endif
   	}

   	// The DeepLLL task
   	*pSwaps = 0;
   	*pReductions = 0;
   	int k = 1;
   	
   	while (k<n) {
      		// Keeps track of C - the projected vector used in Lovasz Condition Checks
      		RR C;
      		C = 0;

      		// Reduce b_k with all previous basis vectors
      		for (int j=k-1; j>=0; j--) {
         		RR closest_integer_RR;
	 		round(closest_integer_RR, ppM[k][j]);
         		ZZ closest_integer_to_mu_kj;
	 		conv(closest_integer_to_mu_kj, closest_integer_RR);
         
	 		if (closest_integer_to_mu_kj != 0) {
            			// Count the number of reductions
            			(*pReductions)++;
            			// Reduce b_k with b_j
            			NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            			for (int i = j-1; i>=0; i--) {
               				// By Exercise 17.4.4 from Galbraith v2
               				sub(ppM[k][i], ppM[k][i], closest_integer_RR * ppM[j][i]);
            			}

            			// mu_kj = mu_kj - [mu_kj]
            			sub(ppM[k][j], ppM[k][j], closest_integer_RR);
         		}
      		}

      		ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);
      
      		// Compute C = ||b_k||^2
      		InnerProduct(C, ppBasisRR[k], ppBasisRR[k]);

      		int p = 0;
      		p = k;

      		// Check Lovasz condition for i=0, i=1,...,i=k-1
      		for (int i=0; i<p; i++) {

         		if (C >= (delta * pB[i])) {

            			C = C - (ppM[k][i] * ppM[k][i] * pB[i]);
            
           		} else {

            			// b_i <- b_k, and vectors b_i, ... , b_k-1 are shifted up one position
            			vec_ZZ pTemp;
            			pTemp = ppBasis[k];

            			for (int j=k-1; j>=i; j--) {
               				ppBasis[j+1] = ppBasis[j];
            			}
            			
				ppBasis[i] = pTemp;

            			(*pSwaps)++;
            
	    			// Update the values of ppM, pB

            			deep_LLL_GSO_update_RR (ppM, pB, k , i, n);
				insertion_positions[i] += 1;

            			k = (1 > (i)) ? 0 : i-1;
            			break;
         		}

      		}

      		k++;
   	}
	for (int i=0; i<n; i++)
	{
	    cout << "Insertions in position " << i << ": " << insertion_positions[i] << endl;
	}
  
 	ppM.kill();
   	ppBasisGSO.kill();
   	ppBasisRR.kill();

   	return ppBasis;
}

// LC-DeepLLL (block restricted insertions)
inline mat_ZZ NTL_blockwise_deepLLL (mat_ZZ ppBasis, int m, int n, RR delta, int epsilon, long long int * pSwaps, long long int * pReductions) 
{
   mat_RR ppBasisGSO, ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);
   
   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   cout << "\n\n\nThe GSO of the basis is:\n" << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i];
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   ppM.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif
   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
         ppM[i][j] = inner / pB[j];
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   // The DeepLLL task
   *pSwaps = 0;
   *pReductions = 0;
   int k = 1;

   while (k<n) {

      // Keeps track of C - the projected vector used in Lovasz Condition Checks
      RR C;
      C = 0;

      // Reduce b_k with all previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
	 round(closest_integer_RR, ppM[k][j]);
         ZZ closest_integer_to_mu_kj;
	 conv(closest_integer_to_mu_kj, closest_integer_RR);
         
	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            for (int i = j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               sub(ppM[k][i], ppM[k][i], closest_integer_RR * ppM[j][i]);
            }

            // mu_kj = mu_kj - [mu_kj]
            sub(ppM[k][j], ppM[k][j], closest_integer_RR);
         }
      }
      ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);
      
      // Compute C = ||b_k||^2
      InnerProduct(C, ppBasisRR[k], ppBasisRR[k]);

      int p = 0;
      if (k <= 2*epsilon) {
         p = k;
         for (int i=0; i<p; i++) {
            if (C >= (delta * pB[i])) {

               C = C - (ppM[k][i] * ppM[k][i] * pB[i]);

               } else {
  
               // b_i <- b_k, and vectors b_i, ... , b_k-1 are shifted up one position
               vec_ZZ pTemp;
               pTemp = ppBasis[k];
 
               for (int j=k-1; j>=i; j--) {
                  ppBasis[j+1] = ppBasis[j];
               }
               ppBasis[i] = pTemp;

               // Update the values of ppM, pB (to be done)
               (*pSwaps)++;

               deep_LLL_GSO_update_RR (ppM, pB, k , i, n);
               k = (1 > (i)) ? 0 : i-1;
               break;
            }
         }

      } else {
	 p = epsilon;
         int flag = 0;
	 for (int i=0; i<p; i++) {
            if (C >= (delta * pB[i])) {

               C = C - (ppM[k][i] * ppM[k][i] * pB[i]);

               } else {

               // b_i <- b_k, and vectors b_i, ... , b_k-1 are shifted up one position
               if (i > epsilon) {	
	          if (k-i > epsilon) {		  
	             cout << endl << "ERROR! INSERTION FOR k = " << k << " and i = " << i << endl;
		     exit(i);
	          }
	       }
	       vec_ZZ pTemp;
               pTemp = ppBasis[k];

               for (int j=k-1; j>=i; j--) {
                  ppBasis[j+1] = ppBasis[j];
               }
               ppBasis[i] = pTemp;

               // Update the values of ppM, pB (to be done)
               (*pSwaps)++;

               deep_LLL_GSO_update_RR (ppM, pB, k , i, n);
               k = (1 > (i)) ? 0 : i-1;
	       flag = 1;
               break;
	    }
	 }
	 if (flag == 0) {

	    for (int i=p; i<=k-p; i++) {
	       C = C - (ppM[k][i] * ppM[k][i] * pB[i]);
	    }
	    for (int i=k-p+1; i<k; i++) {
	       
	       if (C >= (delta * pB[i])) {
	          
	          C = C - (ppM[k][i] * ppM[k][i] * pB[i]);

               } else {

                  // b_i <- b_k, and vectors b_i, ... , b_k-1 are shifted up one position
                  if (i > epsilon) {	
	             if (k-i > epsilon) {		  
		        cout << endl << "ERROR! INSERTION FOR k = " << k << " and i = " << i << endl;
		     }
		  }
		  vec_ZZ pTemp;
                  pTemp = ppBasis[k];

                  for (int j=k-1; j>=i; j--) {
                     ppBasis[j+1] = ppBasis[j];
                  }
                  ppBasis[i] = pTemp;

                  // Update the values of ppM, pB (to be done)
                  (*pSwaps)++;

                  deep_LLL_GSO_update_RR (ppM, pB, k , i, n);
                  k = (1 > (i)) ? 0 : i-1;
                  flag = 1;
                  break;
	       }
	    }
         }
      }
      //Increment k
      k++;
   }
   // Check if there are any pairs (k,i) whose \delta_{k,i} < \delta
   // Compute C = ||b_k||^2
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n); 
   RR C; 
   C = 0;
   for (k=0;k<n;k++) {
     InnerProduct(C, ppBasisRR[k], ppBasisRR[k]);

     // Check Lovasz condition for i=0, i=1,...,i=k-1
     for (int i=0; i<k; i++) { 
        C = C - (ppM[k][i] * ppM[k][i] * pB[i]);
     }
   }
  
   ppM.kill();
   ppBasisGSO.kill();
   
   return ppBasis;
}

// Pot-DeepLLL
inline mat_ZZ NTL_PotLLL (mat_ZZ ppBasis, int m, int n, RR delta, long long int * pSwaps, long long int * pReductions) 
{

   mat_RR ppBasisGSO;
   ppBasisGSO.SetDims(n,m);
   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   cout << ppBasisGSO << endl;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i] << endl;
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   ppM.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif

   mat_RR ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
	 div(ppM[i][j], inner, pB[j]);
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   // The PotLLL task
   *pSwaps = 0;
   *pReductions = 0;
   int k = 1;
   RR C;
   ZZ temp;

   while (k<n) {

      // Reduce b_k with all previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
	 round(closest_integer_RR, ppM[k][j]);
         ZZ closest_integer_to_mu_kj;
	 conv(closest_integer_to_mu_kj, closest_integer_RR);
         
	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            for (int i = j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_RR * ppM[j][i];
            }

            // mu_kj = mu_kj - [mu_kj]
            ppM[k][j] -= closest_integer_RR;
         }
      }

      // Initialise P = 1, P_min = 1, k = 0
      // For j = k-1,...,1,0, we find the index where an insertion of vector b_k in position j reduces the potential the most
      RR P, P_min;
      P = 1;
      P_min = 1;
      int l = 1;

      for (int j=k-1; j>=0; j--) {
         RR sum;
	 sum = 0;
         for (int r=j; r<k; r++) {
            sum += ppM[k][r]*ppM[k][r]*pB[r];
         }

         RR P_mult;
         P_mult = (pB[k] + sum) / pB[j];
         P = P * P_mult;
         if (P<P_min) {
            l = j;
            P_min = P;
         }
      }

      #ifdef __DEBUG__
      printf ("\nP_min = %Lf", P_min);
      #endif

      if (delta > P_min) {
         // b_l <- b_k, and vectors b_l, ... , b_k-1 are shifted up one position
         vec_ZZ pTemp;
         pTemp = ppBasis[k];

         for (int j=k-1; j>=l; j--) {
            ppBasis[j+1] = ppBasis[j];
         }
         ppBasis[l] = pTemp;

         // Update values of ppM[i][j] for j<i, pB[i]
         (*pSwaps)++;

         deep_LLL_GSO_update_RR (ppM, pB, k , l, n);

         k = (1 > (l)) ? 1 : l;
         } else {
         k++;
      }
   }
   ppM.kill();
   ppBasisGSO.kill();
   
   return ppBasis;
}

// SS-DeepLLL
inline mat_ZZ NTL_SS_LLL (mat_ZZ ppBasis, int m, int n, RR eta, long long int * pSwaps, long long int * pReductions) 
{
   mat_RR ppBasisGSO, ppBasisRR;
   ppBasisGSO.SetDims(n,m);
   ppBasisRR.SetDims(n,m);

   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);
   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   cout << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i];
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   ppM.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif
   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
         div(ppM[i][j], inner, pB[j]);
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   // The S^2LLL task
   *pSwaps = 0;
   *pReductions = 0;
   int k = 1;
   RR eta_p;
   eta_p = 0.000001;

   while (k<n) {

      RR squared_sum;
      squared_sum = 0;
      for (int i=0; i<n; i++) {
         squared_sum += pB[i];
      }

      // Reduce b_k with all previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
	 round(closest_integer_RR, ppM[k][j]);  
         ZZ closest_integer_to_mu_kj;
	 conv(closest_integer_to_mu_kj, closest_integer_RR);
         
	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            for (int i = j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_RR * ppM[j][i];
            }

            // mu_kj = mu_kj - [mu_kj]
            sub(ppM[k][j], ppM[k][j], closest_integer_RR);
         }
      }

      // For j = k-1,...,1,0, we find the index i where we insert b_k to reduce the value of SS the most
      RR delta_k_1;
      delta_k_1 = 0.0;
      delta_k_1 = (pB[k] + (ppM[k][k-1] * ppM[k][k-1] * pB[k-1])) / pB[k-1];

      // Initialise D_j = D_{k-1}
      RR D_j;
      D_j = (ppM[k][k-1] * ppM[k][k-1] * pB[k-1]) + pB[k];

      // Initialise S_ik = S_{i-1,k}
      // After initialising S and D, we can update them in the 'for' loop below for general S_jk and D_j for 0 <= j < k-1
      RR S_ik;
      S_ik = ppM[k][k-1] * ppM[k][k-1] * (1 - delta_k_1);

      int i = k-1;

      RR S_lk;
      S_lk = S_ik;

      for (int l=k-2; l>=0; l--) {
         D_j += ppM[k][l] * ppM[k][l] * pB[l];
         S_lk += ppM[k][l] * ppM[k][l] * pB[l] * ((pB[l]/D_j) - 1);

         if (S_lk > S_ik) {
            i = l;
            S_ik = S_lk;
         }
      }

      if (S_ik <= eta_p * squared_sum) {
         k++;
      } else {
         // b_l <- b_k, and vectors b_l, ... , b_k-1 are shifted up one position
         vec_ZZ pTemp;
         pTemp = ppBasis[k];

         for (int j=k-1; j>=i; j--) {
            ppBasis[j+1] = ppBasis[j];
         }
         ppBasis[i] = pTemp;

         // Update values of ppM, pB
         (*pSwaps)++;

         deep_LLL_GSO_update_RR (ppM, pB, k , i, n);
         k = (1 >(i)) ? 1 : i;
      }
   }
   ppM.kill();
   ppBasisGSO.kill();
   ppBasisRR.kill();
   
   return ppBasis;
}

// LC-GG
inline mat_ZZ NTL_blockwise_LC_GG (mat_ZZ ppBasis, int m, int n, RR delta_threshold, int epsilon, long long int * pSwaps, long long int * pReductions) 
{
   mat_RR ppBasisGSO, ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   cout << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i],ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i] << endl;
      #endif
   }

   // This 2-dimensional nxn array ppM (ppDelta) will store all the values of mu_ij (delta_ij)
   mat_RR ppM;
   ppM.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif

   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
         div(ppM[i][j], inner, pB[j]);
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl << ppDelta[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }
   
   // The dynamic LLL task
   *pSwaps = 0;
   *pReductions = 0;
   int pReductionIter = 0;
   RR min_delta, old_min_delta;

   // First reduce every vector by all its previous vectors
   for (int k=1; k<n; k++) {

      // Reduce b_k with all the previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
	 round(closest_integer_RR, ppM[k][j]);
         ZZ closest_integer_to_mu_kj;
	 conv(closest_integer_to_mu_kj, closest_integer_RR);

	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            for (int i=j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_RR * ppM[j][i];
            }

            // mu_kj = mu_kj - [mu_kj]
            ppM[k][j] -= closest_integer_RR;
	 }
      }
   }
   (pReductionIter)++;

   // Find the index where there is a potential swap
   int k = 0;
   int i = 0;
  
   old_min_delta = min_delta; 
   min_delta = RR_find_insertion_indices_blockwise_LC_GG(ppM, ppBasis, pB, n, m, epsilon, &k, &i);

    
   while (min_delta < delta_threshold) {
      vec_ZZ pTemp;
      pTemp = ppBasis[k];

      for (int j=k-1; j>=i; j--) {
         ppBasis[j+1] = ppBasis[j];
      }
      ppBasis[i] = pTemp;

      // Update values of ppM, pB, delta
      deep_LLL_GSO_update_RR (ppM, pB, k , i, n);

      // Count the number of swaps
      (*pSwaps)++;


      // Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
      for (int l=i; l<n; l++)
      {
         // Reduce b_l with all the previous basis vectors
         for (int j=l-1; j>=0; j--)
         {
		 
            RR closest_integer_mu_lj_RR;
	    round(closest_integer_mu_lj_RR, ppM[l][j]); 
            ZZ closest_integer_to_mu_lj;
	    conv(closest_integer_to_mu_lj, closest_integer_mu_lj_RR);
             
	    if (closest_integer_to_mu_lj != 0) {
               // Count the number of reductions
               (*pReductions)++;
               // Reduce b_l with b_j
               NTL_reduce (ppBasis[j], ppBasis[l], m, closest_integer_to_mu_lj);

               for (int i=j-1; i>=0; i--) {
                  // By Exercise 17.4.4 from Galbraith v2
                  ppM[l][i] -= closest_integer_mu_lj_RR * ppM[j][i];
               }

               // mu_lj = mu_lj - [mu_lj]
               ppM[l][j] -= closest_integer_mu_lj_RR;
	    }
         }
      }
      (pReductionIter)++;

      #ifdef __DEBUG__
      printf ("\nThe Basis vectors are:\n");
      cout << ppBasis;
      printf ("\nThe values of mu_ij are:");
      for (int i=0; i<n; i++) {
         for (int j=0; j<i; j++) {
            cout << ppM[i][j] << endl << ppDelta[i][j];
         }
         printf ("\n");
      }
      #endif
      old_min_delta = min_delta;
      // Update k
      min_delta = RR_find_insertion_indices_blockwise_LC_GG(ppM, ppBasis, pB, n, m, epsilon, &k, &i);
   }

   ppM.kill();
   ppBasisGSO.kill();
   
   return ppBasis;
}

// Pot-GG
inline mat_ZZ NTL_Pot_GG (mat_ZZ ppBasis, int m, int n, RR delta_threshold, long long int * pSwaps, long long int * pReductions) 
{
   mat_RR ppBasisGSO, ppBasisRR;
   ppBasisGSO.SetDims(n,m);
   ppBasisRR.SetDims(n,m);

   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);
   ppBasisRR = copyNTLBasisToRR (ppBasis, m, n);
   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   cout << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i] << endl;
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   ppM.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif

   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
         InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
	 div(ppM[i][j], inner, pB[j]);
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   *pSwaps = 0;
   *pReductions = 0;
   RR P_min;

   // First reduce every vector by all its previous vectors
   for (int k=1; k<n; k++) {

      // Reduce b_k with all the previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
	 round(closest_integer_RR, ppM[k][j]);
	 ZZ closest_integer_to_mu_kj;
	 conv(closest_integer_to_mu_kj, closest_integer_RR);

	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            for (int i=j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_RR * ppM[j][i];
            }

            // mu_kj = mu_kj - [mu_kj]
            ppM[k][j] -= closest_integer_RR;
	 }
      }
   }

   // Find the index where there is a potential swap
   int k = 0;
   int i = 0;

   P_min = RR_find_insertion_indices_dynamic_PotLLL (ppM, ppBasis, pB, n, m, &k, &i);
   cout << "k = " << k << "  i = " << i << endl;
   #ifdef __DEBUG__
   printf ("\n\n\nPotLLL_dynamic_delta: P_min=%Lf\n\n\n", P_min);
   #endif

   while (P_min < delta_threshold) {

      #ifdef __DEBUG__
      printf ("Insertion at k=%d, i=%d for delta=%Lf\n", k, i, P_min);
      #endif

      vec_ZZ pTemp;
      pTemp = ppBasis[k];

      for (int j=k-1; j>=i; j--) {
         ppBasis[j+1] = ppBasis[j];
      }
      ppBasis[i] = pTemp;

         // Update values of ppM, pB
      deep_LLL_GSO_update_RR (ppM, pB, k , i, n);


      // Count the number of swaps
      (*pSwaps)++;

      // Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
      for (int l=i; l<n; l++)
      {
         // Reduce b_l with all the previous basis vectors
         for (int j=l-1; j>=0; j--)
         {
            RR closest_int_mu_lj_RR;
	    round(closest_int_mu_lj_RR, ppM[l][j]);
            ZZ closest_integer_to_mu_lj;
	    conv(closest_integer_to_mu_lj, closest_int_mu_lj_RR);

	    if (closest_integer_to_mu_lj != 0) {
               // Count the number of reductions
               (*pReductions)++;
               // Reduce b_l with b_j
               NTL_reduce (ppBasis[j], ppBasis[l], m, closest_integer_to_mu_lj);

               for (int i=j-1; i>=0; i--) {
                  // By Exercise 17.4.4 from Galbraith v2
                  ppM[l][i] -= closest_int_mu_lj_RR * ppM[j][i];
               }

               // mu_lj = mu_lj - [mu_lj]
               ppM[l][j] -= closest_int_mu_lj_RR;
	    }
         }
      }

      #ifdef __DEBUG__
      printf ("\nThe Basis vectors are:\n");
      cout << ppBasis;
      printf ("\nThe values of mu_ij are:");
      for (int i=0; i<n; i++) {
         for (int j=0; j<i; j++) {
            cout << ppM[i][j] << endl; 
         }
         printf ("\n");
      }
      #endif
      // Update k
      P_min = RR_find_insertion_indices_dynamic_PotLLL (ppM, ppBasis, pB, n, m, &k, &i);
      cout << "k = " << k << "  i = " << i << endl;
   }

   ppM.kill();
   ppBasisGSO.kill();
   
   return ppBasis;
}

// SS-GG
mat_ZZ NTL_SS_GG (mat_ZZ ppBasis, int m, int n, RR eta, long long int * pSwaps, long long int * pReductions) 
{
   RR eta_p;
   eta_p = 1 - eta;

   mat_RR ppBasisGSO, ppBasisRR;
   ppBasisGSO.SetDims(n,m);
   ppBasisRR.SetDims(n,m);
   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   cout << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif

   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i], ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i] << endl;
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   ppM.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif
   
   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
         div(ppM[i][j], inner, pB[j]);
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl; 
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   RR S_max, SS_sub_S_max;
   S_max = 0;
   RR squared_sum;
   squared_sum = 0;
   for (int i=0; i<n; i++) {
      squared_sum += pB[i];
   }
   SS_sub_S_max = squared_sum;
  
   cout << "SQUARED SUM ON INPUT: " << squared_sum << endl;

   // First reduce every vector by all its previous vectors
   for (int k=1; k<n; k++) {

      // Reduce b_k with all the previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
         round(closest_integer_RR, ppM[k][j]);
	 ZZ closest_integer_to_mu_kj;
         conv(closest_integer_to_mu_kj, closest_integer_RR);	 

	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);
		
            for (int i=j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_RR * ppM[j][i];
            }

            // mu_kj = mu_kj - [mu_kj]
            ppM[k][j] -= closest_integer_RR;
         }
      }
   }
   // Find the index where there is a potential swap
   int k = 0;
   int i = 0;
   S_max = RR_find_insertion_indices_SS_GG (ppM, ppBasis, pB, n, m, &k, &i);
   cout << "k = " << k << " i = " << i << endl;

   RR SS_recomputed;
   while (S_max > eta_p * squared_sum) {
      vec_ZZ pTemp;
      pTemp = ppBasis[k];

      for (int j=k-1; j>=i; j--) {
         ppBasis[j+1] = ppBasis[j];
      }
      ppBasis[i] = pTemp;

      // Update values of ppM, pB
      deep_LLL_GSO_update_RR (ppM, pB, k , i, n);

      // Update squared_sum
      SS_sub_S_max -= S_max;
      SS_recomputed = 0;
      for (int i=0; i<n; i++) {
         SS_recomputed += pB[i];
      }
      squared_sum = SS_recomputed;
      // Count the number of swaps
      (*pSwaps)++;

      // Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
      for (int l=i; l<n; l++)
      {
         // Reduce b_l with all the previous basis vectors
         for (int j=l-1; j>=0; j--) {
            RR closest_integer_mu_lj_RR;
	    round(closest_integer_mu_lj_RR, ppM[l][j]);
            ZZ closest_integer_to_mu_lj;
            conv(closest_integer_to_mu_lj, closest_integer_mu_lj_RR);

	    if (closest_integer_to_mu_lj != 0) {
               // Count the number of reductions
               (*pReductions)++;
               // Reduce b_l with b_j
               NTL_reduce (ppBasis[j], ppBasis[l], m, closest_integer_to_mu_lj);

               for (int i=j-1; i>=0; i--) {
                  // By Exercise 17.4.4 from Galbraith v2
                  ppM[l][i] -= closest_integer_mu_lj_RR * ppM[j][i];
               }

               // mu_lj = mu_lj - [mu_lj]
               ppM[l][j] -= closest_integer_mu_lj_RR;
	    }
         }
      }

      #ifdef __DEBUG__
      printf ("\nThe Basis vectors are:\n");
      cout << ppBasis << endl;
      printf ("\nThe values of mu_ij are:");
      for (int i=0; i<n; i++) {
         for (int j=0; j<i; j++) {
            cout << ppM[i][j] << endl;
         }
         printf ("\n");
      }
      #endif

      // Update k
      S_max = RR_find_insertion_indices_SS_GG (ppM, ppBasis, pB, n, m, &k, &i);
   cout << "k = " << k << " i = " << i << endl;
   }

   ppM.kill();
   ppBasisGSO.kill();
   ppBasisRR.kill();
   
   return ppBasis;
}

// Pot-Filtered LC-GG
inline mat_ZZ pot_Filtered_LC_GG (mat_ZZ ppBasis, int m, int n, RR delta_threshold, long long int * pSwaps, long long int * pReductions) 
{
   mat_RR ppBasisGSO, ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   cout << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i],ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i] << endl;
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij 
   mat_RR ppM;
   ppM.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif

   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
         div(ppM[i][j], inner, pB[j]);
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }

   *pSwaps = 0;
   *pReductions = 0;
   int pReductionIter = 0;
   RR min_delta, old_min_delta;

   // First reduce every vector by all its previous vectors
   for (int k=1; k<n; k++) {

      // Reduce b_k with all the previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
	 round(closest_integer_RR, ppM[k][j]);
         ZZ closest_integer_to_mu_kj;
	 conv(closest_integer_to_mu_kj, closest_integer_RR);

	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            for (int i=j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_RR * ppM[j][i];
            }

            // mu_kj = mu_kj - [mu_kj]
            ppM[k][j] -= closest_integer_RR;
	 }
      }
   }
   (pReductionIter)++;

   // Find the index where there is a potential swap
   int k = 0;
   int i = 0;
  
   old_min_delta = min_delta; 
   min_delta = find_insertion_indices_Pot_F_GG (ppM, ppBasis, pB, delta_threshold, n, m, &k, &i);
    
   while (min_delta < delta_threshold) {
      vec_ZZ pTemp;
      pTemp = ppBasis[k];

      for (int j=k-1; j>=i; j--) {
         ppBasis[j+1] = ppBasis[j];
      }
      ppBasis[i] = pTemp;

      // Update values of ppM, pB, delta
      deep_LLL_GSO_update_RR (ppM, pB, k , i, n);

      // Count the number of swaps
      (*pSwaps)++;


      // Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
      for (int l=i; l<n; l++)
      {
         // Reduce b_l with all the previous basis vectors
         for (int j=l-1; j>=0; j--)
         {
		 
            RR closest_integer_mu_lj_RR;
	    round(closest_integer_mu_lj_RR, ppM[l][j]); 
            ZZ closest_integer_to_mu_lj;
	    conv(closest_integer_to_mu_lj, closest_integer_mu_lj_RR);
             
	    if (closest_integer_to_mu_lj != 0) {
               // Count the number of reductions
               (*pReductions)++;
               // Reduce b_l with b_j
               NTL_reduce (ppBasis[j], ppBasis[l], m, closest_integer_to_mu_lj);

               for (int i=j-1; i>=0; i--) {
                  // By Exercise 17.4.4 from Galbraith v2
                  ppM[l][i] -= closest_integer_mu_lj_RR * ppM[j][i];
               }

               // mu_lj = mu_lj - [mu_lj]
               ppM[l][j] -= closest_integer_mu_lj_RR;
	    }
         }
      }
      (pReductionIter)++;

      #ifdef __DEBUG__
      printf ("\nThe Basis vectors are:\n");
      cout << ppBasis;
      printf ("\nThe values of mu_ij are:");
      for (int i=0; i<n; i++) {
         for (int j=0; j<i; j++) {
            cout << ppM[i][j] << endl; 
         }
         printf ("\n");
      }
      #endif
      old_min_delta = min_delta;
      // Update k
      min_delta = find_insertion_indices_Pot_F_GG (ppM, ppBasis, pB, delta_threshold, n, m, &k, &i);
   }

   ppM.kill();
   ppBasisGSO.kill();
   
   return ppBasis;
}

// SS-Filtered LC-GG
inline mat_ZZ SS_Filtered_LC_GG (mat_ZZ ppBasis, int m, int n, RR delta_threshold, RR eta, long long int * pSwaps, long long int * pReductions) 
{   
   mat_RR ppBasisGSO, ppBasisRR;
   ppBasisRR = copyNTLBasisToRR(ppBasis, m, n);

   // Find the GSO of the Basis
   ppBasisGSO = GSO_RR (ppBasis, m, n);

   #ifdef __DEBUG__
   printf ("\n\n\nThe GSO of the basis is:\n");
   cout << ppBasisGSO;
   #endif

   // This 1-dimensional array pB will store all the values of ||b_i*||^2
   vec_RR pB;
   pB.SetLength(n);
   #ifdef __DEBUG__
   printf ("\nThe values of ||b_i*||^2 are:\n");
   #endif
   for (int i=0; i<n; i++) {
      InnerProduct(pB[i], ppBasisGSO[i],ppBasisGSO[i]);
      #ifdef __DEBUG__
      cout << pB[i] << endl;
      #endif
   }

   // This 2-dimensional nxn array ppM will store all the values of mu_ij
   mat_RR ppM;
   ppM.SetDims(n,m);
   #ifdef __DEBUG__
   printf ("\n\nThe values of mu_ij are:");
   #endif

   for (int i=1; i<n; i++) {
      for (int j=0; j<i; j++) {
	 RR inner;
	 InnerProduct(inner, ppBasisRR[i], ppBasisGSO[j]);
         div(ppM[i][j], inner, pB[j]);
         #ifdef __DEBUG__
         cout << ppM[i][j] << endl;
         #endif
      }
      #ifdef __DEBUG__
      printf ("\n");
      #endif
   }
   
   *pSwaps = 0;
   *pReductions = 0;
   //*pMinDeltaReductions = 0;
   int pReductionIter = 0;
   RR min_delta, old_min_delta;

   // First reduce every vector by all its previous vectors
   for (int k=1; k<n; k++) {

      // Reduce b_k with all the previous basis vectors
      for (int j=k-1; j>=0; j--) {
         RR closest_integer_RR;
	 round(closest_integer_RR, ppM[k][j]);
         ZZ closest_integer_to_mu_kj;
	 conv(closest_integer_to_mu_kj, closest_integer_RR);

	 if (closest_integer_to_mu_kj != 0) {
            // Count the number of reductions
            (*pReductions)++;
            // Reduce b_k with b_j
            NTL_reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            for (int i=j-1; i>=0; i--) {
               // By Exercise 17.4.4 from Galbraith v2
               ppM[k][i] -= closest_integer_RR * ppM[j][i];
            }

            // mu_kj = mu_kj - [mu_kj]
            ppM[k][j] -= closest_integer_RR;
	 }
      }
   }
   (pReductionIter)++;

   // Find the index where there is a potential swap
   int k = 0;
   int i = 0;
  
   old_min_delta = min_delta; 
   min_delta = find_insertion_indices_SS_F_GG (ppM, ppBasis, pB, delta_threshold, eta, n, m, &k, &i);

   while (min_delta < delta_threshold) {
      vec_ZZ pTemp;
      pTemp = ppBasis[k];

      for (int j=k-1; j>=i; j--) {
         ppBasis[j+1] = ppBasis[j];
      }
      ppBasis[i] = pTemp;

      // Update values of ppM, pB, delta
      deep_LLL_GSO_update_RR (ppM, pB, k , i, n);

      // Count the number of swaps
      (*pSwaps)++;

      // Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
      for (int l=i; l<n; l++)
      {
         // Reduce b_l with all the previous basis vectors
         for (int j=l-1; j>=0; j--)
         {
		 
            RR closest_integer_mu_lj_RR;
	    round(closest_integer_mu_lj_RR, ppM[l][j]); 
            ZZ closest_integer_to_mu_lj;
	    conv(closest_integer_to_mu_lj, closest_integer_mu_lj_RR);
             
	    if (closest_integer_to_mu_lj != 0) {
               // Count the number of reductions
               (*pReductions)++;
               // Reduce b_l with b_j
               NTL_reduce (ppBasis[j], ppBasis[l], m, closest_integer_to_mu_lj);

               for (int i=j-1; i>=0; i--) {
                  // By Exercise 17.4.4 from Galbraith v2
                  ppM[l][i] -= closest_integer_mu_lj_RR * ppM[j][i];
               }

               // mu_lj = mu_lj - [mu_lj]
               ppM[l][j] -= closest_integer_mu_lj_RR;
	    }
         }
      }
      (pReductionIter)++;

      #ifdef __DEBUG__
      printf ("\nThe Basis vectors are:\n");
      cout << ppBasis;
      printf ("\nThe values of mu_ij are:");
      for (int i=0; i<n; i++) {
         for (int j=0; j<i; j++) {
            cout << ppM[i][j] << endl; 
         }
         printf ("\n");
      }
      #endif
      old_min_delta = min_delta;
      // Update k
      min_delta = find_insertion_indices_SS_F_GG (ppM, ppBasis, pB, delta_threshold, eta, n, m, &k, &i);
   }

   ppM.kill();
   ppBasisGSO.kill();
   
   return ppBasis;
}

//----------------------------------------------------------
// ALGORITHMS USING INT AND LONG DOUBLE/TEMPLATE VERSION 
//----------------------------------------------------------

// Pot-GG
inline int ** Pot_GG_std (int ** ppBasis, int m, int n, long double delta_threshold, long long int * pSwaps, long long int * pReductions) 
{
	long double ** ppBasisGSO;
	ppBasisGSO = GSO (ppBasis, m, n);
	long double ** ppBasisDouble;
	ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
	long double ** ppM;
	ppM = (long double **) calloc (n, sizeof(long double *));
   	for (int i=0; i<n; i++) {
   		ppM[i] = (long double *) calloc (m, sizeof(long double));
	}


   	#ifdef __DEBUG__
 	printf ("\n\n\nThe GSO of the basis is:\n");
 	cout << ppBasisGSO;
 	#endif

   	// This 1-dimensional array pB will store all the values of ||b_i*||^2
   	long double * pB;
   	pB = (long double *) calloc (m, sizeof(long double));

   	#ifdef __DEBUG__
   	printf ("\nThe values of ||b_i*||^2 are:\n");
   	#endif
   	
	for (int i=0; i<n; i++) {
     		pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
      		#ifdef __DEBUG__
      		cout << pB[i] << endl;
      		#endif
	}

   	#ifdef __DEBUG__
  	 printf ("\n\nThe values of mu_ij are:");
   	#endif

  	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
	 		long double inner;
         		inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
	 		ppM[i][j] = inner / pB[j];
         		#ifdef __DEBUG__
         		cout << ppM[i][j] << endl;
         		#endif
      		}
      		#ifdef __DEBUG__
      		printf ("\n");
      		#endif
   	}

   	long double P_min = 0;
        int closest_integer = 0;
	long double closest_integer_Lf = 0.0;

   	// First reduce every vector by all its previous vectors
   	for (int k=1; k<n; k++) {

   	   // Reduce b_k with all the previous basis vectors
   		for (int j=k-1; j>=0; j--) {
	 		closest_integer = CLOSEST_INTEGER(ppM[k][j]);
        		closest_integer_Lf = (long double) closest_integer; 

		 	if (closest_integer != 0) {
         			// Count the number of reductions
            			(*pReductions)++;
            			// Reduce b_k with b_j
            			reduce (ppBasis[j], ppBasis[k], m, closest_integer);

            			for (int i=j-1; i>=0; i--) {
               				// By Exercise 17.4.4 from Galbraith v2
               				ppM[k][i] -= closest_integer_Lf * ppM[j][i];
            			}

            			// mu_kj = mu_kj - [mu_kj]
            			ppM[k][j] -= closest_integer_Lf;
	 		}
      		}
   	}

   	// Find the index where there is a potential swap
   	int k = 0;
   	int i = 0;

	int * pTemp;

   	P_min = find_insertion_indices_dynamic_PotLLL (ppM, ppBasis, pB, n, m, &k, &i);
   	#ifdef __DEBUG__
   	printf ("\n\n\nPotLLL_dynamic_delta: P_min=%Lf\n\n\n", P_min);
   	#endif

   	while (P_min < delta_threshold) {

      		#ifdef __DEBUG__
      		printf ("Insertion at k=%d, i=%d for delta=%Lf\n", k, i, P_min);
      		#endif

      		pTemp = ppBasis[k];

      		for (int j=k-1; j>=i; j--) {
         		ppBasis[j+1] = ppBasis[j];
      		}
      		
		ppBasis[i] = pTemp;

         	// Update values of ppM, pB
      		deep_LLL_GSO_update (ppM, pB, k , i, n);


      		// Count the number of swaps
      		(*pSwaps)++;


      		// Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
      		
		for (int l=i; l<n; l++) {
			// Reduce b_l with all the previous basis vectors
			for (int j=l-1; j>=0; j--) {
	    			closest_integer = CLOSEST_INTEGER(ppM[l][j]);
        			closest_integer_Lf = (long double) closest_integer; 

	    			if (closest_integer != 0) {
               				// Count the number of reductions
               				(*pReductions)++;
               				// Reduce b_l with b_j
               				reduce (ppBasis[j], ppBasis[l], m, closest_integer);

               				for (int i=j-1; i>=0; i--) {
                  				// By Exercise 17.4.4 from Galbraith v2
                  				ppM[l][i] -= closest_integer_Lf * ppM[j][i];
               				}

               				// mu_lj = mu_lj - [mu_lj]
               				ppM[l][j] -= closest_integer_Lf;
	    			}
         		}
      		}

      		#ifdef __DEBUG__
      		printf ("\nThe Basis vectors are:\n");
      		cout << ppBasis;
      		printf ("\nThe values of mu_ij are:");
      		for (int i=0; i<n; i++) {
         		for (int j=0; j<i; j++) {
            			cout << ppM[i][j] << endl; 
         		}
         		
			printf ("\n");
      		}
      		#endif
      
		// Update k
      		P_min = find_insertion_indices_dynamic_PotLLL (ppM, ppBasis, pB, n, m, &k, &i);
   	}
	deleteBasisDouble (ppM, n);
	deleteBasisDouble (ppBasisGSO, n);
	deleteBasisDouble (ppBasisDouble, n);
   
   	return ppBasis;
}

// SS-GG
inline int ** SS_GG_std (int ** ppBasis, int m, int n, long double eta, long long int * pSwaps, long long int * pReductions) 
{
   	long double eta_p = 0;
   	eta_p = 1 - eta;

	long double ** ppBasisGSO;
	ppBasisGSO = GSO (ppBasis, m, n);
	long double ** ppBasisDouble;
	ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
	long double ** ppM;
	ppM = (long double **) calloc (n, sizeof(long double *));
   	for (int i=0; i<n; i++) {
   		ppM[i] = (long double *) calloc (m, sizeof(long double));
	}

   	#ifdef __DEBUG__
   	printf ("\n\n\nThe GSO of the basis is:\n");
	printBasisDoubleAsColumns(ppBasisGSO, m, n);
   	#endif

   	// This 1-dimensional array pB will store all the values of ||b_i*||^2
   	long double * pB;
   	pB = (long double *) calloc (m, sizeof(long double));
   	
   	#ifdef __DEBUG__
   	printf ("\nThe values of ||b_i*||^2 are:\n");
   	#endif

   	for (int i=0; i<n; i++) {
     		pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
      		#ifdef __DEBUG__
      		cout << pB[i] << endl;
      		#endif
   	}

   	#ifdef __DEBUG__
   	printf ("\n\nThe values of mu_ij are:");
   	#endif
   
	long double inner = 0;
	int int_inner = 0;

   	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
			inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
         		ppM[i][j] = inner / pB[j];

         		#ifdef __DEBUG__
         		cout << ppM[i][j] << " "; 
         		#endif
      		}
      	#ifdef __DEBUG__
     	printf ("\n");
      	#endif
   	}

   	long double S_max = 0;
        long double SS_sub_S_max = 0;
   	long double squared_sum = 0;
   	
	long double bound = 0;

	for (int i=0; i<n; i++) {
      		squared_sum += pB[i];
   	}

	SS_sub_S_max = squared_sum;
  
	int closest_integer = 0;	
	long double closest_integer_Lf = 0.0;	
	// First reduce every vector by all its previous vectors
   	for (int k=1; k<n; k++) {

      		// Reduce b_k with all the previous basis vectors
      		for (int j=k-1; j>=0; j--) {
         		closest_integer = CLOSEST_INTEGER(ppM[k][j]);
			closest_integer_Lf = (long double) closest_integer;
	 		if (closest_integer != 0) {
            			// Count the number of reductions
            			(*pReductions)++;
            			// Reduce b_k with b_j
            			reduce (ppBasis[j], ppBasis[k], m, closest_integer);
		
            			for (int i=j-1; i>=0; i--) {
               				// By Exercise 17.4.4 from Galbraith v2
               				ppM[k][i] -= closest_integer_Lf * ppM[j][i];
            			}

            			// mu_kj = mu_kj - [mu_kj]
            			ppM[k][j] -= closest_integer_Lf;
         		}
      		}
   	}
   
	// Find the index where there is a potential swap
 	int k = 0;
   	int i = 0;
   
	S_max = find_insertion_indices_SS_GG (ppM, ppBasis, pB, n, m, &k, &i);

	bound = eta_p * squared_sum;
	long double SS_recomputed = 0;

      	int * pTemp;

   	while (S_max > bound) {

      		pTemp = ppBasis[k];

      		for (int j=k-1; j>=i; j--) {
         		ppBasis[j+1] = ppBasis[j];
      		}
      		
		ppBasis[i] = pTemp;

      		// Update values of ppM, pB
      		deep_LLL_GSO_update (ppM, pB, k , i, n);

      		// Update squared_sum
      		SS_sub_S_max -= S_max;
      		SS_recomputed = 0;
      		for (int i=0; i<n; i++) {
        		SS_recomputed += pB[i];
      		}
      		squared_sum = SS_recomputed;
      	
		// Count the number of swaps
      		(*pSwaps)++;
	
      		// Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
      		for (int l=i; l<n; l++) {
        	 // Reduce b_l with all the previous basis vectors
        		for (int j=l-1; j>=0; j--) {
	    			closest_integer = CLOSEST_INTEGER(ppM[l][j]);
				closest_integer_Lf = (long double) closest_integer;
	    			if (closest_integer != 0) {
               				// Count the number of reductions
               				(*pReductions)++;
               				// Reduce b_l with b_j
               				reduce (ppBasis[j], ppBasis[l], m, closest_integer);

               				for (int i=j-1; i>=0; i--) {
                  				// By Exercise 17.4.4 from Galbraith v2
                  				ppM[l][i] -= closest_integer_Lf * ppM[j][i];
               				}
	
        	       			// mu_lj = mu_lj - [mu_lj]
               				ppM[l][j] -= closest_integer_Lf;
	    			}
         		}
      		}
      
		#ifdef __DEBUG__
      		printf ("\nThe Basis vectors are:\n");
      		cout << ppBasis << endl;
      		printf ("\nThe values of mu_ij are:");
      		for (int i=0; i<n; i++) {
         		for (int j=0; j<i; j++) {
            			cout << ppM[i][j] << endl;
         		}
         		printf ("\n");
      		}
      		#endif

      		// Update k
      		S_max = find_insertion_indices_SS_GG (ppM, ppBasis, pB, n, m, &k, &i);
		bound = eta_p * squared_sum;
   	}
   	deleteBasisDouble(ppM, n);
   	deleteBasisDouble(ppBasisGSO, n);
   	deleteBasisDouble(ppBasisDouble, n);
  
   	return ppBasis;
}

// SS-DeepLLL
inline int ** SS_LLL_std (int ** ppBasis, int m, int n, long double eta, long long int * pSwaps, long long int * pReductions) 
{
   	long double eta_p = 0;
   	eta_p = 1 - eta;

	long double ** ppBasisGSO;
	ppBasisGSO = GSO (ppBasis, m, n);
	long double ** ppBasisDouble;
	ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
	long double ** ppM;
	ppM = (long double **) calloc (n, sizeof(long double *));
   	for (int i=0; i<n; i++) {
   		ppM[i] = (long double *) calloc (m, sizeof(long double));
	}

   	#ifdef __DEBUG__
   	printf ("\n\n\nThe GSO of the basis is:\n");
	printBasisDoubleAsColumns(ppBasisGSO, m, n);
   	#endif

   	// This 1-dimensional array pB will store all the values of ||b_i*||^2
   	long double * pB;
   	pB = (long double *) calloc (m, sizeof(long double));
   	
   	#ifdef __DEBUG__
   	printf ("\nThe values of ||b_i*||^2 are:\n");
   	#endif

   	for (int i=0; i<n; i++) {
     		pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
      		#ifdef __DEBUG__
      		cout << pB[i] << endl;
      		#endif
   	}

   	#ifdef __DEBUG__
   	printf ("\n\nThe values of mu_ij are:");
   	#endif
   
	long double inner = 0;
	long long_inner = 0;

   	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
			long_inner = inner_product_long(ppBasis[i], ppBasis[j], m); 
	 		inner = (long double) long_inner; 
         		ppM[i][j] = inner / pB[j];

         		#ifdef __DEBUG__
         		cout << ppM[i][j] << " "; 
         		#endif
      		}
      	#ifdef __DEBUG__
     	printf ("\n");
      	#endif
   	}

	// The S^2LLL task
   	int k = 1;
	int closest_integer = 0;

	long double delta_k_1 = 0.0;
	long double D_j = 0.0;
	long double S_ik = 0.0;
	long double S_lk = 0.0;

	int * pTemp;

   	while (k<n) {
      		long double squared_sum = 0;
      		
		for (int i=0; i<n; i++) {
         		squared_sum += pB[i];
      		}

      		// Reduce b_k with all previous basis vectors
      		for (int j=k-1; j>=0; j--) {
	 		closest_integer = CLOSEST_INTEGER(ppM[k][j]);  
         
	 		if (closest_integer != 0) {
            			// Count the number of reductions
            			(*pReductions)++;
            			// Reduce b_k with b_j
            			reduce (ppBasis[j], ppBasis[k], m, closest_integer);

            			for (int i = j-1; i>=0; i--) {
               				// By Exercise 17.4.4 from Galbraith v2
               				ppM[k][i] -= closest_integer * ppM[j][i];
            			}

            			// mu_kj = mu_kj - [mu_kj]
            			ppM[k][j] -= closest_integer;
         		}
     	 	}

      		// For j = k-1,...,1,0, we find the index i where we insert b_k to reduce the value of SS the most
      		delta_k_1 = (pB[k] + (ppM[k][k-1] * ppM[k][k-1] * pB[k-1])) / pB[k-1];

      		// Initialise D_j = D_{k-1}
      		D_j = (ppM[k][k-1] * ppM[k][k-1] * pB[k-1]) + pB[k];

      		// Initialise S_ik = S_{i-1,k}
      		// After initialising S and D, we can update them in the 'for' loop below for general S_jk and D_j for 0 <= j < k-1
      		S_ik = ppM[k][k-1] * ppM[k][k-1] * (1 - delta_k_1);

      		int i = k-1;

      		S_lk = S_ik;

      		for (int l=k-2; l>=0; l--) {
        		D_j += ppM[k][l] * ppM[k][l] * pB[l];
         		S_lk += ppM[k][l] * ppM[k][l] * pB[l] * ((pB[l]/D_j) - 1);

         		if (S_lk > S_ik) {
            			i = l;
            			S_ik = S_lk;
         		}
      		}

      		if (S_ik <= eta_p * squared_sum) {
         		k++;
      		} else {
         		// b_l <- b_k, and vectors b_l, ... , b_k-1 are shifted up one position
         		pTemp = ppBasis[k];

         		for (int j=k-1; j>=i; j--) {
            			ppBasis[j+1] = ppBasis[j];
         		}
         		
			ppBasis[i] = pTemp;

         		// Update values of ppM, pB
         		(*pSwaps)++;

         		deep_LLL_GSO_update (ppM, pB, k , i, n);
         		k = (1 >(i)) ? 1 : i;
      		}
   	}
  	deleteBasisDouble(ppM, n);
  	deleteBasisDouble(ppBasisGSO, n);
   	deleteBasisDouble(ppBasisDouble, n);
   
   	return ppBasis;
}

// Pot-DeepLLL
inline int ** Pot_LLL_std (int ** ppBasis, int m, int n, long double delta, long long int * pSwaps, long long int * pReductions) 
{
	long double ** ppBasisGSO;
	ppBasisGSO = GSO (ppBasis, m, n);
	long double ** ppBasisDouble;
	ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
	long double ** ppM;
	ppM = (long double **) calloc (n, sizeof(long double *));
   	for (int i=0; i<n; i++) {
   		ppM[i] = (long double *) calloc (m, sizeof(long double));
	}

   	#ifdef __DEBUG__
   	printf ("\n\n\nThe GSO of the basis is:\n");
	printBasisDoubleAsColumns(ppBasisGSO, m, n);
   	#endif

   	// This 1-dimensional array pB will store all the values of ||b_i*||^2
   	long double * pB;
   	pB = (long double *) calloc (m, sizeof(long double));
   	
   	#ifdef __DEBUG__
   	printf ("\nThe values of ||b_i*||^2 are:\n");
   	#endif

   	for (int i=0; i<n; i++) {
     		pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
      		#ifdef __DEBUG__
      		cout << pB[i] << endl;
      		#endif
   	}

   	#ifdef __DEBUG__
   	printf ("\n\nThe values of mu_ij are:");
   	#endif
   
	long double inner = 0;

   	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
	 		inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
         		ppM[i][j] = inner / pB[j];

         		#ifdef __DEBUG__
         		cout << ppM[i][j] << " "; 
         		#endif
      		}
      	#ifdef __DEBUG__
     	printf ("\n");
      	#endif
   	}

   	// The PotLLL task
   	int k = 1;
   	long double C = 0.0;
   	int temp = 0;
	int closest_integer = 0;
	long double closest_integer_Lf = 0.0;

	long double P = 1.0;
	long double P_min = 1.0;
	int l = 1;
	long double sum = 0.0;
	long double P_mult = 0.0;
	int * pTemp;

   	while (k<n) {
      		// Reduce b_k with all previous basis vectors
      		for (int j=k-1; j>=0; j--) {
	 		closest_integer = CLOSEST_INTEGER(ppM[k][j]);
        		closest_integer_Lf = (long double) closest_integer; 
	 		if (closest_integer != 0) {
            			// Count the number of reductions
            			(*pReductions)++;
            			// Reduce b_k with b_j
            			reduce (ppBasis[j], ppBasis[k], m, closest_integer);

            			for (int i = j-1; i>=0; i--) {
               				// By Exercise 17.4.4 from Galbraith v2
               				ppM[k][i] -= closest_integer_Lf * ppM[j][i];
            			}

            			// mu_kj = mu_kj - [mu_kj]
            			ppM[k][j] -= closest_integer_Lf;
         		}
      		}

      		// Initialise P = 1, P_min = 1, k = 0
      		// For j = k-1,...,1,0, we find the index where an insertion of vector b_k in position j reduces the potential the most
      		P = 1;
      		P_min = 1;
      		l = 1;

      		for (int j=k-1; j>=0; j--) {
	 		sum = 0.0;
         		for (int r=j; r<k; r++) {
            			sum += ppM[k][r]*ppM[k][r]*pB[r];
         		}

         		P_mult = (pB[k] + sum) / pB[j];
         		P = P * P_mult;
         		if (P<P_min) {
            			l = j;
            			P_min = P;
         		}
      		}

      		#ifdef __DEBUG__
      		printf ("\nP_min = %Lf", P_min);
      		#endif

      		if (delta > P_min) {
         		// b_l <- b_k, and vectors b_l, ... , b_k-1 are shifted up one position
         		pTemp = ppBasis[k];

         		for (int j=k-1; j>=l; j--) {
            			ppBasis[j+1] = ppBasis[j];
         		}
         		ppBasis[l] = pTemp;

         		// Update values of ppM[i][j] for j<i, pB[i]
         		(*pSwaps)++;

         		deep_LLL_GSO_update (ppM, pB, k , l, n);

         		k = (1 > (l)) ? 1 : l;
         	} else {
         	k++;
      		}
   	}
   	deleteBasisDouble(ppM, n);
   	deleteBasisDouble(ppBasisGSO, n);
   	deleteBasisDouble(ppBasisDouble, n);
   
   	return ppBasis;
}

// LC-GG (block restricted insertions)
inline int ** block_LC_GG_std (int ** ppBasis, int m, int n, long double delta_threshold, int epsilon, long long int * pSwaps, long long int * pReductions, long long int * pGSORecomputations) 
{
	long double ** ppBasisGSO;
        ppBasisGSO = GSO (ppBasis, m, n);
        long double ** ppBasisDouble;
        ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
        long double ** ppM;
        ppM = (long double **) calloc (n, sizeof(long double *));
        for (int i=0; i<n; i++) {
                ppM[i] = (long double *) calloc (m, sizeof(long double));
        }   

	long double * pB;
	pB = (long double *) calloc (n, sizeof(long double));
	
	#ifdef __DEBUG__
  	printf ("\n\n\nThe GSO of the basis is:\n");
   	cout << ppBasisGSO;
   	#endif

	for (int i=0; i<n; i++) {
                pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
                #ifdef __DEBUG__
                cout << pB[i] << endl;
                #endif
        }

        #ifdef __DEBUG__
        printf ("\n\nThe values of mu_ij are:");
        #endif

	long double inner = 0;

   	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
	 		inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
         		ppM[i][j] = inner / pB[j];

         		#ifdef __DEBUG__
         		cout << ppM[i][j] << " "; 
         		#endif
      		}
		
      		#ifdef __DEBUG__
     		printf ("\n");
      		#endif
   	}

   	// The greedy global LLL task
   	*pSwaps = 0;
   	*pReductions = 0;
   	int pReductionIter = 0;
   	long double min_delta, old_min_delta;

	int closest_integer = 0;
        long double closest_integer_Lf = 0.0;
        // First reduce every vector by all its previous vectors
        for (int k=1; k<n; k++) {

                // Reduce b_k with all the previous basis vectors
                for (int j=k-1; j>=0; j--) {
                        closest_integer = CLOSEST_INTEGER(ppM[k][j]);
                        closest_integer_Lf = (long double) closest_integer;
                        if (closest_integer != 0) {
                                // Count the number of reductions
                                (*pReductions)++;
                                // Reduce b_k with b_j
                                reduce (ppBasis[j], ppBasis[k], m, closest_integer);

                                for (int i=j-1; i>=0; i--) {
                                        // By Exercise 17.4.4 from Galbraith v2
                                        ppM[k][i] -= closest_integer_Lf * ppM[j][i];
                                }

                                // mu_kj = mu_kj - [mu_kj]
                                ppM[k][j] -= closest_integer_Lf;
                        }
                }
        }

   	(pReductionIter)++;

   	// Find the index where there is a potential swap
   	int k = 0;
   	int i = 0;
	int old_k = 0;
	int old_i = 0;

	// ppKI is a 2D array, where ppKI[k][i] is the latest iteration number where
	// an insertion of b_k into position i is performed
        int ** ppKI;
        ppKI = (int **) calloc (n+1, sizeof(int *));
        for (int i=0; i<n+1; i++) {
                ppKI[i] = (int *) calloc (n+1, sizeof(int));
        }   

	for (int i=0; i<n+1; i++) {
		for (int j=0; j<n+1; j++) {
			ppKI[i][j] = 0;
		}
	}

	int max_loop = NCR(n,2);

	// pK, pI are 1D arrays being used as stacks
	// They keep track of the indices (k,i) of insertion
	long int * pK;
	long int * pI;
	pK = (long int *) calloc (max_loop, sizeof(long int));
	pI = (long int *) calloc (max_loop, sizeof(long int));

	for (int i=0; i< max_loop; i++) {
		pK[i] = 0;
		pI[i] = 0;
	}

	int iteration_count = 1;
	int j = 0;

	// ik keeps track of the iteration number modulo the maximum possible loop length
	int ik;
	ik = 0;
  
	// loop_flag = 0 if we encounter a loop in the insertions
	int loop_flag = 1;
	// loop_length denotes length of a loop in iterations if one were to exist
	// We use this to iteratively check the stacks for a loop
	// If a loop exists, we recompute the GSO information
	int loop_length = 0;

	// Counts number of GSO recomputations performed
	int loop_count = 0;

   	min_delta = find_insertion_indices_blockwise_LC_GG(ppM, ppBasis, pB, n, m, epsilon, &k, &i);

	pK[ik] = k;
	pI[ik] = i;
	ppKI[k][i] = iteration_count; 
				
	ik = (ik + 1) % max_loop;

   	#ifdef __DEBUG__
   	cout << "\n\n\nLC-GG: min_delta=" << min_delta;
   	#endif

   	int * pTemp;
    
   	while (min_delta < delta_threshold) {
	
      		#ifdef __DEBUG__
      		cout << "Insertion at k=" << k << "i=" << i << "for delta=" << min_delta << endl;
      		#endif

		pTemp = ppBasis[k];

      		for (int j=k-1; j>=i; j--) {
        		ppBasis[j+1] = ppBasis[j];
      		}
      		ppBasis[i] = pTemp;

      		// Update values of ppM, pB, delta
      		deep_LLL_GSO_update (ppM, pB, k , i, n);
	
		// Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
                for (int l=i; l<n; l++) {
                 // Reduce b_l with all the previous basis vectors
                        for (int j=l-1; j>=0; j--) {
                                closest_integer = CLOSEST_INTEGER(ppM[l][j]);
                                closest_integer_Lf = (long double) closest_integer;
                                if (closest_integer != 0) {
                                        // Count the number of reductions
                                        (*pReductions)++;
                                        // Reduce b_l with b_j
                                        reduce (ppBasis[j], ppBasis[l], m, closest_integer);

                                        for (int i=j-1; i>=0; i--) {
                                                // By Exercise 17.4.4 from Galbraith v2
                                                ppM[l][i] -= closest_integer_Lf * ppM[j][i];
                                        }

                                        // mu_lj = mu_lj - [mu_lj]
                                        ppM[l][j] -= closest_integer_Lf;
                                }
                        }
                }

      		(pReductionIter)++;

      		#ifdef __DEBUG__
      		printf ("\nThe Basis vectors are:\n");
      		cout << ppBasis;
      		printf ("\nThe values of mu_ij are:");
      		for (int i=0; i<n; i++) {
         		for (int j=0; j<i; j++) {
            			cout << ppM[i][j] << endl << ppDelta[i][j];
         		}
         		printf ("\n");
      		}
      		#endif
	
		old_k = k;
		old_i = i;
      		// Update k
      		min_delta = find_insertion_indices_blockwise_LC_GG(ppM, ppBasis, pB, n, m, epsilon, &k, &i);

		(iteration_count)++;
		ik = (ik + 1) % max_loop;

		if (iteration_count - ppKI[k][i] <= max_loop) {
			j = ppKI[k][i];
			loop_length = iteration_count - j;
			loop_flag = 1;

			for (int t=1; t<=loop_length; t++) {
				if (pK[(ik - t) % max_loop] != pK[(ik - loop_length - t) % max_loop] || pI[(ik - t) % max_loop] != pI[(ik - loop_length - t) % max_loop]) {
					loop_flag = 0;
					break;

				}
			}

			// If a loop exists, recompute GSO information
			if (loop_flag == 1) {
				(loop_count)++;
                                (*pGSORecomputations)++;
        			ppBasisGSO = GSO (ppBasis, m, n);
				
				for (int i=0; i<n; i++) {
                			pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
        			}

				inner = 0;
				ppBasisDouble = copyBasisToDouble (ppBasis, m, n);

   				for (int i=1; i<n; i++) {
      					for (int j=0; j<i; j++) {
	 					inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
         					ppM[i][j] = inner / pB[j];

      					}
		
   				}
				
				min_delta = find_insertion_indices_blockwise_LC_GG(ppM, ppBasis, pB, n, m, epsilon, &k, &i);
			}
		}
		
		pK[ik] = k;
		pI[ik] = i;
		ppKI[k][i] = iteration_count;
		

   	}

	*pSwaps = iteration_count;

	deleteBasis(ppKI, n+1);
   	free (pK);
   	free (pI);
	deleteBasisDouble(ppM, n);
        deleteBasisDouble(ppBasisGSO, n);
        deleteBasisDouble(ppBasisDouble, n);
   
   	return ppBasis;
}

// LC-GG
inline int ** LC_GG_std (int ** ppBasis, int m, int n, long double delta_threshold, long long int * pSwaps, long long int * pReductions, long long int * pGSORecomputations) 
{	
	long double ** ppBasisGSO;
        ppBasisGSO = GSO (ppBasis, m, n);
        long double ** ppBasisDouble;
        ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
        long double ** ppM;
        ppM = (long double **) calloc (n, sizeof(long double *));
        for (int i=0; i<n; i++) {
                ppM[i] = (long double *) calloc (m, sizeof(long double));
        }   

	long double * pB;
	pB = (long double *) calloc (n, sizeof(long double));
	
	#ifdef __DEBUG__
  	printf ("\n\n\nThe GSO of the basis is:\n");
   	cout << ppBasisGSO;
   	#endif

	for (int i=0; i<n; i++) {
                pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
                #ifdef __DEBUG__
                cout << pB[i] << endl;
                #endif
        }

        #ifdef __DEBUG__
        printf ("\n\nThe values of mu_ij are:");
        #endif

	long double inner = 0;

   	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
	 		inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
         		ppM[i][j] = inner / pB[j];

         		#ifdef __DEBUG__
         		cout << ppM[i][j] << " "; 
         		#endif
      		}
		
      		#ifdef __DEBUG__
     		printf ("\n");
      		#endif
   	}

   	// The greedy global LLL task
   	*pSwaps = 0;
   	*pReductions = 0;
   	int pReductionIter = 0;
   	long double min_delta, old_min_delta;

	int closest_integer = 0;
        long double closest_integer_Lf = 0.0;
        // First reduce every vector by all its previous vectors
        for (int k=1; k<n; k++) {

                // Reduce b_k with all the previous basis vectors
                for (int j=k-1; j>=0; j--) {
                        closest_integer = CLOSEST_INTEGER(ppM[k][j]);
                        closest_integer_Lf = (long double) closest_integer;
                        if (closest_integer != 0) {
                                // Count the number of reductions
                                (*pReductions)++;
                                // Reduce b_k with b_j
                                reduce (ppBasis[j], ppBasis[k], m, closest_integer);

                                for (int i=j-1; i>=0; i--) {
                                        // By Exercise 17.4.4 from Galbraith v2
                                        ppM[k][i] -= closest_integer_Lf * ppM[j][i];
                                }

                                // mu_kj = mu_kj - [mu_kj]
                                ppM[k][j] -= closest_integer_Lf;
                        }
                }
        }

   	(pReductionIter)++;

   	int k = 0;
   	int i = 0;
	// ppKI is a 2D array, where ppKI[k][i] is the latest iteration number where
	// an insertion of b_k into position i is performed
        int ** ppKI;
        ppKI = (int **) calloc (n+1, sizeof(int *));
        for (int i=0; i<n+1; i++) {
                ppKI[i] = (int *) calloc (n+1, sizeof(int));
        }   

	for (int i=0; i<n+1; i++) {
		for (int j=0; j<n+1; j++) {
			ppKI[i][j] = 0;
		}
	}

	int max_loop = NCR(n,2);

	// pK, pI are 1D arrays being used as stacks
	// They keep track of the indices (k,i) of insertion
	long int * pK;
	long int * pI;
	pK = (long int *) calloc (max_loop, sizeof(long int));
	pI = (long int *) calloc (max_loop, sizeof(long int));

	for (int i=0; i< max_loop; i++) {
		pK[i] = 0;
		pI[i] = 0;
	}

	int iteration_count = 1;
	int j = 0;

	// ik keeps track of the iteration number modulo the maximum possible loop length
	int ik;
	ik = 0;
  
	// loop_flag = 0 if we encounter a loop in the insertions
	int loop_flag = 1;
	// loop_length denotes length of a loop in iterations if one were to exist
	// We use this to iteratively check the stacks for a loop
	// If a loop exists, we recompute the GSO information
	int loop_length = 0;

	// Counts number of GSO recomputations performed
	int loop_count = 0;

   	min_delta = find_insertion_indices_LC_GG(ppM, ppBasis, pB, n, m, &k, &i);

	pK[ik] = k;
	pI[ik] = i;
	ppKI[k][i] = iteration_count; 
				
	ik = (ik + 1) % max_loop;

   	#ifdef __DEBUG__
   	cout << "\n\n\nLC-GG: min_delta=" << min_delta;
   	#endif

   	int * pTemp;
    
   	while (min_delta < delta_threshold) {

      		#ifdef __DEBUG__
      		cout << "Insertion at k=" << k << "i=" << i << "for delta=" << min_delta << endl;
      		#endif

		pTemp = ppBasis[k];
	
      		for (int j=k-1; j>=i; j--) {
        		ppBasis[j+1] = ppBasis[j];
      		}
      		
		ppBasis[i] = pTemp;

      		// Update values of ppM, pB, delta
      		deep_LLL_GSO_update (ppM, pB, k , i, n);
	
      		// Count the number of swaps
      		(*pSwaps)++;

		// Reduce every vector b_{l}, i <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
                for (int l=i; l<n; l++) {
                 // Reduce b_l with all the previous basis vectors
                        for (int j=l-1; j>=0; j--) {
                                closest_integer = CLOSEST_INTEGER(ppM[l][j]);
                                closest_integer_Lf = (long double) closest_integer;
                                if (closest_integer != 0) {
                                        // Count the number of reductions
                                        (*pReductions)++;
                                        // Reduce b_l with b_j
                                        reduce (ppBasis[j], ppBasis[l], m, closest_integer);

                                        for (int i=j-1; i>=0; i--) {
                                                // By Exercise 17.4.4 from Galbraith v2
                                                ppM[l][i] -= closest_integer_Lf * ppM[j][i];
                                        }

                                        // mu_lj = mu_lj - [mu_lj]
                                        ppM[l][j] -= closest_integer_Lf;
                                }
                        }
                }

      		(pReductionIter)++;

      		#ifdef __DEBUG__
      		printf ("\nThe Basis vectors are:\n");
      		cout << ppBasis;
      		printf ("\nThe values of mu_ij are:");
      		for (int i=0; i<n; i++) {
         		for (int j=0; j<i; j++) {
            			cout << ppM[i][j] << endl << ppDelta[i][j];
         		}
         		printf ("\n");
      		}
      		#endif
	
      		old_min_delta = min_delta;
      		min_delta = find_insertion_indices_LC_GG(ppM, ppBasis, pB, n, m, &k, &i);

		(iteration_count)++;
		ik = (ik + 1) % max_loop;

		if (iteration_count - ppKI[k][i] <= max_loop) {
			j = ppKI[k][i];
			loop_length = iteration_count - j;
			loop_flag = 1;

			for (int t=1; t<=loop_length; t++) {
				if (pK[(ik - t) % max_loop] != pK[(ik - loop_length - t) % max_loop] || pI[(ik - t) % max_loop] != pI[(ik - loop_length - t) % max_loop]) {
					loop_flag = 0;
					break;

				}
			}

			// If a loop exists, recompute GSO information
			if (loop_flag == 1) {
				(loop_count)++;
                                (*pGSORecomputations)++;
        			ppBasisGSO = GSO (ppBasis, m, n);
				
				for (int i=0; i<n; i++) {
                			pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
        			}

				inner = 0;
				ppBasisDouble = copyBasisToDouble (ppBasis, m, n);

   				for (int i=1; i<n; i++) {
      					for (int j=0; j<i; j++) {
	 					inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
         					ppM[i][j] = inner / pB[j];
      					}
   				}
				
				min_delta = find_insertion_indices_LC_GG(ppM, ppBasis, pB, n, m, &k, &i);
			}
		}
		
		pK[ik] = k;
		pI[ik] = i;
		ppKI[k][i] = iteration_count;
		
   	}
	*pSwaps = iteration_count;

	deleteBasis(ppKI, n+1);
   	free (pK);
   	free (pI);

	deleteBasisDouble(ppM, n);
        deleteBasisDouble(ppBasisGSO, n);
        deleteBasisDouble(ppBasisDouble, n);
   
   return ppBasis;
}

// LC-PGG
inline int ** LC_PGG(int ** ppBasis, int m, int n, long double delta_threshold, long long int * pSwaps, long long int * pReductions, long long int * pIterations) 
{	
	long double ** ppBasisGSO;
        ppBasisGSO = GSO (ppBasis, m, n);
        long double ** ppBasisDouble;
        ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
	long double ** ppDelta;
        long double ** ppM;
        ppM = (long double **) calloc (n, sizeof(long double *));
        ppDelta = (long double **) calloc (n, sizeof(long double *));
        for (int i=0; i<n; i++) {
                ppM[i] = (long double *) calloc (m, sizeof(long double));
                ppDelta[i] = (long double *) calloc (m, sizeof(long double));
        }   

	long double * pB;
	pB = (long double *) calloc (n, sizeof(long double));
	
	#ifdef __DEBUG__
  	printf ("\n\n\nThe GSO of the basis is:\n");
   	cout << ppBasisGSO;
   	#endif

	for (int i=0; i<n; i++) {
                pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
                #ifdef __DEBUG__
                cout << pB[i] << endl;
                #endif
        }

        #ifdef __DEBUG__
        printf ("\n\nThe values of mu_ij are:");
        #endif

	long double inner = 0;

   	for (int i=1; i<n; i++) {
      		for (int j=0; j<i; j++) {
	 		inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
         		ppM[i][j] = inner / pB[j];

         		#ifdef __DEBUG__
         		cout << ppM[i][j] << " "; 
         		#endif
      		}
		
      		#ifdef __DEBUG__
     		printf ("\n");
      		#endif
   	}

   	// The greedy global LLL task
	int nInsertions = 0;
   	*pSwaps = 0;
   	*pReductions = 0;
	*pIterations = 0;

	int closest_integer = 0;
        long double closest_integer_Lf = 0.0;

   	int i_min = 1;
	long double C = 0.0;

   	while (i_min < n) { 

		(*pIterations)++;

		// Reduce every vector b_{l}, i_min <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
		// Initially, i_min = 1, corresponding to full size reduction
                for (int l=i_min; l<n; l++) {
                 // Reduce b_l with all the previous basis vectors
                        for (int j=l-1; j>=0; j--) {
                                closest_integer = CLOSEST_INTEGER(ppM[l][j]);
                                closest_integer_Lf = (long double) closest_integer;
                                if (closest_integer != 0) {
                                        // Count the number of reductions
                                        (*pReductions)++;
                                        // Reduce b_l with b_j
                                        reduce (ppBasis[j], ppBasis[l], m, closest_integer);

                                        for (int i=j-1; i>=0; i--) {
                                                // By Exercise 17.4.4 from Galbraith v2
                                                ppM[l][i] -= closest_integer_Lf * ppM[j][i];
                                        }

                                        // mu_lj = mu_lj - [mu_lj]
                                        ppM[l][j] -= closest_integer_Lf;
                                }
                        }
                }

		// Compute ppDelta
   		for (int k=1; k<n; k++) {
			C = (long double) inner_product_template(ppBasis[k], ppBasis[k], m);
      			for (int j=0; j<k; j++) {
         			// Lovasz condition definition
				ppDelta[k][j] = (C / pB[j]);
			
				// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
         			C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
         		}
   		}

		nInsertions = 0;
		i_min = LC_PGG_recursive_search (ppDelta, ppM, pB, ppBasis, 0, n-1, n, delta_threshold, &nInsertions);

      		// Count the number of swaps
      		*pSwaps += nInsertions;
		
   	}

	deleteBasisDouble(ppM, n);
	deleteBasisDouble(ppDelta, n);
        deleteBasisDouble(ppBasisGSO, n);
        deleteBasisDouble(ppBasisDouble, n);
   
 	return ppBasis;
}

// LC-PGG with globally restricted insertions
inline int ** global_block_LC_PGG(int ** ppBasis, int m, int n, long double delta_threshold, int epsilon,  long long int * pSwaps, long long int * pReductions, long long int * pIterations) 
{	
	long double ** ppBasisGSO;
	ppBasisGSO = GSO (ppBasis, m, n);
	long double ** ppBasisDouble;
	ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
	long double ** ppDelta;
	long double ** ppM;
	ppM = (long double **) calloc (n, sizeof(long double *));
	ppDelta = (long double **) calloc (n, sizeof(long double *));
	for (int i=0; i<n; i++) {
		ppM[i] = (long double *) calloc (m, sizeof(long double));
		ppDelta[i] = (long double *) calloc (m, sizeof(long double));
	}   

	long double * pB;
	pB = (long double *) calloc (n, sizeof(long double));
	
	#ifdef __DEBUG__
	printf ("\n\n\nThe GSO of the basis is:\n");
	cout << ppBasisGSO;
	#endif

	for (int i=0; i<n; i++) {
		pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
		#ifdef __DEBUG__
		cout << pB[i] << endl;
		#endif
	}

	#ifdef __DEBUG__
	printf ("\n\nThe values of mu_ij are:");
	#endif

	long double inner = 0;

	for (int i=1; i<n; i++) {
		for (int j=0; j<i; j++) {
			inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
			ppM[i][j] = inner / pB[j];

			#ifdef __DEBUG__
			cout << ppM[i][j] << " "; 
			#endif
		}
		
		#ifdef __DEBUG__
		printf ("\n");
		#endif
	}

	// The greedy global LLL task
	int nInsertions = 0;
	*pSwaps = 0;
	*pReductions = 0;
	*pIterations = 0;

	int closest_integer = 0;
	long double closest_integer_Lf = 0.0;

	int i_min = 1;
	long double C = 0.0;
	int index = 0;

	while (i_min < n) { 

		(*pIterations)++;

		// Reduce every vector b_{l}, i_min <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
		// Initially, i_min = 1, corresponding to full size reduction
		for (int l=i_min; l<n; l++) {
		 // Reduce b_l with all the previous basis vectors
			for (int j=l-1; j>=0; j--) {
				closest_integer = CLOSEST_INTEGER(ppM[l][j]);
				closest_integer_Lf = (long double) closest_integer;
				if (closest_integer != 0) {
					// Count the number of reductions
					(*pReductions)++;
					// Reduce b_l with b_j
					reduce (ppBasis[j], ppBasis[l], m, closest_integer);

					for (int i=j-1; i>=0; i--) {
						// By Exercise 17.4.4 from Galbraith v2
						ppM[l][i] -= closest_integer_Lf * ppM[j][i];
					}

					// mu_lj = mu_lj - [mu_lj]
					ppM[l][j] -= closest_integer_Lf;
				}
			}
		}

   		for (int k=1; k<n; k++) {
			C = (long double) inner_product_template(ppBasis[k], ppBasis[k], m);
      		
			if (k <= 2*epsilon) {
        			for (int j=0; j<k; j++) {
            			// Lovasz condition definition
					ppDelta[k][j] = C / pB[j];
					// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            				C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
				}

			} else {
	 			index = k - epsilon; 
	 			// Check first epsilon indices
         			for (int j=0; j<epsilon; j++) {
	    				// Lovasz condition definition
					ppDelta[k][j] = C / pB[j];
            			
					// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            				C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 			}
	
		 		// Update C for the indices not being checked
         			for (int j=epsilon; j<index; j++) {
					ppDelta[k][j] = 1.0;
	    				C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 			}

	 			// Check indices from k-epsilon+1 to k-1
	 			for (int j=index; j<k; j++) {
            				// Lovasz condition definition
					ppDelta[k][j] = C / pB[j];
            
					// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            				C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 			}
      			}
   		}

		nInsertions = 0;
		i_min = block_LC_PGG_recursive_search_global (ppDelta, ppM, pB, ppBasis, 0, n-1, n, epsilon, delta_threshold, &nInsertions);

		// Count the number of swaps
		*pSwaps += nInsertions;
		
	}

	deleteBasisDouble(ppM, n);
	deleteBasisDouble(ppDelta, n);
	deleteBasisDouble(ppBasisGSO, n);
	deleteBasisDouble(ppBasisDouble, n);
   
	return ppBasis;
}

// LC-PGG with locally restricted insertions
inline int ** local_block_LC_PGG(int ** ppBasis, int m, int n, long double delta_threshold, int epsilon,  long long int * pSwaps, long long int * pReductions, long long int * pIterations) 
{	
	long double ** ppBasisGSO;
	ppBasisGSO = GSO (ppBasis, m, n);
	long double ** ppBasisDouble;
	ppBasisDouble = copyBasisToDouble (ppBasis, m, n);
	long double ** ppDelta;
	long double ** ppM;
	ppM = (long double **) calloc (n, sizeof(long double *));
	ppDelta = (long double **) calloc (n, sizeof(long double *));
	for (int i=0; i<n; i++) {
		ppM[i] = (long double *) calloc (m, sizeof(long double));
		ppDelta[i] = (long double *) calloc (m, sizeof(long double));
	}   

	long double * pB;
	pB = (long double *) calloc (n, sizeof(long double));
	
	#ifdef __DEBUG__
	printf ("\n\n\nThe GSO of the basis is:\n");
	cout << ppBasisGSO;
	#endif

	for (int i=0; i<n; i++) {
		pB[i] = inner_product_template(ppBasisGSO[i], ppBasisGSO[i], m);
		#ifdef __DEBUG__
		cout << pB[i] << endl;
		#endif
	}

	#ifdef __DEBUG__
	printf ("\n\nThe values of mu_ij are:");
	#endif

	long double inner = 0;

	for (int i=1; i<n; i++) {
		for (int j=0; j<i; j++) {
			inner = inner_product_template(ppBasisDouble[i], ppBasisGSO[j], m);
			ppM[i][j] = inner / pB[j];

			#ifdef __DEBUG__
			cout << ppM[i][j] << " "; 
			#endif
		}
		
		#ifdef __DEBUG__
		printf ("\n");
		#endif
	}

	// The greedy global LLL task
	int nInsertions = 0;
	*pSwaps = 0;
	*pReductions = 0;
	*pIterations = 0;

	int closest_integer = 0;
	long double closest_integer_Lf = 0.0;

	int i_min = 1;
	long double C = 0.0;
	int index = 0;

	while (i_min < n) { 

		(*pIterations)++;

		// Reduce every vector b_{l}, i_min <= l < n, by all its previous vectors b_{k-1} ... b_{l-1}
		// Initially, i_min = 1, corresponding to full size reduction
		for (int l=i_min; l<n; l++) {
		 // Reduce b_l with all the previous basis vectors
			for (int j=l-1; j>=0; j--) {
				closest_integer = CLOSEST_INTEGER(ppM[l][j]);
				closest_integer_Lf = (long double) closest_integer;
				if (closest_integer != 0) {
					// Count the number of reductions
					(*pReductions)++;
					// Reduce b_l with b_j
					reduce (ppBasis[j], ppBasis[l], m, closest_integer);

					for (int i=j-1; i>=0; i--) {
						// By Exercise 17.4.4 from Galbraith v2
						ppM[l][i] -= closest_integer_Lf * ppM[j][i];
					}

					// mu_lj = mu_lj - [mu_lj]
					ppM[l][j] -= closest_integer_Lf;
				}
			}
		}

   		for (int k=1; k<n; k++) {
			C = (long double) inner_product_template(ppBasis[k], ppBasis[k], m);
      		
			if (k <= 2*epsilon) {
        			for (int j=0; j<k; j++) {
            			// Lovasz condition definition
					ppDelta[k][j] = C / pB[j];
					// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            				C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
				}

			} else {
	 			index = k - epsilon; 
	 			// Check first epsilon indices
         			for (int j=0; j<epsilon; j++) {
	    				// Lovasz condition definition
					ppDelta[k][j] = C / pB[j];
            			
					// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            				C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 			}
	
		 		// Update C for the indices not being checked
         			for (int j=epsilon; j<index; j++) {
					ppDelta[k][j] = 1.0;
	    				C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 			}

	 			// Check indices from k-epsilon+1 to k-1
	 			for (int j=index; j<k; j++) {
            				// Lovasz condition definition
					ppDelta[k][j] = C / pB[j];
            
					// Definition 2 [Projection of b_k \perp (b_1, b_2, ..., b_{i-1})] of [YY2019]
            				C = C - (ppM[k][j] * ppM[k][j] * pB[j]);
	 			}
      			}
   		}

		nInsertions = 0;
		i_min = block_LC_PGG_recursive_search_local (ppDelta, ppM, pB, ppBasis, 0, n-1, n, epsilon, delta_threshold, &nInsertions);

		// Count the number of swaps
		*pSwaps += nInsertions;
		
	}

	deleteBasisDouble(ppM, n);
	deleteBasisDouble(ppDelta, n);
	deleteBasisDouble(ppBasisGSO, n);
	deleteBasisDouble(ppBasisDouble, n);
   
	return ppBasis;
}

// LC-DeepLLL
inline int ** LC_DeepLLL_std (int ** ppBasis, int m, int n, long double delta, int max_loop, long long int * pSwaps, long long int * pReductions) 
{
	long double ** ppBasisGSO;

	// Find the GSO of the Basis
   	ppBasisGSO = GSO (ppBasis, m, n);

   	#ifdef __DEBUG__
   	printf ("\n\n\nThe GSO of the basis is:\n");
   	printBasisDoubleAsColumns (ppBasisGSO, m, n);
   	#endif

   	// This 1-dimensional array pB will store all the values of ||b_i*||^2
  	long double * pB;
   	pB = (long double *) calloc (n, sizeof(long double));
   	#ifdef __DEBUG__
  	printf ("\nThe values of ||b_i*||^2 are:\n");
   	#endif
   	
	for (int i=0; i<n; i++) {
      		pB[i] = inner_product_template (ppBasisGSO[i], ppBasisGSO[i], m);
     	 	#ifdef __DEBUG__
      		printf ("%Lf, ", pB[i]);
      		#endif
   	}

   	// This 2-dimensional nxn array ppM will store all the values of mu_ij
   	long double ** ppM;
   	ppM = (long double **) calloc (n, sizeof (long double*));
   	#ifdef __DEBUG__
   	printf ("\n\nThe values of mu_ij are:");
   	#endif
   	for (int i=1; i<n; i++) {
      		ppM[i] = (long double *) calloc (n, sizeof (long double));
      		
		for (int j=0; j<i; j++) {
         		ppM[i][j] = (long double)(inner_product_template (ppBasis[i], ppBasisGSO[j], m)) / (long double)pB[j];
      		}
      		
		#ifdef __DEBUG__
      		printf ("\n");
      		#endif
   	}

   	// The DeepLLL task
   	*pSwaps = 0;
   	*pReductions = 0;
   	long long int potIncrease = 0;
        long double prevPot = 0.0;
	long double currentPot = 0.0;
	prevPot = basisPotential(ppBasis, m, n);
   	int k = 1;
      	int p = 0;
	// ppKI is a 2D array, where ppKI[k][i] is the latest iteration number where
	// an insertion of b_k into position i is performed
        int ** ppKI;
        ppKI = (int **) calloc (n+1, sizeof(int *));
        for (int i=0; i<n+1; i++) {
                ppKI[i] = (int *) calloc (n+1, sizeof(int));
        }   

	for (int i=0; i<n+1; i++) {
		for (int j=0; j<n+1; j++) {
			ppKI[i][j] = 0;
		}
	}

	// pK, pI are 1D arrays being used as stacks
	// They keep track of the indices (k,i) of insertion
	long int * pK;
	long int * pI;
	pK = (long int *) calloc (max_loop, sizeof(long int));
	pI = (long int *) calloc (max_loop, sizeof(long int));

	for (int i=0; i< max_loop; i++) {
		pK[i] = 0;
		pI[i] = 0;
	}

	int iteration_count = 0;
	int j = 0;

	// ik keeps track of the iteration number modulo the maximum possible loop length
	int ik;
	ik = 0;
  
	// loop_flag = 0 if we encounter a loop in the insertions
	int loop_flag = 1;
	// loop_length denotes length of a loop in iterations if one were to exist
	// We use this to iteratively check the stacks for a loop
	// If a loop exists, we recompute the GSO information
	int loop_length = 0;

	// Counts number of GSO recomputations performed
	int loop_count = 0;
				
	ik = (ik + 1) % max_loop;

   	while (k<n) {

      		// Keeps track of C - the projected vector used in Lovasz Condition Checks
      		long double C = 0.0;

      		// Reduce b_k with all previous basis vectors
      		for (int j=k-1; j>=0; j--) {
         		int closest_integer_to_mu_kj = CLOSEST_INTEGER (ppM[k][j]);
         
	 		if (closest_integer_to_mu_kj != 0) {
            			// Count the number of reductions
            			(*pReductions)++;
            			// Reduce b_k with b_j
            			reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            			for (int i = j-1; i>=0; i--) {
               				// By Exercise 17.4.4 from Galbraith v2
               				ppM[k][i] -= closest_integer_to_mu_kj * ppM[j][i];
            			}

            			// mu_kj = mu_kj - [mu_kj]
            			ppM[k][j] -= closest_integer_to_mu_kj;
         		}
      		}

      		// Compute C = ||b_k||^2
      		C = inner_product_template (ppBasis[k], ppBasis[k],m);

      		p = k;
      
		// Check Lovasz condition for i=0, i=1,...,i=k-1
      		for (int i=0; i<p; i++) {

         		if (C >= (delta * pB[i])) {

            			C = C - (ppM[k][i] * ppM[k][i] * pB[i]);
            
           		} else {
				(iteration_count)++;
				ik = (ik + 1) % max_loop;

				if (iteration_count - ppKI[k][i] <= max_loop) {
					j = ppKI[k][i];
					loop_length = iteration_count - j;
					loop_flag = 1;
				}
	
				for (int t=1; t<=loop_length; t++) {
					if (pK[(ik - t) % max_loop] != pK[(ik - loop_length - t) % max_loop] || pI[(ik - t) % max_loop] != pI[(ik - loop_length - t) % max_loop]) {
						loop_flag = 0;
						break;
					}
				}

				// If a loop exists, recompute GSO information
				if (loop_flag == 1) {
					(loop_count)++;
					cout << "Loops = " << loop_count << endl;
        				ppBasisGSO = GSO (ppBasis, m, n);

				        for (int i=0; i<n; i++) {
                				pB[i] = inner_product_template (ppBasisGSO[i], ppBasisGSO[i], m);
        				}

        				for (int i=1; i<n; i++) {
                    			for (int j=0; j<i; j++) {
                        				ppM[i][j] = (long double)(inner_product_template (ppBasis[i], ppBasisGSO[j], m)) / (long double)pB[j];
                				}

        				}
		
      					for (int i=0; i<p; i++) {

         					if (C >= (delta * pB[i])) {

            						C = C - (ppM[k][i] * ppM[k][i] * pB[i]);
            
           					} else {
							int * pTemp;
            						// b_i <- b_k, and vectors b_i, ... , b_k-1 are shifted up one position
            						pTemp = ppBasis[k];

            						for (int j=k-1; j>=i; j--) {
               							ppBasis[j+1] = ppBasis[j];
            						}
            				
							ppBasis[i] = pTemp;

            						// Update the values of ppM, pB (to be done)
            						(*pSwaps)++;

            						deep_LLL_GSO_update (ppM, pB, k , i, n);
            						k = (1 > (i)) ? 0 : i-1;
							currentPot = basisPotential(ppBasis, m, n);
							if (currentPot > prevPot) {
	            					potIncrease++;	    
							}
            						break;
						}
					}
					break;

				} else {
					int * pTemp;
					// b_i <- b_k, and vectors b_i, ... , b_k-1 are shifted up one position
                                        pTemp = ppBasis[k];

                                        for (int j=k-1; j>=i; j--) {
                                        	ppBasis[j+1] = ppBasis[j];
                                        }

                                        ppBasis[i] = pTemp;

                                        // Update the values of ppM, pB (to be done)
                                        (*pSwaps)++;

                                        deep_LLL_GSO_update (ppM, pB, k , i, n);
                                        k = (1 > (i)) ? 0 : i-1;
					currentPot = basisPotential(ppBasis, m, n);
					if (currentPot > prevPot) {
	            			potIncrease++;	    
					}
                                        break;
				}
         		}
         
      		}

		k++;
   	}

   	deleteBasisDouble(ppM, n);
   	deleteBasisDouble(ppBasisGSO, n);
	deleteBasis(ppKI, n+1);
   	free (pK);
   	free (pI);
	long double potNotIncrease = *pSwaps - potIncrease;
	cout << "NUMBER OF INSERTIONS = " << *pSwaps << endl; 
	cout << "NUMBER OF POT INCREASES = " << potIncrease << endl; 
	cout << "NUMBER OF POT NON-INCREASES = " << potNotIncrease << endl; 
   	return ppBasis;
}

// LC-DeepLLL (restricted insertions)
inline int ** LC_DeepLLL_std (int ** ppBasis, int m, int n, long double delta, int epsilon, int max_loop, long long int * pSwaps, long long int * pReductions) 
{
	long double ** ppBasisGSO;

	// Find the GSO of the Basis
   	ppBasisGSO = GSO (ppBasis, m, n);

   	#ifdef __DEBUG__
   	printf ("\n\n\nThe GSO of the basis is:\n");
   	printBasisDoubleAsColumns (ppBasisGSO, m, n);
   	#endif

   	// This 1-dimensional array pB will store all the values of ||b_i*||^2
  	long double * pB;
   	pB = (long double *) calloc (n, sizeof(long double));
   	#ifdef __DEBUG__
  	printf ("\nThe values of ||b_i*||^2 are:\n");
   	#endif
   	
	for (int i=0; i<n; i++) {
      		pB[i] = inner_product_template (ppBasisGSO[i], ppBasisGSO[i], m);
     	 	#ifdef __DEBUG__
      		printf ("%Lf, ", pB[i]);
      		#endif
   	}

   	// This 2-dimensional nxn array ppM will store all the values of mu_ij
   	long double ** ppM;
   	ppM = (long double **) calloc (n, sizeof (long double*));
   	#ifdef __DEBUG__
   	printf ("\n\nThe values of mu_ij are:");
   	#endif
   	for (int i=1; i<n; i++) {
      		ppM[i] = (long double *) calloc (n, sizeof (long double));
      		
		for (int j=0; j<i; j++) {
         		ppM[i][j] = (long double)(inner_product_template (ppBasis[i], ppBasisGSO[j], m)) / (long double)pB[j];
      		}
      		
		#ifdef __DEBUG__
      		printf ("\n");
      		#endif
   	}

   	// The DeepLLL task
   	*pSwaps = 0;
   	*pReductions = 0;
   	long long int potIncrease = 0;
        long double prevPot = 0.0;
	long double currentPot = 0.0;
	prevPot = basisPotential(ppBasis, m, n);

	long double delta_min = 0.0;
   	int k = 1;
      	int p = 0;
	// ppKI is a 2D array, where ppKI[k][i] is the latest iteration number where
	// an insertion of b_k into position i is performed
        int ** ppKI;
        ppKI = (int **) calloc (n+1, sizeof(int *));
        for (int i=0; i<n+1; i++) {
                ppKI[i] = (int *) calloc (n+1, sizeof(int));
        }   

	for (int i=0; i<n+1; i++) {
		for (int j=0; j<n+1; j++) {
			ppKI[i][j] = 0;
		}
	}

	// pK, pI are 1D arrays being used as stacks
	// They keep track of the indices (k,i) of insertion
	long int * pK;
	long int * pI;
	pK = (long int *) calloc (max_loop, sizeof(long int));
	pI = (long int *) calloc (max_loop, sizeof(long int));

	for (int i=0; i< max_loop; i++) {
		pK[i] = 0;
		pI[i] = 0;
	}

	int iteration_count = 0;
	int j = 0;
	int i_min = 0;

	// ik keeps track of the iteration number modulo the maximum possible loop length
	int ik;
	ik = 0;
  
	// loop_flag = 0 if we encounter a loop in the insertions
	int loop_flag = 1;
	// loop_length denotes length of a loop in iterations if one were to exist
	// We use this to iteratively check the stacks for a loop
	// If a loop exists, we recompute the GSO information
	int loop_length = 0;

	// Counts number of GSO recomputations performed
	int loop_count = 0;
				
	ik = (ik + 1) % max_loop;

   	while (k<n) {

      		// Keeps track of C - the projected vector used in Lovasz Condition Checks
      		long double C = 0.0;

      		// Reduce b_k with all previous basis vectors
      		for (int j=k-1; j>=0; j--) {
         		int closest_integer_to_mu_kj = CLOSEST_INTEGER (ppM[k][j]);
         
	 		if (closest_integer_to_mu_kj != 0) {
            			// Count the number of reductions
            			(*pReductions)++;
            			// Reduce b_k with b_j
            			reduce (ppBasis[j], ppBasis[k], m, closest_integer_to_mu_kj);

            			for (int i = j-1; i>=0; i--) {
               				// By Exercise 17.4.4 from Galbraith v2
               				ppM[k][i] -= closest_integer_to_mu_kj * ppM[j][i];
            			}

            			// mu_kj = mu_kj - [mu_kj]
            			ppM[k][j] -= closest_integer_to_mu_kj;
         		}
      		}

      		// Compute C = ||b_k||^2
                delta_min = find_insertion_indices_blockwise_LC_DeepLLL (ppM, ppBasis, pB, n, m, epsilon, k, &i_min);

         	if (delta_min < delta) {
			int * pTemp;
            		// b_i <- b_k, and vectors b_i, ... , b_k-1 are shifted up one position
            		pTemp = ppBasis[k];

			for (int j=k-1; j>=i_min; j--) {
              			ppBasis[j+1] = ppBasis[j];
            		}
            				
			ppBasis[i_min] = pTemp;

			// Update the values of ppM, pB (to be done)
            		(*pSwaps)++;

            		deep_LLL_GSO_update (ppM, pB, k , i_min, n);
            		k = (1 > (i_min)) ? 0 : i_min-1;
			
		 currentPot = basisPotential(ppBasis, m, n);
		 if (currentPot > prevPot) {
	        	 potIncrease++;	    
		 }

		}

		k++;
   	}

   	deleteBasisDouble(ppM, n);
   	deleteBasisDouble(ppBasisGSO, n);
	deleteBasis(ppKI, n+1);
   	free (pK);
   	free (pI);
	long double potNotIncrease = *pSwaps - potIncrease;
	cout << "NUMBER OF INSERTIONS = " << *pSwaps << endl; 
	cout << "NUMBER OF POT INCREASES = " << potIncrease << endl; 
	cout << "NUMBER OF POT NON-INCREASES = " << potNotIncrease << endl; 
   	return ppBasis;
}


