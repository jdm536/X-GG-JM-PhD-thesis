/* 
 * *****************************************************************************
 * Filename: common.h
 * Authors : Sanjay Bhattacherjee and Jack Moyler
 * *****************************************************************************
*/

#ifndef COMMON_H
#define COMMON_H

#include <cstring>
#include "test_utils.h"
#include <fplll.h>

using namespace std;
using namespace fplll;

#define __DEBUG__

#define DEFAULT_PRECISION 53
#define DEFAULT_BETA 0
#define DEFAULT_SUBLATTICE_SIZE 0
#define DEFAULT_GSO 0
#define DEFAULT_SEARCH_SWITCH 0
#define DEFAULT_PRECISION_CORRECTION_LOOPS_ALLOWED 2
#define DEFAULT_PRECISION_CORRECTION_LOOPS_ALLOWED_CFA 1

#define SEARCH_SWITCH_MIN_PRECISION_BINARY 1
#define SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE 2
#define SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE 3

#define PRECISION_STEP_DOWN 5
#define PRECISION_STEP_UP 5
#define PRECISION_JUMP_FACTOR 2

#define ETA 0.51

// -----------------------------------------------------------------------------
// Errors
#define INDEX_SEARCH_FAIL 0
// Note that the reduction logic depends on INDEX_SEARCH_FAIL being set to 0
#define PRECISION_FAIL 1
#define ERROR_COMMANDLINE_ARGUMENT_NUMBER 2
#define ERROR_ALGORITHM_NUMBER 3
#define ERROR_FILE_READ 4
#define ERROR_INVALID_INTEGER_TYPE_IN_BASIS 5
#define ERROR_INVALID_PRECISION_RANGE 6
#define ERROR_INVALID_INNER_PRODUCT 7
// -----------------------------------------------------------------------------

#define TYPE_MPZ_T 1
#define TYPE_LONG_LONG_INT 2
#define TYPE_LONG_INT 3
#define TYPE_INT 4

#define SIZE_LONG_LONG_INT ((4 * sizeof(long long int)) - 1)
#define SIZE_LONG_INT ((4 * sizeof(long int)) - 1)
#define SIZE_INT ((4 * sizeof(int)) - 1)

#define MAX_LONG_LONG_INT ((long long int) ((1 << SIZE_LONG_LONG_INT) - 1))
#define MAX_LONG_INT ((long int) ((1 << SIZE_LONG_INT) - 1))
#define MAX_INT ((int) ((1 << SIZE_INT) - 1))

// Checks equality of two input bases
template<typename Z_T>
inline bool basesEqual (ZZ_mat<Z_T> A, ZZ_mat<Z_T> A_previous, int d) 
{
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
        {
            if (A[i][j] != A_previous[i][j])
            {
                return false;
            }
        }
    }
    return true;
}

// Prints basis in the same way as FPLLL
template<typename Z_T>
inline void printBasis (Z_T ** &ppA, int d) 
{
    cout << "["; 
    for (int i=0; i<d; i++)
    {
        cout << "[";
        for (int j=0; j<d; j++)
        {
            cout << ppA[i][j] << " " ;
        }
        cout << "]" << endl;
    }
    cout << "]" << endl;

    return;
}

// Compute the Gram matrix of inner products <b_i,bj>
inline void compute_Gram_Matrix (
		ZZ_mat<mpz_t> &A,
		ZZ_mat<mpz_t> &A_gram,
		int d
		)
{
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++) 
	{
            A[i].dot_product (A_gram[i][j], A[j]);
	}
    }
    return;
}

// Compute the Gram matrix of inner products <b_i,bj>
template <typename Z_T>
inline void compute_Gram_Matrix (
                Z_T ** &ppA,
                long long int ** &ppA_gram,
		int d
		)
{
    long long int result;

    for (int i=0; i<d; i++) 
    {
	for (int j=0; j<=i; j++)
	{
            result = 0;
	    for (int k=0; k<d; k++) 
	    {
               result += (long long int) ppA[i][k] * (long long int) ppA[j][k];
	    }
	    ppA_gram[i][j] = result;
	}
    }

    return;
}

// Compute row k of the Gram matrix
inline void compute_Gram_Row (
		ZZ_mat<mpz_t> &A,
		ZZ_mat<mpz_t> &A_gram,
		int k,
		int d
		)
{
    for (int j=0; j<=k; j++) 
    {
        A[k].dot_product (A_gram[k][j], A[j]);
    }
    return;
}

// Compute row k of the Gram matrix
template <typename Z_T>
inline void compute_Gram_Row (
                Z_T ** &ppA,
                long long int ** &ppA_gram,
		int k,
		int d
		)
{
    long long int result;

    for (int i=0; i<=k; i++)
    {
        result = 0;
        for (int j=0; j<d; j++) 
        {
            result += (long long int) ppA[k][j] * (long long int) ppA[i][j];
	}
	ppA_gram[k][i] = result;
    }
    
    return;
}

// Computing the CFA from start_index to end_index
// Figure 4 of Nguyen and Stehle 2009
inline void compute_CFA (
                ZZ_mat<mpz_t> &A_gram,
                FP_mat<mpfr_t> &r,
                FP_mat<mpfr_t> &mu,
                FP_mat<mpfr_t> &s,
                int d,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    // Line 1
    for (int i=start_index; i<=end_index; i++)
    {
	// Line 2
        for (int j=0; j<i; j++)
	{
	    // Line 3
            r[i][j].set_z(A_gram[i][j], MPFR_RNDN);

	    // Line 4
	    for (int k=0; k<j; k++)
	    {
                r[i][j] -= mu[j][k] * r[i][k];
	    }

	    // Line 5
	    mu[i][j] = r[i][j] / r[j][j];
	}

	// Line 6
	s[i][0].set_z(A_gram[i][i], MPFR_RNDN);
	
	for (int j=1; j<=i; j++)
	{
	    s[i][j] = s[i][j-1] - (mu[i][j-1] * r[i][j-1]);
	}

	// Line 7
	r[i][i] = s[i][i];
    }

    return;
}

// Computing the CFA from start_index to end_index
// Figure 4 of Nguyen and Stehle 2009
template <typename FP_T>
inline void compute_CFA (
                long long int ** &ppA_gram,
                FP_mat<FP_T> &r,
                FP_mat<FP_T> &mu,
                FP_mat<FP_T> &s,
                int d,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    // Line 1
    for (int i=start_index; i<=end_index; i++)
    {
	// Line 2    
        for (int j=0; j<i; j++)
	{
	    // Line 3
	    r[i][j] = ppA_gram[i][j];

	    // Line 4
	    for (int k=0; k<j; k++)
	    {
                r[i][j] -= mu[j][k] * r[i][k];
	    }

	    // Line 5
	    mu[i][j] = r[i][j] / r[j][j];
	}

	// Line 6
	s[i][0] = ppA_gram[i][i];

	for (int j=1; j<=i; j++)
	{
	    s[i][j] = s[i][j-1] - (mu[i][j-1] * r[i][j-1]);
	}

	// Line 7
	r[i][i] = s[i][i];
    }

    return;
}


// The lazy size reduce algorithm from Nguyen and Stehle 2009
// Takes one index k denoting the vector being size reduced by all earlier vectors
template <typename Z_T>
inline void lazy_size_reduce (
		ZZ_mat<Z_T> &A,
	       	ZZ_mat<Z_T> &A_gram,
	       	FP_mat<mpfr_t> &r,
	       	FP_mat<mpfr_t> &mu,
	       	FP_mat<mpfr_t> &s,
		mpfr_t &nint_mu_f,
                FP_NR<mpfr_t> &nint_mu_FP,
                mpz_t &nint_mu_z,
                vector<Z_NR<Z_T>> &nint_mu_vec,
	       	FP_NR<mpfr_t> eta_bar,
	       	int k,
	       	int d,
		FP_NR<mpfr_t> mu_tmp,
		Z_NR<mpz_t> z_tmp1,
		Z_NR<mpz_t> z_tmp2
		)
{
    bool size_reduction_required = false;

    int flag = 0;
    // Size reducing only the index k
    // NS09 Fig 5, Line 2
    compute_CFA (A_gram, r, mu, s, d, k, k);

    for (int i=0; i<k; i++)
    {
        mu_tmp.abs(mu[k][i]); 
        if (mu_tmp > eta_bar)
        { 
	    size_reduction_required = true;
	    break;
       	}
    }

    // Continue size reducing until max_abs_mu <= eta_bar
    while (size_reduction_required)
    {
	// If flag == 1, multiple rounds of size reduction required
	// Therefore, recompute CFA
	if (flag == 1)
	{
	    cout << "Flag = 1" << endl;
	    compute_Gram_Row (A, A_gram, k, d);
	    compute_CFA (A_gram, r, mu, s, d, k, k);
	}
	// NS09, Fig 5, Line 3 (for loop)
	for (int i=k-1; i>=0; i--)
	{
	    // NS09, Fig 5, Line 4
            mpfr_round(nint_mu_f, mu[k][i].get_data());
            nint_mu_FP = nint_mu_f;
            mpfr_get_z (nint_mu_z, nint_mu_f, MPFR_RNDN);
            nint_mu_vec[i] = nint_mu_z;

	    // NS09, Fig 5, Line 5 (mu update)
	    if (!(nint_mu_vec[i].is_zero()))
	    {
                for (int j=0; j<i; j++) 
	        {
                    mu[k][j] = mu[k][j] - (nint_mu_FP * mu[i][j]);
	        }

	        // Size reducing the basis vector b_k
	        // NS09, Fig 5, Line 6 (i) - Update basis vector
                for (int j=0; j<d; j++)
                {
                    A[k][j].submul(nint_mu_vec[i], A[i][j]);
	        }
	    }
	}

	// Update Gram matrix entries including index k
	// Need to update <b[k], b[j]> for j<k and <b[j],b[k]> for k<j
	// Using NS09, part just before Sec 3.4
	z_tmp2 = 2;

	for (int l=0; l<k; l++)
	{
	    z_tmp1.mul(nint_mu_vec[l], nint_mu_vec[l]);
	    z_tmp1.mul(z_tmp1, A_gram[l][l]);
	    A_gram[k][k].add(A_gram[k][k],  z_tmp1);

	    z_tmp1 = 0;
	    z_tmp1.mul(nint_mu_vec[l], A_gram[k][l]);
	    z_tmp1.mul(z_tmp2, z_tmp1);
	    A_gram[k][k].sub(A_gram[k][k], z_tmp1);
	}

	for (int l=0; l<k; l++)
	{
	    for (int j=0; j<l; j++)
	    {
		z_tmp1 = 0;
	        z_tmp1.mul(z_tmp2, nint_mu_vec[l]);
		z_tmp1.mul(z_tmp1, nint_mu_vec[j]);
		z_tmp1.mul(z_tmp1, A_gram[l][j]);
		A_gram[k][k].add(A_gram[k][k], z_tmp1);
	    }
	}

	// Note that this updates entire Gram matrix
	// Therefore no additional recomputation of Gram required
	for (int l=0; l<k; l++)
	{
	    for (int j=0; j<k; j++)
	    {
		z_tmp1 = 0;
		if (l<j)
		{
	            z_tmp1.mul(nint_mu_vec[j], A_gram[j][l]);
		    A_gram[k][l].sub(A_gram[k][l], z_tmp1);
		} else
		{
		    z_tmp1.mul(nint_mu_vec[j], A_gram[l][j]);
		    A_gram[k][l].sub(A_gram[k][l], z_tmp1);
		}
	    }
	}

	for (int l=k+1; l<d; l++) 
	{
	    for (int j=0; j<k; j++)
	    {
		z_tmp1 = 0;
	        z_tmp1.mul(nint_mu_vec[j], A_gram[l][j]);
		A_gram[l][k].sub(A_gram[l][k], z_tmp1);
	    }
	}
	size_reduction_required = false;
	
	compute_CFA (A_gram, r, mu, s, d, k, k);

        for (int i=0; i<k; i++)
        {
            mu_tmp.abs(mu[k][i]); 
            if (mu_tmp > eta_bar)
            { 
		cout << "mu[" << k << "][" << i << "] = " << mu[k][i] << " > " << eta_bar << endl; 
	        size_reduction_required = true;
	        flag = 1;
	        break;
       	    }
        }
    }
    return;
}

// The lazy size reduce algorithm from NS09
// Takes one index k denoting the vector being size reduced by all earlier vectors
template <typename T>
inline void lazy_size_reduce (
		T ** &ppA,
	       	long long int ** &ppA_gram,
	       	FP_mat<mpfr_t> &r,
	       	FP_mat<mpfr_t> &mu,
	       	FP_mat<mpfr_t> &s,
		mpfr_t &nint_mu_f,
                FP_NR<mpfr_t> &nint_mu_FP,
                T * &pNint_mu,
	       	FP_NR<mpfr_t> eta_bar,
	       	int k,
	       	int d,
		FP_NR<mpfr_t> mu_tmp,
		Z_NR<mpz_t> z_tmp1,
		Z_NR<mpz_t> z_tmp2
		)
{
    bool size_reduction_required = false;
    int flag = 0;

    // Line 1
    // Size reducing only the index k
    compute_CFA (ppA_gram, r, mu, s, d, k, k);

    // Line 3 (Check max |mu|)
    for (int i=0; i<k; i++)
    {
        mu_tmp.abs(mu[k][i]); 
        if (mu_tmp > eta_bar)
        { 
	    size_reduction_required = true;
            break;
       	}
    }

    // Continue size reducing until max_abs_mu < eta_bar
    while (size_reduction_required)
    {
        // Line 3 (for loop)
	if (flag == 1)
	{
            compute_CFA (ppA_gram, r, mu, s, d, k, k);
	}
	for (int i=k-1; i>=0; i--)
	{
	    // Line 4
            mpfr_round(nint_mu_f, mu[k][i].get_data());
	    nint_mu_FP = nint_mu_f;
            pNint_mu[i] = mpfr_get_si (nint_mu_f, MPFR_RNDN);

	    // Line 5 (mu update)
	    for (int j=0; j<i; j++) 
	    {
                mu[k][j] = mu[k][j] - (nint_mu_FP * mu[i][j]);
	    }

	}
	// Line 6 (i) - Update basis vector
	for (int i=0; i<k; i++)
	{
	    for (int j=0; j<d; j++)
	    {
                ppA[k][j] = ppA[k][j] - (pNint_mu[i] * ppA[i][j]);
	    }
	}

	// Update Gram matrix entries including index k
	// Need to update <b[k], b[j]> for j<k and <b[j],b[k]> for k<j
	// Using Nguyen and Stehle 2009, part just before Sec 3.4
	for (int l=0; l<k; l++)
	{
	    ppA_gram[k][k] = ppA_gram[k][k] + (pNint_mu[l] * pNint_mu[l] * ppA_gram[l][l]);
	    ppA_gram[k][k] = ppA_gram[k][k] - (2 * pNint_mu[l] * ppA_gram[k][l]);
	}

	for (int l=0; l<k; l++)
	{
	    for (int j=0; j<l; j++)
	    {
		ppA_gram[k][k] = ppA_gram[k][k] + (2 * pNint_mu[j] * pNint_mu[l] * ppA_gram[l][j]);
	    } 
	}

	for (int l=0; l<k; l++)
	{
	    for (int j=0; j<k; j++)
	    {
		if (l<j)
		{
	            ppA_gram[k][l] = ppA_gram[k][l] - (pNint_mu[j] * ppA_gram[j][l]);
		} else
		{
		    ppA_gram[k][l] = ppA_gram[k][l] - (pNint_mu[j] * ppA_gram[l][j]); 
		}
	    }
	}
	
	for (int l=k+1; l<d; l++)
	{
	    for (int j=0; j<k; j++)
	    {
	        ppA_gram[l][k] = ppA_gram[l][k] - (pNint_mu[j] * ppA_gram[l][j]);
	    }
	}
	size_reduction_required = false;
	// End of Gram matrix updates
	
	compute_CFA (ppA_gram, r, mu, s, d, k, k);
        for (int i=0; i<k; i++)
        {
            mu_tmp.abs(mu[k][i]); 
            if (mu_tmp > eta_bar)
            { 
	        flag = 1;
	        size_reduction_required = true;
	        break;
       	    }
        }
    }

    return;
}

// Size reduction vector b_i with b_j
template <typename Z_T> 
inline void size_reduce (ZZ_mat<Z_T> &A, int d, int i, int j, Z_NR<Z_T> &nint_mu)
{
    nint_mu.neg(nint_mu); 
    for (int l=0; l<d; l++)
    {
        A[i][l].addmul(nint_mu, A[j][l]);
    }
    return;
}

// Size reduction vector b_i with b_j
template <typename T>
inline void size_reduce (T ** &ppA, int d, int i, int j, long int nint_mu)
{
    for (int l=0; l<d; l++)
    {
        ppA[i][l] -= ((T) nint_mu * ppA[j][l]);
    }
    return;
}

// Size reduction vector b_i with b_j
template <typename Z_T> 
inline void size_reduce_wrapper (
                ZZ_mat<Z_T> &A,
                FP_mat<mpfr_t> &mu,
                int d,
                int i,
                int j,
                mpfr_t &nint_mu_f, 
                FP_NR<mpfr_t> &nint_mu_FP,
                mpz_t &nint_mu_z,
                Z_NR<mpz_t> &max_nint_mu,
                Z_NR<mpz_t> &nint_mu,
                int start_index = 0
                )

{        
    mpfr_round(nint_mu_f, mu[i][j].get_data());
    mpfr_get_z (nint_mu_z, nint_mu_f, MPFR_RNDN);
    nint_mu = nint_mu_z;
    max_nint_mu = (nint_mu > max_nint_mu) ? nint_mu : max_nint_mu;

    size_reduce (A, d, i, j, nint_mu);

    // Update mu_ij
    nint_mu_FP = nint_mu_f;
    for (int l=j-1; l>=start_index; l--) 
    {
        mu[i][l] = mu[i][l] - nint_mu_FP * mu[j][l];
    }
    mu[i][j] = mu[i][j] - nint_mu_FP;

    return;
}

// Size reduction vector b_i with b_j
template <typename T> 
inline void size_reduce_wrapper (
                T ** &ppA,
                FP_mat<mpfr_t> &mu,
                int d,
                int i,
                int j,
                mpfr_t &nint_mu_f, 
                FP_NR<mpfr_t> &nint_mu_FP,
                long int &max_nint_mu,
                long int &nint_mu,
                int start_index = 0
                )
{
    mpfr_round(nint_mu_f, mu[i][j].get_data());
    nint_mu = mpfr_get_si (nint_mu_f, MPFR_RNDN);
    max_nint_mu = (nint_mu > max_nint_mu) ? nint_mu : max_nint_mu;

    size_reduce (ppA, d, i, j, nint_mu);
    nint_mu_FP = nint_mu_f;
            
    // Update mu_ij
    for (int l=j-1; l>=start_index; l--) 
    {
        mu[i][l] = mu[i][l] - nint_mu_FP * mu[j][l];
    }
    mu[i][j] = mu[i][j] - nint_mu_FP;

    return;
}

// Finds the integer datatype for lattice reduction of preprocessed bases
inline int getIntegerType (ZZ_mat<mpz_t> A, int d)
{
    int int_type = TYPE_MPZ_T;
    Z_NR<mpz_t> max_element;
    max_element.abs(A[0][0]);
    Z_NR<mpz_t> mpz_temp;

    if (d > 450)
        return int_type;

    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
        {
            mpz_temp.abs(A[i][j]);
            if (mpz_temp > max_element)
            {
                max_element = mpz_temp; 
            }
        }
    }

    if (max_element <= MAX_LONG_LONG_INT)
    {
        int_type = TYPE_LONG_LONG_INT;
    }

    if (max_element <= MAX_LONG_INT)
    {
        int_type = TYPE_LONG_INT;
    }

    if (max_element <= MAX_INT)
    {
        int_type = TYPE_INT;
    }
    cout << "int_type = " << int_type << endl;
    return int_type;
}

// Copy basis from ZZ_mat<mpz_t> to FP_mat<mpfr_t>
inline void copy_basis_mpz_mpfr (ZZ_mat<mpz_t> &A, FP_mat<mpfr_t> &A_mpfr, FP_mat<mpfr_t> &gso, FP_mat<mpfr_t> &mu, int d) 
{
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++) 
        {
            A_mpfr(i,j).set_z(A(i,j));
            gso(i,j).set_z(A(i,j));
            mu[i][j] = 0.0;
        }
    }
    return; 
}

// Copy basis from Z_T ** to FP_mat<mpfr_t>
template <typename Z_T>
inline void copy_basis_Z_T_mpfr (Z_T ** &ppA, FP_mat<mpfr_t> &A_mpfr, FP_mat<mpfr_t> &gso, FP_mat<mpfr_t> &mu, int d) 
{
    Z_NR<mpz_t> temp;
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++) 
        {
            temp = ppA[i][j];
            A_mpfr(i,j).set_z(temp);
            gso(i,j).set_z(temp);
            mu[i][j] = 0.0;
        }
    }
    return; 
}

// Copy basis to another variable name of the same datatype
template<typename T>
inline void copy_basis (T &A_destination, T &A_source, int d) 
{
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++) 
        {
            A_destination[i][j] = A_source[i][j];
        }
    }
    return; 
}

// Copy basis from ZZ_mat<mpz_t> to Z_T **
template <typename Z_T>
void copy_basis_mpz_Z_T (Z_T ** &ppA, ZZ_mat<mpz_t> &A_source, int d) 
{
    long temp;
    ppA = (Z_T **) calloc (d, sizeof(Z_T *));
    for (int i=0; i<d; i++)
    {
        ppA[i] = (Z_T *) calloc (d, sizeof(Z_T));
    }
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++) 
        {
            temp = A_source[i][j].get_si();
            ppA[i][j] = (Z_T) temp;
        }
    }
    return; 
}

// Computing the initial GSO, mu matrices, and B vector
template <typename FP_T>
inline void compute_GSO (
                ZZ_mat<mpz_t> &A,
                FP_mat<FP_T> &gso,
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B ,
                int d,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    // Copy mpz_t basis to mpfr_t basis A_mpfr, GSO
    // Also set each mu[i][j]=0
    FP_mat<FP_T> A_FP(d,d);
    copy_basis_mpz_mpfr (A, A_FP, gso, mu, d);

    FP_NR<FP_T> mu_num, mu_neg;
    for (int i=start_index; i<=end_index; i++) 
    {
        for (int j=i-1; j>=start_index; j--) 
        {
            // Compute numerator and denominator of mu
            A_FP[i].dot_product (mu_num, gso[j]);

            // Compute and store mu_ij in mu matrix
            mu[i][j] = mu_num / B[j];         

            // Update GSO as gso[i] - mu*gso[j]
            mu_neg = -mu[i][j];
            for (int l=0; l<d; l++) 
            {
                gso[i][l].addmul(mu_neg, gso[j][l]);
            }
        }
        // Once GSO updated for all j, store <b_j*,b_j*> in B
        // Then use B[j] instead of denominator computed above? 
        gso[i].dot_product (B[i], gso[i]);
    }
    return;
}

// Computing the initial GSO, mu matrices, and B vector
template <typename Z_T, typename FP_T>
inline void compute_GSO (
                Z_T ** &ppA,
                FP_mat<FP_T> &gso,
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B ,
                int d,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    // Copy mpz_t basis to mpfr_t bases A_mpfr and gso
    // Set all mu[i][j] = 0.0
    FP_mat<FP_T> A_FP(d,d);
    copy_basis_Z_T_mpfr (ppA, A_FP, gso, mu, d);

    FP_NR<FP_T> mu_num, mu_neg;

    for (int i=start_index; i<=end_index; i++) 
    {
        for (int j=i-1; j>=start_index; j--) 
        {
            // Compute numerator and denominator of mu
            A_FP[i].dot_product (mu_num, gso[j]);

            // Compute and store mu_ij in mu matrix
            mu[i][j] = mu_num / B[j];         

            // Update GSO as gso[i] - mu*gso[j]
            mu_neg = -mu[i][j];
            for (int l=0; l<d; l++) 
            {
                gso[i][l].addmul(mu_neg, gso[j][l]);
            }
        }
        // Once GSO updated for all j, store <b_j*,b_j*> in B
        gso[i].dot_product (B[i], gso[i]);
    }
    return;
}

// Compute the log of the potential of a basis
template <typename Z_T>
inline FP_NR<mpfr_t> basisPotential (Z_T ** ppBasis, int d) {
    FP_mat<mpfr_t> gso(d,d);
    FP_mat<mpfr_t> mu(d,d);
    vector<FP_NR<mpfr_t>> B;
    B.resize(d);
    FP_NR<mpfr_t> gsoVectorNormSquared, potTerm;
    double power;

    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
	{
	    gso[i][j] = 0.0;
	    mu[i][j] = 0.0;
	}
	B[i] = 0.0;
    }
    // Compute GSO
    compute_GSO (ppBasis, gso, mu, B, d, 0, 0);

   // Initialise log of potential = 0
   FP_NR<mpfr_t> potential;
   potential = 0.0;
   for (int i=0; i<d; i++) {
      // The potential term for b_i is ||b_i*||^(2(n-i)|| (since indexing from 0)
      // This is the same as (||b_i*||^2)^(n-i)
      // Log of the potential term is (n-i) * log(||b_i*||^2)
      gso[i].dot_product(gsoVectorNormSquared, gso[i]);
      power = d-i;
      potTerm = log(gsoVectorNormSquared);
      potTerm.mul(potTerm, power);  
      potential.add(potential, potTerm);
   }
   return potential;
}


// Using the algorithm from YY18 for GSO update
template <typename FP_T>
inline int deep_LLL_GSO_update (
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B,
                int k,
                int i,
                int d,
                int start_index = 0,
                int end_index = 0
                ) 
{
    end_index = (!end_index) ? d-1 : end_index;

    // Initialise and set P, S, D vectors
    vector<FP_NR<FP_T>> P, S, D;
    P.resize (d, 0.0);
    S.resize (d, 0.0);
    D.resize (d, 0.0);

    // Initialise and set constants
    FP_NR<FP_T> T = 0.0;
    FP_NR<FP_T> fptmp = 0.0;

    P[k] = B[k];
    D[k] = B[k];

    for (int j=k-1; j>=i; j--) 
    {
        P[j] = mu[k][j] * B[j];
        D[j] = D[j+1] + (mu[k][j] * P[j]);
    }

    for (int j=k; j>i; j--) 
    {
        T = mu[k][j-1] / D[j];

        for (int l=end_index; l>k; l--)
        {
            S[l] = S[l] + (mu[l][j] * P[j]);
            mu[l][j] = mu[l][j-1] - (T * S[l]);            
        }

        for (int l=k; l>=j+2; l--)
        {
            S[l] = S[l] + (mu[l-1][j] * P[j]);
            mu[l][j] = mu[l-1][j-1] - (T * S[l]);    
        }

        if (j != k)
        {
            S[j+1] = P[j];
            mu[j+1][j] = mu[j][j-1] - (T * S[j+1]);
        }
        
    }

    T = 1.0 / D[i];

    for (int l=end_index; l>k; l--)
    {
        mu[l][i] = T * (S[l] + mu[l][i] * P[i]);
    }

    for (int l=k; l>=i+2; l--) 
    {
        mu[l][i] = T * (S[l] + mu[l-1][i] * P[i]);
    }

    mu[i+1][i] = T * P[i];

    for (int j=start_index; j<i; j++) 
    {
        fptmp = mu[k][j];

        for (int l=k; l>i; l--)
        {
            mu[l][j] = mu[l-1][j];
        }

        mu[i][j] = fptmp;
    }

    for (int j=k; j>i; j--)
    {
        B[j] = (D[j] * B[j-1]) / D[j-1];
    }

    B[i] = D[i];

    return 0;
}

template <typename FP_T>
inline FP_NR<FP_T> compute_SS (
                FP_NR<FP_T> &squared_sum,
                vector<FP_NR<FP_T>> &B,
                int d,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    squared_sum = 0.0; 
    for (int i=start_index; i<=end_index; i++) 
    {
        squared_sum += B[i];
    } 
    return squared_sum;
}

template <typename FP_T>
inline bool index_search_SS (
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B,
                int d,
                FP_NR<FP_T> squared_sum,
                FP_NR<FP_T> delta_p,
                int &k_max,
                int &i_max,
                int start_index = 0,
                int end_index = 0
                ) 
{
    end_index = (!end_index) ? d-1 : end_index;

    k_max = INDEX_SEARCH_FAIL;
    i_max = start_index; 

    FP_NR<FP_T> SS_change = 0.0;
    FP_NR<FP_T> SS_change_max = delta_p * squared_sum;
    FP_NR<FP_T> projection_l_k = 0.0;
    FP_NR<FP_T> ratio_projection_l_k = 0.0;
    FP_NR<FP_T> inv_ratio_projection_l_k = 0.0;

    for (int k=start_index+1; k<=end_index; k++) 
    {
        // projection_l_k = \pi_{k-1}(\vb_{k}) (initially)
        // D^(k)_{k-1} = (projection_l_k)^2
        projection_l_k = B[k] + (mu[k][k-1] * mu[k][k-1] * B[k-1]);
        // ratio_projection_l_k = projection_l_k / projection_l_l
        ratio_projection_l_k = projection_l_k / B[k-1];

        SS_change = mu[k][k-1] * mu[k][k-1] * (1.0 - ratio_projection_l_k);  

        if (SS_change > SS_change_max) 
        {
            SS_change_max = SS_change;
            k_max = k;
            i_max = k-1;
        }

        for (int i = k-2; i >= start_index; i--) 
        {
            // projection_l_k = \pi_{i}(\vb_{k}) 
            // D^(k)_{i} = (projection_l_k)^2
            projection_l_k += (mu[k][i] * mu[k][i] * B[i]);
            // ratio_projection_l_k = projection_l_k / projection_l_l
            inv_ratio_projection_l_k = B[i] / projection_l_k;
            SS_change += mu[k][i] * mu[k][i] * B[i] * (inv_ratio_projection_l_k - 1.0);  
            
            if (SS_change > SS_change_max) 
            {
                SS_change_max = SS_change;
                k_max = k;
                i_max = i;
            }
        }
    }
    
    if (k_max != INDEX_SEARCH_FAIL)
    {
        return true;
    } 
    return false;
}

// Checks if the basis is SS-DeepLLL reduced
template <typename FP_T>
void ssgg_reduced_check_internal (
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B,
                FP_NR<FP_T> delta_p,
                int precision,
                int d,
                FP_NR<FP_T> eta,
                int &k,
                int &i,
                FP_NR<FP_T> &squared_sum,
                bool &size_reduction_required,
                bool &deep_insertion_required,
                bool &index_search_required,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    size_reduction_required = true;
    deep_insertion_required = true;
    index_search_required = true;
    
    for (k=start_index+1; k<=end_index; k++) 
    {
        for (int j=k-1; j>=start_index; j--)
        {
            if (mu[k][j] < -eta || mu[k][j] > eta)
            {
                cout << "ERROR: BASIS NOT SIZE REDUCED" << endl;
                cout << "mu[k][j] = " << mu[k][j] << endl;
                return;
            }
        }
    }
    
    size_reduction_required = false;
    k = INDEX_SEARCH_FAIL;
    i = -1;

    squared_sum = compute_SS (squared_sum, B, d, start_index, end_index);
    deep_insertion_required = index_search_SS (mu, B, d, squared_sum, delta_p, k, i, start_index, end_index);
    index_search_required = false;

    if (deep_insertion_required) 
    {
        cout << "ERROR: BASIS SIZE REDUCED BUT NOT SS REDUCED" << endl;
    }

    return;
}

// Checks if the basis is Pot-DeepLLL reduced
template <typename FP_T>
void potgg_reduced_check_internal (
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B,
                FP_NR<FP_T> delta,
                int precision,
                int d,
                FP_NR<FP_T> eta,
                int &k,
                int &i,
                bool &size_reduction_required,
                bool &deep_insertion_required,
                bool &index_search_required,
		int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    	
    size_reduction_required = true;
    deep_insertion_required = true;
    index_search_required = true;
    
    for (k=start_index+1; k<=end_index; k++)
    {
        for (int j=k-1; j>=start_index; j--)
        {
            if (mu[k][j] < -eta || mu[k][j] > eta)
            {
                cout << "ERROR: BASIS NOT SIZE REDUCED" << endl;
                cout << "mu[k][j] = " << mu[k][j] << endl;
                return;
            }
        }
    }
    size_reduction_required = false;
    k = INDEX_SEARCH_FAIL;
    i = -1;

    deep_insertion_required = index_search_Pot (mu, B, d, delta, k, i, start_index, end_index);
    index_search_required = false;

    if (deep_insertion_required) 
    {
        cout << "ERROR: BASIS SIZE REDUCED BUT NOT POT REDUCED" << endl;
    }

    return;
}

// Compute the Root Hermite Factor of a Lattice Basis
template <typename FP_T>
FP_NR<FP_T> compute_RHF (
                vector<FP_NR<FP_T>> &B, 
                int d,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;
   
    // Compute volume of lattice as product of GSO lengths
    // Do this by taking inner product of GSO vector with itself and square root
    FP_NR<FP_T> volumeSquared = 1.0;
    FP_NR<FP_T> volume;
   
    for (int i=start_index; i<=end_index; i++) {
       volumeSquared *= B[i];
    }

    volume.sqrt(volumeSquared);

    // Compute the dth root of the volume using root function (from fplll)
    FP_NR<FP_T> dthRootVolume;
    dthRootVolume.root(volume, (end_index - start_index + 1));

    // Compute the norm of b1 using the square root of inner product
    FP_NR<FP_T> b1Norm; 
    b1Norm.sqrt(B[start_index]);

    // Compute the Hermite Factor as gamma = ||b_1|| / Vol(L)^(1/n)
    FP_NR<FP_T> hermiteFactor;
    hermiteFactor = (b1Norm / dthRootVolume);

    // Compute the Root Hermite Factor (RHF) as gamma^(1/n)
    FP_NR<FP_T> RHF;
    RHF.root(hermiteFactor, (end_index - start_index + 1));

    return RHF;
}

// Checks if basis is LC-DeepLLL reduced
template <typename FP_T>
void lcgg_reduced_check_internal (
		ZZ_mat<mpz_t> &A,
		FP_mat<FP_T> &mu,
		vector<FP_NR<FP_T>> &B,
		FP_NR<FP_T> delta,
		FP_NR<FP_T> &delta_kj,
		FP_NR<FP_T> &delta_min,
		int precision,
		int d,
		FP_NR<FP_T> eta,
		int beta, 
		Z_NR<mpz_t> &C_int,
		vector<FP_NR<FP_T>> &C,
		int &k,
		int &i,
		bool &size_reduction_required,
		bool &deep_insertion_required,
		bool &index_search_required,
		int start_index = 0,
		int end_index = 0
		)
{
    end_index = (!end_index) ? d-1 : end_index;

    size_reduction_required = true;
    deep_insertion_required = true;
    index_search_required = true;
    
    for (int l=start_index+1; l<=end_index; l++) 
    {
	for (int j=l-1; j>=start_index; j--)
	{
	    if (mu[l][j] < -eta || mu[l][j] > eta)
	    {
		cout << "ERROR: BASIS NOT SIZE REDUCED" << endl;
		cout << "k = " << l << "   i = " << j << "   mu[k][j] = " << mu[l][j] << endl;
		return;
	    }
	}
    }
    size_reduction_required = false;
    
    for (int l=0; l<d; l++)
    {
        A[l].dot_product (C_int, A[l]);
	C_int.abs(C_int);
	C[l].set_z (C_int, MPFR_RNDN);
    }

    deep_insertion_required = index_search_LC_GG (mu, B, d, C, delta, delta_kj, delta_min, k, i, beta, start_index, end_index);
    if (deep_insertion_required) 
    {
        cout << "ERROR: BASIS SIZE REDUCED BUT NOT LC-DEEP REDUCED" << endl;
        cout << "k = " << k << "\ti = " << i << endl;
        return;
    }
    
    index_search_required = false;

    return;
}

// Checks if basis is LC-DeepLLL reduced
template <typename FP_T, typename Z_T>
void lcgg_reduced_check_internal (
		Z_T ** &ppA,
		FP_mat<FP_T> &mu,
		vector<FP_NR<FP_T>> &B,
		FP_NR<FP_T> delta,
		FP_NR<FP_T> &delta_kj,
		FP_NR<FP_T> &delta_min,
		int precision,
		int d,
		FP_NR<FP_T> eta,
		int beta, 
		long long int &C_int,
		vector<FP_NR<FP_T>> &C,
		int &k,
		int &i,
		bool &size_reduction_required,
		bool &deep_insertion_required,
		bool &index_search_required,
		int start_index = 0,
		int end_index = 0
		)
{
    end_index = (!end_index) ? d-1 : end_index;

    size_reduction_required = true;
    deep_insertion_required = true;
    index_search_required = true;

    Z_T inner;
    inner = 0;
    Z_NR<mpz_t> inner_ZNR;

    for (int l=start_index+1; l<=end_index; l++) 
    {
	for (int j=l-1; j>=start_index; j--)
	{
	    if (mu[l][j] < -eta || mu[l][j] > eta)
	    {
		cout << "ERROR: BASIS NOT SIZE REDUCED" << endl;
		cout << "k = " << l << "   i = " << j << "   mu[k][j] = " << mu[l][j] << endl;
		return;
	    }
	}
    }
    size_reduction_required = false;
    
    for (int l=0; l<d; l++)
    {
        inner = 0;
	// Computing inner product using standard C integer datatype
	for (int j=0; j<d; j++)
	{
	    inner += ppA[l][j] * ppA[l][j];
	}
	if (inner < 0)
	{ 
	    cout << "inner before = " << inner << endl;
	    inner = -inner;
	    cout << "inner after = " << inner << endl;
	}
	// Storing inner in inner_ZNR, the Z_NR<mpz_t> datatype
	// This allows for converting to FP_NR<mpfr_t> using set_z FPLLL function
	inner_ZNR = inner;
	C[l].set_z (inner_ZNR, MPFR_RNDN);
    }

    deep_insertion_required = index_search_LC_GG (mu, B, d, C, delta, delta_kj, delta_min, k, i, beta, start_index, end_index);

    if (deep_insertion_required) 
    {
        cout << "ERROR: BASIS SIZE REDUCED BUT NOT LC-DEEP REDUCED" << endl;
        cout << "k = " << k << "\ti = " << i << endl;
        return;
    }
    return;
}

#endif
