/* 
 * *****************************************************************************
 * Filename: pgg.h
 * Authors : Sanjay Bhattacherjee and Jack Moyler
 * *****************************************************************************
*/

#include "common.h"

#ifndef PGG
#define PGG

#define LC_DELTA 0.99
// Maximum recursion depth for LC-PGG
// Set to 500 so no restriction up until d=1000
#define MAX_RECURSION_DEPTH 500

// Recursive search for LC-PGG
int pgg_recursive_search (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		vector<FP_NR<mpfr_t>> &B,
		ZZ_mat<mpz_t> &A,
		int i,
		int k,
		int d,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    int i_min = k+1;
    int i_left = k+1;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta_mat[i+1][i];
    int i1 = i;
    int k1 = i+1;

    // Check all delta[l][j] in range [i,k] for minDelta
    for (int l=i+2; l<=k; l++) 
    {
    	for (int j=i; j<l; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (minDelta > delta)
    {
        return i_min;
    }

    nInsertions++;
    // Deep insert b_{k1} in position i1 
    A.rotate_right(i1, k1);

    // Update the GSO     
    deep_LLL_GSO_update (mu, B, k1, i1, d, start_index, end_index); 

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 

    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search (delta_mat, mu, B, A, i, i1-1, d, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search (delta_mat, mu, B, A, k1+1, k, d, delta, nInsertions, recursionDepth);

    return i_min;
}

// Recursive search for LC-PGG
template<typename Z_T>
int pgg_recursive_search (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		vector<FP_NR<mpfr_t>> &B,
		Z_T ** ppA,
		int i,
		int k,
		int d,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    int i_min = k+1;
    int i_left = k+1;
    Z_T * pTemp;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta_mat[i+1][i];
    int i1 = i;
    int k1 = i+1;
    // Check all delta[l][j] in range [i,k] for minDelta
    for (int l=i+2; l<=k; l++) 
    {
    	for (int j=i; j<l; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (minDelta > delta)
    {
        return i_min;
    }

    nInsertions++;
    // Deep insert b_{k1} in position i1 
    pTemp = ppA[k1];
    for (int j=k1-1; j>=i1; j--)
    {
        ppA[j+1] = ppA[j];
    }
    ppA[i1] = pTemp;

    // Update the GSO     
    deep_LLL_GSO_update (mu, B, k1, i1, d, start_index, end_index); 

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 
    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search (delta_mat, mu, B, ppA, i, i1-1, d, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search (delta_mat, mu, B, ppA, k1+1, k, d, delta, nInsertions, recursionDepth);

    return i_min;
}

// Recursive search for LC-PGG
// CFA and LSR
int pgg_recursive_search_cfa (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		FP_mat<mpfr_t> &r,
		FP_mat<mpfr_t> &s,
		ZZ_mat<mpz_t> &A,
		ZZ_mat<mpz_t> &A_gram,
		int i,
		int k,
		int d,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    Z_NR<mpz_t> Gram_elem_tmp;
    int i_min = k+1;
    int i_left = k+1;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta_mat[i+1][i];
    int i1 = i;
    int k1 = i+1;
    // Check all delta[l][j] in range [i,k] for minDelta
    for (int l=i+2; l<=k; l++) 
    {
    	for (int j=i; j<l; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (minDelta > delta)
    {
        return i_min;
    }

    nInsertions++;
    // Deep insert b_{k1} in position i1 
    A.rotate_right(i1, k1);

    for (int j=0; j<i1; j++)
    {
        mu[i1][j] = mu[k1][j];
	r[i1][j] = r[k1][j];
    }
    r[i1][i1] = s[k1][i1];

    // Update Gram matrix
    for (int l=k1; l>i1; l--)
    {
        A_gram[l-1][l] = A_gram[l-1][l-1];
	A_gram[l-1][l-1] = A_gram[l][l-1];
	A_gram[l][l-1] = A_gram[l][l];
	A_gram[l][l] = 0;
	A_gram.swap_rows(l-1,l);

	Gram_elem_tmp = A_gram[l][l-1];
        for (int j=l-1; j>i1; j--)
	{
	    A_gram[l][j] = A_gram[l][j-1];
	}
	A_gram[l][i1] = Gram_elem_tmp;
    }

    for (int l=d-1; l>k1; l--)
    {
	Gram_elem_tmp = A_gram[l][k1];
               for (int j=k1; j>i1; j--)
	{
	    A_gram[l][j] = A_gram[l][j-1];
	}
	A_gram[l][i1] = Gram_elem_tmp;
    }

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 

    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search_cfa (delta_mat, mu, r, s, A, A_gram, i, i1-1, d, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search_cfa (delta_mat, mu, r, s, A, A_gram, k1+1, k, d, delta, nInsertions, recursionDepth);

    return i_min;
}

// Recursive search for LC-PGG
// CFA and LSR
template<typename Z_T>
int pgg_recursive_search_cfa (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		FP_mat<mpfr_t> &r,
		FP_mat<mpfr_t> &s,
		Z_T ** ppA,
		long long int ** ppA_gram,
		int i,
		int k,
		int d,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    long long int Gram_elem_tmp;
    Z_T * pTemp;
    long long int * pGramTemp;
    int i_min = k+1;
    int i_left = k+1;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta_mat[i+1][i];
    int i1 = i;
    int k1 = i+1;
    // Check all delta[l][j] in range [i,k] for minDelta
    for (int l=i+2; l<=k; l++) 
    {
    	for (int j=i; j<l; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (minDelta > delta)
    {
        return i_min;
    }

    nInsertions++;
    // Deep insert b_{k1} in position i1 
    pTemp = ppA[k1];
    for (int j=k1-1; j>=i1; j--)
    {
        ppA[j+1] = ppA[j];
    }
    ppA[i1] = pTemp;

    for (int j=0; j<i1; j++)
    {
        mu[i1][j] = mu[k1][j];
        r[i1][j] = r[k1][j];
    }
    r[i1][i1] = s[k1][i1];
    // Update Gram matrix
    for (int l=k1; l>i1; l--)
    {
        ppA_gram[l-1][l] = ppA_gram[l-1][l-1];
        ppA_gram[l-1][l-1] = ppA_gram[l][l-1];
        ppA_gram[l][l-1] = ppA_gram[l][l];
        ppA_gram[l][l] = 0;

        pGramTemp = ppA_gram[l];
        ppA_gram[l] = ppA_gram[l-1];
        ppA_gram[l-1] = pGramTemp;
	Gram_elem_tmp = ppA_gram[l][l-1];

 	for (int j=l-1; j>i1; j--)
	{
	    ppA_gram[l][j] = ppA_gram[l][j-1];
	}
	ppA_gram[l][i1] = Gram_elem_tmp;
    }

    for (int l=d-1; l>k1; l--)
    {
	Gram_elem_tmp = ppA_gram[l][k1];
        for (int j=k1; j>i1; j--)
	{
	    ppA_gram[l][j] = ppA_gram[l][j-1];
	}
	ppA_gram[l][i1] = Gram_elem_tmp;
    }

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 

    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search_cfa (delta_mat, mu, r, s, ppA, ppA_gram, i, i1-1, d, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search_cfa (delta_mat, mu, r, s, ppA, ppA_gram, k1+1, k, d, delta, nInsertions, recursionDepth);

    return i_min;
}

// Recursive search for LC-PGG
// Restricted insertions within block [i,k] (AKA local, two sided)
int pgg_recursive_search_local (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		vector<FP_NR<mpfr_t>> &B,
		ZZ_mat<mpz_t> &A,
		int i,
		int k,
		int d,
		int beta,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    int i_min = k+1;
    int i_left = k+1;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta;
    minDelta += 1;
    int ub = 0;
    int lb = 0;
    int i1 = k+1;
    int k1 = k+1;

    // Check all delta[l][j] in range [i,i+beta-1] and [max(beta,k-beta), k-1]
    for (int l=i+1; l<=k; l++) 
    {
	lb = max(i+beta-1, l-beta);
	ub = min(l,i+beta);
    	for (int j=i; j<ub; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}

	for (int j=lb; j<l; j++)
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (minDelta > delta)
    {
        return i_min;
    }

    nInsertions++;
    // Deep insert b_{k1} in position i1 
    A.rotate_right(i1, k1);

    // Update the GSO     
    deep_LLL_GSO_update (mu, B, k1, i1, d, start_index, end_index); 

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 

    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search_local (delta_mat, mu, B, A, i, i1-1, d, beta, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search_local (delta_mat, mu, B, A, k1+1, k, d, beta, delta, nInsertions, recursionDepth);

    return i_min;
}

// Recursive search for LC-PGG
// Restricted insertions within block [i,k] (local, two sided)
template<typename Z_T>
int pgg_recursive_search_local (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		vector<FP_NR<mpfr_t>> &B,
		Z_T ** ppA,
		int i,
		int k,
		int d,
		int beta,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    int i_min = k+1;
    int i_left = k+1;
    Z_T * pTemp;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta;
    minDelta += 1;
    int ub = 0;
    int lb = 0;
    int i1 = k+1;
    int k1 = k+1;

    // Check all delta[l][j] in range [i,i+beta-1] and [max(beta,k-beta), k-1]
    for (int l=i+1; l<=k; l++) 
    {
	lb = max(i+beta-1, l-beta);
	ub = min(l, i+beta);
    	for (int j=i; j<ub; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    	for (int j=lb; j<l; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (minDelta > delta)
    {
        return i_min;
    }

    nInsertions++;
    // Deep insert b_{k1} in position i1 
    pTemp = ppA[k1];
    for (int j=k1-1; j>=i1; j--)
    {
        ppA[j+1] = ppA[j];
    }
    ppA[i1] = pTemp;

    // Update the GSO     
    deep_LLL_GSO_update (mu, B, k1, i1, d, start_index, end_index); 

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 
    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search_local (delta_mat, mu, B, ppA, i, i1-1, d, beta, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search_local (delta_mat, mu, B, ppA, k1+1, k, d, beta, delta, nInsertions, recursionDepth);

    return i_min;
}

// Recursive search for LC-PGG
// CFA and LSR
// Restricted insertions within block [i,k] (local, two sided)
int pgg_recursive_search_local_cfa (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		FP_mat<mpfr_t> &r,
		FP_mat<mpfr_t> &s,
		ZZ_mat<mpz_t> &A,
		ZZ_mat<mpz_t> &A_gram,
		int i,
		int k,
		int d,
		int beta,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    Z_NR<mpz_t> Gram_elem_tmp;
    int i_min = k+1;
    int i_left = k+1;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta;
    minDelta += 1;
    int ub = 0;
    int lb = 0;
    int i1 = i;
    int k1 = i+1;

    // Check all delta[l][j] in range [i,i+beta-1] and [max(beta,k-beta), k-1]
    for (int l=i+1; l<=k; l++) 
    {
	lb = max(i+beta-1, l-beta);
	ub = min(l, i+beta);
    	for (int j=i; j<ub; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    	for (int j=lb; j<l; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (minDelta > delta)
    {
        return i_min;
    }

    nInsertions++;
    // Deep insert b_{k1} in position i1 
    A.rotate_right(i1, k1);

    // Update the Gram matrix     
    for (int j=0; j<i1; j++)
    {
        mu[i1][j] = mu[k1][j];
	r[i1][j] = r[k1][j];
    }
    r[i1][i1] = s[k1][i1];

    // Update Gram matrix
    for (int l=k1; l>i1; l--)
    {
        A_gram[l-1][l] = A_gram[l-1][l-1];
	A_gram[l-1][l-1] = A_gram[l][l-1];
	A_gram[l][l-1] = A_gram[l][l];
	A_gram[l][l] = 0;
	A_gram.swap_rows(l-1,l);

	Gram_elem_tmp = A_gram[l][l-1];
        for (int j=l-1; j>i1; j--)
	{
	    A_gram[l][j] = A_gram[l][j-1];
	}
	A_gram[l][i1] = Gram_elem_tmp;
    }

    for (int l=d-1; l>k1; l--)
    {
	Gram_elem_tmp = A_gram[l][k1];
               for (int j=k1; j>i1; j--)
	{
	    A_gram[l][j] = A_gram[l][j-1];
	}
	A_gram[l][i1] = Gram_elem_tmp;
    }

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 

    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search_local_cfa (delta_mat, mu, r, s, A, A_gram, i, i1-1, d, beta, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search_local_cfa (delta_mat, mu, r, s, A, A_gram, k1+1, k, d, beta, delta, nInsertions, recursionDepth);

    return i_min;
}

// Recursive search for LC-PGG
// CFA and LSR
// Restricted insertions within block [i,k] (local, two sided)
template<typename Z_T>
int pgg_recursive_search_local_cfa (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		FP_mat<mpfr_t> &r,
		FP_mat<mpfr_t> &s,
		Z_T ** ppA,
		long long int ** ppA_gram,
		int i,
		int k,
		int d,
		int beta,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    long long int Gram_elem_tmp;
    Z_T * pTemp;
    long long int * pGramTemp;
    int ub = 0;
    int lb = 0;
    int i_min = k+1;
    int i_left = k+1;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta;
    minDelta += 1;
    int i1 = k+1;
    int k1 = k+1;
    // Check all delta[l][j] in range [i,i+beta-1] and [max(beta,k-beta), k-1]
    for (int l=i+1; l<=k; l++) 
    {
	lb = max(i+beta-1, l-beta);
	ub = min(l, i+beta);
    	for (int j=i; j<ub; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    	//for (int j=lb; j<k; j++) 
    	for (int j=lb; j<l; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (minDelta > delta)
    {
        return i_min;
    }
    // cout << "minDelta = " << minDelta << "   k = " << k1 << "   i = " << i1 << endl;

    // Deep insert b_{k1} in position i1 
    nInsertions++;
    // Deep insert b_k in position i
    // rotate_right function from fplll
    // (v[first],...,v[last]) becomes (v[last],v[first],...,v[last-1])
    pTemp = ppA[k1];
    for (int j=k1-1; j>=i1; j--)
    {
        ppA[j+1] = ppA[j];
    }
    ppA[i1] = pTemp;

    for (int j=0; j<i1; j++)
    {
        mu[i1][j] = mu[k1][j];
        r[i1][j] = r[k1][j];
    }
    r[i1][i1] = s[k1][i1];
    // Update Gram matrix
    for (int l=k1; l>i1; l--)
    {
        ppA_gram[l-1][l] = ppA_gram[l-1][l-1];
        ppA_gram[l-1][l-1] = ppA_gram[l][l-1];
        ppA_gram[l][l-1] = ppA_gram[l][l];
        ppA_gram[l][l] = 0;
        // SWAP ROWS l and l-1
        pGramTemp = ppA_gram[l];
        ppA_gram[l] = ppA_gram[l-1];
        ppA_gram[l-1] = pGramTemp;
	Gram_elem_tmp = ppA_gram[l][l-1];

 	for (int j=l-1; j>i1; j--)
	{
	    ppA_gram[l][j] = ppA_gram[l][j-1];
	}
	ppA_gram[l][i1] = Gram_elem_tmp;
    }

    for (int l=d-1; l>k1; l--)
    {
	Gram_elem_tmp = ppA_gram[l][k1];
        for (int j=k1; j>i1; j--)
	{
	    ppA_gram[l][j] = ppA_gram[l][j-1];
	}
	ppA_gram[l][i1] = Gram_elem_tmp;
    }

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 

    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search_local_cfa (delta_mat, mu, r, s, ppA, ppA_gram, i, i1-1, d, beta, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search_local_cfa (delta_mat, mu, r, s, ppA, ppA_gram, k1+1, k, d, beta, delta, nInsertions, recursionDepth);

    return i_min;
}

// Recursive search for LC-PGG
// Restricted insertions within blocks [1,beta] and [k-beta,k] (global, two sided)
int pgg_recursive_search_global (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		vector<FP_NR<mpfr_t>> &B,
		ZZ_mat<mpz_t> &A,
		int i,
		int k,
		int d,
		int beta,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    int i_min = k+1;
    int i_left = k+1;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta;
    minDelta += 1;
    int ub = 0;
    int lb = 0;
    int i1 = k+1;
    int k1 = k+1;

    // Check all ppDelta in range [1,beta] and [k-beta, k-1]
    // Find min delta
    for (int l=i+1; l<=k; l++) 
    {
	lb = max(beta, l-beta);
	ub = min(l, beta);
    	for (int j=i; j<ub; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}

	for (int j=lb; j<l; j++)
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (i1 == k+1 && k1 == k+1)
    {
        return i_min;
    }
    if (minDelta > delta)
    {
        return i_min;
    }

    nInsertions++;
    // Deep insert b_{k1} in position i1 
    A.rotate_right(i1, k1);

    // Update the GSO     
    deep_LLL_GSO_update (mu, B, k1, i1, d, start_index, end_index); 

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 

    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search_global (delta_mat, mu, B, A, i, i1-1, d, beta, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search_global (delta_mat, mu, B, A, k1+1, k, d, beta, delta, nInsertions, recursionDepth);

    return i_min;
}

// Recursive search for LC-PGG
// Restricted insertions within block [i,k] (global, two sided)
template<typename Z_T>
int pgg_recursive_search_global (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		vector<FP_NR<mpfr_t>> &B,
		Z_T ** ppA,
		int i,
		int k,
		int d,
		int beta,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    int i_min = k+1;
    int i_left = k+1;
    Z_T * pTemp;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta;
    minDelta += 1;
    int ub = 0;
    int lb = 0;
    int i1 = k+1;
    int k1 = k+1;
    // Check all ppDelta in range [1,beta] and [k-beta, k-1]
    for (int l=i+1; l<=k; l++) 
    {
	lb = max(beta, l-beta);
	ub = min(l, beta);
    	for (int j=i; j<ub; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    	for (int j=lb; j<l; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (i1 == k+1 && k1 == k+1)
    {
        return i_min;
    }
    if (minDelta > delta)
    {
     
        return i_min;
    }

    nInsertions++;
    // Deep insert b_{k1} in position i1 
    pTemp = ppA[k1];
    for (int j=k1-1; j>=i1; j--)
    {
        ppA[j+1] = ppA[j];
    }
    ppA[i1] = pTemp;

    // Update the GSO     
    deep_LLL_GSO_update (mu, B, k1, i1, d, start_index, end_index); 

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 
    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search_global (delta_mat, mu, B, ppA, i, i1-1, d, beta, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search_global (delta_mat, mu, B, ppA, k1+1, k, d, beta, delta, nInsertions, recursionDepth);

    return i_min;
}

// Recursive search for LC-PGG
// CFA and LSR
// Restricted insertions within block [i,k] (global, two sided)
int pgg_recursive_search_global_cfa (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		FP_mat<mpfr_t> &r,
		FP_mat<mpfr_t> &s,
		ZZ_mat<mpz_t> &A,
		ZZ_mat<mpz_t> &A_gram,
		int i,
		int k,
		int d,
		int beta,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    Z_NR<mpz_t> Gram_elem_tmp;
    int i_min = k+1;
    int i_left = k+1;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta;
    minDelta += 1;
    int ub = 0;
    int lb = 0;
    int i1 = i;
    int k1 = i+1;
    // Check all ppDelta in range [1,beta] and [k-beta, k-1]
    for (int l=i+1; l<=k; l++) 
    {
	lb = max(beta, l-beta);
	ub = min(l, beta);
    	for (int j=i; j<ub; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    	for (int j=lb; j<l; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (i1 == k+1 && k1 == k+1) 
    {
        return i_min;
    }
    if (minDelta > delta)
    {
        return i_min;
    }

    nInsertions++;
    // Deep insert b_{k1} in position i1 
    A.rotate_right(i1, k1);

    for (int j=0; j<i1; j++)
    {
        mu[i1][j] = mu[k1][j];
	r[i1][j] = r[k1][j];
    }
    r[i1][i1] = s[k1][i1];

    // Update Gram matrix
    for (int l=k1; l>i1; l--)
    {
        A_gram[l-1][l] = A_gram[l-1][l-1];
	A_gram[l-1][l-1] = A_gram[l][l-1];
	A_gram[l][l-1] = A_gram[l][l];
	A_gram[l][l] = 0;
	A_gram.swap_rows(l-1,l);

	Gram_elem_tmp = A_gram[l][l-1];
        for (int j=l-1; j>i1; j--)
	{
	    A_gram[l][j] = A_gram[l][j-1];
	}
	A_gram[l][i1] = Gram_elem_tmp;
    }

    for (int l=d-1; l>k1; l--)
    {
	Gram_elem_tmp = A_gram[l][k1];
               for (int j=k1; j>i1; j--)
	{
	    A_gram[l][j] = A_gram[l][j-1];
	}
	A_gram[l][i1] = Gram_elem_tmp;
    }

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 

    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search_global_cfa (delta_mat, mu, r, s, A, A_gram, i, i1-1, d, beta, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search_global_cfa (delta_mat, mu, r, s, A, A_gram, k1+1, k, d, beta, delta, nInsertions, recursionDepth);

    return i_min;
}

// Recursive search for LC-PGG
// Z_T ** basis
// Restricted insertions within block [i,k] (global, two sided)
template<typename Z_T>
int pgg_recursive_search_global_cfa (
		FP_mat<mpfr_t> delta_mat,
		FP_mat<mpfr_t> &mu,
		FP_mat<mpfr_t> &r,
		FP_mat<mpfr_t> &s,
		Z_T ** ppA,
		long long int ** ppA_gram,
		int i,
		int k,
		int d,
		int beta,
		FP_NR<mpfr_t> delta,
		int &nInsertions,
		int &recursionDepth,
		int start_index = 0,
		int end_index = 0
		)
{
    long long int Gram_elem_tmp;
    Z_T * pTemp;
    long long int * pGramTemp;
    int ub = 0;
    int lb = 0;
    int i_min = k+1;
    int i_left = k+1;

    if (i>=k)
    {
        return i_min;
    }

    FP_NR<mpfr_t> minDelta;
    minDelta = delta;
    minDelta = delta + 1.0;
    int i1 = k+1;
    int k1 = k+1;
    // Check all ppDelta in range [1,beta] and [k-beta, k-1]
    for (int l=i+1; l<=k; l++) 
    {
	lb = max(beta, l-beta);
	ub = min(l, beta);
    	for (int j=i; j<ub; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    	for (int j=lb; j<l; j++) 
	{
       	    if (delta_mat[l][j] < minDelta) 
	    {
	        minDelta = delta_mat[l][j];
		i1 = j;
		k1 = l;
	    }
	}
    }
    if (i1 == k+1 && k1 == k+1) 
    {
        return i_min;
    }
    if (minDelta > delta)
    {
        return i_min;
    }

    nInsertions++;
    // Deep insert b_{k1} in position i1 
    pTemp = ppA[k1];
    for (int j=k1-1; j>=i1; j--)
    {
        ppA[j+1] = ppA[j];
    }
    ppA[i1] = pTemp;

    for (int j=0; j<i1; j++)
    {
        mu[i1][j] = mu[k1][j];
        r[i1][j] = r[k1][j];
    }
    r[i1][i1] = s[k1][i1];
    // Update Gram matrix
    for (int l=k1; l>i1; l--)
    {
        ppA_gram[l-1][l] = ppA_gram[l-1][l-1];
        ppA_gram[l-1][l-1] = ppA_gram[l][l-1];
        ppA_gram[l][l-1] = ppA_gram[l][l];
        ppA_gram[l][l] = 0;

        pGramTemp = ppA_gram[l];
        ppA_gram[l] = ppA_gram[l-1];
        ppA_gram[l-1] = pGramTemp;
	Gram_elem_tmp = ppA_gram[l][l-1];

 	for (int j=l-1; j>i1; j--)
	{
	    ppA_gram[l][j] = ppA_gram[l][j-1];
	}
	ppA_gram[l][i1] = Gram_elem_tmp;
    }

    for (int l=d-1; l>k1; l--)
    {
	Gram_elem_tmp = ppA_gram[l][k1];
        for (int j=k1; j>i1; j--)
	{
	    ppA_gram[l][j] = ppA_gram[l][j-1];
	}
	ppA_gram[l][i1] = Gram_elem_tmp;
    }

    recursionDepth++;
    if (recursionDepth >= MAX_RECURSION_DEPTH)
    {
	cout << "Max recursion depth met" << endl;
        return i_min;
    } 

    // Index search (i,i1-1), (k1+1,k)
    i_left = pgg_recursive_search_global_cfa (delta_mat, mu, r, s, ppA, ppA_gram, i, i1-1, d, beta, delta, nInsertions, recursionDepth);
    // Check here if the i_left returned by the call above is less than the minimum i obtained thus far
    // If so, we set i_min = i_left (returned by function above)
    if (i_left < i_min) {
        i_min = i_left;
    }

    pgg_recursive_search_global_cfa (delta_mat, mu, r, s, ppA, ppA_gram, k1+1, k, d, beta, delta, nInsertions, recursionDepth);

    return i_min;
}

// LC-PGG reduction for ZZ_mat<mpz_t> type basis
int lcpgg (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    
    mpz_t nint_mu_z;
    mpz_init (nint_mu_z);
    mpfr_t nint_mu_f;
    mpfr_init2 (nint_mu_f, mpfr_prec_t (precision));
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    Z_NR<mpz_t> max_nint_mu, nint_mu, MAX_nint_mu;
    nint_mu = 0;
    max_nint_mu = 0; 
    MAX_nint_mu = MAX_INT; 

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d,0.0);
    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;
    FP_mat<mpfr_t> delta_mat(d,d);

    // GSO - initial
    compute_GSO (A, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
                for (int j=l-1; j>=start_index; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper (A, mu, d, l, j, 
                        nint_mu_f, 
                        nint_mu_FP,
                        nint_mu_z,
                        max_nint_mu,
                        nint_mu,
                        start_index
                        );
                    }
                }
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
	    if (size_reduction_required)
	    {
                for (int l=i_min; l<=end_index; l++)
                {
                    for (int j=l-1; j>=start_index; j--)
                    {
                        // Takes mu[i][j]
                        // Computes nearest mpz_t to mu_[i][j]
                        if (mu[l][j] < -eta || mu[l][j] > eta)
                        {
                            size_reduce_wrapper (A, mu, d, l, j, 
                            nint_mu_f, 
                            nint_mu_FP,
                            nint_mu_z,
                            max_nint_mu,
                            nint_mu,
                            start_index
                            );
                        }
                    }
                }
	    }
	
	    // Compute initial values for C,  C[i]= <b[i],b[i]>
	    for (int l=i_min; l<d; l++)
	    {
                A[l].dot_product (C_int, A[l]);
	        C_int.abs(C_int);
	        C[l].set_z (C_int, MPFR_RNDN);
	    }

	    // Compute delta matrix
   	    for (int k=i_min; k<d; k++) 
	    {
      		for (int j=0; j<k; j++) 
		{
		    delta_mat[k][j] = (C[k] / B[j]);
         	    C[k] = C[k] - (mu[k][j] * mu[k][j] * B[j]);
         	}
   	    }

	    // Call pgg_recursive_search
	    recursionDepth = 0;
	    i_min = pgg_recursive_search (delta_mat, mu, B, A, 0, d-1, d, delta, nInsertions, recursionDepth);
	    size_reduction_required = true;
        }

        // Check the LC_GG_Reducedness
        compute_GSO (A, gso, mu, B, d, start_index, end_index);

        lcgg_reduced_check_internal(
			A,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			C_int,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );
        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    Z_NR<mpz_t> nint_mu_difference;
    nint_mu_difference.sub(MAX_nint_mu, max_nint_mu);
    cout << "MAX_nint_mu - max_nint_mu = " << nint_mu_difference << endl;

    return RED_SUCCESS;
}

// LC-GG reduction for Z_T ** type basis
template <typename Z_T>
int lcpgg (
                Z_T ** &ppA,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    mpfr_t nint_mu_f;
    mpfr_init2(nint_mu_f, mpfr_prec_t (precision));
    long int max_nint_mu = 0;
    long int nint_mu = 0;
    
    Z_T * pTemp;
    long long int inner;
    inner = 0;
    Z_NR<mpz_t> inner_ZNR;
    inner_ZNR = 0;

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d,0.0);

    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;
    FP_mat<mpfr_t> delta_mat(d,d);

    // GSO - initial
    compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
                for (int j=l-1; j>=start_index; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper (ppA, mu, d, l, j, 
                        nint_mu_f, 
                        nint_mu_FP,
                        max_nint_mu,
                        nint_mu,
                        start_index
                        );
                    }
                }
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
	    if (size_reduction_required)
	    {
                for (int l=i_min; l<=end_index; l++)
                {
                    for (int j=l-1; j>=start_index; j--)
                    {
                        // Takes mu[i][j]
                        // Computes nearest mpz_t to mu_[i][j]
                        if (mu[l][j] < -eta || mu[l][j] > eta)
                        {
                            size_reduce_wrapper (ppA, mu, d, l, j, 
                            nint_mu_f, 
                            nint_mu_FP,
                            max_nint_mu,
                            nint_mu,
                            start_index
                            );
                        }
                    }
                }
	    }
	
	    // Compute initial values for C,  C[i]= <b[i],b[i]>
	    for (int l=0; l<d; l++)
	    {
                inner = 0;
	        for (int j=0; j<d; j++)
	        {
	            inner += (ppA[l][j] * ppA[l][j]);
	        }
	        if (inner < 0)
	        {
	            inner = -inner;
	        }
	        inner_ZNR = inner;
	        C[l].set_z (inner_ZNR, MPFR_RNDN);
	    }

	    // Compute delta matrix
   	    for (int k=i_min; k<d; k++) 
	    {
      		for (int j=0; j<k; j++) 
		{
		    delta_mat[k][j] = (C[k] / B[j]);
         	    C[k] = C[k] - (mu[k][j] * mu[k][j] * B[j]);
         	}
   	    }

	    // Call pgg_recursive_search
	    recursionDepth = 0;
	    i_min = pgg_recursive_search (delta_mat, mu, B, ppA, 0, d-1, d, delta, nInsertions, recursionDepth);
	    size_reduction_required = true;
        }

        // Check the LC_GG_Reducedness
        compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

        lcgg_reduced_check_internal(
			ppA,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			inner,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );

        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    cout << "max_nint_mu = " << max_nint_mu << endl;
    cout << "MAX_INT - max_nint_mu = " << MAX_INT - max_nint_mu << endl;
    cout << "MAX_LONG_INT - max_nint_mu = " << MAX_LONG_INT - max_nint_mu << endl;

    return RED_SUCCESS;
}


// LC-GG reduction for ZZ_mat<mpz_t> type basis
// Using CFA and lazy size reduction
int lcpgg_CFA_LSR (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    mpz_t nint_mu_z;
    mpz_init (nint_mu_z);
    mpfr_t nint_mu_f;
    mpfr_init2 (nint_mu_f, mpfr_prec_t (precision));
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    FP_NR<mpfr_t> mu_tmp;
    Z_NR<mpz_t> z_tmp1;
    Z_NR<mpz_t> z_tmp2;
    Z_NR<mpz_t> max_nint_mu, nint_mu, MAX_nint_mu, Gram_elem_tmp;
    nint_mu = 0;
    max_nint_mu = 0; 
    MAX_nint_mu = MAX_INT;
    Gram_elem_tmp = 0; 
    FP_NR<mpfr_t> delta_bar;
    delta_bar = (delta + 1.0) / 2;
    FP_NR<mpfr_t> eta_bar;
    eta_bar = (eta + 0.5) / 2;

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d, 0.0);
    vector<Z_NR<mpz_t>> nint_mu_vec;
    nint_mu_vec.resize(d);
    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;
    FP_mat<mpfr_t> delta_mat(d,d);

    ZZ_mat<mpz_t> A_gram(d,d);
    FP_mat<mpfr_t> r(d,d);
    FP_mat<mpfr_t> s(d,d);

    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
	{
	    A_gram[i][j] = 0;
	    r[i][j] = 0.0;
	    s[i][j] = 0.0;  
	    mu[i][j] = 0.0;
  	}
    }
    // Compute Gram matrix initially
    compute_Gram_Matrix (A, A_gram, d);

    // Initial CFA
    // This updates r, s, mu required for future reduction
    compute_CFA (A_gram, r, mu, s, d, 0, 0); 

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
                lazy_size_reduce(A, A_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_z, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
            if (size_reduction_required)
            {
                for (int l=i_min; l<=end_index; l++)
                {
                    lazy_size_reduce(A, A_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_z, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
                }

                size_reduction_required = false;
                deep_insertion_required = true;
                index_search_required = true;
            }
	    // Compute delta matrix
   	    for (int k=0; k<d; k++) 
	    {
      		for (int j=0; j<k; j++) 
		{
		    delta_mat[k][j] = (s[k][j] / r[j][j]);
         	}
   	    }
    
	    recursionDepth = 0;
            i_min = pgg_recursive_search_cfa (delta_mat, mu, r, s, A, A_gram, 0, d-1, d, delta_bar, nInsertions, recursionDepth);

	    size_reduction_required = true;
        }

        // Check the LC_GG_Reducedness
        compute_GSO (A, gso, mu, B, d, start_index, end_index);

        lcgg_reduced_check_internal(
			A,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			C_int,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );
        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    Z_NR<mpz_t> nint_mu_difference;
    nint_mu_difference.sub(MAX_nint_mu, max_nint_mu);
    cout << "MAX_nint_mu - max_nint_mu = " << nint_mu_difference << endl;

    return RED_SUCCESS;
}

// LC-GG reduction for Z_T ** type basis
// Uses CFA and LSR
template <typename Z_T>
int lcpgg_CFA_LSR (
                Z_T ** &ppA,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )

{
    end_index = (!end_index) ? d-1 : end_index;
    FP_NR<mpfr_t> delta_bar;
    FP_NR<mpfr_t> eta_bar;
    eta_bar = (eta + 0.5) / 2;
    FP_NR<mpfr_t> mu_tmp;
    Z_NR<mpz_t> z_tmp1;
    Z_NR<mpz_t> z_tmp2;
    delta_bar = (delta + 1.0) / 2;

    mpz_t nint_mu_z;
    mpz_init (nint_mu_z);
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    mpfr_t nint_mu_f;
    mpfr_init2(nint_mu_f, mpfr_prec_t (precision));
    long int max_nint_mu = 0;
    long int nint_mu = 0;
    Z_T * nint_mu_vec;
    nint_mu_vec = (Z_T *) calloc (d, sizeof(Z_T));
    
    Z_T * pTemp;

    long long int ** ppA_gram;
    long long int * pGram_tmp;
    long long int Gram_elem_tmp;
    long long int inner;
    ppA_gram = (long long int **) calloc (d, sizeof(long long int *));
    for (int i=0; i<d; i++)
    {
        ppA_gram[i] = (long long int *) calloc (d, sizeof(long long int));
    }
    pGram_tmp = (long long int *) calloc (d, sizeof(long long int));

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d,0.0);

    FP_mat<mpfr_t> delta_mat(d,d);
    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;

    FP_mat<mpfr_t> r(d,d);
    FP_mat<mpfr_t> s(d,d);

    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
	{
            r[i][j] = 0.0;
            s[i][j] = 0.0;
            mu[i][j] = 0.0;
	}
    }

    // Compute Gram matrix initially
    compute_Gram_Matrix (ppA, ppA_gram, d);

    // Initial CFA
    // This updates r, s, mu required for future reduction
    compute_CFA (ppA_gram, r, mu, s, d, 0, 1); 


    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
	        lazy_size_reduce(ppA, ppA_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
	    if (size_reduction_required)
	    {
                for (int l=i_min; l<=end_index; l++)
                {
	            lazy_size_reduce(ppA, ppA_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
                }
		size_reduction_required = false;
		deep_insertion_required = true;
		index_search_required = true;
	    }

	    // Compute delta matrix
	    for (int k=i_min; k<d; k++)
	    {
	        for (int j=0; j<k; j++)
		{
		    delta_mat[k][j] = (s[k][j] / r[j][j]);
		}
	    }
            recursionDepth = 0;
	    i_min = pgg_recursive_search_cfa (delta_mat, mu, r, s, ppA, ppA_gram, 0, d-1, d, delta_bar, nInsertions, recursionDepth);
	    size_reduction_required = true;

        }

        // Check the LC_GG_Reducedness
        compute_GSO (ppA, gso, mu, B, d, start_index, end_index);
	cout << "nInsertions = " << nInsertions << endl;

        lcgg_reduced_check_internal(
			ppA,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			inner,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );
        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    cout << "max_nint_mu = " << max_nint_mu << endl;
    cout << "MAX_INT - max_nint_mu = " << MAX_INT - max_nint_mu << endl;
    cout << "MAX_LONG_INT - max_nint_mu = " << MAX_LONG_INT - max_nint_mu << endl;

    return RED_SUCCESS;
}

// LC-PGG reduction for ZZ_mat<mpz_t> type basis
// Locally restricted insertions
int lcpgg_locally_restricted (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    
    mpz_t nint_mu_z;
    mpz_init (nint_mu_z);
    mpfr_t nint_mu_f;
    mpfr_init2 (nint_mu_f, mpfr_prec_t (precision));
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    Z_NR<mpz_t> max_nint_mu, nint_mu, MAX_nint_mu;
    nint_mu = 0;
    max_nint_mu = 0; 
    MAX_nint_mu = MAX_INT; 

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d,0.0);
    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;
    FP_mat<mpfr_t> delta_mat(d,d);

    // GSO - initial
    compute_GSO (A, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
                for (int j=l-1; j>=start_index; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper (A, mu, d, l, j, 
                        nint_mu_f, 
                        nint_mu_FP,
                        nint_mu_z,
                        max_nint_mu,
                        nint_mu,
                        start_index
                        );
                    }
                }
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
	    if (size_reduction_required)
	    {
                for (int l=i_min; l<=end_index; l++)
                {
                    for (int j=l-1; j>=start_index; j--)
                    {
                        // Takes mu[i][j]
                        // Computes nearest mpz_t to mu_[i][j]
                        if (mu[l][j] < -eta || mu[l][j] > eta)
                        {
                            size_reduce_wrapper (A, mu, d, l, j, 
                            nint_mu_f, 
                            nint_mu_FP,
                            nint_mu_z,
                            max_nint_mu,
                            nint_mu,
                            start_index
                            );
                        }
                    }
                }
	    }
	
	    // Compute initial values for C,  C[i]= <b[i],b[i]>
	    for (int l=i_min; l<d; l++)
	    {
                A[l].dot_product (C_int, A[l]);
	        C_int.abs(C_int);
	        C[l].set_z (C_int, MPFR_RNDN);
	    }

	    // Compute delta matrix
   	    for (int k=i_min; k<d; k++) 
	    {
      		for (int j=0; j<k; j++) 
		{
		    delta_mat[k][j] = (C[k] / B[j]);
         	    C[k] = C[k] - (mu[k][j] * mu[k][j] * B[j]);
         	}
   	    }

	    // Call pgg_recursive_search
	    recursionDepth = 0;
	    i_min = pgg_recursive_search_local (delta_mat, mu, B, A, 0, d-1, d, beta, delta, nInsertions, recursionDepth);
	    size_reduction_required = true;
        }

        // Check the LC_GG_Reducedness
        compute_GSO (A, gso, mu, B, d, start_index, end_index);

        lcgg_reduced_check_internal(
			A,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			C_int,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );
        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    Z_NR<mpz_t> nint_mu_difference;
    nint_mu_difference.sub(MAX_nint_mu, max_nint_mu);
    cout << "MAX_nint_mu - max_nint_mu = " << nint_mu_difference << endl;

    return RED_SUCCESS;
}

// LC-GG reduction for Z_T ** type basis
template <typename Z_T>
int lcpgg_locally_restricted (
                Z_T ** &ppA,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    mpfr_t nint_mu_f;
    mpfr_init2(nint_mu_f, mpfr_prec_t (precision));
    long int max_nint_mu = 0;
    long int nint_mu = 0;
    
    Z_T * pTemp;
    long long int inner;
    inner = 0;
    Z_NR<mpz_t> inner_ZNR;
    inner_ZNR = 0;

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d,0.0);

    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;
    FP_mat<mpfr_t> delta_mat(d,d);

    // GSO - initial
    compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
                for (int j=l-1; j>=start_index; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper (ppA, mu, d, l, j, 
                        nint_mu_f, 
                        nint_mu_FP,
                        max_nint_mu,
                        nint_mu,
                        start_index
                        );
                    }
                }
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
	    if (size_reduction_required)
	    {
                for (int l=i_min; l<=end_index; l++)
                {
                    for (int j=l-1; j>=start_index; j--)
                    {
                        // Takes mu[i][j]
                        // Computes nearest mpz_t to mu_[i][j]
                        if (mu[l][j] < -eta || mu[l][j] > eta)
                        {
                            size_reduce_wrapper (ppA, mu, d, l, j, 
                            nint_mu_f, 
                            nint_mu_FP,
                            max_nint_mu,
                            nint_mu,
                            start_index
                            );
                        }
                    }
                }
	    }
	
	    // Compute initial values for C,  C[i]= <b[i],b[i]>
	    for (int l=0; l<d; l++)
	    {
                inner = 0;
	        for (int j=0; j<d; j++)
	        {
	            inner += (ppA[l][j] * ppA[l][j]);
	        }
	        if (inner < 0)
	        {
	            inner = -inner;
	        }
	        inner_ZNR = inner;
	        C[l].set_z (inner_ZNR, MPFR_RNDN);
	    }

	    // Compute delta matrix
   	    for (int k=i_min; k<d; k++) 
	    {
      		for (int j=0; j<k; j++) 
		{
		    delta_mat[k][j] = (C[k] / B[j]);
         	    C[k] = C[k] - (mu[k][j] * mu[k][j] * B[j]);
         	}
   	    }

	    // Call pgg_recursive_search
	    recursionDepth = 0;
	    i_min = pgg_recursive_search_local (delta_mat, mu, B, ppA, 0, d-1, d, beta, delta, nInsertions, recursionDepth);
	    size_reduction_required = true;
        }

        // Check the LC_GG_Reducedness
        compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

        lcgg_reduced_check_internal(
			ppA,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			inner,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );

        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    cout << "max_nint_mu = " << max_nint_mu << endl;
    cout << "MAX_INT - max_nint_mu = " << MAX_INT - max_nint_mu << endl;
    cout << "MAX_LONG_INT - max_nint_mu = " << MAX_LONG_INT - max_nint_mu << endl;

    return RED_SUCCESS;
}


// LC-GG reduction for ZZ_mat<mpz_t> type basis
// Using CFA and lazy size reduction
int lcpgg_locally_restricted_CFA_LSR (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    mpz_t nint_mu_z;
    mpz_init (nint_mu_z);
    mpfr_t nint_mu_f;
    mpfr_init2 (nint_mu_f, mpfr_prec_t (precision));
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    FP_NR<mpfr_t> mu_tmp;
    Z_NR<mpz_t> z_tmp1;
    Z_NR<mpz_t> z_tmp2;
    Z_NR<mpz_t> max_nint_mu, nint_mu, MAX_nint_mu, Gram_elem_tmp;
    nint_mu = 0;
    max_nint_mu = 0; 
    MAX_nint_mu = MAX_INT;
    Gram_elem_tmp = 0; 
    FP_NR<mpfr_t> delta_bar;
    delta_bar = (delta + 1.0) / 2;
    FP_NR<mpfr_t> eta_bar;
    eta_bar = (eta + 0.5) / 2;

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d, 0.0);
    vector<Z_NR<mpz_t>> nint_mu_vec;
    nint_mu_vec.resize(d);
    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;
    FP_mat<mpfr_t> delta_mat(d,d);

    ZZ_mat<mpz_t> A_gram(d,d);
    FP_mat<mpfr_t> r(d,d);
    FP_mat<mpfr_t> s(d,d);

    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
	{
	    A_gram[i][j] = 0;
	    r[i][j] = 0.0;
	    s[i][j] = 0.0;  
	    mu[i][j] = 0.0;
  	}
    }
    // Compute Gram matrix initially
    compute_Gram_Matrix (A, A_gram, d);

    // Initial CFA
    // This updates r, s, mu required for future reduction
    compute_CFA (A_gram, r, mu, s, d, 0, 0); 

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
                lazy_size_reduce(A, A_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_z, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
            if (size_reduction_required)
            {
                for (int l=i_min; l<=end_index; l++)
                {
                    lazy_size_reduce(A, A_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_z, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
                }

                size_reduction_required = false;
                deep_insertion_required = true;
                index_search_required = true;
            }
	    // Compute delta matrix
   	    for (int k=0; k<d; k++) 
	    {
      		for (int j=0; j<k; j++) 
		{
		    delta_mat[k][j] = (s[k][j] / r[j][j]);
         	}
   	    }
    
	    recursionDepth = 0;
            i_min = pgg_recursive_search_local_cfa (delta_mat, mu, r, s, A, A_gram, 0, d-1, d, beta, delta_bar, nInsertions, recursionDepth);

	    size_reduction_required = true;
        }

        // Check the LC_GG_Reducedness
        compute_GSO (A, gso, mu, B, d, start_index, end_index);

        lcgg_reduced_check_internal(
			A,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			C_int,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );
        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    Z_NR<mpz_t> nint_mu_difference;
    nint_mu_difference.sub(MAX_nint_mu, max_nint_mu);
    cout << "MAX_nint_mu - max_nint_mu = " << nint_mu_difference << endl;

    return RED_SUCCESS;
}

// LC-GG reduction for Z_T ** type basis
// Uses CFA and LSR
template <typename Z_T>
int lcpgg_locally_restricted_CFA_LSR (
                Z_T ** &ppA,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )

{
    end_index = (!end_index) ? d-1 : end_index;
    FP_NR<mpfr_t> delta_bar;
    delta_bar = (delta + 1.0) / 2;
    FP_NR<mpfr_t> eta_bar;
    eta_bar = (eta + 0.5) / 2;

    mpz_t nint_mu_z;
    mpz_init (nint_mu_z);
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    FP_NR<mpfr_t> mu_tmp;
    Z_NR<mpz_t> z_tmp1;
    Z_NR<mpz_t> z_tmp2;
    mpfr_t nint_mu_f;
    mpfr_init2(nint_mu_f, mpfr_prec_t (precision));
    long int max_nint_mu = 0;
    long int nint_mu = 0;
    Z_T * nint_mu_vec;
    nint_mu_vec = (Z_T *) calloc (d, sizeof(Z_T));
    
    Z_T * pTemp;

    long long int ** ppA_gram;
    long long int * pGram_tmp;
    long long int Gram_elem_tmp;
    long long int inner;
    ppA_gram = (long long int **) calloc (d, sizeof(long long int *));
    for (int i=0; i<d; i++)
    {
        ppA_gram[i] = (long long int *) calloc (d, sizeof(long long int));
    }
    pGram_tmp = (long long int *) calloc (d, sizeof(long long int));

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d,0.0);

    FP_mat<mpfr_t> delta_mat(d,d);
    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;

    FP_mat<mpfr_t> r(d,d);
    FP_mat<mpfr_t> s(d,d);

    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
	{
            r[i][j] = 0.0;
            s[i][j] = 0.0;
            mu[i][j] = 0.0;
	}
    }

    // Compute Gram matrix initially
    compute_Gram_Matrix (ppA, ppA_gram, d);

    // Initial CFA
    // This updates r, s, mu required for future reduction
    compute_CFA (ppA_gram, r, mu, s, d, 0, 1); 


    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
	        lazy_size_reduce(ppA, ppA_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
	    if (size_reduction_required)
	    {
                for (int l=i_min; l<=end_index; l++)
                {
	            lazy_size_reduce(ppA, ppA_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
                }
		size_reduction_required = false;
		deep_insertion_required = true;
		index_search_required = true;
	    }

	    // Compute delta matrix
	    for (int k=i_min; k<d; k++)
	    {
	        for (int j=0; j<k; j++)
		{
		    delta_mat[k][j] = (s[k][j] / r[j][j]);
		}
	    }
            recursionDepth = 0;
            i_min = pgg_recursive_search_local_cfa (delta_mat, mu, r, s, ppA, ppA_gram, 0, d-1, d, beta, delta_bar, nInsertions, recursionDepth);
	    size_reduction_required = true;

        }

        // Check the LC_GG_Reducedness
        compute_GSO (ppA, gso, mu, B, d, start_index, end_index);
	cout << "nInsertions = " << nInsertions << endl;

        lcgg_reduced_check_internal(
			ppA,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			inner,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );
        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    cout << "max_nint_mu = " << max_nint_mu << endl;
    cout << "MAX_INT - max_nint_mu = " << MAX_INT - max_nint_mu << endl;
    cout << "MAX_LONG_INT - max_nint_mu = " << MAX_LONG_INT - max_nint_mu << endl;

    return RED_SUCCESS;
}

// LC-PGG reduction for ZZ_mat<mpz_t> type basis
// Globally restricted insertions
int lcpgg_globally_restricted (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    
    mpz_t nint_mu_z;
    mpz_init (nint_mu_z);
    mpfr_t nint_mu_f;
    mpfr_init2 (nint_mu_f, mpfr_prec_t (precision));
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    Z_NR<mpz_t> max_nint_mu, nint_mu, MAX_nint_mu;
    nint_mu = 0;
    max_nint_mu = 0; 
    MAX_nint_mu = MAX_INT; 

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d,0.0);
    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;
    FP_mat<mpfr_t> delta_mat(d,d);

    // GSO - initial
    compute_GSO (A, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
                for (int j=l-1; j>=start_index; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper (A, mu, d, l, j, 
                        nint_mu_f, 
                        nint_mu_FP,
                        nint_mu_z,
                        max_nint_mu,
                        nint_mu,
                        start_index
                        );
                    }
                }
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
	    if (size_reduction_required)
	    {
                for (int l=i_min; l<=end_index; l++)
                {
                    for (int j=l-1; j>=start_index; j--)
                    {
                        // Takes mu[i][j]
                        // Computes nearest mpz_t to mu_[i][j]
                        if (mu[l][j] < -eta || mu[l][j] > eta)
                        {
                            size_reduce_wrapper (A, mu, d, l, j, 
                            nint_mu_f, 
                            nint_mu_FP,
                            nint_mu_z,
                            max_nint_mu,
                            nint_mu,
                            start_index
                            );
                        }
                    }
                }
	    }
	
	    // Compute initial values for C,  C[i]= <b[i],b[i]>
	    for (int l=i_min; l<d; l++)
	    {
                A[l].dot_product (C_int, A[l]);
	        C_int.abs(C_int);
	        C[l].set_z (C_int, MPFR_RNDN);
	    }

	    // Compute delta matrix
   	    for (int k=i_min; k<d; k++) 
	    {
      		for (int j=0; j<k; j++) 
		{
		    delta_mat[k][j] = (C[k] / B[j]);
         	    C[k] = C[k] - (mu[k][j] * mu[k][j] * B[j]);
         	}
   	    }

	    // Call pgg_recursive_search
	    recursionDepth = 0;
	    i_min = pgg_recursive_search_global (delta_mat, mu, B, A, 0, d-1, d, beta, delta, nInsertions, recursionDepth);
	    size_reduction_required = true;
        }

        // Check the LC_GG_Reducedness
        compute_GSO (A, gso, mu, B, d, start_index, end_index);

        lcgg_reduced_check_internal(
			A,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			C_int,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );
        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    Z_NR<mpz_t> nint_mu_difference;
    nint_mu_difference.sub(MAX_nint_mu, max_nint_mu);
    cout << "MAX_nint_mu - max_nint_mu = " << nint_mu_difference << endl;

    return RED_SUCCESS;
}

// LC-GG reduction for Z_T ** type basis
template <typename Z_T>
int lcpgg_globally_restricted (
                Z_T ** &ppA,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    mpfr_t nint_mu_f;
    mpfr_init2(nint_mu_f, mpfr_prec_t (precision));
    long int max_nint_mu = 0;
    long int nint_mu = 0;
    
    Z_T * pTemp;
    long long int inner;
    inner = 0;
    Z_NR<mpz_t> inner_ZNR;
    inner_ZNR = 0;

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d,0.0);

    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;
    FP_mat<mpfr_t> delta_mat(d,d);

    // GSO - initial
    compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
                for (int j=l-1; j>=start_index; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper (ppA, mu, d, l, j, 
                        nint_mu_f, 
                        nint_mu_FP,
                        max_nint_mu,
                        nint_mu,
                        start_index
                        );
                    }
                }
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
	    if (size_reduction_required)
	    {
                for (int l=i_min; l<=end_index; l++)
                {
                    for (int j=l-1; j>=start_index; j--)
                    {
                        // Takes mu[i][j]
                        // Computes nearest mpz_t to mu_[i][j]
                        if (mu[l][j] < -eta || mu[l][j] > eta)
                        {
                            size_reduce_wrapper (ppA, mu, d, l, j, 
                            nint_mu_f, 
                            nint_mu_FP,
                            max_nint_mu,
                            nint_mu,
                            start_index
                            );
                        }
                    }
                }
	    }
	
	    // Compute initial values for C,  C[i]= <b[i],b[i]>
	    for (int l=0; l<d; l++)
	    {
                inner = 0;
	        for (int j=0; j<d; j++)
	        {
	            inner += (ppA[l][j] * ppA[l][j]);
	        }
	        if (inner < 0)
	        {
	            inner = -inner;
	        }
	        inner_ZNR = inner;
	        C[l].set_z (inner_ZNR, MPFR_RNDN);
	    }

	    // Compute delta matrix
   	    for (int k=1; k<d; k++) 
	    {
      		for (int j=0; j<k; j++) 
		{
		    delta_mat[k][j] = (C[k] / B[j]);
         	    C[k] = C[k] - (mu[k][j] * mu[k][j] * B[j]);
         	}
   	    }

	    // Call pgg_recursive_search
	    recursionDepth = 0;
	    i_min = pgg_recursive_search_global (delta_mat, mu, B, ppA, 0, d-1, d, beta, delta, nInsertions, recursionDepth);
	    size_reduction_required = true;
        }

        // Check the LC_GG_Reducedness
        compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

        lcgg_reduced_check_internal(
			ppA,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			inner,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );

        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    cout << "max_nint_mu = " << max_nint_mu << endl;
    cout << "MAX_INT - max_nint_mu = " << MAX_INT - max_nint_mu << endl;
    cout << "MAX_LONG_INT - max_nint_mu = " << MAX_LONG_INT - max_nint_mu << endl;

    return RED_SUCCESS;
}


// LC-GG reduction for ZZ_mat<mpz_t> type basis
// Using CFA and lazy size reduction
int lcpgg_globally_restricted_CFA_LSR (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )
{
    end_index = (!end_index) ? d-1 : end_index;

    mpz_t nint_mu_z;
    mpz_init (nint_mu_z);
    mpfr_t nint_mu_f;
    mpfr_init2 (nint_mu_f, mpfr_prec_t (precision));
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    FP_NR<mpfr_t> mu_tmp;
    Z_NR<mpz_t> z_tmp1;
    Z_NR<mpz_t> z_tmp2;
    Z_NR<mpz_t> max_nint_mu, nint_mu, MAX_nint_mu, Gram_elem_tmp;
    nint_mu = 0;
    max_nint_mu = 0; 
    MAX_nint_mu = MAX_INT;
    Gram_elem_tmp = 0; 
    FP_NR<mpfr_t> delta_bar;
    delta_bar = (delta + 1.0) / 2;
    FP_NR<mpfr_t> eta_bar;
    eta_bar = (eta + 0.5) / 2;

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d, 0.0);
    vector<Z_NR<mpz_t>> nint_mu_vec;
    nint_mu_vec.resize(d);
    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;
    FP_mat<mpfr_t> delta_mat(d,d);

    ZZ_mat<mpz_t> A_gram(d,d);
    FP_mat<mpfr_t> r(d,d);
    FP_mat<mpfr_t> s(d,d);

    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
	{
	    A_gram[i][j] = 0;
	    r[i][j] = 0.0;
	    s[i][j] = 0.0;  
	    mu[i][j] = 0.0;
  	}
    }
    // Compute Gram matrix initially
    compute_Gram_Matrix (A, A_gram, d);

    // Initial CFA
    // This updates r, s, mu required for future reduction
    compute_CFA (A_gram, r, mu, s, d, 0, 0); 

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
                lazy_size_reduce(A, A_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_z, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
            if (size_reduction_required)
            {
                for (int l=i_min; l<=end_index; l++)
                {
                    lazy_size_reduce(A, A_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_z, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
                }

                size_reduction_required = false;
                deep_insertion_required = true;
                index_search_required = true;
            }
	    // Compute delta matrix
   	    for (int k=0; k<d; k++) 
	    {
      		for (int j=0; j<k; j++) 
		{
		    delta_mat[k][j] = (s[k][j] / r[j][j]);
         	}
   	    }
    
	    recursionDepth = 0;
            i_min = pgg_recursive_search_global_cfa (delta_mat, mu, r, s, A, A_gram, 0, d-1, d, beta, delta_bar, nInsertions, recursionDepth);

	    size_reduction_required = true;
        }

        // Check the LC_GG_Reducedness
        compute_GSO (A, gso, mu, B, d, start_index, end_index);

        lcgg_reduced_check_internal(
			A,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			C_int,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );
        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    Z_NR<mpz_t> nint_mu_difference;
    nint_mu_difference.sub(MAX_nint_mu, max_nint_mu);
    cout << "MAX_nint_mu - max_nint_mu = " << nint_mu_difference << endl;

    return RED_SUCCESS;
}

// LC-GG reduction for Z_T ** type basis
// Uses CFA and LSR
template <typename Z_T>
int lcpgg_globally_restricted_CFA_LSR (
                Z_T ** &ppA,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
		int beta,
                FP_mat<mpfr_t> &gso,
                FP_mat<mpfr_t> &mu,
                vector<FP_NR<mpfr_t>> &B,
                int start_index = 0,
                int end_index = 0
                )

{
    end_index = (!end_index) ? d-1 : end_index;
    FP_NR<mpfr_t> delta_bar;
    delta_bar = (delta + 1.0) / 2;
    FP_NR<mpfr_t> eta_bar;
    eta_bar = (eta + 0.5) / 2;

    mpz_t nint_mu_z;
    mpz_init (nint_mu_z);
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    FP_NR<mpfr_t> mu_tmp;
    Z_NR<mpz_t> z_tmp1;
    Z_NR<mpz_t> z_tmp2;
    mpfr_t nint_mu_f;
    mpfr_init2(nint_mu_f, mpfr_prec_t (precision));
    long int max_nint_mu = 0;
    long int nint_mu = 0;
    Z_T * nint_mu_vec;
    nint_mu_vec = (Z_T *) calloc (d, sizeof(Z_T));
    
    Z_T * pTemp;

    long long int ** ppA_gram;
    long long int * pGram_tmp;
    long long int Gram_elem_tmp;
    long long int inner;
    ppA_gram = (long long int **) calloc (d, sizeof(long long int *));
    for (int i=0; i<d; i++)
    {
        ppA_gram[i] = (long long int *) calloc (d, sizeof(long long int));
    }
    pGram_tmp = (long long int *) calloc (d, sizeof(long long int));

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d,0.0);

    FP_mat<mpfr_t> delta_mat(d,d);
    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;

    FP_mat<mpfr_t> r(d,d);
    FP_mat<mpfr_t> s(d,d);

    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
	{
            r[i][j] = 0.0;
            s[i][j] = 0.0;
            mu[i][j] = 0.0;
	}
    }

    // Compute Gram matrix initially
    compute_Gram_Matrix (ppA, ppA_gram, d);

    // Initial CFA
    // This updates r, s, mu required for future reduction
    compute_CFA (ppA_gram, r, mu, s, d, 0, 1); 


    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int i_min = start_index;
    int nInsertions = 0;
    int recursionDepth = 0;
    int nLoops = 1;

    long int totalInsertionDepth = 0;

    while (size_reduction_required || deep_insertion_required)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        if (size_reduction_required)
        {
            // Size reduction - initial
            for (int l=k; l<=end_index; l++)
            {
	        lazy_size_reduce(ppA, ppA_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        // Loop
        while (i_min < d) 
        {            
	    if (size_reduction_required)
	    {
                for (int l=0; l<=end_index; l++)
                {
	            lazy_size_reduce(ppA, ppA_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
                }
		size_reduction_required = false;
		deep_insertion_required = true;
		index_search_required = true;
	    }

	    // Compute delta matrix
	    for (int k=i_min; k<d; k++)
	    {
	        for (int j=0; j<k; j++)
		{
		    delta_mat[k][j] = (s[k][j] / r[j][j]);
		}
	    }
            recursionDepth = 0;
            i_min = pgg_recursive_search_global_cfa (delta_mat, mu, r, s, ppA, ppA_gram, 0, d-1, d, beta, delta_bar, nInsertions, recursionDepth);
	    size_reduction_required = true;

        }

        // Check the LC_GG_Reducedness
        compute_GSO (ppA, gso, mu, B, d, start_index, end_index);
	cout << "nInsertions = " << nInsertions << endl;

        lcgg_reduced_check_internal(
			ppA,
                        mu,
                        B,
                        delta,
		        delta_kj,
		        delta_min,
                        precision,
                        d,
                        eta,
			beta,
			inner,
			C,
                        k,
                        i,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );
        nLoops++;
    }
    long double avgInsertionDepth = (long double) totalInsertionDepth /(long double) nInsertions; 
    cout << "Average Insertion Depth = " << avgInsertionDepth << endl;

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    cout << "max_nint_mu = " << max_nint_mu << endl;
    cout << "MAX_INT - max_nint_mu = " << MAX_INT - max_nint_mu << endl;
    cout << "MAX_LONG_INT - max_nint_mu = " << MAX_LONG_INT - max_nint_mu << endl;

    return RED_SUCCESS;
}


// Wrapper for LC-PGG reduction
int lcpgg_wrapper (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta_LC,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
                long double &time_taken,
                FP_NR<mpfr_t> &RHF,
                FP_NR<mpfr_t> &last_GSO_length,
                int start_index = 0,
                int end_index = 0,
                int sublattice_size = DEFAULT_SUBLATTICE_SIZE,
		int beta = DEFAULT_BETA
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    sublattice_size = (sublattice_size <= 0) ? d : sublattice_size;

    cout << "Precision = " << precision << endl;
    cout << "Using standard GSO" << endl;
    FP_NR<mpfr_t>::set_prec(precision);
    if (beta == 0) 
    {
        cout << "No insertion restriction." << endl;
	beta = d;
    }
    else
    {
        cout << "This option does not allow for insertion restrictions." << endl;
        cout << "If restricted insertions required, choose LC-PGG with local or global restrictions." << endl;
        cout << "Running with no insertion restriction." << endl;
    }
    
    FP_mat<mpfr_t> gso(d,d);
    FP_mat<mpfr_t> mu(d,d);
    vector<FP_NR<mpfr_t>> B;
    B.resize(d,0.0);

    int prec_status; 

    int ** ppIntA = NULL;
    long int ** ppLongIntA = NULL;
    long long int ** ppLongLongIntA = NULL;

    clock_t begin_time;
    clock_t end_time;

    int int_type=0;
    int int_type_new=0;

    // *************************************************************************
    begin_time = clock();

    int_type = getIntegerType(A, d);
    cout << "Int type: " << int_type << endl;
    switch (int_type)
    {
        case TYPE_MPZ_T:
            break;
        case TYPE_LONG_LONG_INT:
            copy_basis_mpz_Z_T (ppLongLongIntA, A, d); 
            break;
        case TYPE_LONG_INT:
            copy_basis_mpz_Z_T (ppLongIntA, A, d); 
            break;
        case TYPE_INT:
            copy_basis_mpz_Z_T (ppIntA, A, d); 
            break;
        default:
            cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
            exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
    };

    int sublattice_start, sublattice_end;
    sublattice_start = start_index;
    sublattice_end = 0;
    while (sublattice_start < end_index)
    {
        sublattice_end += sublattice_size;
        sublattice_end = (sublattice_end > end_index) ? end_index : sublattice_end;
        cout << endl;

        switch (int_type)
        {
            case TYPE_MPZ_T:
		cout << "in switch case mpz_t..." << endl;
                prec_status = lcpgg (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
		cout << "in switch case long long int..." << endl;
                prec_status = lcpgg (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
		cout << "in switch case long int..." << endl;
                prec_status = lcpgg (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
		cout << "in switch case int..." << endl;
                prec_status = lcpgg (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
        sublattice_start += sublattice_size;
    }

    int_type_new = getIntegerType(A, d);

    if (sublattice_size != d)
    {
        cout << "Preprocessing complete. Now running LCGG globally." << endl;
        switch (int_type)
        {
            case TYPE_MPZ_T:
                prec_status = lcpgg (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = lcpgg (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = lcpgg (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = lcpgg (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
    }
    end_time = clock();
    // *************************************************************************

    time_taken = (long double) (end_time - begin_time) / CLOCKS_PER_SEC;
    cout << "Time=" << time_taken;
    RHF = compute_RHF (B, d, start_index, end_index);
    cout << "\tRHF=" << RHF << endl;
    last_GSO_length = B[d-1];
    cout << "Final GSO Length = " << B[d-1] << endl << endl;
            
    return prec_status;
}

// Wrapper for LC-GG reduction (CFA and LSR)
int lcpgg_cfa_wrapper (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta_LC,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
                long double &time_taken,
                FP_NR<mpfr_t> &RHF,
                FP_NR<mpfr_t> &last_GSO_length,
                int start_index = 0,
                int end_index = 0,
                int sublattice_size = DEFAULT_SUBLATTICE_SIZE,
		int beta = DEFAULT_BETA
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    sublattice_size = (sublattice_size <= 0) ? d : sublattice_size;

    cout << "Precision = " << precision << endl;
    cout << "Using CFA and Lazy Size Reduction" << endl;
    FP_NR<mpfr_t>::set_prec(precision);
    if (beta == 0) 
    {
        cout << "No beta on insertions." << endl;
	beta = d;
    }
    else
    {
        cout << "Beta = " << beta << endl;
    }
    
    FP_mat<mpfr_t> gso(d,d);
    FP_mat<mpfr_t> mu(d,d);
    vector<FP_NR<mpfr_t>> B;
    B.resize(d,0.0);

    int prec_status; 

    int ** ppIntA = NULL;
    long int ** ppLongIntA = NULL;
    long long int ** ppLongLongIntA = NULL;

    clock_t begin_time;
    clock_t end_time;

    int int_type=0;
    int int_type_new=0;

    // *************************************************************************
    begin_time = clock();

    int_type = getIntegerType(A, d);
    cout << "Int type: " << int_type << endl;
    switch (int_type)
    {
        case TYPE_MPZ_T:
            break;
        case TYPE_LONG_LONG_INT:
            copy_basis_mpz_Z_T (ppLongLongIntA, A, d); 
            break;
        case TYPE_LONG_INT:
            copy_basis_mpz_Z_T (ppLongIntA, A, d); 
            break;
        case TYPE_INT:
            copy_basis_mpz_Z_T (ppIntA, A, d); 
            break;
        default:
            cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
            exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
    };

    int sublattice_start, sublattice_end;
    sublattice_start = start_index;
    sublattice_end = 0;
    while (sublattice_start < end_index)
    {
        sublattice_end += sublattice_size;
        sublattice_end = (sublattice_end > end_index) ? end_index : sublattice_end;
        cout << endl;

        switch (int_type)
        {
            case TYPE_MPZ_T:
		cout << "in switch case mpz_t..." << endl;
                prec_status = lcpgg_CFA_LSR (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
		cout << "in switch case long long int..." << endl;
                prec_status = lcpgg_CFA_LSR (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
		cout << "in switch case long int..." << endl;
                prec_status = lcpgg_CFA_LSR (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
		cout << "in switch case int..." << endl;
                prec_status = lcpgg_CFA_LSR (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
        sublattice_start += sublattice_size;
    }

    int_type_new = getIntegerType(A, d);

    if (sublattice_size != d)
    {
        cout << "Preprocessing complete. Now running LCGG globally." << endl;
        switch (int_type)
        {
            case TYPE_MPZ_T:
                prec_status = lcpgg_CFA_LSR (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = lcpgg_CFA_LSR (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = lcpgg_CFA_LSR (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = lcpgg_CFA_LSR (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
    }
    end_time = clock();
    // *************************************************************************

    time_taken = (long double) (end_time - begin_time) / CLOCKS_PER_SEC;
    cout << "Time=" << time_taken;
    RHF = compute_RHF (B, d, start_index, end_index);
    cout << "\tRHF=" << RHF << endl;
    last_GSO_length = B[d-1];
    cout << "Final GSO Length = " << B[d-1] << endl << endl;
            
    return prec_status;
}

// Wrapper for LC-PGG reduction
// Locally restricted insertions
int lcpgg_local_wrapper (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta_LC,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
                long double &time_taken,
                FP_NR<mpfr_t> &RHF,
                FP_NR<mpfr_t> &last_GSO_length,
                int start_index = 0,
                int end_index = 0,
                int sublattice_size = DEFAULT_SUBLATTICE_SIZE,
		int beta = DEFAULT_BETA
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    sublattice_size = (sublattice_size <= 0) ? d : sublattice_size;

    cout << "Precision = " << precision << endl;
    cout << "Using standard GSO" << endl;
    FP_NR<mpfr_t>::set_prec(precision);
    if (beta == 0) 
    {
        cout << "No beta on insertions." << endl;
        cout << "Running standard LC-PGG" << endl;
	beta = d;
    }
    else
    {
        cout << "Beta = " << beta << endl;
    }
    
    FP_mat<mpfr_t> gso(d,d);
    FP_mat<mpfr_t> mu(d,d);
    vector<FP_NR<mpfr_t>> B;
    B.resize(d,0.0);

    int prec_status; 

    int ** ppIntA = NULL;
    long int ** ppLongIntA = NULL;
    long long int ** ppLongLongIntA = NULL;

    clock_t begin_time;
    clock_t end_time;

    int int_type=0;
    int int_type_new=0;

    // *************************************************************************
    begin_time = clock();

    int_type = getIntegerType(A, d);
    cout << "Int type: " << int_type << endl;
    switch (int_type)
    {
        case TYPE_MPZ_T:
            break;
        case TYPE_LONG_LONG_INT:
            copy_basis_mpz_Z_T (ppLongLongIntA, A, d); 
            break;
        case TYPE_LONG_INT:
            copy_basis_mpz_Z_T (ppLongIntA, A, d); 
            break;
        case TYPE_INT:
            copy_basis_mpz_Z_T (ppIntA, A, d); 
            break;
        default:
            cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
            exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
    };

    int sublattice_start, sublattice_end;
    sublattice_start = start_index;
    sublattice_end = 0;
    while (sublattice_start < end_index)
    {
        sublattice_end += sublattice_size;
        sublattice_end = (sublattice_end > end_index) ? end_index : sublattice_end;
        cout << endl;

        switch (int_type)
        {
            case TYPE_MPZ_T:
		cout << "in switch case mpz_t..." << endl;
                prec_status = lcpgg_locally_restricted (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
		cout << "in switch case long long int..." << endl;
                prec_status = lcpgg_locally_restricted (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
		cout << "in switch case long int..." << endl;
                prec_status = lcpgg_locally_restricted (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
		cout << "in switch case int..." << endl;
                prec_status = lcpgg_locally_restricted (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
        sublattice_start += sublattice_size;
    }

    int_type_new = getIntegerType(A, d);

    if (sublattice_size != d)
    {
        cout << "Preprocessing complete. Now running LCGG globally." << endl;
        switch (int_type)
        {
            case TYPE_MPZ_T:
                prec_status = lcpgg_locally_restricted (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = lcpgg_locally_restricted (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = lcpgg_locally_restricted (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = lcpgg_locally_restricted (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
    }
    end_time = clock();
    // *************************************************************************

    time_taken = (long double) (end_time - begin_time) / CLOCKS_PER_SEC;
    cout << "Time=" << time_taken;
    RHF = compute_RHF (B, d, start_index, end_index);
    cout << "\tRHF=" << RHF << endl;
    last_GSO_length = B[d-1];
    cout << "Final GSO Length = " << B[d-1] << endl << endl;
            
    return prec_status;
}

// Wrapper for LC-GG reduction (CFA and LSR)
int lcpgg_local_cfa_wrapper (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta_LC,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
                long double &time_taken,
                FP_NR<mpfr_t> &RHF,
                FP_NR<mpfr_t> &last_GSO_length,
                int start_index = 0,
                int end_index = 0,
                int sublattice_size = DEFAULT_SUBLATTICE_SIZE,
		int beta = DEFAULT_BETA
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    sublattice_size = (sublattice_size <= 0) ? d : sublattice_size;

    cout << "Precision = " << precision << endl;
    cout << "Using CFA and Lazy Size Reduction" << endl;
    FP_NR<mpfr_t>::set_prec(precision);
    if (beta == 0) 
    {
        cout << "No beta on insertions." << endl;
	beta = d;
    }
    else
    {
        cout << "Beta = " << beta << endl;
    }
    
    FP_mat<mpfr_t> gso(d,d);
    FP_mat<mpfr_t> mu(d,d);
    vector<FP_NR<mpfr_t>> B;
    B.resize(d,0.0);

    int prec_status; 

    int ** ppIntA = NULL;
    long int ** ppLongIntA = NULL;
    long long int ** ppLongLongIntA = NULL;

    clock_t begin_time;
    clock_t end_time;

    int int_type=0;
    int int_type_new=0;

    // *************************************************************************
    begin_time = clock();

    int_type = getIntegerType(A, d);
    cout << "Int type: " << int_type << endl;
    switch (int_type)
    {
        case TYPE_MPZ_T:
            break;
        case TYPE_LONG_LONG_INT:
            copy_basis_mpz_Z_T (ppLongLongIntA, A, d); 
            break;
        case TYPE_LONG_INT:
            copy_basis_mpz_Z_T (ppLongIntA, A, d); 
            break;
        case TYPE_INT:
            copy_basis_mpz_Z_T (ppIntA, A, d); 
            break;
        default:
            cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
            exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
    };

    int sublattice_start, sublattice_end;
    sublattice_start = start_index;
    sublattice_end = 0;
    while (sublattice_start < end_index)
    {
        sublattice_end += sublattice_size;
        sublattice_end = (sublattice_end > end_index) ? end_index : sublattice_end;
        cout << endl;

        switch (int_type)
        {
            case TYPE_MPZ_T:
		cout << "in switch case mpz_t..." << endl;
                prec_status = lcpgg_locally_restricted_CFA_LSR (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
		cout << "in switch case long long int..." << endl;
                prec_status = lcpgg_locally_restricted_CFA_LSR (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
		cout << "in switch case long int..." << endl;
                prec_status = lcpgg_locally_restricted_CFA_LSR (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
		cout << "in switch case int..." << endl;
                prec_status = lcpgg_locally_restricted_CFA_LSR (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
        sublattice_start += sublattice_size;
    }

    int_type_new = getIntegerType(A, d);

    if (sublattice_size != d)
    {
        cout << "Preprocessing complete. Now running LCGG globally." << endl;
        switch (int_type)
        {
            case TYPE_MPZ_T:
                prec_status = lcpgg_locally_restricted_CFA_LSR (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = lcpgg_locally_restricted_CFA_LSR (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = lcpgg_locally_restricted_CFA_LSR (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = lcpgg_locally_restricted_CFA_LSR (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
    }
    end_time = clock();
    // *************************************************************************

    time_taken = (long double) (end_time - begin_time) / CLOCKS_PER_SEC;
    cout << "Time=" << time_taken;
    RHF = compute_RHF (B, d, start_index, end_index);
    cout << "\tRHF=" << RHF << endl;
    last_GSO_length = B[d-1];
    cout << "Final GSO Length = " << B[d-1] << endl << endl;
            
    return prec_status;
}

// Wrapper for LC-PGG reduction
// Locally restricted insertions
int lcpgg_global_wrapper (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta_LC,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
                long double &time_taken,
                FP_NR<mpfr_t> &RHF,
                FP_NR<mpfr_t> &last_GSO_length,
                int start_index = 0,
                int end_index = 0,
                int sublattice_size = DEFAULT_SUBLATTICE_SIZE,
		int beta = DEFAULT_BETA
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    sublattice_size = (sublattice_size <= 0) ? d : sublattice_size;

    cout << "Precision = " << precision << endl;
    cout << "Using standard GSO" << endl;
    FP_NR<mpfr_t>::set_prec(precision);
    if (beta == 0) 
    {
        cout << "No beta on insertions." << endl;
        cout << "Running standard LC-PGG" << endl;
	beta = d;
    }
    else
    {
        cout << "Beta = " << beta << endl;
    }
    
    FP_mat<mpfr_t> gso(d,d);
    FP_mat<mpfr_t> mu(d,d);
    vector<FP_NR<mpfr_t>> B;
    B.resize(d,0.0);

    int prec_status; 

    int ** ppIntA = NULL;
    long int ** ppLongIntA = NULL;
    long long int ** ppLongLongIntA = NULL;

    clock_t begin_time;
    clock_t end_time;

    int int_type=0;
    int int_type_new=0;

    // *************************************************************************
    begin_time = clock();

    int_type = getIntegerType(A, d);
    cout << "Int type: " << int_type << endl;
    switch (int_type)
    {
        case TYPE_MPZ_T:
            break;
        case TYPE_LONG_LONG_INT:
            copy_basis_mpz_Z_T (ppLongLongIntA, A, d); 
            break;
        case TYPE_LONG_INT:
            copy_basis_mpz_Z_T (ppLongIntA, A, d); 
            break;
        case TYPE_INT:
            copy_basis_mpz_Z_T (ppIntA, A, d); 
            break;
        default:
            cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
            exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
    };

    int sublattice_start, sublattice_end;
    sublattice_start = start_index;
    sublattice_end = 0;
    while (sublattice_start < end_index)
    {
        sublattice_end += sublattice_size;
        sublattice_end = (sublattice_end > end_index) ? end_index : sublattice_end;
        cout << endl;

        switch (int_type)
        {
            case TYPE_MPZ_T:
		cout << "in switch case mpz_t..." << endl;
                prec_status = lcpgg_globally_restricted (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
		cout << "in switch case long long int..." << endl;
                prec_status = lcpgg_globally_restricted (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
		cout << "in switch case long int..." << endl;
                prec_status = lcpgg_globally_restricted (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
		cout << "in switch case int..." << endl;
                prec_status = lcpgg_globally_restricted (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
        sublattice_start += sublattice_size;
    }

    int_type_new = getIntegerType(A, d);

    if (sublattice_size != d)
    {
        cout << "Preprocessing complete. Now running LCGG globally." << endl;
        switch (int_type)
        {
            case TYPE_MPZ_T:
                prec_status = lcpgg_globally_restricted (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = lcpgg_globally_restricted (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = lcpgg_globally_restricted (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = lcpgg_globally_restricted (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
    }
    end_time = clock();
    // *************************************************************************

    time_taken = (long double) (end_time - begin_time) / CLOCKS_PER_SEC;
    cout << "Time=" << time_taken;
    RHF = compute_RHF (B, d, start_index, end_index);
    cout << "\tRHF=" << RHF << endl;
    last_GSO_length = B[d-1];
    cout << "Final GSO Length = " << B[d-1] << endl << endl;
            
    return prec_status;
}

// Wrapper for LC-GG reduction (CFA and LSR)
int lcpgg_global_cfa_wrapper (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta_LC,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
                long double &time_taken,
                FP_NR<mpfr_t> &RHF,
                FP_NR<mpfr_t> &last_GSO_length,
                int start_index = 0,
                int end_index = 0,
                int sublattice_size = DEFAULT_SUBLATTICE_SIZE,
		int beta = DEFAULT_BETA
                )
{
    end_index = (!end_index) ? d-1 : end_index;
    sublattice_size = (sublattice_size <= 0) ? d : sublattice_size;

    cout << "Precision = " << precision << endl;
    cout << "Using CFA and Lazy Size Reduction" << endl;
    FP_NR<mpfr_t>::set_prec(precision);
    if (beta == 0) 
    {
        cout << "No beta on insertions." << endl;
	beta = d;
    }
    else
    {
        cout << "Beta = " << beta << endl;
    }
    
    FP_mat<mpfr_t> gso(d,d);
    FP_mat<mpfr_t> mu(d,d);
    vector<FP_NR<mpfr_t>> B;
    B.resize(d,0.0);

    int prec_status; 

    int ** ppIntA = NULL;
    long int ** ppLongIntA = NULL;
    long long int ** ppLongLongIntA = NULL;

    clock_t begin_time;
    clock_t end_time;

    int int_type=0;
    int int_type_new=0;

    // *************************************************************************
    begin_time = clock();

    int_type = getIntegerType(A, d);
    cout << "Int type: " << int_type << endl;
    switch (int_type)
    {
        case TYPE_MPZ_T:
            break;
        case TYPE_LONG_LONG_INT:
            copy_basis_mpz_Z_T (ppLongLongIntA, A, d); 
            break;
        case TYPE_LONG_INT:
            copy_basis_mpz_Z_T (ppLongIntA, A, d); 
            break;
        case TYPE_INT:
            copy_basis_mpz_Z_T (ppIntA, A, d); 
            break;
        default:
            cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
            exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
    };

    int sublattice_start, sublattice_end;
    sublattice_start = start_index;
    sublattice_end = 0;
    while (sublattice_start < end_index)
    {
        sublattice_end += sublattice_size;
        sublattice_end = (sublattice_end > end_index) ? end_index : sublattice_end;
        cout << endl;

        switch (int_type)
        {
            case TYPE_MPZ_T:
		cout << "in switch case mpz_t..." << endl;
                prec_status = lcpgg_globally_restricted_CFA_LSR (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
		cout << "in switch case long long int..." << endl;
                prec_status = lcpgg_globally_restricted_CFA_LSR (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
		cout << "in switch case long int..." << endl;
                prec_status = lcpgg_globally_restricted_CFA_LSR (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
		cout << "in switch case int..." << endl;
                prec_status = lcpgg_globally_restricted_CFA_LSR (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
        sublattice_start += sublattice_size;
    }

    int_type_new = getIntegerType(A, d);

    if (sublattice_size != d)
    {
        cout << "Preprocessing complete. Now running LCGG globally." << endl;
        switch (int_type)
        {
            case TYPE_MPZ_T:
                prec_status = lcpgg_globally_restricted_CFA_LSR (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = lcpgg_globally_restricted_CFA_LSR (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = lcpgg_globally_restricted_CFA_LSR (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = lcpgg_globally_restricted_CFA_LSR (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            default:
                cout << "ERROR: BASIS INTEGER DATATYPE NOT SUPPORTED" << endl;
                exit(ERROR_INVALID_INTEGER_TYPE_IN_BASIS);
        };
    }
    end_time = clock();
    // *************************************************************************

    time_taken = (long double) (end_time - begin_time) / CLOCKS_PER_SEC;
    cout << "Time=" << time_taken;
    RHF = compute_RHF (B, d, start_index, end_index);
    cout << "\tRHF=" << RHF << endl;
    last_GSO_length = B[d-1];
    cout << "Final GSO Length = " << B[d-1] << endl << endl;
            
    return prec_status;
}
#endif
