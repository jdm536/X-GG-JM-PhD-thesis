/* 
 * *****************************************************************************
 * Filename: lcgg.h
 * Authors : Sanjay Bhattacherjee and Jack Moyler
 * *****************************************************************************
*/

#include "common.h"

#ifndef LCGG
#define LCGG

#define LC_DELTA 0.99

// Index search for LC-DeepLLL
template <typename FP_T>
inline bool index_search_LC_GG (
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B,
                int d,
                vector<FP_NR<FP_T>> &C,
                FP_NR<FP_T> delta_T,
                FP_NR<FP_T> &delta_kj,
                FP_NR<FP_T> &delta_min,
                int &k_max,
                int &i_max,
		int beta,
                int start_index = 0,
                int end_index = 0
                ) 
{
    end_index = (!end_index) ? d-1 : end_index;

    // Set delta_min to be the initial comparison (b[1] in position 0 in general)
    // Check that the minimum delta is allowed insertion at end by comparing with delta_T
    k_max = 1;
    i_max = 0; 
    delta_min = C[1] / B[0];

    for (int k=2; k<=end_index; k++) 
    {
	// Block restricted insertions not required if k <= 2*beta
        if (k <= 2*beta)
        { 
            for (int j=0; j<k; j++) 
            {
		delta_kj = (C[k] / B[j]); 
    	        if (delta_kj < delta_min)
	        {
		    k_max = k;
		    i_max = j;
	            delta_min = delta_kj;
	        } 
	        C[k] = C[k] - (mu[k][j] * mu[k][j] * B[j]);
            }
        }
        else
        {
    	    // First try insertion in the first beta vectors
            int index = k - beta;
            // for (int j=start_index; j<beta; j++) 
            for (int j=0; j<beta; j++) 
            {
		delta_kj = (C[k] / B[j]); 
    	        if (delta_kj < delta_min)
	        {
		    k_max = k;
		    i_max = j;
	            delta_min = delta_kj;
	        } 
	        C[k] = C[k] - (mu[k][j] * mu[k][j] * B[j]);
            }

	    for (int j=beta; j<index; j++)
	    {
	        C[k] = C[k] - (mu[k][j] * mu[k][j] * B[j]);
	    }

	    // Then try inserting in the beta vectors preceding b_k
	    for (int j=index; j<k; j++)
            {
		delta_kj = (C[k] / B[j]); 
    	        if (delta_kj < delta_min)
	        {
	            delta_min = delta_kj;
		    k_max = k;
		    i_max = j;
	        } 
	        C[k] = C[k] - (mu[k][j] * mu[k][j] * B[j]);
            }
        }
    }
    if (delta_min < delta_T)
    {
        return true;
    }

    return false;
}

// Index search for LC-DeepLLL (using CFA)
template <typename FP_T>
inline bool index_search_LC_GG_cfa (
                FP_mat<FP_T> &r,
                FP_mat<FP_T> &s,
                int d,
                FP_NR<FP_T> delta_T,
                FP_NR<FP_T> &delta_kj,
                FP_NR<FP_T> &delta_min,
                int &k_max,
                int &i_max,
		int beta,
                int start_index = 0,
                int end_index = 0
                ) 
{
    end_index = (!end_index) ? d-1 : end_index;

    // Set delta_min to be the initial comparison (b[1] in position 0 in general)
    // Check that the minimum delta is allowed insertion at end by comparing with delta_T
    k_max = 1;
    i_max = 0; 
    delta_min = s[1][0] / r[0][0];

    for (int k=2; k<=end_index; k++) 
    {
	// Block restricted insertions not required if k <= 2*beta
        if (k <= 2*beta)
        { 
            for (int j=0; j<k; j++) 
            {
		delta_kj = (s[k][j] / r[j][j]); 
    	        if (delta_kj < delta_min)
	        {
		    k_max = k;
		    i_max = j;
	            delta_min = delta_kj;
	        } 
            }
        }
        else
        {
    	    // First try insertion in the first beta vectors
            int index = k - beta;
            // for (int j=start_index; j<beta; j++) 
            for (int j=0; j<beta; j++) 
            {
		delta_kj = (s[k][j] / r[j][j]); 
    	        if (delta_kj < delta_min)
	        {
		    k_max = k;
		    i_max = j;
	            delta_min = delta_kj;
	        } 
            }

	    // Then try inserting in the beta vectors preceding b_k
	    for (int j=index; j<k; j++)
            {
		delta_kj = (s[k][j] / r[j][j]); 
    	        if (delta_kj < delta_min)
	        {
	            delta_min = delta_kj;
		    k_max = k;
		    i_max = j;
	        } 
            }
        }
    }
    if (delta_min < delta_T)
    {
        return true;
    }

    return false;
}

// LC-GG reduction for ZZ_mat<mpz_t> type basis
int lcgg (
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

    // GSO - initial
    compute_GSO (A, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int nInsertions = 0;
    int FirstPosLastInsertion = 0;
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
	// Compute initial values for C,  C[i]= <b[i],b[i]>
	for (int l=0; l<d; l++)
	{
            A[l].dot_product (C_int, A[l]);
	    C_int.abs(C_int);
	    C[l].set_z (C_int, MPFR_RNDN);
	}

        if (index_search_required)
        {
            // Index search - initial
            deep_insertion_required = index_search_LC_GG (mu, B, d, C, delta, delta_kj, delta_min, k, i, beta, start_index, end_index);

	    // Recording iteration number of insertions into position 0
            if (i == 0)
	    {
	        FirstPosLastInsertion = nInsertions;
	    }
        }

        // Main loop
        while (deep_insertion_required) 
        {            
            nInsertions++;
            totalInsertionDepth += (k-i);
            // Deep insert b_k in position i
            A.rotate_right(i, k);

            // Update the GSO     
            deep_LLL_GSO_update (mu, B, k, i, d, start_index, end_index); 

            // Size reduce from i to end_index
            for (int l=i; l<=end_index; l++)
            {
                for (int j=l-1; j>=start_index; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper(A, mu, d, l, j, 
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
	    for (int l=0; l<d; l++)
	    {
                A[l].dot_product (C_int, A[l]);
	        C_int.abs(C_int);
	        C[l].set_z (C_int, MPFR_RNDN);
	    }
	    // Index search
            deep_insertion_required = index_search_LC_GG (mu, B, d, C, delta, delta_kj, delta_min, k, i, beta, start_index, end_index);

	    // Recording iteration number of insertions into position 0
            if (i == 0)
	    {
	        FirstPosLastInsertion = nInsertions;
	    }
        }

        // Check the LC-GG-Reducedness
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
    cout << "Insertion Number of final insertion into Pos 0 = " << FirstPosLastInsertion << endl;
    Z_NR<mpz_t> nint_mu_difference;
    nint_mu_difference.sub(MAX_nint_mu, max_nint_mu);
    cout << "MAX_nint_mu - max_nint_mu = " << nint_mu_difference << endl;

    return RED_SUCCESS;
}

// LC-GG reduction for Z_T ** type basis
template <typename Z_T>
int lcgg (
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
    Z_T tmp;
    inner = 0;
    tmp = 0;
    Z_NR<mpz_t> inner_ZNR;
    inner_ZNR = 0;

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d,0.0);

    FP_NR<mpfr_t> delta_min = 0.0;
    FP_NR<mpfr_t> delta_kj = 0.0;

    // GSO - initial
    compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int nInsertions = 0;
    int FirstPosLastInsertion = 0;
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
	// Compute initial values for C,  C[i]= <b[i],b[i]>
	for (int l=0; l<d; l++)
	{
            inner = 0;
	    tmp = 0;
	    for (int j=0; j<d; j++)
	    {
		tmp = (ppA[l][j] * ppA[l][j]);
		if (tmp < 0)
		{
	            cout << "ERROR: Inner product term ppA[" << l << "][" << j <<"] < 0" << endl;
		    return ERROR_INVALID_INNER_PRODUCT; 
		}
	        inner += tmp;
	    }

	    if (inner < 0)
	    {
		cout << "ERROR: INNER PRODUCT < 0" << endl;
		return ERROR_INVALID_INNER_PRODUCT; 
	    }
	    inner_ZNR = inner;
	    C[l].set_z (inner_ZNR, MPFR_RNDN);

	    if (C[l] < 0)
	    {
		cout << "ERROR: C[l] < 0" << endl;
		return ERROR_INVALID_INNER_PRODUCT; 
	    }
	}

        if (index_search_required)
        {
            // Index search - initial
            deep_insertion_required = index_search_LC_GG (mu, B, d, C, delta, delta_kj, delta_min, k, i, beta, start_index, end_index);
            if (i == 0)
	    {
	        FirstPosLastInsertion = nInsertions;
	    }
        }

        // Loop
        while (deep_insertion_required) 
        {            
            nInsertions++;
            totalInsertionDepth += (k-i);
            // Deep insert b_k in position i
            pTemp = ppA[k];
            for (int j=k-1; j>=i; j--)
            {
                ppA[j+1] = ppA[j];
            }
            ppA[i] = pTemp;

            // Update the GSO     
            deep_LLL_GSO_update (mu, B, k, i, d, start_index, end_index); 

            // Size reduce from i to end_index
            for (int l=i; l<d; l++)
            {
                for (int j=l-1; j>=0; j--)
                {
                    // Takes mu[i][j]
                    // Computes nearest mpz_t to mu_[i][j]
                    if (mu[l][j] < -eta || mu[l][j] > eta)
                    {
                        size_reduce_wrapper(ppA, mu, d, l, j, 
                        nint_mu_f, 
                        nint_mu_FP,
                        max_nint_mu,
                        nint_mu,
                        start_index
                        );
                    }
                }
            }
	    for (int l=0; l<d; l++)
	    { 
		// A[l] is size reduced now, so update C[l] = <A[l],A[l]>
                inner = 0;
		inner_ZNR = 0;
	        for (int j=0; j<d; j++)
	        {
	            inner += (ppA[l][j] * ppA[l][j]);
	        }
	        inner_ZNR = inner;
	        C[l].set_z (inner_ZNR, MPFR_RNDN);
		if (C[l] <= 0)
		{
		    cout << "C[l] < 0 BEFORE INDEX SEARCH" << endl;
		}
	    }

            // Index search
            deep_insertion_required = index_search_LC_GG (mu, B, d, C, delta, delta_kj, delta_min, k, i, beta, start_index, end_index);
            if (i == 0)
	    {
	        FirstPosLastInsertion = nInsertions;
	    }
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
    cout << "Insertion Number of final insertion into Pos 0 = " << FirstPosLastInsertion << endl;
    cout << "max_nint_mu = " << max_nint_mu << endl;
    cout << "MAX_INT - max_nint_mu = " << MAX_INT - max_nint_mu << endl;
    cout << "MAX_LONG_INT - max_nint_mu = " << MAX_LONG_INT - max_nint_mu << endl;

    return RED_SUCCESS;
}

// LC-GG reduction for ZZ_mat<mpz_t> type basis
// Using CFA and lazy size reduction
int lcgg_CFA_LSR (
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
    compute_CFA (A_gram, r, mu, s, d, 0, 1); 

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int nInsertions = 0;
    int FirstPosLastInsertion = 0;
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
                lazy_size_reduce(A, A_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_z, nint_mu_vec, eta_bar, k, d, mu_tmp, z_tmp1, z_tmp2);
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        if (index_search_required)
        {
            // Index search - initial
            deep_insertion_required = index_search_LC_GG_cfa (r, s, d, delta_bar, delta_kj, delta_min, k, i, beta, start_index, end_index);
            if (i == 0)
	    {
	        FirstPosLastInsertion = nInsertions;
	    }
        }

        // Main loop 
        while (deep_insertion_required) 
        {            
            nInsertions++;
            totalInsertionDepth += (k-i);
            // Deep insert b_k in position i
            A.rotate_right(i, k);

	    for (int j=0; j<i; j++)
	    {
	        mu[i][j] = mu[k][j];
		r[i][j] = r[k][j];
	    }
	    r[i][i] = s[k][i];

	    // Update Gram matrix
	    for (int l=k; l>i; l--)
	    {
	        A_gram[l-1][l] = A_gram[l-1][l-1];
		A_gram[l-1][l-1] = A_gram[l][l-1];
		A_gram[l][l-1] = A_gram[l][l];
		A_gram[l][l] = 0;
		A_gram.swap_rows(l-1,l);

		Gram_elem_tmp = A_gram[l][l-1];
	        for (int j=l-1; j>i; j--)
		{
		    A_gram[l][j] = A_gram[l][j-1];
		}
		A_gram[l][i] = Gram_elem_tmp;
	    }

	    for (int l=d-1; l>k; l--)
	    {
		Gram_elem_tmp = A_gram[l][k];
                for (int j=k; j>i; j--)
		{
		    A_gram[l][j] = A_gram[l][j-1];
		}
		A_gram[l][i] = Gram_elem_tmp;
	    }

            // Size reduce from i to end_index
            for (int l=i; l<=end_index; l++)
            {
                lazy_size_reduce(A, A_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_z, nint_mu_vec, eta_bar, k, d, mu_tmp, z_tmp1, z_tmp2);
            }

	    // Index search
            deep_insertion_required = index_search_LC_GG_cfa (r, s, d, delta_bar, delta_kj, delta_min, k, i, beta, start_index, end_index);
            if (i == 0)
	    {
	        FirstPosLastInsertion = nInsertions;
	    }
        }
	// END OF MAIN LOOP

        compute_GSO (A, gso, mu, B, d, start_index, end_index);

        // Check the LC_GG_Reducedness
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
    cout << "Insertion Number of final insertion into Pos 0 = " << FirstPosLastInsertion << endl;
    Z_NR<mpz_t> nint_mu_difference;
    nint_mu_difference.sub(MAX_nint_mu, max_nint_mu);
    cout << "MAX_nint_mu - max_nint_mu = " << nint_mu_difference << endl;

    return RED_SUCCESS;
}

// LC-GG reduction for Z_T ** type basis
// Uses CFA and LSR
template <typename Z_T>
int lcgg_CFA_LSR (
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

    FP_NR<mpfr_t> mu_tmp;
    Z_NR<mpz_t> z_tmp1;
    Z_NR<mpz_t> z_tmp2;
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
    long long int inner;

    long long int ** ppA_gram;
    long long int * pGram_tmp;
    long long int Gram_elem_tmp;
    ppA_gram = (long long int **) calloc (d, sizeof(Z_T *));
    for (int i=0; i<d; i++)
    {
        ppA_gram[i] = (long long int *) calloc (d, sizeof(long long int));
    }
    pGram_tmp = (long long int *) calloc (d, sizeof(long long int));

    Z_NR<mpz_t> C_int;
    C_int = 0;
    vector<FP_NR<mpfr_t>> C;
    C.resize(d,0.0);

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
    compute_CFA (ppA_gram, r, mu, s, d, 0, 0); 


    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int nInsertions = 0;
    int FirstPosLastInsertion = 0;
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
	        lazy_size_reduce(ppA, ppA_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_vec, eta_bar, k, d, mu_tmp, z_tmp1, z_tmp2);
            }

            size_reduction_required = false;
            deep_insertion_required = true;
            index_search_required = true;
        }

        if (index_search_required)
        {
            // Index search - initial
            deep_insertion_required = index_search_LC_GG_cfa (r, s, d, delta_bar, delta_kj, delta_min, k, i, beta, start_index, end_index);
            if (i == 0)
	    {
	        FirstPosLastInsertion = nInsertions;
	    }
        }

        // Main loop
        while (deep_insertion_required) 
        {            
            nInsertions++;
            totalInsertionDepth += (k-i);
            // Deep insert b_k in position i
            pTemp = ppA[k];
            for (int j=k-1; j>=i; j--)
            {
                ppA[j+1] = ppA[j];
            }
            ppA[i] = pTemp;

	    for (int j=0; j<i; j++)
	    {
	        mu[i][j] = mu[k][j];
	        r[i][j] = r[k][j];
	    }
	    r[i][i] = s[k][i];

	    // Update Gram matrix
            for (int l=k; l>i; l--)
	    {
	        ppA_gram[l-1][l] = ppA_gram[l-1][l-1];
	        ppA_gram[l-1][l-1] = ppA_gram[l][l-1];
	        ppA_gram[l][l-1] = ppA_gram[l][l];
	        ppA_gram[l][l] = 0;

	        pGram_tmp = ppA_gram[l];
	        ppA_gram[l] = ppA_gram[l-1];
	        ppA_gram[l-1] = pGram_tmp;

		Gram_elem_tmp = ppA_gram[l][l-1];
	        for (int j=l-1; j>i; j--)
		{
		    ppA_gram[l][j] = ppA_gram[l][j-1];
		}
		ppA_gram[l][i] = Gram_elem_tmp;
	    }

	    for (int l=d-1; l>k; l--)
	    {
		Gram_elem_tmp = ppA_gram[l][k];
                for (int j=k; j>i; j--)
		{
		    ppA_gram[l][j] = ppA_gram[l][j-1];
		}
		ppA_gram[l][i] = Gram_elem_tmp;
	    }

            // Size reduce from i to end_index
            for (int l=i; l<=end_index; l++)
            {
	        lazy_size_reduce(ppA, ppA_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_vec, eta_bar, l, d, mu_tmp, z_tmp1, z_tmp2);
            }

            // Index search
            deep_insertion_required = index_search_LC_GG_cfa (r, s, d, delta_bar, delta_kj, delta_min, k, i, beta, start_index, end_index);
            if (i == 0)
	    {
	        FirstPosLastInsertion = nInsertions;
	    }
        }

        compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

        // Check the LC_GG_Reducedness
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
    cout << "Insertion Number of final insertion into Pos 0 = " << FirstPosLastInsertion << endl;
    cout << "max_nint_mu = " << max_nint_mu << endl;
    cout << "MAX_INT - max_nint_mu = " << MAX_INT - max_nint_mu << endl;
    cout << "MAX_LONG_INT - max_nint_mu = " << MAX_LONG_INT - max_nint_mu << endl;

    return RED_SUCCESS;
}

// Wrapper for LC-GG reduction
int lcgg_wrapper (
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
	    cout << "COPYING TO LONG LONG INT" << endl;
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
                prec_status = lcgg (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
		cout << "in switch case long long int..." << endl;
                prec_status = lcgg (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
		cout << "in switch case long int..." << endl;
                prec_status = lcgg (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
		cout << "in switch case int..." << endl;
                prec_status = lcgg (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
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
                prec_status = lcgg (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = lcgg (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = lcgg (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = lcgg (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
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
int lcgg_cfa_wrapper (
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
                prec_status = lcgg_CFA_LSR (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
		cout << "in switch case long long int..." << endl;
                prec_status = lcgg_CFA_LSR (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
		cout << "in switch case long int..." << endl;
                prec_status = lcgg_CFA_LSR (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
		cout << "in switch case int..." << endl;
                prec_status = lcgg_CFA_LSR (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
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
                prec_status = lcgg_CFA_LSR (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = lcgg_CFA_LSR (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = lcgg_CFA_LSR (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = lcgg_CFA_LSR (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
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
