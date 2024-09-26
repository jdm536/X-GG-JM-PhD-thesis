/* 
 * *****************************************************************************
 * Filename: lcdeep.h
 * Authors : Sanjay Bhattacherjee and Jack Moyler
 * *****************************************************************************
*/

#include "common.h"

#ifndef LCDEEP
#define LCDEEP

#define LC_DELTA 0.99

// Index search for LC-DeepLLL
template <typename FP_T>
inline bool index_search_LC_Deep (
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B,
                int d,
                FP_NR<FP_T> &C,
                FP_NR<FP_T> delta,
                int k,
                int &i,
		int beta
                ) 
{
    // If this case is true, then check all indices < k
    if (k <= 2*beta)
    { 
        for (int j=0; j<k; j++) 
        {
    	    if (C >= delta * B[j])
	    {
		if (C < 0) 
		{
		    cout << "Error: C < 0" << endl;
		    cout << "k = " << k << "\ti = " << j << endl;
		}
	        C = C - (mu[k][j] * mu[k][j] * B[j]);
	    } 
	    else
	    {
	        i = j;
	        return true;
	    }
        }
    }
    // Else, check for insertion in indices (1,...,beta) and (beta-k,...,k-1)
    else
    {
	// First try insertion in the first beta vectors
        int index = k - beta;
	for (int j=0; j<beta; j++)
	{
    	    if (C >= delta * B[j])
	    {
	        C = C - (mu[k][j] * mu[k][j] * B[j]);
	    } 
	    else
	    {
	        i = j;
	        return true;
	    }
	}

	for (int j=beta; j<index; j++)
	{
	    C = C - (mu[k][j] * mu[k][j] * B[j]);
	}

	// Then try inserting in the beta vectors preceding b_k
	for (int j=index; j<k; j++)
	{
    	    if (C >= delta * B[j])
	    {
	        C = C - (mu[k][j] * mu[k][j] * B[j]);
	    } 
	    else
	    {
	        i = j;
	        return true;
	    }
	}
    }
    
    return false;
}

// Index search for LC-DeepLLL
// CFA and LSR
template <typename FP_T>
inline bool index_search_LC_Deep_CFA (
                FP_mat<FP_T> &r,
                FP_mat<FP_T> &s,
                FP_NR<FP_T> delta_bar,
                int k,
                int &i,
		int beta
                ) 
{
    // If this case is true, then check all indices < k
    if (k <= 2*beta)
    { 
        for (int j=0; j<k; j++) 
        {
    	    if (delta_bar * r[j][j] > s[k][j])
	    {
                i = j;
		return true;
	    } 
	}
    }
    // Else, check for insertion in indices (1,...,beta) and (beta-k,...,k-1)
    else
    {
	// First try insertion in the first beta vectors
        int index = k - beta;
	for (int j=0; j<beta; j++)
	{
    	    if (delta_bar * r[j][j] > s[k][j])
	    {
	        i = j;
	        return true;
	    } 
	}

	// Then try inserting in the beta vectors preceding b_k
	for (int j=index; j<k; j++)
	{
    	    if (delta_bar * r[j][j] > s[k][j])
	    {
	        i = j;
	        return true;
	    } 
	}
    }
    
    return false;
}

template <typename FP_T>
void lcdeep_reduced_check_internal (
		ZZ_mat<mpz_t> &A,
		FP_mat<FP_T> &mu,
		vector<FP_NR<FP_T>> &B,
		FP_NR<FP_T> delta,
		int precision,
		int d,
		FP_NR<FP_T> eta,
		int beta, 
		Z_NR<mpz_t> &C_int,
		FP_NR<mpfr_t> &C,
		int &k,
		int &i,
		bool &deep_insertion_required,
		bool &basis_reduced,
		int start_index = 0,
		int end_index = 0
		)
{
    end_index = (!end_index) ? d-1 : end_index;

    deep_insertion_required = false;
    basis_reduced = false;
    
    // Check the basis is size reduced
    for (int l=start_index+1; l<=end_index; l++) 
    {
	for (int j=l-1; j>=start_index; j--)
	{
	    if (mu[l][j] < -eta || mu[l][j] > eta)
	    {
		cout << "ERROR: BASIS NOT SIZE REDUCED" << endl;
		cout << "mu[k][j] = " << mu[l][j] << endl;
		return;
	    }
	}
    }
    
    // Checks if more deep insertions are possible
    for (int l=1; l<d; l++)
    {
        A[l].dot_product (C_int, A[l]);
	C.set_z (C_int, MPFR_RNDN);

        deep_insertion_required = index_search_LC_Deep (mu, B, d, C, delta, l, i, beta);

	if (deep_insertion_required) 
	{
	    cout << "ERROR: BASIS SIZE REDUCED BUT NOT LC-DEEP REDUCED" << endl;
	    cout << "k = " << l << "\ti = " << i << endl;
	    return;
	}
    }

    basis_reduced = true;
    return;
}

template <typename FP_T, typename Z_T>
void lcdeep_reduced_check_internal (
		Z_T ** &ppA,
		FP_mat<FP_T> &mu,
		vector<FP_NR<FP_T>> &B,
		FP_NR<FP_T> delta,
		int precision,
		int d,
		FP_NR<FP_T> eta,
		int beta,
		Z_T &C_int,
		FP_NR<mpfr_t> &C,
		int &k,
		int &i,
		bool &deep_insertion_required,
		bool &basis_reduced,
		int start_index = 0,
		int end_index = 0
		)
{
    Z_T inner;
    inner = 0;
    Z_NR<mpz_t> inner_ZNR;
    end_index = (!end_index) ? d-1 : end_index;

    deep_insertion_required = false;
    basis_reduced = false;
    
    // Checks size reducedness
    for (int l=start_index+1; l<=end_index; l++) 
    {
	for (int j=l-1; j>=start_index; j--)
	{
	    if (mu[l][j] < -eta || mu[l][j] > eta)
	    {
		cout << "ERROR: BASIS NOT SIZE REDUCED" << endl;
		cout << "mu[l][j] = " << mu[l][j] << endl;
		return;
	    }
	}
    }
    
    // Checks if any more deep insertions are possible
    for (int l=1; l<d; l++)
    {
        inner = 0;
	for (int j=0; j<d; j++)
	{
	    inner += ppA[l][j] * ppA[l][j];
	}
	inner_ZNR = inner;
	C.set_z (inner_ZNR, MPFR_RNDN);

        deep_insertion_required = index_search_LC_Deep (mu, B, d, C, delta, l, i, beta);

	if (deep_insertion_required) 
	{
	    cout << "ERROR: BASIS SIZE REDUCED BUT NOT LC-DEEP REDUCED" << endl;
	    cout << "k = " << l << "\ti = " << i << endl;
	    return;
	}
    }

    basis_reduced = true;
    return;
}

// LC-DeepLLL reduction for ZZ_mat<mpz_t> type basis
int lcdeep (
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
    FP_NR<mpfr_t> C; 
    C = 0.0;

    // GSO - initial
    compute_GSO (A, gso, mu, B, d, start_index, end_index);

    bool deep_insertion_required = true;
    bool index_search_required = true;
    bool basis_reduced = false;

    int k = start_index+1;
    int i = start_index;
    int nInsertions = 0;
    int nLoops = 1;

    // Allows for recomputation of the GSO and additional reduction nLoops times
    while (basis_reduced = false)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        while (k<d)
        {
            // Size reduction of b_k
            for (int j=k-1; j>=start_index; j--)
       	    {
                // Takes mu[i][j]
                // Computes nearest mpz_t to mu_[i][j]
                if (mu[k][j] < -eta || mu[k][j] > eta)
                {
                    size_reduce_wrapper (A, mu, d, k, j, 
                    nint_mu_f, 
                    nint_mu_FP,
                    nint_mu_z,
                    max_nint_mu,
                    nint_mu,
                    start_index
                    );
                }
            }

            A[k].dot_product (C_int, A[k]);
	    C.set_z (C_int, MPFR_RNDN);
            // Index search
            deep_insertion_required = index_search_LC_Deep (mu, B, d, C, delta, k, i, beta);

	    // Perform deep insertion and GSO update
            if (deep_insertion_required) 
            {            
                nInsertions++;
                // Deep insert b_k in position i
                A.rotate_right(i, k);

                // Update the GSO     
                deep_LLL_GSO_update (mu, B, k, i, d, start_index, end_index); 

		k=i-1;
            }

	    k++; 
        }
        // Check the LC-Deep-Reducedness
        compute_GSO (A, gso, mu, B, d, start_index, end_index);

        lcdeep_reduced_check_internal(
			A,
                        mu,
                        B,
                        delta,
                        precision,
                        d,
                        eta,
			beta,
			C_int,
			C,
                        k,
                        i,
                        deep_insertion_required,
		        basis_reduced,
                        start_index,
                        end_index
                        );
        nLoops++;
    }

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    Z_NR<mpz_t> nint_mu_difference;
    nint_mu_difference.sub(MAX_nint_mu, max_nint_mu);
    cout << "MAX_nint_mu - max_nint_mu = " << nint_mu_difference << endl;

    return RED_SUCCESS;
}

// LC_Deep reduction for Z_T ** type basis
template <typename Z_T>
int lcdeep (
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
    
    Z_T C_int = 0;
    FP_NR<mpfr_t> C; 
    C = 0.0;

    Z_T * pTemp;
    Z_T inner;
    inner = 0;
    Z_NR<mpz_t> inner_ZNR;

    // GSO - initial
    compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

    bool deep_insertion_required = true;
    bool index_search_required = true;
    bool basis_reduced = false;

    int k = start_index+1;
    int i = start_index;
    int nInsertions = 0;
    int nLoops = 1;

    while (basis_reduced == false)
    {
        if (nLoops > precision_correction_loops_allowed)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

	while (k<d)
	{
            for (int j=k-1; j>=start_index; j--)
       	    {
                if (mu[k][j] < -eta || mu[k][j] > eta)
                {
                    size_reduce_wrapper(ppA, mu, d, k, j, 
                    nint_mu_f, 
                    nint_mu_FP,
                    max_nint_mu,
                    nint_mu,
                    start_index
                    );
                }
	    }

            inner = 0;
	    for (int j=0; j<d; j++)
	    {
	        inner += ppA[k][j] * ppA[k][j];
	    }
	    inner_ZNR = inner;
	    C.set_z (inner_ZNR, MPFR_RNDN);

            // Index search
            deep_insertion_required = index_search_LC_Deep (mu, B, d, C, delta, k, i, beta);

	    // Perform deep insertion and GSO update
            if (deep_insertion_required) 
            {            
                nInsertions++;
                // Deep insert b_k in position i
                pTemp = ppA[k];
                for (int j=k-1; j>=i; j--)
                {
                    ppA[j+1] = ppA[j];
                }
                ppA[i] = pTemp;

                // Update the GSO     
                deep_LLL_GSO_update (mu, B, k, i, d, start_index, end_index); 

    	        k = i-1;
            }
            k++;
	}

        // Check the LC-Deep_Reducedness
        compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

        lcdeep_reduced_check_internal(
			ppA,
                        mu,
                        B,
                        delta,
                        precision,
                        d,
                        eta,
			beta,
			C_int,
			C,
                        k,
                        i,
                        deep_insertion_required,
		        basis_reduced,
                        start_index,
                        end_index
                        );
	cout << "Reduced check done. " << endl;
        nLoops++;
    }

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    cout << "max_nint_mu = " << max_nint_mu << endl;
    cout << "MAX_INT - max_nint_mu = " << MAX_INT - max_nint_mu << endl;
    cout << "MAX_LONG_INT - max_nint_mu = " << MAX_LONG_INT - max_nint_mu << endl;

    return RED_SUCCESS;
}

// LC-Deep reduction for ZZ_mat<mpz_t> type basis
// Uses CFA and lazy size reduction (LSR)
int lcdeep_CFA_LSR (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
                int beta,
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
    Z_NR<mpz_t> max_nint_mu, MAX_nint_mu, Gram_elem_tmp;
    vector<Z_NR<mpz_t>> nint_mu_vec;
    nint_mu_vec.resize(d);
    max_nint_mu = 0;
    MAX_nint_mu = MAX_INT;
    Gram_elem_tmp = 0; 
    FP_NR<mpfr_t> delta_bar;
    delta_bar = (delta + 1.0) / 2;
    FP_NR<mpfr_t> eta_bar;
    eta_bar = (eta + 0.5) / 2;
    FP_NR<mpfr_t> mu_tmp;
    Z_NR<mpz_t> z_tmp1;
    Z_NR<mpz_t> z_tmp2;

    bool insertion_required = false;
    
    int * insertion_positions;
    insertion_positions = (int *) calloc (d, sizeof(int));
    for (int i=0; i<d; i++)
    {
        insertion_positions[i] = 0;
    }

    int nInsertions = 0;

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

    // Initial CFA on index 0
    // This updates the r, mu, s required for future reduction
    compute_CFA (A_gram, r, mu, s, d, 0, 1);

    int k = 1;
    int i = 0;

    while (k < d)
    {

	// Size reduce index k and update CFA
        lazy_size_reduce(A, A_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_z, nint_mu_vec, eta_bar, k, d, mu_tmp, z_tmp1, z_tmp2);
        insertion_required = index_search_LC_Deep_CFA (r, s, delta_bar, k, i, beta);
	cout << "Insertion required = " << insertion_required << endl;

	if (insertion_required)
	{
	    nInsertions++;

	    // Update basis - Insert b_k in position i
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
	    insertion_positions[i] += 1;

	    // Reduce k
	    k = (1 > (i)) ? 0 : i-1;   
	}   
	k++;
	insertion_required = false;
    }

    mpfr_clear(nint_mu_f);
    for (int i=0; i<d; i++) 
    {
        B[i] = r[i][i];
    }

    return RED_SUCCESS;
}

// LC-DeepLLL using standard integer datatypes
template <typename Z_T>
int lcdeep_CFA_LSR (
                Z_T ** &ppA,
                FP_NR<mpfr_t> delta,
                int d,
                int precision,
                int precision_correction_loops_allowed,
                FP_NR<mpfr_t> eta,
                int beta,
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
    Z_NR<mpz_t> max_nint_mu, MAX_nint_mu;
    Z_T * nint_mu_vec;
    nint_mu_vec = (Z_T *) calloc (d, sizeof(Z_T));
    max_nint_mu = 0;
    MAX_nint_mu = MAX_INT;
    FP_NR<mpfr_t> delta_bar;
    delta_bar = (delta + 1.0) / 2;
    FP_NR<mpfr_t> eta_bar;
    eta_bar = (eta + 0.5) / 2;
    FP_NR<mpfr_t> mu_tmp;
    Z_NR<mpz_t> z_tmp1;
    Z_NR<mpz_t> z_tmp2;
    
    Z_T * insertion_positions;
    insertion_positions = (Z_T *) calloc (d, sizeof(Z_T));
    for (int i=0; i<d; i++)
    {
        insertion_positions[i] = 0;
    }

    int nInsertions = 0;
    bool insertion_required = false;

    Z_T * pTemp;

    long long int ** ppA_gram;
    long long int * pGram_tmp;
    long long int Gram_elem_tmp;
    Gram_elem_tmp = 0;

    long potIncrease = 0;
    FP_NR<mpfr_t> currentPot = 0.0;
    FP_NR<mpfr_t> prevPot;
    prevPot = basisPotential (ppA, d);

    
    ppA_gram = (long long int **) calloc (d, sizeof(long long int *));
    for (int i=0; i<d; i++) 
    {
        ppA_gram[i] = (long long int *) calloc (d, sizeof(long long int));
    }

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

    // Compute Gram matrix
    compute_Gram_Matrix (ppA, ppA_gram, d);

    // Initial CFA on index 0
    compute_CFA (ppA_gram, r, mu, s, d, 0, 1);

    int k = 1;
    int i = 1;

    while (k < d)
    {
	// Size reduce index k
	lazy_size_reduce(ppA, ppA_gram, r, mu, s, nint_mu_f, nint_mu_FP, nint_mu_vec, eta_bar, k, d, mu_tmp, z_tmp1, z_tmp2);

	// Lovasz condition check for i=1,...,k-1
	for (int i=0; i<k; i++) 
	{ 
	    while (r[i][i] < 0 || s[k][i] < 0) 
	    {
		cout << "i = " << i << "\tk = " << k << endl;
                compute_Gram_Matrix (ppA, ppA_gram, d);
                compute_CFA (ppA_gram, r, mu, s, d, d, 0);
		cout << "RECOMPUTING GRAM" << endl;
	    }
	}
        insertion_required = index_search_LC_Deep_CFA (r, s, delta_bar, k, i, beta);

	if (insertion_required)
	{
            nInsertions++;

	    // Deep insert b_k in position i
	    pTemp = ppA[k];
            for (int j=k; j>i; j--)
            {
                ppA[j] = ppA[j-1];
            }
            ppA[i] = pTemp;
	
    	    // Update mu and r matrices
	    for (int j=0; j<i; j++) 
	    {
	        mu[i][j] = mu[k][j];
	        r[i][j] = r[k][j];
	    }
	    r[i][i] = s[k][i];
                
	    // Gram Matrix updates 
            for (int l=k; l>i; l--)
	    {
	        ppA_gram[l-1][l] = ppA_gram[l-1][l-1];
	        ppA_gram[l-1][l-1] = ppA_gram[l][l-1];
	        ppA_gram[l][l-1] = ppA_gram[l][l];
	        ppA_gram[l][l] = 0;

	        // SWAP ROWS l and l-1
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
	    insertion_positions[i] += 1;
	    
	    k = (1 > (i)) ? 0 : i-1;   
	}   
	k++;
	insertion_required = false;
    }

    mpfr_clear(nint_mu_f);
    for (int i=0; i<d; i++) 
    {
        B[i] = r[i][i];
    }

    return RED_SUCCESS;
}

// Wrapper for LC-DeepLLL reduction
int lcdeep_wrapper (
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
                prec_status = lcdeep (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
		cout << "in switch case long long int..." << endl;
                prec_status = lcdeep (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
		cout << "in switch case long int..." << endl;
                prec_status = lcdeep (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
		cout << "in switch case int..." << endl;
                prec_status = lcdeep (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
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
        cout << "Preprocessing complete. Now running SSGG globally." << endl;
        switch (int_type)
        {
            case TYPE_MPZ_T:
                prec_status = lcdeep (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = lcdeep (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = lcdeep (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = lcdeep (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
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

// Wrapper for LC-DeepLLL reduction
int lcdeep_cfa_wrapper (
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
                prec_status = lcdeep_CFA_LSR (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
		cout << "in switch case long long int..." << endl;
                prec_status = lcdeep_CFA_LSR (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
		cout << "in switch case long int..." << endl;
                prec_status = lcdeep_CFA_LSR (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
		cout << "in switch case int..." << endl;
                prec_status = lcdeep_CFA_LSR (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, mu, B, sublattice_start, sublattice_end);
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
        cout << "Preprocessing complete. Now running SSGG globally." << endl;
        switch (int_type)
        {
            case TYPE_MPZ_T:
                prec_status = lcdeep_CFA_LSR (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = lcdeep_CFA_LSR (ppLongLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = lcdeep_CFA_LSR (ppLongIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = lcdeep_CFA_LSR (ppIntA, delta_LC, d, precision, precision_correction_loops_allowed, eta, beta, mu, B, start_index, end_index);
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
    last_GSO_length = B[d-1];
    cout << "\tRHF=" << RHF << endl;
    cout << "Final GSO Length = " << B[d-1] << endl << endl;
            
    return prec_status;
}
#endif
