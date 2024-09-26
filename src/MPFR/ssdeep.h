/* 
 * *****************************************************************************
 * Filename: ssdeep.h
 * Authors : Sanjay Bhattacherjee and Jack Moyler
 * *****************************************************************************
*/

#include "common.h"

#ifndef SSDEEP
#define SSDEEP

#define SS_DEEP_DELTA 0.999999

// SS-Deep index search
template <typename FP_T>
inline bool index_search_SSDeep (
                FP_mat<FP_T> &mu,
                vector<FP_NR<FP_T>> &B,
                int d,
                FP_NR<FP_T> &squared_sum,
                FP_NR<FP_T> delta_p,
                int k,
                int &i,
                int start_index = 0,
                int end_index = 0
                ) 
{
    end_index = (!end_index) ? d-1 : end_index;

    FP_NR<FP_T> SS_change = 0.0;
    FP_NR<FP_T> SS_change_max = 0.0;
    FP_NR<FP_T> projection_l_k = 0.0;
    FP_NR<FP_T> ratio_projection_l_k = 0.0;
    FP_NR<FP_T> inv_ratio_projection_l_k = 0.0;

    // D^(k)_{k-1} = (projection_l_k)^2
    projection_l_k = B[k] + (mu[k][k-1] * mu[k][k-1] * B[k-1]);
    // ratio_projection_l_k = projection_l_k / projection_l_l
    ratio_projection_l_k = projection_l_k / B[k-1];
    // SS_change = Old_SS - New_SS
    SS_change = mu[k][k-1] * mu[k][k-1] * (1.0 - ratio_projection_l_k);  

    i = k-1;
    SS_change_max = SS_change;

    for (int l=k-2; l>=0; l--) 
    {
        // D^(k)_{i} = (projection_l_k)^2
        projection_l_k += (mu[k][l] * mu[k][l] * B[l]);
        // ratio_projection_l_k = projection_l_k / projection_l_l
        inv_ratio_projection_l_k = B[l] / projection_l_k;
        // SS_change = Old_SS - New_SS
        SS_change += mu[k][l] * mu[k][l] * B[l] * (inv_ratio_projection_l_k - 1.0);  
            
        if (SS_change > SS_change_max) 
        {
            SS_change_max = SS_change;
	    i = l;
        }
    }

    if (SS_change_max > delta_p * squared_sum)
    {
        return true;
    } 
    return false;
}


// SS-Deep reduction for ZZ_mat<mpz_t> type basis
int ssdeep (
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
    FP_NR<mpfr_t> delta_p;
    delta_p = 1.0 - delta;
    
    FP_NR<mpfr_t> squared_sum;
    FP_NR<mpfr_t> previous_squared_sum;
    
    mpz_t nint_mu_z;
    mpz_init (nint_mu_z);
    mpfr_t nint_mu_f;
    mpfr_init2 (nint_mu_f, mpfr_prec_t (precision));
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    Z_NR<mpz_t> max_nint_mu, nint_mu, MAX_nint_mu;
    nint_mu = 0;
    max_nint_mu = 0; 
    MAX_nint_mu = MAX_INT; 

    // GSO - initial
    compute_GSO (A, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true;
    bool deep_insertion_required = true;
    bool index_search_required = true;

    int k = start_index+1;
    int i = start_index;
    int nInsertions = 0;
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

            // Index search
            squared_sum = compute_SS (squared_sum, B, d, start_index, end_index);
            deep_insertion_required = index_search_SS (mu, B, d, squared_sum, delta_p, k, i, start_index, end_index);

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
        // Check the SS-Deep_Reducedness
        compute_GSO (A, gso, mu, B, d, start_index, end_index);

        ssgg_reduced_check_internal(
                        mu,
                        B,
                        delta_p,
                        precision,
                        d,
                        eta,
                        k,
                        i,
                        squared_sum,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );

        if (squared_sum >= previous_squared_sum)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "previous_squared_sum " << previous_squared_sum << endl;
            cout << "squared_sum " << squared_sum << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        nLoops++;
    }

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    Z_NR<mpz_t> nint_mu_difference;
    nint_mu_difference.sub(MAX_nint_mu, max_nint_mu);
    cout << "MAX_nint_mu - max_nint_mu = " << nint_mu_difference << endl;

    return RED_SUCCESS;
}

// SS_Deep reduction for Z_T ** type basis
template <typename Z_T>
int ssdeep (
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

    FP_NR<mpfr_t> delta_p;
    delta_p = 1.0 - delta;
    
    FP_NR<mpfr_t> squared_sum;
    FP_NR<mpfr_t> previous_squared_sum;
    
    FP_NR<mpfr_t> nint_mu_FP = 0.0;
    mpfr_t nint_mu_f;
    mpfr_init2(nint_mu_f, mpfr_prec_t (precision));
    long int max_nint_mu = 0;
    long int nint_mu = 0;

    Z_T * pTemp;

    // GSO - initial
    compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

    bool size_reduction_required = true; 
    bool deep_insertion_required = true; 
    bool index_search_required = true; 

    int k = start_index+1;
    int i = start_index;
    int nInsertions = 0;
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
        previous_squared_sum = compute_SS (previous_squared_sum, B, d, start_index, end_index);

	while (k<d)
	{
            squared_sum = 0.0;
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

            // Index search
            squared_sum = compute_SS (squared_sum, B, d, start_index, end_index);
            deep_insertion_required = index_search_SSDeep (mu, B, d, squared_sum, delta_p, k, i, start_index, end_index);

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

    	        k = (1>(i)) ? 1 : i-1;
            }
            k++;
	}

        // Check the LC-Deep_Reducedness
        compute_GSO (ppA, gso, mu, B, d, start_index, end_index);

        ssgg_reduced_check_internal (
                        mu,
                        B,
                        delta_p,
                        precision,
                        d,
                        eta,
                        k,
                        i,
                        squared_sum,
                        size_reduction_required,
                        deep_insertion_required,
                        index_search_required,
                        start_index,
                        end_index
                        );

        if (squared_sum >= previous_squared_sum)
        {
            cout << "nLoops " << nLoops << endl;
            cout << "Precision " << precision << " insufficient." << endl;
            return PRECISION_FAIL;
        }

        nLoops++;
    }

    mpfr_clear(nint_mu_f);
    cout << "nInsertions = " << nInsertions << endl;
    cout << "max_nint_mu = " << max_nint_mu << endl;
    cout << "MAX_INT - max_nint_mu = " << MAX_INT - max_nint_mu << endl;
    cout << "MAX_LONG_INT - max_nint_mu = " << MAX_LONG_INT - max_nint_mu << endl;

    return RED_SUCCESS;
}

// Wrapper for LC-DeepLLL reduction
int ssdeep_wrapper (
                ZZ_mat<mpz_t> &A,
                FP_NR<mpfr_t> delta_SS,
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
                prec_status = ssdeep (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_LONG_INT:
		cout << "in switch case long long int..." << endl;
                prec_status = ssdeep (ppLongLongIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_LONG_INT:
		cout << "in switch case long int..." << endl;
                prec_status = ssdeep (ppLongIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
                break;
            case TYPE_INT:
		cout << "in switch case int..." << endl;
                prec_status = ssdeep (ppIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, sublattice_start, sublattice_end);
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
                prec_status = ssdeep (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_LONG_INT:
                prec_status = ssdeep (ppLongLongIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_LONG_INT:
                prec_status = ssdeep (ppLongIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
                break;
            case TYPE_INT:
                prec_status = ssdeep (ppIntA, delta_SS, d, precision, precision_correction_loops_allowed, eta, beta, gso, mu, B, start_index, end_index);
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
