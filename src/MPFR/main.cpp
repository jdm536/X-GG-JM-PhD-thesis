/* 
 * *****************************************************************************
 * Filename: main.cpp
 * Authors : Sanjay Bhattacherjee and Jack Moyler
 * *****************************************************************************
*/
#include <iostream>     
#include <iomanip>
#include "ssdeep.h"
#include "ssgg.h"
#include "potdeep.h"
#include "potgg.h"
#include "lcdeep.h"
#include "lcgg.h"
#include "pgg.h"

int main (int argc, char * argv[])
{
    if (argc<3 || argc>9) 
    {
        cout << "Usage: <executable> <input-file-name> <algorithm> [precision = 53] [beta = 0] [orthogonalisation-switch = 0] [search-switch = 0] [sublattice_size = 0] [precision-correction-loops-allowed = 2 (std GS0) = 1 (CFA)] [eta = 0.51]" << endl;
        cout << "Algorithm:" << endl;
        cout << " 1 :  SS-DeepLLL" << endl;
        cout << " 2 :  SS-GGLLL" << endl;
        cout << " 3 :  Pot-DeepLLL" << endl;
        cout << " 4 :  Pot-GGLLL" << endl;
        cout << " 5 :  LC-DeepLLL" << endl;
        cout << " 6 :  LC-GG" << endl;
        cout << " 7 :  LC-PGG" << endl;
        cout << " 8 :  LC-PGG (locally restricted insertions)" << endl;
        cout << " 9 :  LC-PGG (globally restricted insertions)" << endl;
        cout << " 10 :  No reduction" << endl;
        cout << "  NOTE: Option 10 does not perform reduction. It computes the RHF, V1 length etc. for the basis taken as input." << endl;
        cout << "\tIf the input basis is not LLL-preprocessed, large precision may be required for correct computation." << endl;
        cout << "precision:" << endl;
        cout << "\tthe number of bits of significand" << endl;
        cout << "beta:" << endl;
        cout << "\tblock size for restricted insertions (insertions in first beta or last beta positions)" << endl;
        cout << "orthogonalisation-switch:" << endl;
        cout << "\t0: standard GSO (default)" << endl;
        cout << "\t1: CFA and Lazy Size Reduction" << endl;
        cout << "search-switch: " << endl;
        cout << "\t0: no search, using only the specified precision" << endl;
        cout << "\t1: (min precision) binary search" << endl;
        cout << "\t2: (min precision) exhaustive search starting from the specified precision" << endl;
        cout << "\t3: (precision for min RHF + min time) exhaustive search starting from the specified precision" << endl;
        cout << "sublattice_size:" << endl;
        cout << "\tnon-zero if X-GG sublattice preprocessing is desired" << endl;
        cout << "precision-correction-loops-allowed: " << endl;
        cout << "\tk: GSO to be recomputed k-1 times" << endl;
        cout << "eta:" << endl;
        cout << "\tfloating point size reduction relaxation" << endl;
        exit (ERROR_COMMANDLINE_ARGUMENT_NUMBER);
    }
    std::setprecision(12);

    int precision;
    FP_NR<mpfr_t> RHF;
    FP_NR<mpfr_t> last_GSO_length;
    long double time_taken;
    int sublattice_size;
    int beta;
    int search_switch;
    int orthogonalisation_switch;
    int precision_correction_loops_allowed;
    FP_NR<mpfr_t> eta, delta_SS, delta_Pot, delta_LC;

    int precision_min_RHF, precision_min_runtime;
    FP_NR<mpfr_t> RHF_min, RHF_min_runtime;
    long double runtime_min, runtime_min_RHF;

    int precision_fail, precision_mid, precision_success;
    int precision_jump;
    FP_NR<mpfr_t> prev_RHF;

    string input_filename (argv[1]);
    int algo_num = atoi(argv[2]); 
    precision =
            (argc>=4) ? atoi (argv[3]) : DEFAULT_PRECISION;
    beta =
            (argc>=5) ? atoi (argv[4]) : DEFAULT_BETA;
    orthogonalisation_switch = 
            (argc>=6) ? atoi (argv[5]) : DEFAULT_GSO;
    search_switch = 
            (argc>=7) ? atoi (argv[6]) : DEFAULT_SEARCH_SWITCH;
    sublattice_size =
            (argc>=8) ? atoi (argv[7]) : DEFAULT_SUBLATTICE_SIZE;
    if (orthogonalisation_switch == 1)
    {
        precision_correction_loops_allowed = 
            (argc>=9) ? atoi (argv[8]) : DEFAULT_PRECISION_CORRECTION_LOOPS_ALLOWED_CFA;
    } else
    {
        precision_correction_loops_allowed = 
            (argc>=9) ? atoi (argv[8]) : DEFAULT_PRECISION_CORRECTION_LOOPS_ALLOWED;
    }
    eta = 
            (argc>=10) ? atoi (argv[9]) : ETA;
    FP_NR<mpfr_t>::set_prec(precision);
    delta_SS = SS_GG_DELTA;
    delta_Pot = POT_GG_DELTA;
    delta_LC = LC_DELTA;

    ZZ_mat<mpz_t> A;

    int status = 0;
    status |= read_file(A, input_filename.c_str());
    if (status) exit (ERROR_FILE_READ);

    int d = A.get_cols();

    /*
     * Comment/Uncomment line below to toggle input basis printing
    */	
    //cout << A << endl;
    int prec_status;

    if (algo_num == 1)
    {
	cout << endl << "ALGORITHM: SS-Deep" << endl << "------------------------------" << endl;
	FP_NR<mpfr_t> check;
        FP_NR<mpfr_t>::set_prec(precision);
        cout << endl << "Starting with precision = " << precision << endl;

	// delta = 1 - 10 ^{-6}
        delta_SS = 0.999999;
	
        cout << "Delta = " << std::setprecision(12) << delta_SS.get_d() << endl;

        prec_status = ssdeep_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

        cout << endl;

	// Binary search for minimum precision which ensures correct reduction
        if (search_switch == SEARCH_SWITCH_MIN_PRECISION_BINARY)
        {
            // Find an insufficient precison
            precision_fail = precision;
            precision_jump = PRECISION_STEP_DOWN;
            bool min_check = false;

            while (prec_status == RED_SUCCESS)
            {
                min_check = true;
		precision_success = precision_fail;
                cout << endl;
                cout << "Succeeded with precision_fail = " << precision_fail << endl;
                precision_fail -= precision_jump;
		precision_fail = (precision_fail <= 0) ? 1 : precision_fail;
                precision_jump *= PRECISION_JUMP_FACTOR;
                cout << "Checking with precision_fail = " << precision_fail << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssdeep_wrapper (A, delta_SS, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }

            // Find a sufficient precison
            if (!min_check)
            {
                precision_jump = PRECISION_STEP_UP;
                precision_success = precision_fail + precision_jump;

                cout << endl;
                cout << "Checking with precision_success = " << precision_success << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssdeep_wrapper (A, delta_SS, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

                while (prec_status != RED_SUCCESS)
                {
                    cout << "Failed with precision_success = " << precision_success << endl;
                    precision_fail = precision_success;
                    precision_jump *= PRECISION_JUMP_FACTOR;
                    precision_success += precision_jump;
                    cout << "Checking with precision_success = " << precision_success << endl;

                    status |= read_file(A, input_filename.c_str());
                    if (status) exit (ERROR_FILE_READ);
                    prec_status = ssdeep_wrapper (A, delta_SS, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
                }
            }

            if (precision_fail > precision_success)
            {
                cout << "Error: In finding precision range; precision_fail = " << precision_fail << " precision_success = " << precision_success << endl;
                exit (ERROR_INVALID_PRECISION_RANGE);
            }

            // Conduct binary search for the minimum precision
            cout << endl << "Starting binary search with" << endl;
            cout << "Inital precision_fail = " << precision_fail << endl;
            cout << "Inital precision_success = " << precision_success << endl << endl;
            while (precision_fail+1 < precision_success)
            {
                precision_mid = precision_fail + precision_success;
                precision_mid = (precision_mid % 2) ? (precision_mid+1)/2 : precision_mid/2; 
                cout << endl << "precision_fail = " << precision_fail << endl;
                cout << "precision_mid = " << precision_mid << endl;
                cout << "precision_success = " << precision_success << endl << endl;
                cout << "Checking with precision_mid = " << precision_mid << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = 
                    ssdeep_wrapper (A, delta_SS, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
                if (prec_status != RED_SUCCESS)
                    precision_fail = precision_mid;
                else
                    precision_success = precision_mid;
            }

            precision = precision_success;
            cout << "Running with Precision " << precision << endl;

            status |= read_file(A, input_filename.c_str());
            if (status) exit (ERROR_FILE_READ);
            prec_status = 
                ssdeep_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
        }
	// Exhaustive search for a precision which ensures correct reduction
	// Starts at the precision passed as input (default = 53), and increments precision by 1 if incorrect reduction 
        else if (search_switch == SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssdeep_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }
        }
        else if (search_switch == SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssdeep_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }
            precision_min_RHF = precision;
            precision_min_runtime = precision;
            RHF_min = RHF;
            RHF_min_runtime = RHF;
            runtime_min = time_taken;
            runtime_min_RHF = time_taken;
            do
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                prev_RHF = RHF;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssdeep_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

                if (prec_status == RED_SUCCESS)
                {
                    if (RHF < RHF_min)
                    {
                        RHF_min = RHF;
                        runtime_min_RHF = time_taken;
                        precision_min_RHF = precision;
                    }

                    if (time_taken < runtime_min)
                    {
                        RHF_min_runtime = RHF;
                        runtime_min = time_taken;
                        precision_min_runtime = precision;
                    }
                }
            }
            while (prev_RHF != RHF || prec_status != RED_SUCCESS);

            cout << "Min time:\ttime=" << runtime_min;
            cout << "\tRHF=" << RHF_min_runtime;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_runtime;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

            cout << "Min RHF:\ttime=" << runtime_min_RHF;
            cout << "\tRHF=" << RHF_min;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_RHF;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

        }
        else
        {
            // The default case has been handled already
        }
    } 
    // algo_num == 2 - SS-GGLLL
    else if (algo_num == 2)
    {
	cout << endl << "ALGORITHM: SS-GG" << endl << "------------------------------";
	FP_NR<mpfr_t> check;
        FP_NR<mpfr_t>::set_prec(precision);
        cout << endl;
        cout << endl << "Starting with precision = " << precision << endl;

	// delta = 1 - 10 ^{-6}
        delta_SS = 0.999999;
	
        cout << "Delta = " << std::setprecision(12) << delta_SS.get_d() << endl;

        prec_status = ssgg_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

        cout << endl;

	// Binary search for minimum precision which ensures correct reduction
        if (search_switch == SEARCH_SWITCH_MIN_PRECISION_BINARY)
        {
            // Find an insufficient precison
            precision_fail = precision;
            precision_jump = PRECISION_STEP_DOWN;
            bool min_check = false;

            while (prec_status == RED_SUCCESS)
            {
                min_check = true;
		precision_success = precision_fail;
                cout << endl;
                cout << "Succeeded with precision_fail = " << precision_fail << endl;
                precision_fail -= precision_jump;
		precision_fail = (precision_fail <= 0) ? 1 : precision_fail;
                precision_jump *= PRECISION_JUMP_FACTOR;
                cout << "Checking with precision_fail = " << precision_fail << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssgg_wrapper (A, delta_SS, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }

            // Find a sufficient precison
            if (!min_check)
            {
                precision_jump = PRECISION_STEP_UP;
                precision_success = precision_fail + precision_jump;

                cout << endl;
                cout << "Checking with precision_success = " << precision_success << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssgg_wrapper (A, delta_SS, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

                while (prec_status != RED_SUCCESS)
                {
                    cout << "Failed with precision_success = " << precision_success << endl;
                    precision_fail = precision_success;
                    precision_jump *= PRECISION_JUMP_FACTOR;
                    precision_success += precision_jump;
                    cout << "Checking with precision_success = " << precision_success << endl;

                    status |= read_file(A, input_filename.c_str());
                    if (status) exit (ERROR_FILE_READ);
                    prec_status = ssgg_wrapper (A, delta_SS, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
                }
            }

            if (precision_fail > precision_success)
            {
                cout << "Error: In finding precision range; precision_fail = " << precision_fail << " precision_success = " << precision_success << endl;
                exit (ERROR_INVALID_PRECISION_RANGE);
            }

            // Conduct binary search for the minimum precision
            cout << endl << "Starting binary search with" << endl;
            cout << "Inital precision_fail = " << precision_fail << endl;
            cout << "Inital precision_success = " << precision_success << endl << endl;
            while (precision_fail+1 < precision_success)
            {
                precision_mid = precision_fail + precision_success;
                precision_mid = (precision_mid % 2) ? (precision_mid+1)/2 : precision_mid/2; 
                cout << endl << "precision_fail = " << precision_fail << endl;
                cout << "precision_mid = " << precision_mid << endl;
                cout << "precision_success = " << precision_success << endl << endl;
                cout << "Checking with precision_mid = " << precision_mid << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = 
                    ssgg_wrapper (A, delta_SS, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
                if (prec_status != RED_SUCCESS)
                    precision_fail = precision_mid;
                else
                    precision_success = precision_mid;
            }

            precision = precision_success;
            cout << "Running with Precision " << precision << endl;

            status |= read_file(A, input_filename.c_str());
            if (status) exit (ERROR_FILE_READ);
            prec_status = 
                ssgg_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
        }
	// Exhaustive search for a precision which ensures correct reduction
	// Starts at the precision passed as input (default = 53), and increments precision by 1 if incorrect reduction 
        else if (search_switch == SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssgg_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }
        }
        else if (search_switch == SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssgg_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }
            precision_min_RHF = precision;
            precision_min_runtime = precision;
            RHF_min = RHF;
            RHF_min_runtime = RHF;
            runtime_min = time_taken;
            runtime_min_RHF = time_taken;
            do
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                prev_RHF = RHF;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = ssgg_wrapper (A, delta_SS, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

                if (prec_status == RED_SUCCESS)
                {
                    if (RHF < RHF_min)
                    {
                        RHF_min = RHF;
                        runtime_min_RHF = time_taken;
                        precision_min_RHF = precision;
                    }

                    if (time_taken < runtime_min)
                    {
                        RHF_min_runtime = RHF;
                        runtime_min = time_taken;
                        precision_min_runtime = precision;
                    }
                }
            }
            while (prev_RHF != RHF || prec_status != RED_SUCCESS);

            cout << "Min time:\ttime=" << runtime_min;
            cout << "\tRHF=" << RHF_min_runtime;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_runtime;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

            cout << "Min RHF:\ttime=" << runtime_min_RHF;
            cout << "\tRHF=" << RHF_min;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_RHF;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

        }
        else
        {
            // The default case has been handled already
        }
    } 
    // algo_num == 3 - Pot-DeepLLL
    else if (algo_num == 3) 
    {
	cout << endl << "ALGORITHM: POT-Deep" << endl << "------------------------------";
        cout << endl;
        cout << endl << "Starting with precision = " << precision << endl;

        prec_status = potdeep_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

        cout << endl;

	// Binary search for minimum precision which ensures correct reduction
        if (search_switch == SEARCH_SWITCH_MIN_PRECISION_BINARY)
        {
            // Find an insufficient precison
            precision_fail = precision;
            precision_jump = PRECISION_STEP_DOWN;
            bool min_check = false;

            while (prec_status == RED_SUCCESS)
            {
                min_check = true;
		precision_success = precision_fail;
                cout << endl;
                cout << "Succeeded with precision_fail = " << precision_fail << endl;
                precision_fail -= precision_jump;
                precision_jump *= PRECISION_JUMP_FACTOR;
                cout << "Checking with precision_fail = " << precision_fail << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potdeep_wrapper (A, delta_Pot, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }

            // Find a sufficient precison
            if (!min_check)
            {
                precision_jump = PRECISION_STEP_UP;
                precision_success = precision_fail + precision_jump;

                cout << endl;
                cout << "Checking with precision_success = " << precision_success << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potdeep_wrapper (A, delta_Pot, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

                while (prec_status != RED_SUCCESS)
                {
                    cout << "Failed with precision_success = " << precision_success << endl;
                    precision_fail = precision_success;
                    precision_jump *= PRECISION_JUMP_FACTOR;
                    precision_success += precision_jump;
                    cout << "Checking with precision_success = " << precision_success << endl;

                    status |= read_file(A, input_filename.c_str());
                    if (status) exit (ERROR_FILE_READ);
                    prec_status = potdeep_wrapper (A, delta_Pot, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
                }
            }

            if (precision_fail > precision_success)
            {
                cout << "Error: In finding precision range; precision_fail = " << precision_fail << " precision_success = " << precision_success << endl;
                exit (ERROR_INVALID_PRECISION_RANGE);
            }

            // Conduct binary search for the minimum precision
            cout << endl << "Starting binary search with" << endl;
            cout << "Inital precision_fail = " << precision_fail << endl;
            cout << "Inital precision_success = " << precision_success << endl << endl;
            while (precision_fail+1 < precision_success)
            {
                precision_mid = precision_fail + precision_success;
                precision_mid = (precision_mid % 2) ? (precision_mid+1)/2 : precision_mid/2; 
                cout << endl << "precision_fail = " << precision_fail << endl;
                cout << "precision_mid = " << precision_mid << endl;
                cout << "precision_success = " << precision_success << endl << endl;
                cout << "Checking with precision_mid = " << precision_mid << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = 
                    potdeep_wrapper (A, delta_Pot, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
                if (prec_status != RED_SUCCESS)
                    precision_fail = precision_mid;
                else
                    precision_success = precision_mid;
            }

            precision = precision_success;
            cout << "Running with Precision " << precision << endl;

            status |= read_file(A, input_filename.c_str());
            if (status) exit (ERROR_FILE_READ);
            prec_status = 
                potdeep_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
        }
	// Exhaustive search for a precision which ensures correct reduction
	// Starts at the precision passed as input (default = 53), and increments precision by 1 if incorrect reduction 
        else if (search_switch == SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potdeep_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }
        }
        else if (search_switch == SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potdeep_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }
            precision_min_RHF = precision;
            precision_min_runtime = precision;
            RHF_min = RHF;
            RHF_min_runtime = RHF;
            runtime_min = time_taken;
            runtime_min_RHF = time_taken;
            do
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                prev_RHF = RHF;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potdeep_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

                if (prec_status == RED_SUCCESS)
                {
                    if (RHF < RHF_min)
                    {
                        RHF_min = RHF;
                        runtime_min_RHF = time_taken;
                        precision_min_RHF = precision;
                    }

                    if (time_taken < runtime_min)
                    {
                        RHF_min_runtime = RHF;
                        runtime_min = time_taken;
                        precision_min_runtime = precision;
                    }
                }
            }
            while (prev_RHF != RHF || prec_status != RED_SUCCESS);

            cout << "Min time:\ttime=" << runtime_min;
            cout << "\tRHF=" << RHF_min_runtime;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_runtime;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

            cout << "Min RHF:\ttime=" << runtime_min_RHF;
            cout << "\tRHF=" << RHF_min;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_RHF;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

        }
        else
        {
            // The default case has been handled already
        }
    } 
    // algo_num == 4 - Pot-GGLLL
    else if (algo_num == 4) 
    {
	cout << endl << "ALGORITHM: POT-GG" << endl << "------------------------------";
        cout << endl;
        cout << endl << "Starting with precision = " << precision << endl;

        prec_status = potgg_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

        cout << endl;

	// Binary search for minimum precision which ensures correct reduction
        if (search_switch == SEARCH_SWITCH_MIN_PRECISION_BINARY)
        {
            // Find an insufficient precison
            precision_fail = precision;
            precision_jump = PRECISION_STEP_DOWN;
            bool min_check = false;

            while (prec_status == RED_SUCCESS)
            {
                min_check = true;
		precision_success = precision_fail;
                cout << endl;
                cout << "Succeeded with precision_fail = " << precision_fail << endl;
                precision_fail -= precision_jump;
                precision_jump *= PRECISION_JUMP_FACTOR;
                cout << "Checking with precision_fail = " << precision_fail << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potgg_wrapper (A, delta_Pot, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }

            // Find a sufficient precison
            if (!min_check)
            {
                precision_jump = PRECISION_STEP_UP;
                precision_success = precision_fail + precision_jump;

                cout << endl;
                cout << "Checking with precision_success = " << precision_success << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potgg_wrapper (A, delta_Pot, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

                while (prec_status != RED_SUCCESS)
                {
                    cout << "Failed with precision_success = " << precision_success << endl;
                    precision_fail = precision_success;
                    precision_jump *= PRECISION_JUMP_FACTOR;
                    precision_success += precision_jump;
                    cout << "Checking with precision_success = " << precision_success << endl;

                    status |= read_file(A, input_filename.c_str());
                    if (status) exit (ERROR_FILE_READ);
                    prec_status = potgg_wrapper (A, delta_Pot, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
                }
            }

            if (precision_fail > precision_success)
            {
                cout << "Error: In finding precision range; precision_fail = " << precision_fail << " precision_success = " << precision_success << endl;
                exit (ERROR_INVALID_PRECISION_RANGE);
            }

            // Conduct binary search for the minimum precision
            cout << endl << "Starting binary search with" << endl;
            cout << "Inital precision_fail = " << precision_fail << endl;
            cout << "Inital precision_success = " << precision_success << endl << endl;
            while (precision_fail+1 < precision_success)
            {
                precision_mid = precision_fail + precision_success;
                precision_mid = (precision_mid % 2) ? (precision_mid+1)/2 : precision_mid/2; 
                cout << endl << "precision_fail = " << precision_fail << endl;
                cout << "precision_mid = " << precision_mid << endl;
                cout << "precision_success = " << precision_success << endl << endl;
                cout << "Checking with precision_mid = " << precision_mid << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                    prec_status = potgg_wrapper (A, delta_Pot, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
                if (prec_status != RED_SUCCESS)
                    precision_fail = precision_mid;
                else
                    precision_success = precision_mid;
            }

            precision = precision_success;
            cout << "Running with Precision " << precision << endl;

            status |= read_file(A, input_filename.c_str());
            if (status) exit (ERROR_FILE_READ);
            prec_status = potgg_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
        }
	// Exhaustive search for a precision which ensures correct reduction
	// Starts at the precision passed as input (default = 53), and increments precision by 1 if incorrect reduction 
        else if (search_switch == SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potgg_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }
        }
        else if (search_switch == SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potgg_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);
            }
            precision_min_RHF = precision;
            precision_min_runtime = precision;
            RHF_min = RHF;
            RHF_min_runtime = RHF;
            runtime_min = time_taken;
            runtime_min_RHF = time_taken;
            do
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                prev_RHF = RHF;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
                prec_status = potgg_wrapper (A, delta_Pot, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size);

                if (prec_status == RED_SUCCESS)
                {
                    if (RHF < RHF_min)
                    {
                        RHF_min = RHF;
                        runtime_min_RHF = time_taken;
                        precision_min_RHF = precision;
                    }

                    if (time_taken < runtime_min)
                    {
                        RHF_min_runtime = RHF;
                        runtime_min = time_taken;
                        precision_min_runtime = precision;
                    }
                }
            }
            while (prev_RHF != RHF || prec_status != RED_SUCCESS);

            cout << "Min time:\ttime=" << runtime_min;
            cout << "\tRHF=" << RHF_min_runtime;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_runtime;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

            cout << "Min RHF:\ttime=" << runtime_min_RHF;
            cout << "\tRHF=" << RHF_min;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_RHF;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

        }
        else
        {
            // The default case has been handled already
        }
    } 
    // algo_num == 5 - LC-DeepLLL
    else if (algo_num == 5) 
    {
	cout << endl << "ALGORITHM: LC-DEEP" << endl << "------------------------------";
        cout << endl;
        cout << endl << "Starting with precision = " << precision << endl;

	if (orthogonalisation_switch == DEFAULT_GSO)
	{
            prec_status = lcdeep_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	} else
	{
            prec_status = lcdeep_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	}

        cout << endl;

	// Binary search for minimum precision which ensures correct reduction
        if (search_switch == SEARCH_SWITCH_MIN_PRECISION_BINARY)
        {
            // Find an insufficient precison
            precision_fail = precision;
            precision_jump = PRECISION_STEP_DOWN;
            bool min_check = false;

            while (prec_status == RED_SUCCESS)
            {
                min_check = true;
		precision_success = precision_fail;
                cout << endl;
                cout << "Succeeded with precision_fail = " << precision_fail << endl;
                precision_fail -= precision_jump;
		precision_fail = (precision_fail <= 0) ? 1 : precision_fail;
                precision_jump *= PRECISION_JUMP_FACTOR;
                cout << "Checking with precision_fail = " << precision_fail << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);

	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcdeep_wrapper (A, delta_LC, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcdeep_cfa_wrapper (A, delta_LC, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        }
            }

            // Find a sufficient precison
            if (!min_check)
            {
                precision_jump = PRECISION_STEP_UP;
                precision_success = precision_fail + precision_jump;

                cout << endl;
                cout << "Checking with precision_success = " << precision_success << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);

	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcdeep_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcdeep_cfa_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
		} 

                while (prec_status != RED_SUCCESS)
                {
                    cout << "Failed with precision_success = " << precision_success << endl;
                    precision_fail = precision_success;
                    precision_jump *= PRECISION_JUMP_FACTOR;
                    precision_success += precision_jump;
                    cout << "Checking with precision_success = " << precision_success << endl;

                    status |= read_file(A, input_filename.c_str());
                    if (status) exit (ERROR_FILE_READ);
	            if (orthogonalisation_switch == DEFAULT_GSO)
	            {
                        prec_status = lcdeep_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	            } else
	            {
                        prec_status = lcdeep_cfa_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
		    } 
                }
            }

            if (precision_fail > precision_success)
            {
                cout << "Error: In finding precision range; precision_fail = " << precision_fail << " precision_success = " << precision_success << endl;
                exit (ERROR_INVALID_PRECISION_RANGE);
            }

            // Conduct binary search for the minimum precision
            cout << endl << "Starting binary search with" << endl;
            cout << "Inital precision_fail = " << precision_fail << endl;
            cout << "Inital precision_success = " << precision_success << endl << endl;
            while (precision_fail+1 < precision_success)
            {
                precision_mid = precision_fail + precision_success;
                precision_mid = (precision_mid % 2) ? (precision_mid+1)/2 : precision_mid/2; 
                cout << endl << "precision_fail = " << precision_fail << endl;
                cout << "precision_mid = " << precision_mid << endl;
                cout << "precision_success = " << precision_success << endl << endl;
                cout << "Checking with precision_mid = " << precision_mid << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcdeep_wrapper (A, delta_LC, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcdeep_cfa_wrapper (A, delta_LC, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
		} 

                if (prec_status != RED_SUCCESS)
                    precision_fail = precision_mid;
                else
                    precision_success = precision_mid;
            }

            precision = precision_success;
            cout << "Running with Precision " << precision << endl;

            status |= read_file(A, input_filename.c_str());
            if (status) exit (ERROR_FILE_READ);
	    if (orthogonalisation_switch == DEFAULT_GSO)
	    {
                prec_status = lcdeep_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	    } else
	    {
                prec_status = lcdeep_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	    } 
        }
	// Exhaustive search for a precision which ensures correct reduction
	// Starts at the precision passed as input (default = 53), and increments precision by 1 if incorrect reduction 
        else if (search_switch == SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcdeep_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcdeep_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }
        }
        else if (search_switch == SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcdeep_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcdeep_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }
            precision_min_RHF = precision;
            precision_min_runtime = precision;
            RHF_min = RHF;
            RHF_min_runtime = RHF;
            runtime_min = time_taken;
            runtime_min_RHF = time_taken;
            do
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                prev_RHF = RHF;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcdeep_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcdeep_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 

                if (prec_status == RED_SUCCESS)
                {
                    if (RHF < RHF_min)
                    {
                        RHF_min = RHF;
                        runtime_min_RHF = time_taken;
                        precision_min_RHF = precision;
                    }

                    if (time_taken < runtime_min)
                    {
                        RHF_min_runtime = RHF;
                        runtime_min = time_taken;
                        precision_min_runtime = precision;
                    }
                }
            }
            while (prev_RHF != RHF || prec_status != RED_SUCCESS);

            cout << "Min time:\ttime=" << runtime_min;
            cout << "\tRHF=" << RHF_min_runtime;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_runtime;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

            cout << "Min RHF:\ttime=" << runtime_min_RHF;
            cout << "\tRHF=" << RHF_min;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_RHF;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;
	}
        
        else
        {
            // The default case has been handled already
        }
    }
    // algo_num == 6 - LC-GG
    else if (algo_num == 6) 
    {
	cout << endl << "ALGORITHM: LC-GG" << endl << "------------------------------";
        cout << endl;
        cout << endl << "Starting with precision = " << precision << endl;
	if (orthogonalisation_switch == DEFAULT_GSO)
	{
            prec_status = lcgg_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	} else
	{
            prec_status = lcgg_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	} 

        cout << endl;

	// Binary search for minimum precision which ensures correct reduction
        if (search_switch == SEARCH_SWITCH_MIN_PRECISION_BINARY)
        {
            // Find an insufficient precison
            precision_fail = precision;
            precision_jump = PRECISION_STEP_DOWN;
            bool min_check = false;

            while (prec_status == RED_SUCCESS)
            {
                min_check = true;
		precision_success = precision_fail;
                cout << endl;
                cout << "Succeeded with precision_fail = " << precision_fail << endl;
                precision_fail -= precision_jump;
		precision_fail = (precision_fail <= 0) ? 1 : precision_fail;
                precision_jump *= PRECISION_JUMP_FACTOR;
                cout << "Checking with precision_fail = " << precision_fail << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcgg_wrapper (A, delta_LC, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcgg_cfa_wrapper (A, delta_LC, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }

            // Find a sufficient precison
            if (!min_check)
            {
                precision_jump = PRECISION_STEP_UP;
                precision_success = precision_fail + precision_jump;

                cout << endl;
                cout << "Checking with precision_success = " << precision_success << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	
		if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcgg_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcgg_cfa_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 

                while (prec_status != RED_SUCCESS)
                {
                    cout << "Failed with precision_success = " << precision_success << endl;
                    precision_fail = precision_success;
                    precision_jump *= PRECISION_JUMP_FACTOR;
                    precision_success += precision_jump;
                    cout << "Checking with precision_success = " << precision_success << endl;

                    status |= read_file(A, input_filename.c_str());
                    if (status) exit (ERROR_FILE_READ);
	            if (orthogonalisation_switch == DEFAULT_GSO)
	            {
                        prec_status = lcgg_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	            } else
	            {
                        prec_status = lcgg_cfa_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	            } 
                }
            }

            if (precision_fail > precision_success)
            {
                cout << "Error: In finding precision range; precision_fail = " << precision_fail << " precision_success = " << precision_success << endl;
                exit (ERROR_INVALID_PRECISION_RANGE);
            }

            // Conduct binary search for the minimum precision
            cout << endl << "Starting binary search with" << endl;
            cout << "Inital precision_fail = " << precision_fail << endl;
            cout << "Inital precision_success = " << precision_success << endl << endl;
            while (precision_fail+1 < precision_success)
            {
                precision_mid = precision_fail + precision_success;
                precision_mid = (precision_mid % 2) ? (precision_mid+1)/2 : precision_mid/2; 
                cout << endl << "precision_fail = " << precision_fail << endl;
                cout << "precision_mid = " << precision_mid << endl;
                cout << "precision_success = " << precision_success << endl << endl;
                cout << "Checking with precision_mid = " << precision_mid << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcgg_wrapper (A, delta_LC, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcgg_cfa_wrapper (A, delta_LC, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
                if (prec_status != RED_SUCCESS)
                    precision_fail = precision_mid;
                else
                    precision_success = precision_mid;
            }

            precision = precision_success;
            cout << "Running with Precision " << precision << endl;

            status |= read_file(A, input_filename.c_str());
            if (status) exit (ERROR_FILE_READ);
	    if (orthogonalisation_switch == DEFAULT_GSO)
	    {
                prec_status = lcgg_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
    	    } else
	    {
                prec_status = lcgg_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	    } 
        }
	// Exhaustive search for a precision which ensures correct reduction
	// Starts at the precision passed as input (default = 53), and increments precision by 1 if incorrect reduction 
        else if (search_switch == SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcgg_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcgg_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }
        }
        else if (search_switch == SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcgg_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcgg_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }
            precision_min_RHF = precision;
            precision_min_runtime = precision;
            RHF_min = RHF;
            RHF_min_runtime = RHF;
            runtime_min = time_taken;
            runtime_min_RHF = time_taken;
            do
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                prev_RHF = RHF;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
                    prec_status = lcgg_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
                    prec_status = lcgg_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 

                if (prec_status == RED_SUCCESS)
                {
                    if (RHF < RHF_min)
                    {
                        RHF_min = RHF;
                        runtime_min_RHF = time_taken;
                        precision_min_RHF = precision;
                    }

                    if (time_taken < runtime_min)
                    {
                        RHF_min_runtime = RHF;
                        runtime_min = time_taken;
                        precision_min_runtime = precision;
                    }
                }
            }
            while (prev_RHF != RHF || prec_status != RED_SUCCESS);

            cout << "Min time:\ttime=" << runtime_min;
            cout << "\tRHF=" << RHF_min_runtime;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_runtime;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

            cout << "Min RHF:\ttime=" << runtime_min_RHF;
            cout << "\tRHF=" << RHF_min;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_RHF;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;
	}
        
        else
        {
            // The default case has been handled already
        }
    }
    // algo_num == 7 - LC-PGG
    else if (algo_num == 7) 
    {
	cout << endl << "------------------------------" << endl << "ALGORITHM: LC-PGG" << endl << "------------------------------";
        cout << endl;
        cout << endl << "Starting with precision = " << precision << endl;
	if (orthogonalisation_switch == DEFAULT_GSO)
	{
	    cout << "Default GSO" << endl << "------------------------------" << endl;
            prec_status = lcpgg_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	} else
	{
	    cout << "CFA and Lazy Size Reduce" << endl << "------------------------------" << endl;
            prec_status = lcpgg_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	} 

        cout << endl;

	// Binary search for minimum precision which ensures correct reduction
        if (search_switch == SEARCH_SWITCH_MIN_PRECISION_BINARY)
        {
            // Find an insufficient precison
            precision_fail = precision;
            precision_jump = PRECISION_STEP_DOWN;
            bool min_check = false;

            while (prec_status == RED_SUCCESS)
            {
                min_check = true;
		precision_success = precision_fail;
                cout << endl;
                cout << "Succeeded with precision_fail = " << precision_fail << endl;
                precision_fail -= precision_jump;
		precision_fail = (precision_fail <= 0) ? 1 : precision_fail;
                precision_jump *= PRECISION_JUMP_FACTOR;
                cout << "Checking with precision_fail = " << precision_fail << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_wrapper (A, delta_LC, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_cfa_wrapper (A, delta_LC, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }

            // Find a sufficient precison
            if (!min_check)
            {
                precision_jump = PRECISION_STEP_UP;
                precision_success = precision_fail + precision_jump;

                cout << endl;
                cout << "Checking with precision_success = " << precision_success << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	
		if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_cfa_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 

                while (prec_status != RED_SUCCESS)
                {
                    cout << "Failed with precision_success = " << precision_success << endl;
                    precision_fail = precision_success;
                    precision_jump *= PRECISION_JUMP_FACTOR;
                    precision_success += precision_jump;
                    cout << "Checking with precision_success = " << precision_success << endl;

                    status |= read_file(A, input_filename.c_str());
                    if (status) exit (ERROR_FILE_READ);
	            if (orthogonalisation_switch == DEFAULT_GSO)
	            {
	                cout << "Default GSO" << endl << "------------------------------" << endl;
                        prec_status = lcpgg_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	            } else
	            {
	                cout << "CFA and LSR" << endl << "------------------------------" << endl;
                        prec_status = lcpgg_cfa_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	            } 
                }
            }

            if (precision_fail > precision_success)
            {
                cout << "Error: In finding precision range; precision_fail = " << precision_fail << " precision_success = " << precision_success << endl;
                exit (ERROR_INVALID_PRECISION_RANGE);
            }

            // Conduct binary search for the minimum precision
            cout << endl << "Starting binary search with" << endl;
            cout << "Inital precision_fail = " << precision_fail << endl;
            cout << "Inital precision_success = " << precision_success << endl << endl;
            while (precision_fail+1 < precision_success)
            {
                precision_mid = precision_fail + precision_success;
                precision_mid = (precision_mid % 2) ? (precision_mid+1)/2 : precision_mid/2; 
                cout << endl << "precision_fail = " << precision_fail << endl;
                cout << "precision_mid = " << precision_mid << endl;
                cout << "precision_success = " << precision_success << endl << endl;
                cout << "Checking with precision_mid = " << precision_mid << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_wrapper (A, delta_LC, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_cfa_wrapper (A, delta_LC, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
                if (prec_status != RED_SUCCESS)
                    precision_fail = precision_mid;
                else
                    precision_success = precision_mid;
            }

            precision = precision_success;
            cout << "Running with Precision " << precision << endl;

            status |= read_file(A, input_filename.c_str());
            if (status) exit (ERROR_FILE_READ);
	    if (orthogonalisation_switch == DEFAULT_GSO)
	    {
	        cout << "Default GSO" << endl << "------------------------------" << endl;
                prec_status = lcpgg_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
    	    } else
	    {
	        cout << "CFA and LSR" << endl << "------------------------------" << endl;
                prec_status = lcpgg_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	    } 
        }
	// Exhaustive search for a precision which ensures correct reduction
	// Starts at the precision passed as input (default = 53), and increments precision by 1 if incorrect reduction 
        else if (search_switch == SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }
        }
        else if (search_switch == SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }
            precision_min_RHF = precision;
            precision_min_runtime = precision;
            RHF_min = RHF;
            RHF_min_runtime = RHF;
            runtime_min = time_taken;
            runtime_min_RHF = time_taken;
            do
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                prev_RHF = RHF;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 

                if (prec_status == RED_SUCCESS)
                {
                    if (RHF < RHF_min)
                    {
                        RHF_min = RHF;
                        runtime_min_RHF = time_taken;
                        precision_min_RHF = precision;
                    }

                    if (time_taken < runtime_min)
                    {
                        RHF_min_runtime = RHF;
                        runtime_min = time_taken;
                        precision_min_runtime = precision;
                    }
                }
            }
            while (prev_RHF != RHF || prec_status != RED_SUCCESS);

            cout << "Min time:\ttime=" << runtime_min;
            cout << "\tRHF=" << RHF_min_runtime;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_runtime;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

            cout << "Min RHF:\ttime=" << runtime_min_RHF;
            cout << "\tRHF=" << RHF_min;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_RHF;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;
	}
        
        else
        {
            // The default case has been handled already
        }
    }
    // algo_num == 8 - LC-PGG (locally restricted insertions)
    else if (algo_num == 8) 
    {
	cout << endl << "------------------------------" << endl << "ALGORITHM: LC-PGG (local restricted insertions)" << endl << "------------------------------";
        cout << endl;
        cout << endl << "Starting with precision = " << precision << endl;
	if (orthogonalisation_switch == DEFAULT_GSO)
	{
	    cout << "Default GSO" << endl << "------------------------------" << endl;
            prec_status = lcpgg_local_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	} else
	{
	    cout << "CFA and Lazy Size Reduce" << endl << "------------------------------" << endl;
            prec_status = lcpgg_local_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	} 

        cout << endl;

	// Binary search for minimum precision which ensures correct reduction
        if (search_switch == SEARCH_SWITCH_MIN_PRECISION_BINARY)
        {
            // Find an insufficient precison
            precision_fail = precision;
            precision_jump = PRECISION_STEP_DOWN;
            bool min_check = false;

            while (prec_status == RED_SUCCESS)
            {
                min_check = true;
		precision_success = precision_fail;
                cout << endl;
                cout << "Succeeded with precision_fail = " << precision_fail << endl;
                precision_fail -= precision_jump;
		precision_fail = (precision_fail <= 0) ? 1 : precision_fail;
                precision_jump *= PRECISION_JUMP_FACTOR;
                cout << "Checking with precision_fail = " << precision_fail << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_wrapper (A, delta_LC, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_cfa_wrapper (A, delta_LC, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }

            // Find a sufficient precison
            if (!min_check)
            {
                precision_jump = PRECISION_STEP_UP;
                precision_success = precision_fail + precision_jump;

                cout << endl;
                cout << "Checking with precision_success = " << precision_success << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	
		if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_cfa_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 

                while (prec_status != RED_SUCCESS)
                {
                    cout << "Failed with precision_success = " << precision_success << endl;
                    precision_fail = precision_success;
                    precision_jump *= PRECISION_JUMP_FACTOR;
                    precision_success += precision_jump;
                    cout << "Checking with precision_success = " << precision_success << endl;

                    status |= read_file(A, input_filename.c_str());
                    if (status) exit (ERROR_FILE_READ);
	            if (orthogonalisation_switch == DEFAULT_GSO)
	            {
	                cout << "Default GSO" << endl << "------------------------------" << endl;
                        prec_status = lcpgg_local_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	            } else
	            {
	                cout << "CFA and LSR" << endl << "------------------------------" << endl;
                        prec_status = lcpgg_local_cfa_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	            } 
                }
            }

            if (precision_fail > precision_success)
            {
                cout << "Error: In finding precision range; precision_fail = " << precision_fail << " precision_success = " << precision_success << endl;
                exit (ERROR_INVALID_PRECISION_RANGE);
            }

            // Conduct binary search for the minimum precision
            cout << endl << "Starting binary search with" << endl;
            cout << "Inital precision_fail = " << precision_fail << endl;
            cout << "Inital precision_success = " << precision_success << endl << endl;
            while (precision_fail+1 < precision_success)
            {
                precision_mid = precision_fail + precision_success;
                precision_mid = (precision_mid % 2) ? (precision_mid+1)/2 : precision_mid/2; 
                cout << endl << "precision_fail = " << precision_fail << endl;
                cout << "precision_mid = " << precision_mid << endl;
                cout << "precision_success = " << precision_success << endl << endl;
                cout << "Checking with precision_mid = " << precision_mid << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_wrapper (A, delta_LC, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_cfa_wrapper (A, delta_LC, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
                if (prec_status != RED_SUCCESS)
                    precision_fail = precision_mid;
                else
                    precision_success = precision_mid;
            }

            precision = precision_success;
            cout << "Running with Precision " << precision << endl;

            status |= read_file(A, input_filename.c_str());
            if (status) exit (ERROR_FILE_READ);
	    if (orthogonalisation_switch == DEFAULT_GSO)
	    {
	        cout << "Default GSO" << endl << "------------------------------" << endl;
                prec_status = lcpgg_local_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
    	    } else
	    {
	        cout << "CFA and LSR" << endl << "------------------------------" << endl;
                prec_status = lcpgg_local_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	    } 
        }
	// Exhaustive search for a precision which ensures correct reduction
	// Starts at the precision passed as input (default = 53), and increments precision by 1 if incorrect reduction 
        else if (search_switch == SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }
        }
        else if (search_switch == SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }
            precision_min_RHF = precision;
            precision_min_runtime = precision;
            RHF_min = RHF;
            RHF_min_runtime = RHF;
            runtime_min = time_taken;
            runtime_min_RHF = time_taken;
            do
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                prev_RHF = RHF;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_local_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 

                if (prec_status == RED_SUCCESS)
                {
                    if (RHF < RHF_min)
                    {
                        RHF_min = RHF;
                        runtime_min_RHF = time_taken;
                        precision_min_RHF = precision;
                    }

                    if (time_taken < runtime_min)
                    {
                        RHF_min_runtime = RHF;
                        runtime_min = time_taken;
                        precision_min_runtime = precision;
                    }
                }
            }
            while (prev_RHF != RHF || prec_status != RED_SUCCESS);

            cout << "Min time:\ttime=" << runtime_min;
            cout << "\tRHF=" << RHF_min_runtime;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_runtime;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

            cout << "Min RHF:\ttime=" << runtime_min_RHF;
            cout << "\tRHF=" << RHF_min;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_RHF;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;
	}
        
        else
        {
            // The default case has been handled already
        }
    }
    // algo_num == 9 - LC-PGG (globally restricted insertions)
    else if (algo_num == 9) 
    {
	cout << endl << "------------------------------" << endl << "ALGORITHM: LC-PGG (global restricted insertions)" << endl << "------------------------------";
        cout << endl;
        cout << endl << "Starting with precision = " << precision << endl;
	if (orthogonalisation_switch == DEFAULT_GSO)
	{
	    cout << "Default GSO" << endl << "------------------------------" << endl;
            prec_status = lcpgg_global_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	} else
	{
	    cout << "CFA and Lazy Size Reduce" << endl << "------------------------------" << endl;
            prec_status = lcpgg_global_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	} 

        cout << endl;

	// Binary search for minimum precision which ensures correct reduction
        if (search_switch == SEARCH_SWITCH_MIN_PRECISION_BINARY)
        {
            // Find an insufficient precison
            precision_fail = precision;
            precision_jump = PRECISION_STEP_DOWN;
            bool min_check = false;

            while (prec_status == RED_SUCCESS)
            {
                min_check = true;
		precision_success = precision_fail;
                cout << endl;
                cout << "Succeeded with precision_fail = " << precision_fail << endl;
                precision_fail -= precision_jump;
		precision_fail = (precision_fail <= 0) ? 1 : precision_fail;
                precision_jump *= PRECISION_JUMP_FACTOR;
                cout << "Checking with precision_fail = " << precision_fail << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_wrapper (A, delta_LC, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_cfa_wrapper (A, delta_LC, d, precision_fail, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }

            // Find a sufficient precison
            if (!min_check)
            {
                precision_jump = PRECISION_STEP_UP;
                precision_success = precision_fail + precision_jump;

                cout << endl;
                cout << "Checking with precision_success = " << precision_success << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	
		if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_cfa_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 

                while (prec_status != RED_SUCCESS)
                {
                    cout << "Failed with precision_success = " << precision_success << endl;
                    precision_fail = precision_success;
                    precision_jump *= PRECISION_JUMP_FACTOR;
                    precision_success += precision_jump;
                    cout << "Checking with precision_success = " << precision_success << endl;

                    status |= read_file(A, input_filename.c_str());
                    if (status) exit (ERROR_FILE_READ);
	            if (orthogonalisation_switch == DEFAULT_GSO)
	            {
	                cout << "Default GSO" << endl << "------------------------------" << endl;
                        prec_status = lcpgg_global_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	            } else
	            {
	                cout << "CFA and LSR" << endl << "------------------------------" << endl;
                        prec_status = lcpgg_global_cfa_wrapper (A, delta_LC, d, precision_success, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	            } 
                }
            }

            if (precision_fail > precision_success)
            {
                cout << "Error: In finding precision range; precision_fail = " << precision_fail << " precision_success = " << precision_success << endl;
                exit (ERROR_INVALID_PRECISION_RANGE);
            }

            // Conduct binary search for the minimum precision
            cout << endl << "Starting binary search with" << endl;
            cout << "Inital precision_fail = " << precision_fail << endl;
            cout << "Inital precision_success = " << precision_success << endl << endl;
            while (precision_fail+1 < precision_success)
            {
                precision_mid = precision_fail + precision_success;
                precision_mid = (precision_mid % 2) ? (precision_mid+1)/2 : precision_mid/2; 
                cout << endl << "precision_fail = " << precision_fail << endl;
                cout << "precision_mid = " << precision_mid << endl;
                cout << "precision_success = " << precision_success << endl << endl;
                cout << "Checking with precision_mid = " << precision_mid << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_wrapper (A, delta_LC, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_cfa_wrapper (A, delta_LC, d, precision_mid, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
                if (prec_status != RED_SUCCESS)
                    precision_fail = precision_mid;
                else
                    precision_success = precision_mid;
            }

            precision = precision_success;
            cout << "Running with Precision " << precision << endl;

            status |= read_file(A, input_filename.c_str());
            if (status) exit (ERROR_FILE_READ);
	    if (orthogonalisation_switch == DEFAULT_GSO)
	    {
	        cout << "Default GSO" << endl << "------------------------------" << endl;
                prec_status = lcpgg_global_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
    	    } else
	    {
	        cout << "CFA and LSR" << endl << "------------------------------" << endl;
                prec_status = lcpgg_global_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	    } 
        }
	// Exhaustive search for a precision which ensures correct reduction
	// Starts at the precision passed as input (default = 53), and increments precision by 1 if incorrect reduction 
        else if (search_switch == SEARCH_SWITCH_MIN_PRECISION_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }
        }
        else if (search_switch == SEARCH_SWITCH_MIN_RHF_RUNTIME_EXHAUSTIVE)
        {
            while (prec_status != RED_SUCCESS)
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 
            }
            precision_min_RHF = precision;
            precision_min_runtime = precision;
            RHF_min = RHF;
            RHF_min_runtime = RHF;
            runtime_min = time_taken;
            runtime_min_RHF = time_taken;
            do
            {
                precision++;
                cout << "Running with Precision " << precision << endl;

                prev_RHF = RHF;

                status |= read_file(A, input_filename.c_str());
                if (status) exit (ERROR_FILE_READ);
	        if (orthogonalisation_switch == DEFAULT_GSO)
	        {
	            cout << "Default GSO" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } else
	        {
	            cout << "CFA and LSR" << endl << "------------------------------" << endl;
                    prec_status = lcpgg_global_cfa_wrapper (A, delta_LC, d, precision, precision_correction_loops_allowed, eta, time_taken, RHF, last_GSO_length, 0, d-1, sublattice_size, beta);
	        } 

                if (prec_status == RED_SUCCESS)
                {
                    if (RHF < RHF_min)
                    {
                        RHF_min = RHF;
                        runtime_min_RHF = time_taken;
                        precision_min_RHF = precision;
                    }

                    if (time_taken < runtime_min)
                    {
                        RHF_min_runtime = RHF;
                        runtime_min = time_taken;
                        precision_min_runtime = precision;
                    }
                }
            }
            while (prev_RHF != RHF || prec_status != RED_SUCCESS);

            cout << "Min time:\ttime=" << runtime_min;
            cout << "\tRHF=" << RHF_min_runtime;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_runtime;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;

            cout << "Min RHF:\ttime=" << runtime_min_RHF;
            cout << "\tRHF=" << RHF_min;
            cout << "\td=" << d;
            cout << "\talgo_num=" << algo_num;
            cout << "\tprec=" << precision_min_RHF;
            cout << "\tsublattice=" << sublattice_size;
            cout << "\tswitch=" << search_switch;
            cout << "\tloops=" << precision_correction_loops_allowed;
            cout << "\teta=" << eta << endl;
	}
        
        else
        {
            // The default case has been handled already
        }
    }
    // Computes RHF for input (LLL-reduced) basis
    else if (algo_num == 10) 
    {
        FP_NR<mpfr_t>::set_prec(precision);
        FP_mat<mpfr_t> gso(d,d);
        FP_mat<mpfr_t> mu(d,d);
        vector<FP_NR<mpfr_t>> B;
        B.resize(d,0.0);
        compute_GSO (A, gso, mu, B, d);
        RHF = compute_RHF (B, d);
	cout << "RHF on Input: " << RHF << endl;
	return 0;
    }
    else 
    {
        cout << "Algorithm number should be in [1,8]" << endl;
        cout << "Algorithm 1 :  SS-DeepLLL" << endl;
        cout << "Algorithm 2 :  SS-GGLLL" << endl;
        cout << "Algorithm 3 :  Pot-DeepLLL" << endl;
        cout << "Algorithm 4 :  Pot-GGLLL" << endl;
        cout << "Algorithm 5 :  LC-DeepLLL" << endl;
        cout << "Algorithm 6 :  LC-GG" << endl;
        cout << "Algorithm 7 :  LC-PGG" << endl;
        cout << "Algorithm 8 :  LC-PGG (locally restricted insertions)" << endl;
        cout << "Algorithm 9 :  LC-PGG (globally restricted insertions)" << endl;
        cout << "Option 10    :  No Reduction" << endl;
        exit (ERROR_ALGORITHM_NUMBER);
    }

    /*
     * Comment/Uncomment line below to toggle output basis printing
    */	
    //cout << A << endl;

    cout << endl;
    cout << "Final   :\ttime=" << time_taken;
    cout << "\tRHF=" << RHF;
    cout << "\tFinalGSO=" << last_GSO_length;
    cout << "\td=" << d;
    cout << "\talgo_num=" << algo_num;
    cout << "\tprec=" << precision;
    cout << "\tsublattice=" << sublattice_size;
    cout << "\tswitch=" << search_switch;
    cout << "\tloops=" << precision_correction_loops_allowed;
    cout << "\teta=" << eta << endl;

    return 0;
}
