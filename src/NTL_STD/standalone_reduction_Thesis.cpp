/*************************************************
 * Program:     standalone_reduction_Thesis
*************************************************/
#include "commonThesis.h"

using namespace std;

#define EXECUTABLE_NAME "standalone_reduction_Thesis"

#define OUTPUT_FILE_NAME "_out_"

// The main function
int main (int argc, char * argv[]) {

	// argv[0]: executable name
	// argv[1]: m
	// argv[2]: start basis number
	// argv[3]: end basis number
	// argv[4]: Algorithm

	if (argc != 5) {
		cout << "Usage: <executable> <m> <start_basis_#> <end_basis_#> <algorithm>" << endl;
		cout << "1: " << _LLL_ << endl;
		cout << "2: " << _SS_LLL_ << endl;
		cout << "3: " << _SS_GGLLL_ << endl;
		cout << "4: " << _Pot_LLL_ << endl;
		cout << "5: " << _Pot_GGLLL_ << endl;
		cout << "6: " << _LC_DeepLLL_ << endl;
		cout << "7: " << _DeepLLL_5_ << endl;
		cout << "8: " << _DeepLLL_10_ << endl;
		cout << "9: " << _DeepLLL_15_ << endl;
		cout << "10: " << _DeepLLL_20_ << endl;
		cout << "11: " << _LC_GGLLL_5_ << endl;
		cout << "12: " << _LC_GGLLL_10_ << endl;
		cout << "13 " << _LC_GGLLL_15_ << endl;
		cout << "14: " << _LC_GGLLL_20_ << endl;
		cout << "15: " << _Pot_F_LC_GGLLL_ << endl;
		cout << "16: " << _SS_F_LC_GGLLL_ << endl;

		exit(-1);
	}
	
	int m = atoi(argv[1]);
	int n = m;
	int basis_number_start = atoi(argv[2]);
	int basis_number_end = atoi(argv[3]);
	int algo = atoi(argv[4]);

        int nIterations = basis_number_end - basis_number_start + 1;
	RR RRIterations; 
	conv(RRIterations, nIterations);

        mat_ZZ ppOriginalBasis;
	ppOriginalBasis.SetDims(m,m);
	ZZ volume;

	char filename[20000];
	char dir[10000];
	struct stat stats;

	long double delta_int = 0.99;
	long double delta_threshold_int = 0.99;
   	RR delta_threshold, delta;
   	delta_threshold = 0.99;
   	delta = 0.99;
        int a = 99;
	int b = 100;

   	// The value of eta used in S^2: 1 - 10^(-6)
	long double eta_int = 0.999999;
   	RR eta;
   	eta = 0.999999;

   	long long int totalSwaps = 0;
 	long long int totalReductions = 0;
 	long double totalTime = 0.0;
  	RR totalV1Norm; totalV1Norm = 0;
	RR totalSVPFactor; totalSVPFactor = 0;
	RR totalRHF; totalRHF = 0;
 	RR totalOD; totalOD = 0;
 	RR totalPot; totalPot = 0;
   	RR totalSS; totalSS = 0;
   	RR V1NormInput; V1NormInput = 0; 
   	RR V1NormOutput; V1NormOutput = 0;
	RR RRSVPConstant; RRSVPConstant = 0;
	RR SVPInputFactor; SVPInputFactor = 0;
	RR SVPOutputFactor; SVPOutputFactor = 0;
   	RR RHFInput; RHFInput = 0;
   	RR RHFOutput; RHFOutput = 0;
   	RR ODInput; ODInput = 0;
   	RR ODOutput; ODOutput = 0;
   	RR PotInput; PotInput = 0;
   	RR PotOutput; PotOutput = 0;
   	RR SSInput; SSInput = 0;
   	RR SSOutput; SSOutput = 0;

      	long long int nNonReducedPairs = 0;
      	long long int nMinDeltaReductions = 0;
        
	setPrecision(m);

	sprintf (dir,"%s/%s",_PARENT_INPUT_DIR_, _INPUT_DIR_);
	cout << dir << endl;
	stat(dir, &stats);
	if (!S_ISDIR(stats.st_mode)){
		cout << "Input directory doesn't exist." << endl;
		exit(4);
	}

	switch (algo) {
		case 1:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _LLL_);
			break;
		case 2:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _SS_LLL_);
			break;
		case 3:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _SS_GGLLL_);
			break;
		case 4:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _Pot_LLL_);
			break;
		case 5:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _Pot_GGLLL_);
			break;
		case 6:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _DeepLLL_);
			break;
		case 7: 
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _DeepLLL_5_);
			break;
		case 8: 
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _DeepLLL_10_);
			break;
		case 9: 
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _DeepLLL_15_);
			break;
		case 10: 
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _DeepLLL_20_);
			break;
		case 11: 
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _LC_GGLLL_5_);
			break;
		case 12: 
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _LC_GGLLL_10_);
			break;
		case 13: 
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _LC_GGLLL_15_);
			break;
		case 14: 
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _LC_GGLLL_20_);
			break;
		case 15: 
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _Pot_F_LC_GGLLL_);
			break;
		case 16: 
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_, _SS_F_LC_GGLLL_);
			break;
		default:
			cout << "No such basis reduction algorithm." << endl;
			exit(algo);
	}

	stat(dir, &stats);
	if (!S_ISDIR(stats.st_mode)){
		if (!mkdir (dir, S_IRWXO|S_IRWXU|S_IRGRP)) {
			cout << "Directory created successfully." << endl;
		} else {
			cout << "Directory could not be created." << endl;
			exit(5);
		}
	} else {
		cout << "Directory already exists." << endl;
	}

	for (int i=basis_number_start; i<=basis_number_end; i++) {
		
		sprintf (filename, "%s/%s/%s_%d_%06d.txt", _PARENT_INPUT_DIR_ , _INPUT_DIR_, _INPUT_FILE_PREFIX_, m, i);
		cout << filename << endl;

		ifstream inFile (filename);
		inFile >> m; 

		for (int i=0; i<m; i++) {
			for (int j=0; j<m; j++) {
	      			inFile >> ppOriginalBasis[i][j]; 
			}
		}
		
		inFile.close();

		/*
		 * Comment/Uncomment line below to toggle input basis printing
		*/	
		//cout << ppOriginalBasis << endl;

		volume = ppOriginalBasis[0][0];
		RRSVPConstant = SVP_Challenge_Factor(n, volume);

      		long long int nSwapCount = 0;
      		long long int nReductionCount = 0;
      		long long int nMinDeltaReductions = 0;
      		long double time_taken = 0.0;
      		V1NormInput = NTL_vector_norm (ppOriginalBasis[0], m);
   		SVPInputFactor = V1NormInput / RRSVPConstant;
      		RHFInput = NTLrootHermiteFactor (ppOriginalBasis, n, m);
      		ODInput = NTLorthogonalityDefect (ppOriginalBasis, n, m);
      		PotInput = NTLbasisPotential (ppOriginalBasis, n, m);
      		SSInput = NTLsquaredSum (ppOriginalBasis, n, m);

      		clock_t begin_time;
      		clock_t end_time;
		switch (algo) {
			case 1: // _LLL_
      				begin_time = clock();
      				ppOriginalBasis = LLL_mat_ZZ (ppOriginalBasis, n, m, delta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 2: // _SS_LLL_
      				begin_time = clock();
      				ppOriginalBasis = NTL_SS_LLL (ppOriginalBasis, m, n, eta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 3: // _SS_GG_
      				begin_time = clock();
      				ppOriginalBasis = NTL_SS_GG (ppOriginalBasis, m, n, eta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 4: // _Pot_LLL_
      				begin_time = clock();
      				ppOriginalBasis = NTL_PotLLL (ppOriginalBasis, m, n, delta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 5: // _Pot_GG_
      				begin_time = clock();
      				ppOriginalBasis = NTL_Pot_GG (ppOriginalBasis, m, n, delta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 6: // _LC_DeepLLL_
      				begin_time = clock();
      				ppOriginalBasis = LC_DeepLLL_NTL (ppOriginalBasis, m, n, delta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
		        case 7: // _LC_DeepLLL_5_
				begin_time = clock();
				ppOriginalBasis = NTL_blockwise_deepLLL(ppOriginalBasis, m, n, delta_threshold, 5, &nSwapCount, &nReductionCount);  
				end_time = clock();
                                break;
		        case 8: // _LC_DeepLLL_10_
				begin_time = clock();
				ppOriginalBasis = NTL_blockwise_deepLLL(ppOriginalBasis, m, n, delta_threshold, 10, &nSwapCount, &nReductionCount);  
				end_time = clock();
                                break;
		        case 9: // _LC_DeepLLL_15_
				begin_time = clock();
				ppOriginalBasis = NTL_blockwise_deepLLL(ppOriginalBasis, m, n, delta_threshold, 15, &nSwapCount, &nReductionCount);  
				end_time = clock();
                                break;
		        case 10: // _LC_DeepLLL_20_
				begin_time = clock();
				ppOriginalBasis = NTL_blockwise_deepLLL(ppOriginalBasis, m, n, delta_threshold, 20, &nSwapCount, &nReductionCount);  
				end_time = clock();
                                break;
			case 11: // _LC_GG_5_
				begin_time = clock();
				ppOriginalBasis = NTL_blockwise_LC_GG (ppOriginalBasis, m, n, delta_threshold, 5, &nSwapCount, &nReductionCount);
                                end_time = clock();
				break;
			case 12: // _LC_GG_10_
				begin_time = clock();
				ppOriginalBasis = NTL_blockwise_LC_GG (ppOriginalBasis, m, n, delta_threshold, 10, &nSwapCount, &nReductionCount);
                                end_time = clock();
				break;
			case 13: // _LC_GG_15_
				begin_time = clock();
				ppOriginalBasis = NTL_blockwise_LC_GG (ppOriginalBasis, m, n, delta_threshold, 15, &nSwapCount, &nReductionCount);
                                end_time = clock();
				break;
			case 14: // _LC_GG_20_
				begin_time = clock();
				ppOriginalBasis = NTL_blockwise_LC_GG (ppOriginalBasis, m, n, delta_threshold, 20, &nSwapCount, &nReductionCount);
                                end_time = clock();
				break;
			case 15: // _Pot_F_LC_GG_
				begin_time = clock();
                                ppOriginalBasis =  pot_Filtered_LC_GG (ppOriginalBasis, m, n, delta_threshold, &nSwapCount, &nReductionCount);
                                end_time = clock();
				break;
			case 16: // _SS_F_LC_GG_
				begin_time = clock();
                                ppOriginalBasis =  SS_Filtered_LC_GG (ppOriginalBasis, m, n, delta_threshold, eta, &nSwapCount, &nReductionCount);
                                end_time = clock();
				break;
			default:
				cout << "No such basis reduction algorithm." << endl;
				exit(algo);
		}
      		
		time_taken = (long double) (end_time - begin_time) / CLOCKS_PER_SEC;

		/*
		 * Uncomment line below if output basis should be printed
		*/	
		//cout << ppOriginalBasis << endl;

		totalSwaps += nSwapCount;
      		totalReductions += nReductionCount;
      		totalTime += time_taken;
      		V1NormOutput = NTL_vector_norm (ppOriginalBasis[0], m);
      		totalV1Norm += V1NormOutput;
   		SVPOutputFactor = V1NormOutput / RRSVPConstant;
		totalSVPFactor += SVPOutputFactor;
      		RHFOutput = NTLrootHermiteFactor (ppOriginalBasis, m, n);
      		totalRHF += RHFOutput;
      		ODOutput = NTLorthogonalityDefect (ppOriginalBasis, m, n);
      		totalOD += ODOutput;
      		PotOutput = NTLbasisPotential (ppOriginalBasis, m, n);
      		totalPot += PotOutput;
      		SSOutput = NTLsquaredSum (ppOriginalBasis, m, n);
      		totalSS += SSOutput;
	
		sprintf (filename, "%s/%s/%s_%d_%06d.txt", _PARENT_OUTPUT_DIR_, dir, _OUTPUT_FILE_PREFIX_, m, i);
		cout << filename << endl;

		ofstream outFile (filename);

		outFile << nSwapCount << "\t";
 		outFile << nReductionCount << "\t";
   		outFile << time_taken << "\t";
   		outFile << V1NormOutput << "\t";
   		outFile << SVPOutputFactor << "\t";
        	outFile << RHFOutput << "\t";
   		outFile << ODOutput << "\t";
   		outFile << PotOutput << "\t";
   		outFile << SSOutput << "\n";

		outFile.close();


		cout << "\n\nNumber of deep insertions: " << nSwapCount << "\n";
 		cout << "Number of size reductions: " <<nReductionCount << "\n";
   		cout << "Time taken: " << time_taken << "\n";
   		cout << "Length of first vector: " << V1NormOutput << "\n";
   		cout << "SVP Challenge Factor: " << SVPOutputFactor << "\n";
        cout << "RHF: " << RHFOutput << "\n";
   		cout << "Orthogonality Defect: " << ODOutput << "\n";
   		cout << "Log potential: " << PotOutput << "\n";
   		cout << "Squared sum: " << SSOutput << "\n";
	}

	ppOriginalBasis.kill();

	exit(0);
}
