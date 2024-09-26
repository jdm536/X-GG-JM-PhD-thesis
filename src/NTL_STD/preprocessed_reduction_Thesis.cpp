/*************************************************
 * Program:     preprocessed_reduction_Thesis
*************************************************/
#include "commonThesis.h"

using namespace std;

#define EXECUTABLE_NAME "preprocessed_reduction_Thesis"

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
		cout << "1: " << _SS_LLL_ << endl;
		cout << "2: " << _SS_GGLLL_ << endl;
		cout << "3: " << _Pot_LLL_ << endl;
		cout << "4: " << _Pot_GGLLL_ << endl;
		cout << "5: " << _LC_GGLLL_ << endl;
		cout << "6: " << _LC_GGLLL_5_ << endl;
		cout << "7: " << _LC_GGLLL_10_ << endl;
		cout << "8: " << _LC_GGLLL_15_ << endl;
		cout << "9: " << _LC_GGLLL_20_ << endl;
		cout << "10: " << _LC_PGGLLL_ << endl;
		cout << "11: " << _LC_PGGLLL_5G_ << endl;
		cout << "12: " << _LC_PGGLLL_10G_ << endl;
		cout << "13: " << _LC_PGGLLL_15G_ << endl;
		cout << "14: " << _LC_PGGLLL_20G_ << endl;
		cout << "15: " << _LC_PGGLLL_5L_ << endl;
		cout << "16: " << _LC_PGGLLL_10L_ << endl;
		cout << "17: " << _LC_PGGLLL_15L_ << endl;
		cout << "18: " << _LC_PGGLLL_20L_ << endl;
		cout << "19: " << _LC_DeepLLL_ << endl;
		cout << "20: " << _DeepLLL_5_ << endl;
		cout << "21: " << _DeepLLL_10_ << endl;
		cout << "22: " << _DeepLLL_15_ << endl;
		cout << "23: " << _DeepLLL_20_ << endl;
		cout << "24: " << _LLL_ << endl << endl;
		cout << "Opt. 24 computes RHF, V1 Length etc. for the FPLLL reduced basis which would be used as input to other reduction algorithms" << endl;
		cout << "So no reduction takes place for Opt. 20." << endl;

		exit(-1);
	}
	
	int m = atoi(argv[1]);
	int n = m;
	int basis_number_start = atoi(argv[2]);
	int basis_number_end = atoi(argv[3]);
	int algo = atoi(argv[4]);
	int max_loop = 5;

        int nIterations = basis_number_end - basis_number_start + 1;
	RR RRIterations; 
	conv(RRIterations, nIterations);

	int ** ppBasisInt;
   	ppBasisInt = (int **) calloc (n, sizeof(int *));
   	for (int i=0; i<n; i++) {
      		ppBasisInt[i] = (int *) calloc (m, sizeof(int));
	}

        mat_ZZ ppBasisZZ;
	ppBasisZZ.SetDims(m,m);
	ZZ volume;

	char filename[20000];
	char dir[10000];
	struct stat stats;

	long double delta_threshold = 0.99;

   	// The value of eta used in S^2: 1 - 10^(-6)
	long double eta = 0.999999;

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
        
	sprintf (dir,"%s/%s",_PARENT_INPUT_DIR_, _INPUT_DIR_FPLLL_RED_);
	cout << dir << endl;
	stat(dir, &stats);
	if (!S_ISDIR(stats.st_mode)){
		cout << "Input directory doesn't exist." << endl;
		exit(4);
	}

	switch (algo) {
		case 1:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _SS_LLL_);
			break;
		case 2:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _SS_GGLLL_);
			break;
		case 3:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _Pot_LLL_);
			break;
		case 4:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _Pot_GGLLL_);
			break;
		case 5:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_GGLLL_);
			break;
		case 6:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_GGLLL_5_);
			break;
		case 7:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_GGLLL_10_);
			break;
		case 8:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_GGLLL_15_);
			break;
		case 9:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_GGLLL_20_);
			break;
		case 10:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_PGGLLL_);
			break;
		case 11:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_PGGLLL_5G_);
			break;
		case 12:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_PGGLLL_10G_);
			break;
		case 13:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_PGGLLL_15G_);
			break;
		case 14:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_PGGLLL_20G_);
			break;
		case 15:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_PGGLLL_5L_);
			break;
		case 16:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_PGGLLL_10L_);
			break;
		case 17:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_PGGLLL_15L_);
			break;
		case 18:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_PGGLLL_20L_);
			break;
		case 19:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LC_DeepLLL_);
			break;
		case 20:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _DeepLLL_5_);
			break;
		case 21:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _DeepLLL_10_);
			break;
		case 22:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _DeepLLL_15_);
			break;
		case 23:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _DeepLLL_20_);
			break;
		case 24:
			sprintf (dir, "%s/%s",_PARENT_OUTPUT_DIR_PREPROCESSED_, _LLL_);
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
		
		sprintf (filename, "%s/%s/%s_%d_%06d.txt", _PARENT_INPUT_DIR_ , _INPUT_DIR_FPLLL_RED_, _INPUT_FILE_PREFIX_, m, i);

		char c;
		ifstream inFile (filename);
		inFile >> c;
		for (int i=0; i<m; i++) {
			inFile >> c;
			for (int j=0; j<m; j++) {
				inFile >> ppBasisInt[i][j];
			}
			inFile >> c;
		}
		inFile >> c;
		inFile.close();	
		
		ppBasisZZ = copyBasisToNTL (ppBasisInt, m, n);

		/*
		 * Comment/Uncomment line below to toggle input basis printing
		*/	
		//cout << ppBasisZZ << endl;

		volume = abs(determinant(ppBasisZZ, 1));
		RRSVPConstant = SVP_Challenge_Factor(n, volume);

      		long long int nSwapCount = 0;
      		long long int nReductionCount = 0;
      		long long int nGSORecomputations = 0;
      		long long int nIterations = 0;
      		long long int nMinDeltaReductions = 0;
      		long double time_taken = 0.0;
      		V1NormInput = NTL_vector_norm (ppBasisZZ[0], m);
   		SVPInputFactor = V1NormInput / RRSVPConstant;
      		RHFInput = NTLrootHermiteFactor (ppBasisZZ, n, m);
      		ODInput = NTLorthogonalityDefect (ppBasisZZ, n, m);
      		PotInput = NTLbasisPotential (ppBasisZZ, n, m);
      		SSInput = NTLsquaredSum (ppBasisZZ, n, m);

      		clock_t begin_time;
      		clock_t end_time;
		switch (algo) {
			case 1: // _SS_LLL_
      				begin_time = clock();
      				ppBasisInt = SS_LLL_std (ppBasisInt, m, n, eta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 2: // _SS_GG_
      				begin_time = clock();
      				ppBasisInt = SS_GG_std (ppBasisInt, m, n, eta, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 3: // _Pot_LLL_
      				begin_time = clock();
      				ppBasisInt = Pot_LLL_std (ppBasisInt, m, n, delta_threshold, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 4: // _Pot_GG_
      				begin_time = clock();
      				ppBasisInt = Pot_GG_std (ppBasisInt, m, n, delta_threshold, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 5: // _LC_GG_
      				begin_time = clock();
      				ppBasisInt = LC_GG_std (ppBasisInt, m, n, delta_threshold, &nSwapCount, &nReductionCount, &nGSORecomputations);
      				end_time = clock();
				break;
			case 6: // _LC_GG_5_
      				begin_time = clock();
      				ppBasisInt = block_LC_GG_std (ppBasisInt, m, n, delta_threshold, 5, &nSwapCount, &nReductionCount, &nGSORecomputations);
      				end_time = clock();
				break;
			case 7: // _LC_GG_10_
      				begin_time = clock();
      				ppBasisInt = block_LC_GG_std (ppBasisInt, m, n, delta_threshold, 10, &nSwapCount, &nReductionCount, &nGSORecomputations);
      				end_time = clock();
				break;
			case 8: // _LC_GG_15_
      				begin_time = clock();
      				ppBasisInt = block_LC_GG_std (ppBasisInt, m, n, delta_threshold, 15, &nSwapCount, &nReductionCount, &nGSORecomputations);
      				end_time = clock();
				break;
			case 9: // _LC_GG_20_
      				begin_time = clock();
      				ppBasisInt = block_LC_GG_std (ppBasisInt, m, n, delta_threshold, 20,  &nSwapCount, &nReductionCount, &nGSORecomputations);
      				end_time = clock();
				break;
			case 10: // _LC_PGG_
      				begin_time = clock();
				ppBasisInt = LC_PGG(ppBasisInt, m, n, delta_threshold, &nSwapCount, &nReductionCount, &nIterations); 
      				end_time = clock();
				break;
			case 11: // _LC_PGG_5G_
      				begin_time = clock();
      				ppBasisInt = global_block_LC_PGG (ppBasisInt, m, n, delta_threshold, 5, &nSwapCount, &nReductionCount, &nIterations);
      				end_time = clock();
				break;
			case 12: // _LC_PGG_10G_
      				begin_time = clock();
      				ppBasisInt = global_block_LC_PGG (ppBasisInt, m, n, delta_threshold, 10, &nSwapCount, &nReductionCount, &nIterations);
      				end_time = clock();
				break;
			case 13: // _LC_PGG_15G_
      				begin_time = clock();
      				ppBasisInt = global_block_LC_PGG (ppBasisInt, m, n, delta_threshold, 15, &nSwapCount, &nReductionCount, &nIterations);
      				end_time = clock();
				break;
			case 14: // _LC_PGG_20G_
      				begin_time = clock();
      				ppBasisInt = global_block_LC_PGG (ppBasisInt, m, n, delta_threshold, 20,  &nSwapCount, &nReductionCount, &nIterations);
      				end_time = clock();
				break;
			case 15: // _LC_PGG_5L_
      				begin_time = clock();
      				ppBasisInt = local_block_LC_PGG (ppBasisInt, m, n, delta_threshold, 5, &nSwapCount, &nReductionCount, &nIterations);
      				end_time = clock();
				break;
			case 16: // _LC_PGG_10L_
      				begin_time = clock();
      				ppBasisInt = local_block_LC_PGG (ppBasisInt, m, n, delta_threshold, 10, &nSwapCount, &nReductionCount, &nIterations);
      				end_time = clock();
				break;
			case 17: // _LC_PGG_15L_
      				begin_time = clock();
      				ppBasisInt = local_block_LC_PGG (ppBasisInt, m, n, delta_threshold, 15, &nSwapCount, &nReductionCount, &nIterations);
      				end_time = clock();
				break;
			case 18: // _LC_PGG_20L_
      				begin_time = clock();
      				ppBasisInt = local_block_LC_PGG (ppBasisInt, m, n, delta_threshold, 20,  &nSwapCount, &nReductionCount, &nIterations);
      				end_time = clock();
				break;
			case 19: // _LC_DeepLLL_
      				begin_time = clock();
      				ppBasisInt = LC_DeepLLL_std (ppBasisInt, m, n, delta_threshold, n, max_loop, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 20: // _LC_DeepLLL_5_
      				begin_time = clock();
      				ppBasisInt = LC_DeepLLL_std (ppBasisInt, m, n, delta_threshold, 5, max_loop, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 21: // _LC_DeepLLL_10_
      				begin_time = clock();
      				ppBasisInt = LC_DeepLLL_std (ppBasisInt, m, n, delta_threshold, 10, max_loop, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 22: // _LC_DeepLLL_15_
      				begin_time = clock();
      				ppBasisInt = LC_DeepLLL_std (ppBasisInt, m, n, delta_threshold, 15, max_loop, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 23: // _LC_DeepLLL_20_
      				begin_time = clock();
      				ppBasisInt = LC_DeepLLL_std (ppBasisInt, m, n, delta_threshold, 20, max_loop, &nSwapCount, &nReductionCount);
      				end_time = clock();
				break;
			case 24: // _LLL_ (for basis quantities only, not runtime, swaps etc.)
      				begin_time = clock();
      				end_time = clock();
				break;

			default:
				cout << "No such basis reduction algorithm." << endl;
				exit(algo);
		}
      		
		time_taken = (long double) (end_time - begin_time) / CLOCKS_PER_SEC;
		
		ppBasisZZ = copyBasisToNTL (ppBasisInt, m, n);

		/*
		 * Comment/Uncomment line below to toggle output basis printing
		*/	
		//cout << ppBasisZZ << endl;

		totalSwaps += nSwapCount;
      		totalReductions += nReductionCount;
      		totalTime += time_taken;

      		V1NormOutput = NTL_vector_norm (ppBasisZZ[0], m);
      		totalV1Norm += V1NormOutput;
   		SVPOutputFactor = V1NormOutput / RRSVPConstant;
		totalSVPFactor += SVPOutputFactor;
      		RHFOutput = NTLrootHermiteFactor (ppBasisZZ, m, n);
      		totalRHF += RHFOutput;
      		ODOutput = NTLorthogonalityDefect (ppBasisZZ, m, n);
      		totalOD += ODOutput;
      		PotOutput = NTLbasisPotential (ppBasisZZ, m, n);
      		totalPot += PotOutput;
      		SSOutput = NTLsquaredSum (ppBasisZZ, m, n);
      		totalSS += SSOutput;
	
		sprintf (filename, "%s/%s/%s_%d_%06d.txt", _PARENT_OUTPUT_DIR_, dir, _OUTPUT_FILE_PREFIX_, m, i);

		ofstream outFile (filename);

		outFile << nSwapCount << "\t";
 		outFile << nReductionCount << "\t";
   		outFile << time_taken << "\t";
   		outFile << V1NormOutput << "\t";
   		outFile << SVPOutputFactor << "\t";
        	outFile << RHFOutput << "\t";
   		outFile << ODOutput << "\t";
   		outFile << PotOutput << "\t";
   		outFile << SSOutput << "\t";
   		outFile << nGSORecomputations << "\t";
   		outFile << nIterations << "\n";

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
   		cout << "Number of GSO recomputations: " << nGSORecomputations << "\n";
   		cout << "Number of iterations: " << nIterations << "\n";
	}

	ppBasisZZ.kill();
	deleteBasis (ppBasisInt, n);

	exit(0);
}
