#include "read.h"

// Read the phenotype and covariate information

void readpheno(char* filename, struct FAMILY *family, struct HASH *hash, struct DATA_STRUCT data_struct)  {

	int i, j, i_total = 0, fam = 1, cov_check, ret;
	int n_cov = data_struct.n_cov;
	int n_fam = data_struct.n_fam;
	int n_total = data_struct.n_total;
	double pheno_tmp, cov_tmp;
	char ind_id[MAXINDLEN];

	// Initialize hash table
	if (hash->ind2ind != NULL) kh_destroy(str, hash->ind2ind);
	hash->ind2ind = kh_init(str);
	khiter_t k_i;

	// Initialize arrays
	(family+fam)->n_total = n_total;
	(family+fam)->trait_obs = dvector(1, n_total);
	(family+fam)->Z = dvector(1, n_total);
	(family+fam)->cov = dmatrix(1, n_total, 1, n_cov);
	// Store the subject IDs
	(family+fam)->indiv = (char**)malloc((n_total+1) * sizeof(char*));
	for (i = 1; i <= n_total; i++) (family+fam)->indiv[i] = (char *)malloc(MAXINDLEN * sizeof(char));

	// Read and store phenotype and covariates
	FILE *pedfile;
	pedfile = fopen(filename, "r");
	if (pedfile == NULL) {
		printf("ERROR: Can't open phenotype data file!\n");
		exit(1);
	}

	while (!feof(pedfile)){
		if(n_cov > 1) fscanf(pedfile, "%*d %s %*s %*s %*d %lf ", ind_id, &pheno_tmp);
		if (n_cov == 1) fscanf(pedfile, "%*d %s %*s %*s %*d %lf\n", ind_id, &pheno_tmp);

		// Fill hash table
		k_i = kh_put(str,hash->ind2ind, strdup(ind_id), &ret);
		kh_value(hash->ind2ind,k_i) = ++i_total;

		// Store pheno
		(family+fam)->trait_obs[i_total] = pheno_tmp;
		(family+fam)->cov[i_total][1] = 1; // Intercept

		// Store covariates
		if (n_cov > 1) {
			for (i = 2; i < n_cov; i++)
				fscanf(pedfile, "%lf ", &(family + fam)->cov[i_total][i]);
			fscanf(pedfile, "%lf\n", &(family + fam)->cov[i_total][i]);
		}

		// Store ID
		strcpy((family+fam)->indiv[i_total], ind_id);
	}

	// Check that none of covariates are constant
	for(i = 2; i <= n_cov; i++){
		for(cov_check = 0, j = 2; j <= (family+fam)->n_total; j++) {
			if ((family + fam)->cov[j][i] != (family + fam)->cov[1][i]) {
				cov_check = 1;
				break;
			}
		}
		if(cov_check == 0){
			printf("ERROR: Phenotype file has at least one constant covariate!\n");
			exit(1);
		}
	}

	fclose(pedfile);
}

/* Read in the eigen values and eigen vectors (in binary mode)
  Assumes order of individuals is same as in phenotype file
*/
void readeigdecomp(char *filename, struct FAMILY *family) {

	int i, j, fam = 1;
	int n_total = (family + fam)->n_total;

	FILE *eigfile;
	eigfile = fopen(filename, "rb");
	if(eigfile == NULL){
		printf("ERROR: Can't open eigen-decomposition file: %s!\n", filename);
		exit(1);
	}

	// Eigen vectors and values for covariance matrix used in CARAT
	(family + fam)->eigvector = dmatrix(1, n_total, 1, n_total);
	(family + fam)->eigvalue = dvector(1, n_total);

	// Read in decomposition results
	for(i = 1; i <= n_total; i++){
		fread(&(family + fam)->eigvalue[i], sizeof(double), 1, eigfile);
	}
	for(i = 1; i <= n_total; i++){
		for(j = 1; j <= n_total; j++){
			fread(&(family + fam)->eigvector[i][j], sizeof(double), 1, eigfile);
		}
	}
	fclose(eigfile);
}


/* Read in the GRM
Assumes that pairs not included have GRM entry of 0 (or 1 if on diagonal)
*/
void readkin(char *filename, struct FAMILY *family, struct HASH hash) {

	int ind1, ind2, fam = 1, fam_tmp, k_i;
	int n_total = (family + fam)->n_total;
	double kcoef;
	char ind1_id[MAXINDLEN], ind2_id[MAXINDLEN];


	FILE *kinfile;
	kinfile = fopen(filename, "r");
	if (kinfile == NULL) {
		printf("ERROR: Can't open GRM file: %s!\n", filename);
		exit(1);
	}

	// Initialize GRM phi
	(family + fam)->phi = dmatrix(1, n_total, 1, n_total);
	for (ind1 = 1; ind1 <= n_total ; ++ind1) {
		for (ind2 = 1; ind2 < ind1 ; ++ind2)
			(family + fam)->phi[ind1][ind2] = (family + fam)->phi[ind2][ind1] = 0;
		(family + fam)->phi[ind1][ind1] = 1;
	}

	// Store GRM
	while (!feof(kinfile)) {
		if (fscanf(kinfile, "%d %s %s %lf\n", &fam_tmp, ind1_id, ind2_id, &kcoef) != 4) exit(1);
		if(fam_tmp != 1){
			printf("ERROR: Family is not 1 in GRM file!\n", filename);
			exit(1);
		}

		// Get indices from hash table
		k_i = kh_get(str, hash.ind2ind, ind1_id);
		ind1 = kh_value(hash.ind2ind, k_i);

		k_i = kh_get(str, hash.ind2ind, ind2_id);
		ind2 = kh_value(hash.ind2ind, k_i);

		// Check if GRM entry is already filled (if it is, keep the new kcoef value)
		if(ind1 == ind2) {
			if((family + fam)->phi[ind1][ind1]!=1 && (family + fam)->phi[ind1][ind1]!=kcoef)
				printf("WARNING: In GRM file, the diagonal element for individual %s is included twice with different values! The second value, %lf, is used.\n", ind1_id, kcoef);

			(family + fam)->phi[ind1][ind1] = kcoef;
		} else {
			if((family + fam)->phi[ind1][ind2]!=0 && (family + fam)->phi[ind1][ind2]!=kcoef)
				printf("WARNING: In GRM file, the individual pair (%s and %s) is included twice with different values! The second value, %lf, is used.\n", ind1_id, ind2_id, kcoef);

			(family + fam)->phi[ind1][ind2] = kcoef;
			(family + fam)->phi[ind2][ind1] = kcoef;
		}
	}
	fclose(kinfile);
}
