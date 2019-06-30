#include "datacheck.h"

void format_checks(char *pedfilename, struct DATA_STRUCT *data_struct) {

	int col = 0, n_cov, fam, ntotal = 0, i;
	double pheno;
	char car, ind[MAXINDLEN];

	/* check if pedfile is present and can be opened */
	FILE *pedfile;
	pedfile = fopen(pedfilename, "r");
	if (pedfile == NULL) {
		printf("ERROR: Can't open phenotype data file!\n");
		exit(1);
	}

	/* read through the pedfile to count the number of columns, to then find the number of covariates */
	car = 0;
	while (car != EOF && car != '\n') {
		car = getc(pedfile);
		if (car == ' ' || car == '\t') {
			col++;
		}
	}
	col++;
	n_cov = col - 6 + 1; // +1 for the intercept
	if(n_cov < 2) printf("WARNING: No covariates present in phenotype data file!\n");
	fclose(pedfile);

	// Get total number of individuals
	pedfile = fopen(pedfilename, "r");
	while (!feof(pedfile)) {
		if( fscanf(pedfile, "%d %s %*s %*s %*d %lf", &fam, ind, &pheno) != 3) {
			printf("ERROR: Check format of phenotype file!");
			exit(1);
		}
		// check only one family is present
		if (fam != 1) {
			printf("WARNING: Family for individual %s is NOT 1 (only Fam = 1 is allowed in program)!\n", ind);
			exit(1);
		}
//		// check coding for phenotype
//		if(pheno!=0 && pheno!=1){
//			printf("ERROR: Phenotype for individual %s is NOT coded in 0/1 (no missing phenotype is allowed by the program)!\n", ind);
//			exit(1);
//		}
		ntotal++;

		for (col=1; col < (n_cov - 1); col++) fscanf(pedfile, "%*lf ");
		fscanf(pedfile, "%*lf\n");

	}

	fclose(pedfile);
	//Save no. of covariates present
	data_struct->n_cov = n_cov;
	data_struct->n_total = ntotal;
}
