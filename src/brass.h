#ifndef BRASS_H
#define BRASS_H


#ifdef __cplusplus
extern "C"
{
#endif

#define MISSVAL -9.0
#define NUMTOL 1e-6
#define ANALYSIS_STATUS 1
#define MAXFILELEN 2001 //maximum length of file names
#define MAXINDLEN 20 //maximum length of individual ID
#define MAXPHELEN 1000 //maximum length of lines in phenotype file
#define POSDEF_CUTOFF 1e-8
#define NPERM_MAX 1e6 // maximum number of permutations


    struct FAMILY {
        int n_total;			   // total number individuals in family
        char **indiv;			   // Subject IDs
        double *trait_obs;		   // trait[index among phenotyped] */
        double naff; 			   // Number of cases for observed trait
        double **cov;		   	   // cov[i][k] = cov value for ind i at cov k
        double **phi;			   // GRM
        double *eigvalue;		   // Eigenvalues of GRM
        double **eigvector;		   // Eigenvectors of GRM
        double *trait;			   // Vector used for fitting null model
        double *Z;			   // Vector used for fitting null model
        double *beta_hat;		   // Vector used for fitting null model
        double *Dinv;		       // Vector used for fitting null model
        double *delta;			   // Linear transformation of residuals w/ uncorrelated entries
        double *mu_hat_QL;         // Estimate of the mean of Y in null model
        double **pre_mult_mat;     // Scaling matrix used to un-scale permuted vector delta
        double *trait_rep;		   // Trait replicate
    };

    struct DATA_STRUCT {
        int n_fam;		/* total number of families */
        int n_total;		/* total number of individuals (across all families) */
        int n_cov;		/* total number of covariates including intercept */
        double tol;
        int n_perm;     // Number of trait replicates to be generated
    };

    void leigcomp_kin(struct FAMILY *family);
    void eigcomp_kin(struct FAMILY *family);
    void write_eigendecomp(char *outeig_filename, struct FAMILY *family);
    void compute_logistic_coeff(double *yaffec, double **Cov, int totpc, int n_cov, double *beta_init);
    void estimate_gee_parameters(int method, int totpc, double **V, double *L, double *Dinv, struct FAMILY *family, int n_cov, int n_fam, double *xi, double *betahat);
    double target_func_xi(double xi, int totpc, double **V, double *L, double *Dinv, struct FAMILY *family, double *current_beta, int n_fam, int n_cov);
    int compute_gee_betahat(double **V, double *Dinv, double *current_beta, struct FAMILY *family, int n_fam, int n_cov);
    void pre_permute(int get_perm, double xi, int n_cov, struct FAMILY *family, char *filename);
    void get_quant_rep(struct FAMILY *family, struct DATA_STRUCT data_struct);
    void naive_shuffle(double *trait, int n_total);

#ifdef __cplusplus
}
#endif

#endif
