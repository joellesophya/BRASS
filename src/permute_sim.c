#include "permute_sim.h"

void permute(struct FAMILY *family, struct HASH hash, struct DATA_STRUCT data_struct, char *filename) {

    int i_perm, fam = 1;
    int n_total = (family + fam)->n_total;
    int n_perm = data_struct.n_perm;
    clock_t time1, time2;
    double time_out;
    // Initialize array for trait replicates
    (family + fam)->trait_rep = dvector(1, n_total);

    time1 = clock();

    FILE *yrepfile;
    yrepfile = fopen(filename, "w");

    print_rep(0, family, hash, yrepfile);
    for (i_perm = 1; i_perm <= n_perm; i_perm++) {
        get_quant_rep(family, data_struct);
        print_rep(i_perm, family, hash, yrepfile);
    }
    time2 = clock();
    time_out = (double) (time2 - time1) / CLOCKS_PER_SEC;
    printf("Time for BRASS with %d replicates: ", n_perm);
    if(time_out > 60)
        printf("%d min %d secs\n", (int) (time_out / 60), ((int) time_out) % 60);
    else printf("%lf secs\n", time_out);
    fclose(yrepfile);
}

void print_rep(int i_perm, struct FAMILY *family, struct HASH hash, FILE *file_out){

    int i, fam = 1;
    int n_total = (family + fam)->n_total;

    // Print original individual IDs on first line
    for (i = 1; i <= n_total; i++) {

        if (i_perm == 0) {
            fprintf(file_out, "%s%s", (family + fam)->indiv[i], (i == n_total ? "" : " "));
        } else {
            fprintf(file_out, "%.17g%s", (family + fam)->trait_rep[i], (i == n_total ? "" : " "));
        }
    }
    fprintf(file_out, "\n");
    fflush(file_out);

}
