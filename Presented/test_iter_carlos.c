#include "mtwist-1.1/mtwist.c"
#include "game.h"
#include <time.h>
//TYPE 0==COOPERATORS
//TYPE 1==DEFECTORS

int main(int argc, const char *argv[])
{
    mt_seed();
    unsigned int i,j;
    unsigned int total_games=100;
    double b;

    unsigned int **count;
    unsigned int **count_aux;
    count=malloc(total_games*sizeof(unsigned int *));
    count_aux=malloc(total_games*sizeof(unsigned int *));

    unsigned int **coop_cluster;
    coop_cluster=malloc(total_games*sizeof(unsigned int *));
    for (i = 0; i < total_games; i++)
        coop_cluster[i]=NULL;
    double coop_mean;
    double coop_fixed_mean;
    double coop_isolated_mean;
    double coop_cluster_mean;
    double coop_effective_surface_mean;

    unsigned int **def_cluster;
    def_cluster=malloc(total_games*sizeof(unsigned int *));
    for (i = 0; i < total_games; i++)
        def_cluster[i]=NULL;
    double def_mean;
    double def_fixed_mean;
    double def_isolated_mean;
    double def_cluster_mean;
    double def_effective_surface_mean;
    
    unsigned int *iter;
    unsigned int iter0=5000;
    unsigned int iter1;
    unsigned int iter_mean;
    iter=malloc(total_games*sizeof(unsigned int));

    FILE *f;
    f=fopen("data_iter2.csv","w");
    fprintf(f,"b\t");
    fprintf(f,"coop_mean\tcoop_fixed_mean\tcoop_isolated_mean\tcoop_cluster_mean\tcoop_effective_surface_mean\t");
    fprintf(f,"def_mean\tdef_fixed_mean\tdef_isolated_mean\tdef_cluster_mean\tdef_effective_surface_mean\t");
    fprintf(f,"iter_mean\n");
    fclose(f);

    struct game gam;
    struct game *g_aux;
    g_aux=&gam;

    for (b = 1.5; b < 2.0; b+=0.05)
    {
        for (j = 0; j < total_games;) 
        {
            gam=create_prisioner_dilemma_game(b,4000,8000,1);

            for (i = 0; i < iter0; i++)
                evolve_game(g_aux,0);
            iter[j]=iter0;

            iter1=1000;
            do
            {
                count[j]=count_type(gam);
                for (i = 0; i < iter1; i++)
                    evolve_game(g_aux,0);
                iter[j]+=iter1;
                count_aux[j]=count_type(gam);
            }
            while (abs(count[j][0]-count_aux[j][0])>10);

            iter1=10000;
            for (i = 0; i < iter1; i++)
                evolve_game(g_aux,0);
                
            iter1=1000;
            for (i = 0; i < iter1; i++)
                evolve_game(g_aux,1);

            free(count[j]);
            free(coop_cluster[j]);
            free(def_cluster[j]);

            count[j]=count_type(gam);
            coop_cluster[j]=cluster_analize(gam,0);
            def_cluster[j]=cluster_analize(gam,1);
/*
            if(coop_cluster[j][2]==0)
                printf("b: %.3lf game: %u FAILED, COOPERATORS DESTROYED\n",b,j);
            else if(def_cluster[j][2]==0)
                printf("b: %.3lf game: %u FAILED, DEFECTORS DESTROYED\n",b,j);
            else
            {*/
                printf("b: %.3lf game: %u DONE\n",b,j);
                j++;
            //}

            //graph_to_file("grafo.gdf",gam.g,1);

            free_game(g_aux);
        }

        coop_mean=0;
        coop_fixed_mean=0;
        coop_isolated_mean=0;
        coop_cluster_mean=0;
        coop_effective_surface_mean=0;

        def_mean=0;
        def_fixed_mean=0;
        def_isolated_mean=0;
        def_cluster_mean=0;
        def_effective_surface_mean=0;

        iter_mean=0;
        for (i = 0; i < total_games; i++)
        {
            coop_mean+=(double)count[i][0]/total_games;
            coop_fixed_mean+=(double)coop_cluster[i][2]/total_games;
            coop_isolated_mean+=(double)coop_cluster[i][4]/total_games;
            coop_cluster_mean+=(double)(coop_cluster[i][0]-coop_cluster[i][1])/total_games;
            coop_effective_surface_mean+=(double)coop_cluster[i][3]/(coop_cluster[i][2]-coop_cluster[i][4])/total_games;

            def_mean+=(double)count[i][1]/total_games;
            def_fixed_mean+=(double)def_cluster[i][2]/total_games;
            def_isolated_mean+=(double)def_cluster[i][4]/total_games;
            def_cluster_mean+=(double)(def_cluster[i][0]-def_cluster[i][1])/total_games;
            def_effective_surface_mean+=(double)def_cluster[i][3]/(def_cluster[i][2]-def_cluster[i][4])/total_games;

            iter_mean+=iter[i]/total_games;
        }

        coop_mean/=4000.;
        coop_fixed_mean/=4000.;
        coop_isolated_mean/=4000.;

        def_mean/=4000.;
        def_fixed_mean/=4000.;
        def_isolated_mean/=4000.;

        f=fopen("data_iter2.csv","a");
        fprintf(f,"%lf\t",b);
        fprintf(f,"%lf\t%lf\t%lf\t%lf\t%lf\t",coop_mean,coop_fixed_mean,coop_isolated_mean,coop_cluster_mean,coop_effective_surface_mean);
        fprintf(f,"%lf\t%lf\t%lf\t%lf\t%lf\t",def_mean,def_fixed_mean,def_isolated_mean,def_cluster_mean,def_effective_surface_mean);
        fprintf(f,"%u\n",iter_mean);
        fclose(f);
        printf("\n");
    }
    return 0;
}
