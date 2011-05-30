#include "Graphs-adapted/graph.h"
#include <string.h>

struct game
{
    struct graph g;
    double **payoff_matrix;
};

struct game create_prisioner_dilemma_game(double b, unsigned int nods, unsigned int edgs, int graph)
{
    struct game pdg;
    
    if(graph==0)
        pdg.g=create_erdos_renyi_graph(nods,edgs);
    else
        pdg.g=create_barabasi_albert_graph(2,nods-2,2,1);

    pdg.payoff_matrix=malloc(2*sizeof(double *));
    pdg.payoff_matrix[0]=malloc(2*sizeof(double));
    pdg.payoff_matrix[1]=malloc(2*sizeof(double));

    pdg.payoff_matrix[0][0]=1;
    pdg.payoff_matrix[0][1]=b;
    pdg.payoff_matrix[1][0]=0;
    pdg.payoff_matrix[1][1]=0;

    return pdg;
}

struct game create_snowdrift_game(double r, unsigned int nods, unsigned int edgs, int graph)
{
    struct game sdg;
    
    if(graph==0)
        sdg.g=create_erdos_renyi_graph(nods,edgs);
    else
        sdg.g=create_barabasi_albert_graph(2,nods-2,2,1);

    sdg.payoff_matrix=malloc(2*sizeof(double *));
    sdg.payoff_matrix[0]=malloc(2*sizeof(double));
    sdg.payoff_matrix[1]=malloc(2*sizeof(double));

    sdg.payoff_matrix[0][0]=1/(2*r);
    sdg.payoff_matrix[0][1]=(1+r)/(2*r);
    sdg.payoff_matrix[1][0]=(1-r)/(2*r);
    sdg.payoff_matrix[1][1]=0;

    return sdg;
}

void change(struct node *nod1, double b, int passed);
void actualice_game(struct game *gam);
void evolve_game(struct game *gam, int t0_passed)
{
    unsigned int i;
    for (i = 0; i < gam->g.n_edges; i++)
    {
        gam->g.edges[i].a->acumulated_payoff+=gam->payoff_matrix[gam->g.edges[i].b->strategy][gam->g.edges[i].a->strategy];
        gam->g.edges[i].b->acumulated_payoff+=gam->payoff_matrix[gam->g.edges[i].a->strategy][gam->g.edges[i].b->strategy];
    }

    struct node *nod;
    for (i = 0; i < gam->g.n_nodes; i++)
    {
        nod=&gam->g.nodes[i];
        change(nod,gam->payoff_matrix[0][1], t0_passed);
    }

    actualice_game(gam);
}

void change(struct node *nod1, double b, int passed)
{
    if(nod1->n_neighbours==0)
        return;
        
    unsigned int i;
    i=(unsigned int)mt_lrand()%nod1->n_neighbours;

    struct edge_list_node *ed_aux;
    ed_aux=nod1->edges.first;
    while(0<i)
    {
        ed_aux=ed_aux->next;
        i--;
    }

    struct node *nod_aux;
    if(ed_aux->ed->a==nod1)
        nod_aux=ed_aux->ed->b;
    else
        nod_aux=ed_aux->ed->a;

    if(nod1->acumulated_payoff>=nod_aux->acumulated_payoff)
        return;
    else if(nod1->strategy!=nod_aux->strategy)
    {
        float aux;
        aux=nod_aux->acumulated_payoff-nod1->acumulated_payoff;
        
        if(nod1->n_neighbours>=nod_aux->n_neighbours)
            aux/=(nod1->n_neighbours*b);
        else
            aux/=(nod_aux->n_neighbours*b);

        if((float)mt_drand()<=aux)
        {
            nod1->future_strategy=nod_aux->strategy;
            if(passed==1)
                nod1->changed_since_t0=1;
        }
    }
}

void actualice_game(struct game *gam)
{
    unsigned int i;

    for (i = 0; i < gam->g.n_nodes; i++)
    {
        gam->g.nodes[i].strategy=gam->g.nodes[i].future_strategy;
        gam->g.nodes[i].acumulated_payoff=0;
    }
}

void free_game(struct game *gam)
{
    struct graph *g_aux;
    g_aux=&(gam->g);
    free_graph(g_aux);
}

unsigned int * count_type(struct game gam)
{
    unsigned int i;
    unsigned int *aux;
    aux=malloc(2*sizeof(unsigned int));
    memset(aux,0,2*sizeof(unsigned int));

    for (i = 0; i < gam.g.n_nodes; i++)
    {
        if(gam.g.nodes[i].strategy==0)
            aux[0]++;
        else if(gam.g.nodes[i].strategy==1)
            aux[1]++;
    }
    return aux;
}

/* 0 number of clusters
 * 1 number of dynamicaly isolated clusters
 * 2 number of fixed nodes
 * 3 number of frontier nodes
 * 4 number of totally isolated nodes (inside isolated clusters)
 * 5 number of solitary fixed nodes not isolated                */
unsigned int * cluster_analize(struct game gam, int str)
{
    unsigned int i, j;
    unsigned int node_added=0;
    unsigned int front_check=0;
    unsigned int isolated_check_cluster=0;
    unsigned int solitary_isolated=0;
    unsigned int nodes_in_cluster=0;
    unsigned int *aux;
    aux=malloc(6*sizeof(unsigned int));
    memset(aux,0,6*sizeof(unsigned int));

    struct edge_list_node *aux_list_node;
    struct node *aux_node2;

    int *processed;
    processed=malloc(gam.g.n_nodes*sizeof(int));
    memset(processed,0,gam.g.n_nodes*sizeof(unsigned int));

    for (i = 0; i < gam.g.n_nodes; i++)
    {
        if (gam.g.nodes[i].strategy==str && \
                gam.g.nodes[i].changed_since_t0==0 && \
                processed[i]==0)
        {
            processed[i]=1;

            for (j = i; j < gam.g.n_nodes ; j++)
            {
                if(processed[j]==1) 
                {
                    if(gam.g.nodes[j].n_neighbours!=0)
                    {
                        aux_list_node=gam.g.nodes[j].edges.first;

                        do
                        {
                            if(aux_list_node->ed->a->id==j)
                                aux_node2=aux_list_node->ed->b;
                            else
                                aux_node2=aux_list_node->ed->a;

                            if(aux_node2->changed_since_t0 || \
                                aux_node2->strategy!=str)
                            {
                                front_check++;
                                isolated_check_cluster++;
                            }

                            if(aux_node2->strategy==str && \
                                aux_node2->changed_since_t0==0)
                            {
                                solitary_isolated++;
                                if(processed[aux_node2->id]==0)
                                {
                                    processed[aux_node2->id]=1;
                                    node_added++;
                                }
                            }

                            aux_list_node=aux_list_node->next;
                        }
                        while (aux_list_node!=NULL);

                        if(solitary_isolated==0)
                            aux[5]++;   //Solitary isolated
                        solitary_isolated=0;

                        nodes_in_cluster++;
                    }
                    else
                        aux[4]++;   //Isolated node

                    //CHECKS NODES
                    if(front_check!=0)
                    {
                        aux[3]++;   //Node in the frontier
                        front_check=0;
                    }

                    aux[2]++;   //Fixed node

                    processed[j]=2;
                }
                if(j==(gam.g.n_nodes-1)&&node_added!=0)
                {
                    j=i;
                    node_added=0;
                }
            }

            //CHECKS CLUSTERS
            if(isolated_check_cluster==0 && \
                nodes_in_cluster<=gam.g.n_nodes/2.)
            {
                aux[1]++; //Isolated cluster
                aux[4]+=nodes_in_cluster;
            }
            nodes_in_cluster=0;
            isolated_check_cluster=0;

            aux[0]++; //Cluster
        }
    }

    free(processed);
    processed=NULL;
    return aux;
}
