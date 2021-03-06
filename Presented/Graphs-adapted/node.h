#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NODE

#ifndef EDGE
#include "edge.h"
#endif

struct node
{
    int id;
    int strategy;
    int future_strategy;
    int changed_since_t0;
    float acumulated_payoff;
    struct edge_list edges;
    unsigned int n_neighbours;
};

/* Connect two nodes using an edge */
void connect(struct node * nod1, struct node * nod2, struct edge * edg)
{
    struct edge_list *list;
    list=&(nod1->edges);

    if(!check_neighbour(list,nod2))
    {
        edg->a=nod1;
        edg->b=nod2;
        add_edge(list,edg);

        list=&(nod2->edges);
        add_edge(list,edg);

        nod1->n_neighbours++;
        nod2->n_neighbours++;
        return;
    }
    edg->a=NULL;
    edg->b=NULL;
}
