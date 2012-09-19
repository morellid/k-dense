//
//  main.c
//  KDense
//
//  Created by Gabriele Cocco on 8/8/12.
//  Copyright (c) 2012 Gabriele Cocco. All rights reserved.
//

#include <stdio.h>
#include <igraph.h>
#include <string.h>
#include <stdlib.h>

#define DATA_FILE "/Users/Gabriele/Desktop/del/k-dense/topologie/topology_2012.txt"
#define KCORE_OUT_FILE "/Users/Gabriele/Desktop/del/k-dense/topologie/kcore_2012_true.txt"
#define KDENSE_OUT_FILE "/Users/Gabriele/Desktop/del/k-dense/topologie/kdense_2012.txt"
#define INFO_FILE "/Users/Gabriele/Desktop/del/k-dense/topologie/info_2012.txt"
#define MAX_STRING_LENGTH 256

int igraph_kdense(const igraph_t *graph,
                  unsigned int k,
                  int* deleted_count,
                  igraph_neimode_t mode);

int igraph_kcore(const igraph_t *graph,
                 unsigned int k,
                 int* deleted_count,
                 igraph_neimode_t mode);

int igraph_kdense(const igraph_t *graph,
                  unsigned int k,
                  int* deleted_count,
                  igraph_neimode_t mode) {
    
    // List of edges to delete
    //igraph_vector_t edges_to_delete;
    //igraph_vector_init(&edges_to_delete, 0);
    
    // Neightbors of first and second nodes of each edge
    igraph_vector_t first_neighbors;
    igraph_vector_t second_neighbors;
    igraph_vector_init(&first_neighbors, 0);
    igraph_vector_init(&second_neighbors, 0);
    
    int delete_at_least_one_edge = 1;
    while(delete_at_least_one_edge) {
        delete_at_least_one_edge = 0;
    
        // Do kcore(k)
        igraph_kcore(graph, k, deleted_count, mode);
        
        // Foreach edge...
        //Avoid iterator: sometimes failes in deleting vector of edges
        //igraph_eit_t iterator;
        //igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &iterator);
        //while(!IGRAPH_EIT_END(iterator)) {mode);
        
        long ecount = igraph_ecount(graph);
        for(long edge = ecount - 1; edge >= 0; edge--) {
            // Get vertices connected by the edge
            int first, second;
            igraph_edge(graph, (int)edge, &first, &second);
            
            // Get respective neighbors
            igraph_neighbors(graph, &first_neighbors, first, mode);
            igraph_neighbors(graph, &second_neighbors, second, mode);
            
            // Count common neighbors
            int common = 0;
            for(unsigned int i = 0; i < igraph_vector_size(&first_neighbors); i++) {
                if(igraph_vector_contains(&second_neighbors, VECTOR(first_neighbors)[i]))
                    common++;
            }
            
            // If less than k - 2 shared neighbors push back edge to be deleted
            if(common < k - 2) {
                //igraph_vector_push_back(&edges_to_delete, edge);
                IGRAPH_CHECK(igraph_delete_edges(graph, igraph_ess_1((int)edge)));
                delete_at_least_one_edge = 1;
            }
            
            // Clear neighbors vectors
            igraph_vector_clear(&first_neighbors);
            igraph_vector_clear(&second_neighbors);
            
            // Next edge
            //IGRAPH_EIT_NEXT(iterator);
        }
        
        // Delete edges from graph        
        //igraph_delete_edges(graph, igraph_ess_vector(&edges_to_delete));
        
        // Release iterator
        //igraph_eit_destroy(&iterator);
        
    }
     // Destroy vectors
    //igraph_vector_destroy(&edges_to_delete);
    igraph_vector_destroy(&first_neighbors);
    igraph_vector_destroy(&second_neighbors);
    
    return 0;
}

// NOTE: graph will be changed by this function
int igraph_kcore(const igraph_t *graph,
                 unsigned int k,
                 int* deleted_count,
                 igraph_neimode_t mode) {
    
    int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t degrees;
    //igraph_vector_t vertices_to_delete;
    igraph_neimode_t omode;
    
    if (mode != IGRAPH_ALL && mode != IGRAPH_OUT && mode != IGRAPH_IN) {
        IGRAPH_ERROR("Invalid mode in k-cores", IGRAPH_EINVAL);
    }
    if (!igraph_is_directed(graph) || mode==IGRAPH_ALL) {
        mode=omode=IGRAPH_ALL;
    } else if (mode==IGRAPH_IN) {
        omode=IGRAPH_OUT;
    } else {
        omode=IGRAPH_IN;
    }
    
    // Init degrees vector
    igraph_vector_init(&degrees, 0);
    
    int delete_at_least_one_vertex = 1;
    while(delete_at_least_one_vertex) {
        delete_at_least_one_vertex = 0;
        no_of_nodes = igraph_vcount(graph);
        
        // Get node degrees
        IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), mode, IGRAPH_LOOPS));
        
        // Delete vertices from the last to the first (cause delete changes vertices ids)
        for(long i = no_of_nodes - 1; i >= 0; i--) {
            if(VECTOR(degrees)[i] < k - 1) {
                // Store original index in deleted vector
                (*deleted_count)++;
                // Remember I've to go on 
                delete_at_least_one_vertex = 1;
                // Actually delete vertex
                IGRAPH_CHECK(igraph_delete_vertices(graph, igraph_vss_1((int)i)));
            }
        }
        
        // Clear degree and delete vectors
        igraph_vector_clear(&degrees);
    }
    
    // Destroy vectors
    igraph_vector_destroy(&degrees);
    
    return 0;
}

int main(int argc, const char * argv[])
{
    igraph_t graph;
    igraph_vector_t kcore_per_vertex;
    
    printf("Initializing...");
    // Init vector
    igraph_vector_init(&kcore_per_vertex, 0);
    
    // Create empty graph
    igraph_empty(&graph, 0, 0);
    printf("DONE\n");
    
    printf("Loading file and normalizing IDs...");
    // Open topology data file
    FILE* file = fopen(DATA_FILE, "r");
    
    // Read file and normalize IDs
    int hash[1 << 20];
    memset(hash, -1, (1 << 20) * sizeof(int));
    int current_id = 0;
    
    char first_node[MAX_STRING_LENGTH];
    char second_node[MAX_STRING_LENGTH];
    while(!feof(file)){
        fscanf(file, "%s", first_node);
        fscanf(file, "%s", second_node);
        
        int f_node = atoi(first_node);
        int s_node = atoi(second_node);
        
        if(hash[f_node] == -1) {
            // f_node id not inserted: associate normalized id and insert into graph
            printf("Found node %8d: normalized to %8d\n", f_node, current_id);
            hash[f_node] = current_id;
            current_id++;
            IGRAPH_CHECK(igraph_add_vertices(&graph, 1, 0));
        }
        if(hash[s_node] == -1) {
            // s_node id not inserted: associate normalized id and insert into graph
            printf("Found node %8d: normalized to %8d\n", s_node, current_id);
            hash[s_node] = current_id;
            current_id++;
            IGRAPH_CHECK(igraph_add_vertices(&graph, 1, 0));
        }
        igraph_add_edge(&graph, hash[f_node], hash[s_node]);
    }
    
    printf("DONE\n");
    
    // Get number of vertices
    long vcount = igraph_vcount(&graph);
  
    // Get k-core
    IGRAPH_CHECK(igraph_coreness((const igraph_t*)&graph, &kcore_per_vertex, IGRAPH_ALL));
     /*
    unsigned int* vertices_per_kcore;
    // Output kcore: number of vertex with core K (for each K)
    vertices_per_kcore = igraph_Calloc(igraph_vector_size(&kcore_per_vertex), unsigned int);
    for(unsigned int i = 0; i < igraph_vector_size(&kcore_per_vertex); i++)
    {
        int core = VECTOR(kcore_per_vertex)[i];
        vertices_per_kcore[core]++;
    }    
    for(unsigned int i = 0; i < igraph_vector_size(&kcore_per_vertex); i++)
    {
        int count = vertices_per_kcore[i];
        printf("Numero nodi con core %d = %d\n", i, count);
    }
    FILE* kcore_file = fopen(KCORE_OUT_FILE, "w");
    for(unsigned int i = 0; i < igraph_vector_size(&kcore_per_vertex); i++)
        fprintf(kcore_file, "%d %d\n", i, vertices_per_kcore[i]);
    
    printf("\n");
    */
    // Get k-core alternative
    // Store original ids
    igraph_t kcore;
    int deleted_count = 0;
    unsigned int* vertices_per_kcore_alt;
    int kcore_k_max = 0;
    
    vertices_per_kcore_alt = igraph_Calloc(vcount, unsigned int);
    
    igraph_copy(&kcore, &graph);
    
    // For k = 1 .. N compute kcore and check which vertex are deleted
    int can_stop = 0;
    for(int k = 2; !can_stop && k < vcount; k++) {
        deleted_count = 0;
        
        igraph_kcore(&kcore, k, &deleted_count, IGRAPH_ALL);
        
        vertices_per_kcore_alt[k - 1] = deleted_count;        
        if(igraph_vcount(&kcore) == 0) {
            can_stop = 1;
            kcore_k_max = k - 1;
        }
    }
    
    for(unsigned int i = 0; i < vcount; i++)
    {
        int count = vertices_per_kcore_alt[i];
        printf("Numero nodi con core %d = %d\n", i, count);
    }
    
    // Get k-dense
    FILE* kcore_file = fopen(KCORE_OUT_FILE, "w");
    for(unsigned int i = 0; i < vcount; i++)
        fprintf(kcore_file, "%d %d\n", i, vertices_per_kcore_alt[i]);

    // Get k-core alternative end
    
    // Get k-dense
    igraph_t kdense;
    igraph_copy(&kdense, &graph);
    int kdense_k_max = 0;
    
    // Clear vertices list
    memset(vertices_per_kcore_alt, 0, vcount * sizeof(unsigned int));
    
    // For k = 1 .. N compute kcore and check which vertex are deleted
    can_stop = 0;
    for(int k = 3; !can_stop && k < igraph_vcount(&graph); k++) {
        deleted_count = 0;
        
        igraph_kdense(&kdense, k, &deleted_count, IGRAPH_ALL);
        
        vertices_per_kcore_alt[k - 1] = deleted_count;        
        if(igraph_vcount(&kdense) == 0) {
            can_stop = 1;
            kdense_k_max = k - 1;
        }
    }
    
    for(unsigned int i = 0; i < vcount; i++)
    {
        int count = vertices_per_kcore_alt[i];
        printf("Numero nodi con core %d = %d\n", i, count);
    }
    // Get k-dense
    FILE* kdense_file = fopen(KDENSE_OUT_FILE, "w");
    for(unsigned int i = 0; i < vcount; i++)
        fprintf(kdense_file, "%d %d\n", i, vertices_per_kcore_alt[i]);

    // Calcolare indici prestazione
    FILE* info_file = fopen(INFO_FILE, "w");
    igraph_t graphs[2];
    igraph_copy(&kcore, &graph);
    igraph_kcore(&kcore, kcore_k_max, &deleted_count, IGRAPH_ALL);
    igraph_copy(&kdense, &graph);
    igraph_kdense(&kdense, kdense_k_max, &deleted_count, IGRAPH_ALL);
    graphs[0] = kcore;
    graphs[1] = kdense;
    
    for(unsigned int i = 0; i < 2; i++) {
        fprintf(info_file, "Method %d\n", i + 1);
        
        igraph_integer_t diameter;
        igraph_real_t unconnected, assortativity;
        igraph_diameter(&graphs[i], &diameter, NULL, NULL, NULL, 0, 1);
        
        fprintf(info_file, "Diameter;%d\n", diameter);
        
        igraph_assortativity_degree(&graphs[i], &assortativity, 0);
        fprintf(info_file, "Assortativity;%f\n", assortativity);
        
        fprintf(info_file, "Node;Betwenness;Closeness;Transitivity\n");
        
        igraph_vector_t betwenness, transitivity, closeness;
        igraph_vector_init(&betwenness, 0);
        igraph_vector_init(&closeness, 0);
        igraph_vector_init(&transitivity, 0);
        
        igraph_betweenness(&graphs[i], &betwenness, igraph_vss_all(), 0, NULL, 1);
        igraph_closeness(&graphs[i], &closeness, igraph_vss_all(), IGRAPH_ALL, NULL);
        igraph_transitivity_local_undirected(&graphs[i], &transitivity, igraph_vss_all(), IGRAPH_TRANSITIVITY_ZERO);
        
        /* Classificazione ? */
        for(unsigned int v = 0; v < igraph_vcount(&graphs[i]); v++) {
            fprintf(info_file, "%d;%f;%f;%f\n", v, VECTOR(betwenness)[v], VECTOR(closeness)[v], VECTOR(transitivity)[v]);
        }
        fprintf(info_file, "\n");
        
        igraph_vector_destroy(&betwenness);
        igraph_vector_destroy(&closeness);
        igraph_vector_destroy(&transitivity);
    }    
    
    // Release and free
    igraph_destroy(&graph);
    igraph_destroy(&kcore);
    igraph_destroy(&kdense);
    
    igraph_vector_destroy(&kcore_per_vertex);
    
    fclose(file);
    fclose(kcore_file);
    fclose(kdense_file);
    fclose(info_file);
    //igraph_free(vertices_per_kcore);
    free(vertices_per_kcore_alt);
    
    return 0;
}

