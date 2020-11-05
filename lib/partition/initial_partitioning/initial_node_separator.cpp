//
// Author: Christian Schulz <christian.schulz.phone@gmail.com>
// 

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "initial_node_separator.h"
#include "graph_partitioner.h"
#include "tools/quality_metrics.h"
#include "tools/random_functions.h"
#include "partition/uncoarsening/separator/vertex_separator_algorithm.h"
#include "partition/uncoarsening/refinement/node_separators/fm_ns_local_search.h"

initial_node_separator::initial_node_separator() {
                
}

initial_node_separator::~initial_node_separator() {
                
}

NodeWeight initial_node_separator::single_run( const PartitionConfig & config, graph_access & G) {

        std::streambuf* backup = std::cout.rdbuf();
        std::ofstream ofs;
        ofs.open("/dev/null");
        std::cout.rdbuf(ofs.rdbuf()); 

        graph_partitioner partitioner;
        PartitionConfig partition_config         = config;
        partition_config.mode_node_separators    = false;
        partition_config.global_cycle_iterations = 1;
        partition_config.repetitions             = 1;

        //computing a partition
        partitioner.perform_partitioning(partition_config, G);

        complete_boundary boundary(&G);
        boundary.build();

        ofs.close();
        std::cout.rdbuf(backup);

        vertex_separator_algorithm vsa; std::vector<NodeID> separator;
        //create a very simple separator from that partition
        if( partition_config.sep_full_boundary_ip ) {
                vsa.compute_vertex_separator_simpler(partition_config, G, boundary, separator);
        } else {
                vsa.compute_vertex_separator_simple(partition_config, G, boundary, separator);
        }
 
        std::vector<NodeID> output_separator;
        //improve the separator using flow based techniques
        vsa.improve_vertex_separator(partition_config, G, separator, output_separator);

        quality_metrics qm;
        return qm.separator_weight(G);

}

void initial_node_separator::compute_node_separator( const PartitionConfig & config, graph_access & G) {
        if(config.graph_allready_partitioned) return;

        std::vector< NodeID > best_separator(G.number_of_nodes(),0);
        NodeWeight best_separator_size = std::numeric_limits< NodeWeight >::max();

	int unsucc_counter = 0;

        for( int i  = 0; i  < config.max_initial_ns_tries; i++) {
                NodeWeight cur_separator_size = single_run(config, G);
                if(cur_separator_size < best_separator_size) {
                        forall_nodes(G, node) {
                                best_separator[node] = G.getPartitionIndex(node);
                        } endfor
                        best_separator_size = cur_separator_size;
                
                        PRINT(std::cout <<  "improved initial separator size " <<  cur_separator_size  << std::endl;)
			unsucc_counter = 0;
                } else {
			unsucc_counter++;
		}

		if( config.faster_ns && unsucc_counter >= 5) {
			break;
		}
        }

        forall_nodes(G, node) {
                G.setPartitionIndex(node, best_separator[node]);
        } endfor
        
        std::ofstream myfile;
        myfile.open ("example");
        myfile << G.number_of_nodes() << " " << G.number_of_edges()/2 << "\n";
        
        forall_nodes(G, node) {
                EdgeID first = G.get_first_edge(node);
                EdgeID end = G.get_first_invalid_edge(node);
                for (EdgeID i = first; i < end; i++) {
                    myfile << G.getEdgeTarget(i)+1 << " ";
                }
                myfile << "\n";
        } endfor
        myfile.close();
    
        system("python utils.py --path example");
    
        std::ifstream infile("coarsestvsp");
        std::string line;
        int counter = 0;
        while (std::getline(infile, line))
        {
            std::istringstream iss(line);
            
            int a;
            
            if (!(iss >> a)) { break; } // error
            //std::cout << "node " << counter << " in part " << a << std::endl;
            G.setPartitionIndex(counter, a);
            counter++;

            // process separator
        }

        
}
