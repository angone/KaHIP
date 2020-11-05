import re, os, sys, copy
import numpy as np
import argparse
import random
import math
import time
import operator
from dwave_qbsolv import QBSolv
import pickle
import logging
import networkx as nx

def read_metis_graph(file_location):
    
    '''Read graph from Metis file format'''
    ### INPUT:
    # file_location: file location of the graph
    ### OUTPUT:
    # G: networkx graph with n nodes, node index from 1 to n

    with open(file_location, "r") as input_data_file:
        lines = input_data_file.read()
    lines = lines.split("\n")
    info = lines[0].split(" ")
    G = nx.Graph()
    
    if len(info) == 2:
        n = int(info[0])
        m = int(info[1])
        logging.info(f"Graph has {n} nodes and {m} edges")
        for i in range(1, n+1):
            line = lines[i].split(" ")
            G.add_node(i, weight = 1)
            if line[-1] == "":
                d = len(line) - 1
            else:
                d = len(line)
            for j in range(d):
                G.add_edge(i, int(line[j]), weight = 1)
            assert(G.degree[i] == d)
        assert(G.number_of_nodes() == n)
        assert(G.number_of_edges() == m)

    elif len(info) == 3:
        n = int(info[0])
        m = int(info[1])
        f = int(info[2])
        if f == 1:
            logging.info(f"Graph has {n} nodes and {m} edges, graph has edge weights")
            logging.error("Not implemented")
        elif f == 10:
            logging.info(f"Graph has {n} nodes and {m} edges, graph has node weights")
            logging.error("Not implemented")
        elif f == 11:
            logging.info(f"Graph has {n} nodes and {m} edges, graph has edge and node weights")
            for i in range(1, n+1):
                line = lines[i].split(" ")
                G.add_node(i, weight = int(line[0]))
                d = int((len(line) - 1) / 2)
                for j in range(d):
                    G.add_edge(i, int(line[2*j+1]), weight = int(line[2*j+2]))
                assert(G.degree[i] == d)
            assert(G.number_of_nodes() == n)
            assert(G.number_of_edges() == m)
            
        else:
            logging.error("Error: f value")
    else:
        logging.error("Error: file format not supported")

    return G
    
    
def vsp_qubo_qbsolv(G, A, B, C):
    
    '''Formulate QUBO of VSP to qbsolv format'''
    ### INPUT:
    # G: networkx graph object
    ### OUTPUT:
    # Q: dictionary of QUBO in qbsolv format
    
    Q = {}
    
    n = G.number_of_nodes()
    
    for i in range(n):
        Q[(i, i)] = - G.nodes[i+1]["weight"]
        Q[(i+n, i+n)] = - G.nodes[i+1]["weight"]
    for i, j in G.edges():
        Q[(i-1, j+n-1)] = G[i][j]["weight"] * A
        Q[(j-1, i+n-1)] = G[i][j]["weight"] * A
    for i in range(n):
        Q[(i, i+n)] = B
    for i in range(n):
        Q[(i, i)] += G.nodes[i+1]["weight"] ** 2 * C
        Q[(i+n, i+n)] += G.nodes[i+1]["weight"] ** 2 * C
    for i in range(n):
        for j in range(i+1, n):
            Q[(i, j)] = 2 * G.nodes[i+1]["weight"] * G.nodes[j+1]["weight"] * C
            Q[(i+n, j+n)] = 2 * G.nodes[i+1]["weight"] * G.nodes[j+1]["weight"] * C
    for i in range(n):
        for j in range(n):
            if (i, j+n) in Q:
                Q[(i, j+n)] -= 2 * G.nodes[i+1]["weight"] * G.nodes[j+1]["weight"] * C
            else:
                Q[(i, j+n)] = - 2 * G.nodes[i+1]["weight"] * G.nodes[j+1]["weight"] * C
    
    return Q


def get_separator_qbsolv(G, response):
    
    '''Formulate QUBO of VSP to qbsolv format'''
    ### INPUT:
    # G: networkx graph object
    # response: response from qbsolv
    ### OUTPUT:
    # solution: dictionary of solution, node index from 1 to n
    solution = {}
    n = G.number_of_nodes()
    for i in range(n):
        if response[i] == 1 and response[i+n] == 0:
            solution[i+1] = 0
        elif response[i] == 0 and response[i+n] == 1:
            solution[i+1] = 1
        elif response[i] == 0 and response[i+n] == 0:
            solution[i+1] = 2
        else:
            logging.error(f"Error: node {i} does not belong to any partition")
    part1 = set()
    part2 = set()
    separator = set()
    for i in range(n):
        if solution[i+1] == 0:
            part1.add(i+1)
        elif solution[i+1] == 2:
            separator.add(i+1)
        else:
            part2.add(i+1)
    for i, j in G.edges():
        if solution[i] + solution[j] == 1:
            logging.error(f"Error: node {i} and node {j} not separated correctly")
    
    f = open("coarsestvsp", "w")
    for i in range(n):
        f.write(f"{solution[i+1]}\n")
    f.close()
    
    return part1, part2, separator
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type = str,
                                  default = "example",
                                  help = "path to graph")
    parser.add_argument("--log", type = str,
                                 help = "logging level")
    parser.add_argument("--A", type = int,
                               default = 100,
                               help = "large postive number to penalize violation of connection constraint")
    parser.add_argument("--B", type = int,
                               default = 100,
                               help = "large postive number to penalize violation of uniqueness constraint")
    parser.add_argument("--C", type = int,
                               default = 100,
                               help = "large postive number to penalize imbalance of two partitions")
    
    args = parser.parse_args()
    if args.log == "INFO":
        logging.basicConfig(level = logging.INFO)
    elif args.log == "DEBUG":
        logging.basicConfig(level = logging.DEBUG)
    
    G = read_metis_graph(args.path)
    Q = vsp_qubo_qbsolv(G, args.A, args.B, args.C)
    response = QBSolv().sample_qubo(Q)
    response = response.samples()[0]
    get_separator_qbsolv(G, response)
    