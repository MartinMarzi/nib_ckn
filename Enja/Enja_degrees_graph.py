import networkx as nx
import numpy as np
import csv
from operator import itemgetter
import pandas as pd

pd.options.display.max_rows = 20
from collections import defaultdict
import matplotlib.pyplot as plt

plt.rcParams.update({'figure.figsize': (10, 10)})

dt = "../data/LKN_ath_v3B_2016-08-30.txt"
f = open(dt, "r")
df = pd.read_table(dt, sep="\t", header=None, names=["entity1", "entity2", "relation", "r", "PIS"])


# print(df.head())

# Vid's function to import data
def parse_ckn_csv(fname):
    g = nx.DiGraph()
    fields = ['from', 'to', 'type']
    with open(fname, newline='') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(2048))
        csvfile.seek(0)
        reader = csv.DictReader(csvfile, fieldnames=fields, dialect=dialect, restkey='rest', )
        for row in reader:
            g.add_edge(row['from'], row['to'], type=row['type'])
            if row['type'] == 'binding':
                g.add_edge(row['to'], row['from'], type=row['type'])
    return g


# out/in/degree
def get_degree(G):
    out_degree = defaultdict(int)
    in_degree = defaultdict(int)
    for e1, e2, t in G.edges.data():
        out_degree[e1] += 1
        in_degree[e2] += 1
    degree_ = defaultdict(int)
    for i in [out_degree, in_degree]:
        for ii in i:
            degree_[ii] += i[ii]
    return degree_, out_degree, in_degree


def dict_to_list(x):
    return [(xx, x[xx]) for xx in x]


def sort_descending(s):
    return sorted(s, key=itemgetter(1), reverse=True)


def print_line(x):
    for i, ii in x:
        print(i, ii)


# edge types
def get_edges_type(G):
    "edges by types"
    edges_by_type = defaultdict(list)
    for e1, e2, t in G.edges.data():
        tt = t["type"]
        edges_by_type[tt].append((e1, e2))
    edges_by_type_len = {x: len(edges_by_type[x]) for x in edges_by_type}

    # edges all {edge: n}
    edges_all_n = defaultdict(int)
    for type_ in edges_by_type:
        # print("\n" + str(type_))
        edges_by_type_freq = defaultdict(int)
        all_ = 0
        for x, y in edges_by_type[type_]:
            all_ += 2
            edges_by_type_freq[x] += 1
            edges_by_type_freq[y] += 1
        # print(all_, sum([y for _, y in dict_to_list(edges_by_type_freq)]))

        for x in edges_by_type_freq:
            edges_all_n[x] += edges_by_type_freq[x]

    return edges_by_type, edges_by_type_len, edges_all_n


def subset_to_graph(edges_all_n, top_n=10):
    sorted_all_ = sort_descending(dict_to_list(edges_all_n))[:top_n*10]
    sorted_all = sort_descending(dict_to_list(edges_all_n))[:top_n]
    sorted_all_nodes_ = [x for x, _ in sorted_all_]
    sorted_all_nodes = [x for x, _ in sorted_all]
    print("\n")
    print_line(sorted_all)  # ubistvu isto kot da bi vzela degree()[:10]

    # subset za graf
    to_add = []
    edge_colors_ = defaultdict(list)
    for e1, e2, t in G.edges.data():
        tt = t["type"]
        if {e1, e2}.issubset(set(sorted_all_nodes_)):
            # if e1 == e2:
            #     print(e1, e2, t)
            if e1 != e2:
                to_add.append((e1, e2, t))
                edge_colors_[tt].append((e1, e2))
    return to_add, edge_colors_


colors_ = {"binding": "red",
           "unk_TF": "blue",
           "activation": "green",
           "act_TF": "orange",
           "inh_TF": "black",
           "inhibition": "grey",
           "synthesis": "yellow",
           "inh_ncRNA": "purple"}


if __name__ == "__main__":
    G = parse_ckn_csv(dt)

    # get degree (in_out)
    degree_, out_degree, in_degree = get_degree(G)
    print("degree_legit", sort_descending(G.degree)[:20])
    print("degree_mine", sort_descending(dict_to_list(degree_))[:20])
    print("out_degree_legit", sort_descending(G.out_degree)[:20])
    print("out_degree_mine", sort_descending(dict_to_list(out_degree))[:20])
    print("in_degree_legit", sort_descending(G.in_degree)[:20])
    print("in_degree_mine", sort_descending(dict_to_list(in_degree))[:20])

    # po tipu edge-ov
    edges_by_type, edges_by_type_len, edges_all_n = get_edges_type(G)
    print("\n")
    print_line(sorted(dict_to_list(edges_by_type_len), key=itemgetter(1), reverse=True))

    # get subset za graf
    top_n = 10
    to_add, edge_colors_ = subset_to_graph(edges_all_n, top_n=top_n)
    GG = nx.DiGraph()
    GG.add_edges_from(to_add)
    pos = nx.spring_layout(GG)
    # pos = nx.graphviz_layout(GG)
    # nx.draw_networkx(GG, pos=pos, arrowsize=5, with_labels=False, node_size=100)
    # plt.show()

    nx.draw_networkx_nodes(GG, pos=pos, node_size=100, label=True)
    nx.draw_networkx_labels(GG, pos=pos, font_size=10)
    for cl in edge_colors_:
        nx.draw_networkx_edges(GG, pos, edgelist=edge_colors_[cl], edge_color=colors_[cl])  # style="dashed")
    plt.savefig(fname="graph_subset_n={}.png".format(top_n))
    plt.show()
