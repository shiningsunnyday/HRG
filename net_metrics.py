__author__ = ['Salvador Aguinaga', 'Rodrigo Palacios', 'David Chaing', 'Tim Weninger']

import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
import collections
from collections import Counter
from random import sample
import os


def hops(all_succs, start, level=0, debug=False):
    if debug: print("level:", level)

    succs = all_succs[start] if start in all_succs else []
    if debug: print("succs:", succs)

    lensuccs = len(succs)
    if debug: print("lensuccs:", lensuccs)
    if debug: print()
    if not succs:
        yield level, 0
    else:
        yield level, lensuccs

    for succ in succs:
        # print("succ:", succ)
        for h in hops(all_succs, succ, level + 1):
            yield h


def get_graph_hops(graph, num_samples):
    c = Counter()
    for i in range(0, num_samples):
        node = sample(graph.nodes(), 1)[0]
        b = nx.bfs_successors(graph, node)

        for l, h in hops(b, node):
            c[l] += h

    hopper = Counter()
    for l in c:
        hopper[l] = float(c[l]) / float(num_samples)
    return hopper


def bfs_eff_diam(G, NTestNodes, P):
    EffDiam = -1
    FullDiam = -1
    AvgSPL = -1

    DistToCntH = {}

    NodeIdV = list(nx.nodes(G))
    random.shuffle(NodeIdV)

    for tries in range(0, min(NTestNodes, nx.number_of_nodes(G))):
        NId = NodeIdV[tries]
        b = nx.bfs_successors(G, NId)
        for l, h in hops(b, NId):
            if h is 0: continue
            if not l + 1 in DistToCntH:
                DistToCntH[l + 1] = h
            else:
                DistToCntH[l + 1] += h

    DistNbrsPdfV = {}
    SumPathL = 0.0
    PathCnt = 0.0
    for i in DistToCntH.keys():
        DistNbrsPdfV[i] = DistToCntH[i]
        SumPathL += i * DistToCntH[i]
        PathCnt += DistToCntH[i]

    oDistNbrsPdfV = collections.OrderedDict(sorted(DistNbrsPdfV.items()))

    CdfV = oDistNbrsPdfV
    for i in range(1, len(CdfV)):
        if not i + 1 in CdfV:
            CdfV[i + 1] = 0
        CdfV[i + 1] = CdfV[i] + CdfV[i + 1]

    try:
        EffPairs = P * CdfV[next(reversed(CdfV))]
    except:
        return 0

    for ValN in CdfV.keys():
        if CdfV[ValN] > EffPairs: break

    if ValN >= len(CdfV): return next(reversed(CdfV))
    if ValN is 0: return 1
    # interpolate
    DeltaNbrs = CdfV[ValN] - CdfV[ValN - 1];
    if DeltaNbrs is 0: return ValN;
    return ValN - 1 + (EffPairs - CdfV[ValN - 1]) / DeltaNbrs


def draw_diam_plot(orig_g, mG):
    df = pd.DataFrame(mG)
    gD = bfs_eff_diam(orig_g, 20, .9)
    ori_degree_seq = []
    for i in range(0, len(max(mG))):
        ori_degree_seq.append(gD)
    fig, ax = plt.subplots()
    ax.fill_between(df.columns, df.mean() - df.sem(), df.mean() + df.sem(), color='blue', alpha=0.2, label="se")
    h, = ax.plot(df.mean(), color='blue', aa=True, linewidth=4, ls='--', label="H*")
    orig, = ax.plot(ori_degree_seq, color='black', linewidth=2, ls='-', label="H")

    ax.set_title('Diameter Plot')
    ax.set_ylabel('Diameter')
    ax.set_xlabel('Growth')

    ax.xaxis.set_tick_params(
        # axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom='off',  # ticks along the bottom edge are off
        top='off',  # ticks along the top edge are off
        labelbottom='off')  # labels along the bottom edge are off
    ax.legend([orig, h], ['$H$', 'HRG $H^*$'], loc=4)
    # ax = plt.gcf()
    # ax.set_size_inches(5, 4, forward=True)
    fig.savefig(os.path.join(os.getcwd(), 'diam_plot.png'))



def draw_graphlet_plot(orig_g, mG):
    df = pd.DataFrame(mG)
    width=.25

    N=11
    dforig = pd.DataFrame(orig_g)
    means = (dforig.mean()['e0'], dforig.mean()['e1'], dforig.mean()['e2'], dforig.mean()['e2c'], dforig.mean()['tri'], dforig.mean()['p3'], dforig.mean()['star'], dforig.mean()['tritail'], dforig.mean()['square'], dforig.mean()['squarediag'], dforig.mean()['k4'])
    sem = (dforig.sem()['e0'], dforig.sem()['e1'], dforig.sem()['e2'], dforig.sem()['e2c'], dforig.sem()['tri'], dforig.sem()['p3'], dforig.sem()['star'], dforig.sem()['tritail'], dforig.sem()['square'], dforig.sem()['squarediag'], dforig.sem()['k4'])
    ind = np.arange(N)
    fig,ax = plt.subplots()
    rects = ax.bar(ind+.02, means, width-.02, color = 'k', yerr=sem)

    means = (df.mean()['e0'], df.mean()['e1'], df.mean()['e2'], df.mean()['e2c'], df.mean()['tri'], df.mean()['p3'], df.mean()['star'], df.mean()['tritail'], df.mean()['square'], df.mean()['squarediag'], df.mean()['k4'])
    sem = (df.sem()['e0'], df.sem()['e1'], df.sem()['e2'], df.sem()['e2c'], df.sem()['tri'], df.sem()['p3'], df.sem()['star'], df.sem()['tritail'], df.sem()['square'], df.sem()['squarediag'], df.sem()['k4'])
    rects = ax.bar(ind+width+.02, means, width-.02, color = 'b', yerr=sem)

    plt.ylim(ymin=0)
    #fig = plt.gcf()
    #fig.set_size_inches(5, 3, forward=True)
    fig.savefig(os.path.join(os.getcwd(), 'graphlet_plot.png'))


def draw_degree_rank_plot(orig_g, mG):
    ori_degree_seq = sorted(dict(nx.degree(orig_g)).values(), reverse=True)  # degree sequence
    deg_seqs = []
    for newg in mG:
        deg_seqs.append(sorted(dict(nx.degree(newg)).values(), reverse=True))  # degree sequence
    df = pd.DataFrame(deg_seqs)
    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.fill_between(df.columns, df.mean() - df.sem(), df.mean() + df.sem(), color='blue', alpha=0.2, label="se")
    h, = ax.plot(df.mean(), color='blue', aa=True, linewidth=4, ls='--', label="H*")
    orig, = ax.plot(ori_degree_seq, color='black', linewidth=4, ls='-', label="H")

    ax.set_title('Degree Distribution')
    ax.set_ylabel('Degree')
    ax.set_ylabel('Ordered Vertices')

    ax.xaxis.set_tick_params(
        # axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom='off',  # ticks along the bottom edge are off
        top='off',  # ticks along the top edge are off
        labelbottom='off')  # labels along the bottom edge are off

    ax.legend([orig, h], ['$H$', 'HRG $H^*$'], loc=3)
    # fig = ax.gcf()
    # fig.set_size_inches(5, 4, forward=True)
    fig.savefig(os.path.join(os.getcwd(), 'degree_rank_plot.png'))


def draw_network_value(orig_g, mG):
    """
    Network values: The distribution of eigenvector components (indicators of "network value")
    associated to the largest eigenvalue of the graph adjacency matrix has also been found to be
    skewed (Chakrabarti et al., 2004).
    """
    eig_cents = [nx.eigenvector_centrality(g) for g in mG]  # nodes with eigencentrality

    # srt_eig_cents = sorted(eig_cents, reverse=True)
    net_vals = []
    for cntr in eig_cents:
        net_vals.append(sorted(cntr.values(), reverse=True))
    df = pd.DataFrame(net_vals)

    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.fill_between(df.columns, df.mean() - df.sem(), df.mean() + df.sem(), color='blue', alpha=0.2, label="se")

    h, = ax.plot(df.mean(), color='blue', aa=True, linewidth=4, ls='--', label="H*")
    orig, = ax.plot(sorted(nx.eigenvector_centrality(orig_g).values(), reverse=True), color='black', linewidth=4,
                     ls='-', label="H")

    ax.set_title('Principle Eigenvector Distribution')
    ax.set_ylabel('Principle Eigenvector')
    ax.xaxis.set_tick_params(
        # axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom='off',  # ticks along the bottom edge are off
        top='off',  # ticks along the top edge are off
        labelbottom='off')  # labels along the bottom edge are off

    ax.legend([orig, h], ['$H$', 'HRG $H^*$'], loc=3)
    # fig = plt.gcf()
    # fig.set_size_inches(5, 4, forward=True)
    fig.savefig(os.path.join(os.getcwd(), 'eigen_plot.png'))


def draw_hop_plot(orig_g, mG):
    m_hops_ar = []
    for g in mG:
        c = get_graph_hops(g, 20)
        d = dict(c)
        m_hops_ar.append(d.values())
        print("H* hops finished")

    df = pd.DataFrame(m_hops_ar)

    ## original plot
    c = get_graph_hops(orig_g, 20)
    dorig = dict(c)
    fig, ax = plt.subplots()
    ax.fill_between(df.columns, df.mean() - df.sem(), df.mean() + df.sem(), color='blue', alpha=0.2, label="se")
    h, = ax.plot(df.mean(), color='blue', aa=True, linewidth=4, ls='--', label="H*")
    orig, = ax.plot(dorig.values(), color='black', linewidth=4, ls='-', label="H")
    ax.set_title('Hop Plot')
    ax.set_ylabel('Reachable Pairs')
    ax.set_xlabel('Number of Hops')
    # ax.ylim(ymax=max(dorig.values()) + max(dorig.values()) * .10)

    ax.legend([orig, h, ], ['$H$', 'HRG $H^*$'], loc=1)
    #fig = plt.gcf()
    #fig.set_size_inches(5, 4, forward=True)
    fig.savefig(os.path.join(os.getcwd(), 'hop_plot.png'))
    