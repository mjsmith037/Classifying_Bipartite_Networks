# This repository is a work-in-progress, please have patience as we continue working to improve this documentation

# Classifying_Bipartite_Networks

This repository contains code and data to replicate and expand upon results
published in XXXXXX. Want to jump right in? [Click here to skip to a quick demo.](a-demonstration)

In that work [Stefano Allesina](http://allesinalab.uchicago.edu/) and myself
attemted to classify ecological bipartite networks using size- and
connectance-independent measures of network structure. While we found that such
classification is in principal possible and in practice relatively easy for
non-biological data sets, we were unable to achieve similar success with the
large collection of ecological network data that has been published over the
past hundred years.



As such, we cast this




# Repository Organization

There are four main folders in this repository:

### [`/Data`](Data)

This folder contains the processed data collected thus far organized into the original data files in
[`/Data/edgelists`](Data/edgelists). Subdivided into gross categories, each data file is a plain text,
`.csv` file containing three columns: the first two form an edgelist, each row signifying a link
between the node indicated in the first column with the node indicated in the second column. The
third column indicates the strength of this link. This will be either a real number (if provided) or
1 in the case of unweighted raw data.

For each empirical data file, we attempted to produce two randomizations. First, we produced a
constrained Erdős–Rényi randomized which preserved the number of nodes and links in the original
network as well as the fact that the network forms a single connected component. These
randomizations are stored in [`/Data/edgelists_er`](Data/edgelists_er). Second, we attempted a
configuration model randomization of each network, producing a random graph which preserves the
degree distribution of the original network, while again confirming a single connected component.
Being a more restrictive randomization, we were less successful at successfully generating these
graphs, in particular for the `movies` subcategory. Those we were able to generate are stored in
[`/Data/edgelists_cm`](Data/edgelists_cm).

Finally, we recorded various meta-data for each network collected. These are stored in a master file
entitled [`/Data/Metadata.csv`](Data/Metadata.csv). The `name` and `type` columns in this file act as the link between each row of metadata and its corresponding edgelist stored in `/Data/edgelist/[type]/[name].csv`.

### [`/Code`](Code)

This Folder contains the code required to run the analysis, subdivided into three stages. The first, [`/Code/Processing`](Code/Processing), will be unnecessary unless you are adding additional data prior to analysis. The code in this folder was used to generate the randomized versions of each empirical network.

The second, [`/Code/Analysis`](Code/Analysis) contains the code needed to calculate various structural metrics for each of the data files as well as some code for ancillary analyses reported in XXXXXX. CONTINUE

Finally, [`/Code/Plotting`](Code/Plotting) contains the code needed to produce the figures presented in XXXXXX. CONTINUE

### [`/Figures`](Figures)

This folder, now empty, will be the default output folder for figures generated by the plotting scripts mentioned above.

### [`/Results`](Results)

Running

# A Demonstration
