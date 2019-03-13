# Telling ecological networks apart by their structure: a computational challenge

This repository contains code and data to replicate and expand upon results
published in "Telling ecological networks apart by their structure: a
computational challenge" (In Revision). Want to jump right in? [Click here to
skip to a quick demo.](#a-demonstration)

In the aformentioned work [Stefano Allesina](http://allesinalab.uchicago.edu/)
and myself attempted to classify ecological bipartite networks using size- and
connectance-independent measures of network structure. While we found that such
classification is in principal possible and in practice relatively easy for
non-biological data sets, we were unable to achieve similar success with a large
collection of ecological network data that has been published over the past
hundred years.

Naturally, there are a number of possible explanations for this lack of success,
ranging from an unsurpassable level of variation in ecological data to a failure
of metric/method selection. Having attempted to solve this ourselves for several
years, we now turn to the scientific community writ large&mdash;from the computer
scientists with the most state-of-the-art techniques to the field ecologists
with the most insight into biological relevance&mdash;to pick up the gauntlet.

We have sought to clearly enunciate the problem as well as the relevant aspects
of a successful solution and provide here a framework for the easy application
of new metrics and methods on a database of unprecedented scale and conformity.

# Repository Organization

There are four main folders in this repository:

### [`/Data`](Data)

This folder contains the processed data collected thus far organized into the
original data files in [`/Data/edgelists`](Data/edgelists). Subdivided into
gross categories, each data file is a plain text, `.csv` file containing three
columns: the first two form an edgelist, each row signifying a link between the
node indicated in the first column with the node indicated in the second column.
The third column indicates the strength of this link. This will be either a real
number (if provided) or 1 in the case of unweighted raw data.

For each empirical data file, we attempted to produce two randomizations. First,
we produced a constrained Erdős–Rényi randomized which preserved the number of
nodes and links in the original network as well as the fact that the network
forms a single connected component. These randomizations are stored in
[`/Data/edgelists_er`](Data/edgelists_er). Second, we attempted a configuration
model randomization of each network, producing a random graph which preserves
the degree distribution of the original network, while again confirming a single
connected component. Being a more restrictive randomization, we were less
successful at successfully generating these graphs, in particular for the
`movies` subcategory. Those we were able to generate are stored in
[`/Data/edgelists_cm`](Data/edgelists_cm).

Finally, we recorded various meta-data for each network collected. These are
stored in a master file entitled [`/Data/Metadata.csv`](Data/Metadata.csv). The
`name` and `type` columns in this file act as the link between each row of
metadata and its corresponding edgelist stored in
`/Data/edgelist/[type]/[name].csv`.

### [`/Code`](Code)

This Folder contains the code required to run the analysis, subdivided into
three stages. The first, [`/Code/Processing`](Code/Processing), will be
unnecessary unless you are adding additional data prior to analysis. The code in
this folder was used to generate the randomized versions of each empirical
network.

The second, [`/Code/Analysis`](Code/Analysis) contains the code needed to
calculate various structural metrics for each of the data files as well as some
code for ancillary analyses reported in the associated publication.

We have included all of the metrics we calculated in the [`/Results`](Results)
folder already, but if you wanted to re-run these calculations or wanted to add
new metrics to the list, the key scripts for this are
[`full_analysis_of_one_file.R`](/Code/Analysis/full_analysis_of_one_file.R) and
[`run_analysis_over_all_files.R`](/Code/Analysis/run_analysis_over_all_files.R).
The latter is a wrapper for the former, running it for each data file in a
directory. The former takes in an edgelist file (such as one of the files in the
[`/Data`](Data) directory) and calculates a number of network-structural
metrics, returning a corresponding output file in the [`/Results`](Results)
directory. This file will consist of a single row table with each column holding
the value for a different metric.

Finally, [`/Code/Plotting`](Code/Plotting) contains the code needed to produce
the figures presented in the associated publication.
<!-- TODO: continue -->

### [`/Figures`](Figures)

This folder, now empty except for the [`full_pca_simple_demonstration.pdf`](/Figures/full_pca_simple_demonstration.pdf), which results from the simple demonstration outlined below, will be the default output folder for figures generated
by the plotting scripts mentioned above.
<!-- TODO: continue -->

### [`/Results`](Results)

Running the analysis as instructed results in the generation of a one-line
`.csv` result file for each network in the data folder. As mentioned above,
this folder is pre-filled with the metrics we calculated for the paper.

# A Demonstration

The data files and their randomizations are already provided in the /Data
directory. To get a feel for the analysis, we will start with a simple
demonstration using the metrics reported in the associated publication.

To start, we will calculate a few structural metrics for each network: the two
rightmost (largest real part) eigenvalues (λ₁ and λ₂) and
three analytical estimates for these values (a configuration model λ₁, an Erdős–Rényi λ₁,
and a Marchenko–Pastur λ₂).

The code
[`run_simple_demonstration_analysis.R`](/Code/Analysis/run_simple_demonstration_analysis.R)
performs these calculations and produces an output file entitled
[`SimpleDemo_results.csv`](/Results/SimpleDemo_results.csv) in the /Results
folder.

Our hope is that these metrics will allow us to partition networks by type, and
we can visualize the relative success of this using principal component analysis
(PCA). For this, we turn to the code in the [`/Code/Plotting`](Code/Plotting)
directory.

The key figure for this can be generated using
[`plot_pca_simple_demonstration.R`](/Code/Plotting/plot_pca_simple_demonstration.R)
which will save a new `.pdf` file into the [`/Figures`](Figures) folder entitled
`pca_simple_demonstration.pdf`.

 Other figures can be generated using the
 [`plot_pca_results.R`](/Code/Plotting/plot_pca_results.R) script, which
 constructs the pca space more generally and then calls the relevant
 figure-specific code to generate each of the first three figures of the main
 text, as well as their analogues in the supplementary information. The data and
 figures associated with the Nestedness/Modularity component of the work can be
 generated using the
 [`nest_mod_comparison.R`](/Code/Plotting/nest_mod_comparison.R) and
 [`pca_only_mod_nest.R`](/Code/Plotting/pca_only_mod_nest.R) scripts, which will
 be being cleaned and annotated over the coming days.
