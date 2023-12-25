# Estimating the Number of Morphological Rate Partitions in a Phylogenetic Tree

Professor Steve Wang
Student Researchers: Joyce Ben, Zachary Potthoff, Pradip Sharma Poudel, Anhad Singh

## Abstract

In a phylogenetic tree, each branch can be thought of as having a true underlying rate of evolution. These rates are unknown, but they may be estimated from the observed number of character changes and the estimated duration of each branch. Given a tree with NB branches, there may be as few as one rate of evolution (if all branches have the same rate) or as many as NB rates (if each branch has a different rate).

Here, our goal is to estimate the number of distinct rates of evolution in the tree when the observed data are discrete counts of character changes. First, we describe an exhaustive search algorithm that examines each possible partition of k rates assigned to contiguous non-overlapping regions of the tree, where k ranges from 1 to NB. For each possible partition, we calculate how well it fits the observed data using AIC. The partition with the smallest AIC provides an estimate for the number of distinct rates, their magnitudes, and their corresponding regions of the tree. This exhaustive search algorithm is guaranteed to find the best-fitting partition, but it is impractically slow for trees with large numbers of tips (i.e., approximately 20 or more).

We next describe two fast algorithms that can be applied to large trees: a forward (splitting) algorithm that starts by assuming that all branches have the same rate and then checks whether adding additional rates improves the fit, and a backwards (merging) algorithm that starts by assuming each branch has a different rate and checks whether merging contiguous branches improves the fit. We compare these fast algorithms with the exhaustive search algorithm and assess their performance on simulated datasets. Finally, we apply our methods to a dataset of lungfish fossils to better understand their evolutionary dynamics.
