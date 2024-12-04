In this file I'll try to keep traces of my ideas concerning the parallel genetic algorithm
1. add the min size clusters, that should be settled to maxSierCluster/2
2. add other subpopulations with other encoding schemes such as VAE, FC, PMEB
3. do migrate a set of solutions instead of one

## Assignment Algorithm
* Step 1 : Create clusters fro each median vertex
* Step 2 : For each non median vertices, check if it has a link with a cluster/clusters 
(direct link to a median vertex/with already assinged vertices), then choose among them (clusters) the one having the greatest weight
By summing the weights of all the links with each cluster.
* Step 3 : The remaining vertices are assigned according to the following situations :
* * It has neighbors already assigned in the step 2, in this situation we follow the same principle as followed in step 3
* * all its neighbors are assinged and their corresponding clusters are all full, here we chose an non-full cluster and we assign the vertex to it
* * it has non assigned neighbors, we follow the same process as the previous point

