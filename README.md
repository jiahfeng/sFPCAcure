# sFPCAcure

Simulation code of implementing supervised functional principal component analysis (sFPCA) under the mixture cure rate model

The repository includes the following files:

`fun.R`: The R source code for fitting the sFPCA for imaging data, which contains the following key functions.

`sfpca_img = function(type = "bernstein", Y , train_dat.id, theta, lambda, npc, V.est, Tr.est, d.est, r, Z, ncr, tau)`

`simu_curerate.R`: A demo script for a simulation study of the sFPCA for imaging data.

-   Tr.est: coordinates of triangles used to set up the simulation study.

-   V.est: vertices of the triangles.

-   Z: expanded grid for the rectangular box containing the object with an irregular boundary.

## Contact

Jiahui Feng  <<jhfeng731@gmail.com>>

## Reference

Feng, J., Shi, H., Ma, D., Beg, M. F., and Cao, J., Supervised functional principal component analysis under the mixture cure rate model: an application to Alzheimer's disease.
