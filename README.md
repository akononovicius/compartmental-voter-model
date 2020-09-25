# Compartmental voter model

Here you can find C source code implementing the model proposed in
[1]. Idea behind the
model is quite simple: agents of different types switch between compartments,
which are of finite capacity. Switching rate is influenced by the number of
agents with the same type in the destination compartment. This makes the model
similar to the well-known
[voter model](http://rf.mokslasplius.lt/rinkejo-modelis/) [2].

Model source code is stored in the "model" folder, while in the root folder you
can find a simple example, which runs the simulation and prints out the time
series generated by the model. Example is based upon Leicester example for the
paper.

Anyone may use the source code for any purpose as long as the paper [1]
is appropriately referenced.

## References

1. A. Kononovicius. *Compartmental voter model*.
Journal of Statistical Mechanics 2019: 103402 (2019).
doi: [10.1088/1742-5468/ab409b](https://dx.doi.org/10.1088/1742-5468/ab409b).
[arXiv: 1906.01842 [physics.soc-ph]](https://arxiv.org/abs/1906.01842).
2. P. Clifford, A. Sudbury. *A model for spatial conflict*. Biometrika **60**:
581-588 (1973). doi: [10.1093/biomet/60.3.581](https://dx.doi.org/10.1093/biomet/60.3.581).
