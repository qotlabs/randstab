# randstab
Generate random stabilizer quantum states explicitly.

The proof of the implemented algorithm is presented in the article G.I. Struchalin _et al._ "Experimental Estimation of Quantum State Properties from Classical Shadows", [PRX Quantum 2, 010307 (2021)](https://link.aps.org/doi/10.1103/PRXQuantum.2.010307), in Appendix A.

## Dependencies
* [Numpy](https://numpy.org/) – the fundamental package for scientific computing with Python.

## Usage

To start sampling n-qubit random stabilizer states simply import the module and call the function `random_stabilizer_state(n)`. It will produce a Numpy array with 2<sup>n</sup> complex elements. For example:
```
import randstab as rs
rs.random_stabilizer_state(3)
```
