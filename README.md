# randstab
Generate random stabilizer quantum states explicitly.

## Dependencies
* [Numpy](https://numpy.org/) – the fundamental package for scientific computing with Python.

## Usage

To start sampling n-qubit random stabilizer states simply import the module and call the function `random_stabilizer_state(n)`. It will produce a Numpy array with 2<sup>n</sup> complex elements. For example:
```
import randstab as rs
rs.random_stabilizer_state(3)
```
