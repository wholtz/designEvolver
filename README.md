# Design Evolver

Design Evolver is an optimizer that runs on top of the [j5](https://j5.jbei.org/) DNA assembly design software.  j5 (without Design Evolver) does not consider many of the possible ways to design a DNA assembly. Design Evolver gives j5 more degrees of freedom to find assemblies that are more likely to assemble successfully.  This workflow is particularly useful for creating difficult to assemble DNA such as, sequences with repeats or DNA coding for RNA with strong secondary structures.

## Depedencies

* [Python](https://www.python.org/) v3.3.x
* [Numpy](https://github.com/numpy/numpy) v1.7.1
* [DEAP](https://github.com/DEAP/deap/) v0.91
* [BioPython](https://github.com/biopython/biopython) v1.62
* [j5](https://j5.jbei.org/) v2.2.5 (not required if using a remote j5 server)

Newer versions of these dependencies may work but have not yet been tested by the maintainers.

## Installation:
1. If you plan to run j5 locally, then follow the installation instructions for [j5](https://j5.jbei.org/).
1. Clone this git repository to your local machine.
1. In the directory containing this README.md file, run: `pip3 install -r requirements.txt`

## Running Design Evolver
1. See example invocation in the [example_usage](/example_usage) file.
   * You must edit this file before using it as you must supply a j5 username and password.
1. See the [j5 manual](https://j5usermanual.lbl.gov/) for definitions of the input and output files.

## Limitations
Design Evolver does not support all j5 inputs. The following currently cannot be optimized by Design Evolver:
1. Combinatorial designs

## Bug reports and comments:
Please use the [associated issue tracker](https://github.com/wholtz/designEvolver/issues) on GitHub.

## Commercial Licensing

As noted in the [LICENSE](/LICENSE) file, commercial usage is not permitted without an additional license. To obtain a license for commercial usage contact the [University of California Tech Transfer Office](https://techtransfer.universityofcalifornia.edu/NCD/23597.html).
