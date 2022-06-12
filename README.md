# Engineering Attribute-Based Encryption Schemes Code Base
This repository contains the code for the authors' thesis project.

## Installation Instructions
  * Clone project ``git clone https://github.com/Osolat/thesis_code.git``

### Dependencies
The current version of the code requires Relic (https://github.com/relic-toolkit/relic) to be installed in /usr/lib, although that could be easily changed. Use latest commit, not the latest release. We have only tested on curve bls12-381. Therefore, follow and install Relic using preset ``x64-pbc-bls12-381.sh``

#### Prequisites for installing Relic
Follow guidelines from Relic, but here's an reiteration of important dependencies:
  * Install GMP
    * Download .lz file from https://gmplib.org and unzip
    * Install through https://gmplib.org/manual/Installing-GMP
    * May require m4 and flex.``apt update && apt upgrade`` then ``sudo apt install libfl-dev``to install flex. ``sudo apt-get install m4`` to install m4
  * Install gmp dev lib for some gmpxx.h header file needed in relic
    * ``sudo apt-get install libgmp-dev``       

## Running
For large number of attributes the policy tree will need a lot of stack space due to space allocation of Relic

``ulimit -s unlimited && objects/main 10000``

## Helpful notes about the policy tree implementation
- Recover coefficients reverses the order of children. A flat tree will have leaf n at pos 0, leaf n-1 at pos 1 etc. 
- In the n-standard tree the recursion goes down to bottom of tree first.
- The n-standard tree goes to the right
- Use leaf_index to index the coefficients properly.


## Authors
 * Benjamin B. Hansen
 * Viktor H. Miltersen
 * Jonas H. Salomonsson

## Supervisor
  * Diego F. Aranha
