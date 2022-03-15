# Engineering Attribute-Based Encryption Schemes Code Base
This repository contains the code for the authors' thesis project.

## Licensing!!
This is a work in progress project and is only visible for the sake of ease between group members. Proper use of licensing is yet to be added to the code base, particularly concerning legacy code that is directly taken from https://github.com/abecryptools/abe_squared, which in turn is taken from https://github.com/zeutro/openabe/.

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


## Authors
 * Benjamin B. Hansen
 * Viktor H. Miltersen
 * Jonas H. Salomonsson

## Supervisor
  * Diego F. Aranha
