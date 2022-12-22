# Probe Points

This repository consists of a small program which is helpful for analysing where the electron density is most and least common from various PDB models. This work was initially started to improve the probe points that were in Nautilus (Cowtan 2014).

### Development

#### Prerequisites

You must have
- Clipper
- MMDB2

installed, which come included in the CCP4. To ensure they are in your path, you must source the appropriate script. To do this run:

    source /opt/xtal/ccp4-X.X/bin/ccp4.setup-sh 
where X.X is your CCP4 version.

#### Development

To compile this, just simply run:

    make
and the executable 'probe' should be created in the root directory of the project.

To run the executable run:

    ./probe
    
