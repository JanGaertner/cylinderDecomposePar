# Cylinder Decompose

Decompose the domain in cylinder coordinates by providing the 
number of axial, radial, and circumferential splits. 

## Usage

 In the system/decomposeParDict set:
 ```
    numberOfSubdomains 64;

    method         cylinder;

    cylinderCoeffs
    {
        n       (4 4 4);    // axial, radial, circumferential
        rotAxis z;          // rotation axis
    }
```

Provide the rotation axis either as word e.g., x,y, or z or alternative
provide a vector of the rotation axis
```
    ...
    rotAxis     (1 1 0);
```

## OpenFOAM Versions

Only works for OpenFOAM v2206 and above!


## Install
Source OpenFOAM and execute `wmake` in the src folder. 

