# 3D-to-2D Slice Projection Tool for Nek5000 Data

This MATLAB script extracts a 2D slice from a 3D flow field file (e.g., `WING0.f00001`) and maps it onto the mesh of a 2D dummy file (e.g., `BASvortex_dir0.f00001`).  
The slice is taken in the **(Y, Z)** plane and assumed to be **extruded along the X-axis**, enabling simplified visualization or comparative analysis of 3D results on a 2D geometry.

- Reads 3D Nek5000 field files (`.f00001`)  
- Extracts a specified 2D slice (Y–Z plane)  
- Projects the slice onto a 2D dummy mesh  
- Outputs a Nek-compatible 2D field file  

## Dependencies
This script relies on MATLAB functions from the open-source repository  
[**eX-Mech/nekmatlab**](https://github.com/eX-Mech/nekmatlab)

Specifically, the following functions are used:
- `readnek.m` – for reading Nek5000 field files  
- `writenek.m` – for writing modified Nek5000 field files  

Full credit for these functions goes to the original authors of **nekmatlab**.

