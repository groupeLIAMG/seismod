# seismod
C codes for seismic modeling

*Warning* These codes are in alpha form and have not been thoroughly tested


## 2D codes

1. e_iso
   * Elastic wave propagation in 2D isotropic media, on a cartesian staggered grid
   * FDTD scheme
   * split pml
   * OpenMP

2. e_cyl_iso
   * Elastic wave propagation in 2.5D isotropic media, on a staggered grid in
     cylindrical coordinates with vertical symmetry axis
   * FDTD scheme
   * CPML
   * OpenMP

3. ve_vti_pml
   * ViscoElastic wave propagation in 2D VTI anisotropic media, for L
     attenuation quasi-dilatational mechanism and 1 shear relaxation mechanism,
     on a cartesian staggered grid
   * PSTD scheme, 4th order Runge-Kutta time stepping
   * CPML
   * OpenMP enabled fftw

4. ve_vti_sh_pml
   * ViscoElastic SH-wave propagation in 2D VTI anisotropic media, with 1 shear
     relaxation mechanism, on a cartesian staggered grid
   * PSTD scheme, 4th order Runge-Kutta time stepping
   * CPML
   * OpenMP enabled fftw

5. pve_iso_pml
   * PoroViscoElastic wave propagation in 2D isotropic media, for one attenuation
     mechanism, on a cartesian staggered grid
   * PSTD scheme, 4th order Runge-Kutta time stepping
   * CPML
   * OpenMP enabled fftw

6. pve_vti_pml
   * PoroViscoElastic wave propagation in 2D VTI anisotropic media, for one
     attenuation mechanism, on a cartesian staggered grid
   * PSTD scheme, 4th order Runge-Kutta time stepping
   * CPML
   * OpenMP enabled fftw



## 3D code

1. a_iso_3d
   * Acoustic wave propagation in 3D isotropic media, on a cartesian staggered grid
   * FDTD scheme
   * CPML
   * OpenMP
