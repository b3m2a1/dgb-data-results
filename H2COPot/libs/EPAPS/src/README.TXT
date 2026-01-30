=====================README.TXT ===========================

Authors: A. Yachmenev, S.N. Yurchenko, P. Jensen, and W. Thiel

Title:  A new `spectroscopic' potential energy surface for formaldehyde in its ground electronic state

Total No. of Files: 10

File Types: 	ASCII, Excel

Description: 'Spectroscpoic' PES of H2CO (H2CO-2011) (ground electronic state) 
              and the experimental energies of H2CO comprising our fitting set used in the PES refinements

1. README.TXT.............. This file.

2. abinitpes_h2co.f90.......Fortran 90 routine for calculating the ab initio (CCSD(T)/aug-cc-pVQZ) 
                            potential energy values in combination with the input files 
                            H2CO_ABINIT_PES and h2co_sampleinp.

3. refpes_h2co.f90..........Fortran 90 routine for calculating the refined potential energy values 
                            in combination with the input files H2CO_REF_PES_2011 and h2co_sampleinp.

4. H2CO_ABINIT_PES..........Input file for abinitpes_h2co.f90 containing the potential parameters 
                            defining the ab initio (CCSD(T)/aug-cc-pVQZ) PES of H2CO.

5. H2CO_REF_PES_2011........Input file for refpes_h2co.f90 containing the potential parameters 
                            defining H2CO-2011 (see Eq. 1 of the paper).

6. h2co_sampleinp...........Sample input file for abinitpes_h2co.f90 and refpes_h2co.f90 
                            containing the grid points.

7. h2co_sampleout_abinit....Output produced by abinitpes_h2co.f90 containing the ab initio 
                            potential energy values of the grid.

8. h2co_sampleout_ref.......Output produced by refpes_h2co.f90 containing the refined 
                            potential energy values of the grid.

9. Fitting_set_H2CO.xls ... Fitting set: The set of experimental energies (cm-1) of H2CO for J=0,1,2,3,4, 
                            and 5 selected for our refinements. Fitting_set_H2CO.xls also containes 
                            the corresponing calculated term values as well as the Obs.-Calc. deviations.

10. Predicted_energies.txt. A set of the predicted energies (cm-1) of H2CO for J=0,1,2,3,4,5, and 6 
                            below 8000.0 cm-1 (relative to the zero point energy) using the H2CO-2011 PES. 


Contacts: 
       W. Thiel thiel@mpi-muelheim.mpg.de  
       and
       A. Yachmenev, andrey@mpi-muelheim.mpg.de 
       S.N. Yurchenko, s.yurchenko@chemie.tu-dresden.de
       P. Jensen jensen@uni-wuppertal.de

=====================README.TXT ===========================
