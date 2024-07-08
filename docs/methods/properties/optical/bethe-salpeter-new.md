# Bethe-Salpeter Equation

$\text{\color{red}Todo: Referenzen alphabetisch}$

In this section, we discuss the basics for computing optical properties of molecules using
the Bethe-Salpeter equation (BSE) in CP2K [Graml2024b](#Graml2024b). The BSE enables the computation of 
electronic excitation energies and optical absorption spectra, for a review, see [[Blase2018](#Blase2018), [Blase2020](#Blase2020), [Bruneval2015](#Bruneval2015), [Sander2015](#Sander2015)]. 
In this howto, we describe in Sec. [1](#theory-and-implementation-of-bse) the Theory and implementation of BSE, in Sec. [2](#bse-input) the BSE keywords and in Sec. [3](#minimal-example-for-a-bse-calculation) the input and the output of a BSE calculation.



## 1. Theory and implementation of BSE

A central goal of a BSE calculation is to compute electronic excitation energies $\Omega^{(n)}, n=1,2,\ldots$ (cf. Refs. [[Blase2018](#Blase2018), [Blase2020](#Blase2020), [Graml2024b](#Graml2024b)] for more usecases and details).

The following ingredients are necessary for computing $\Omega^{(n)}$:
- Occupied Kohn-Sham (KS) orbitals $\varphi_i(\mathbf{r})$ and empy KS orbitals $\varphi_a(\mathbf{r})$ from a DFT calculation, where $i=1,N_{occ}$ and $a=N_{occ}+1,N_{occ}+N_{empty}$,
- $GW$ eigenvalues $\varepsilon_i^{GW}$ and $\varepsilon_a^{GW}$ of corresponding KS orbitals

In CP2K, we use $G_0W_0$ eigenvalues, see details in [GW] and in Ref. [[Golze2019](#Golze2019)], i.e. we perform BSE@$G_0W_0$@DFT (see full input in Sec. [3.1](#input-file)).

We obtain optical properties in BSE from the generalized diagonalization of a block-matrix $ABBA$:

$$\left( \begin{array}{cc}A &  B\\B &  A\end{array} \right)\left( \begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array} \right) = \Omega^{(n)}\left(\begin{array}{cc}1&0\\0&-1\end{array}\right)\left(\begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array}\right) \quad .$$

We abbreviate $A$ and $B$ as matrices with index $A_{ia,jb}$, i.e. they have $N_{occ}N_{empty}$ rows and $N_{occ}N_{empty}$ columns. The entries of A and B are given by [[Blase2018](#Blase2018)]

$$ \begin{align}
    A_{ia,jb} &= (\varepsilon_a^GW-\varepsilon_i^GW)\delta_{ij}\delta_{ab} + \alpha^\mathrm{S/T}
    v_{ia,jb} - W_{ij,ab}(\omega=0) \quad ,\\
    B_{ia,jb} &= \alpha^\mathrm{(S/T)} v_{ia,bj} - W_{ib,aj}(\omega=0) \quad .
\end{align}$$
where $\delta_{ij}$ is the Kronecker delta. 
The user sets $\alpha^S=2$ for computing singlet excitations and $\alpha^T=0$ for computing triplet excitations. 
$v_{pq,rs}$ is the bare Coulomb interaction  and $W_{pq,rs}(\omega=0)$ the
statically ($\omega=0$) screened Coulomb interaction , where $p,q,r,s \in [ 1, N_{occ}+N_{empty}]$ are KS orbital indices.
$(\mathbf{X}^{(n)},\mathbf{Y}^{(n)})$ with elements $X_{ia}^{(n)}$ and $Y_{ia}^{(n)}$ are the eigenvectors which relate to the wavefunction of the electronic excitation [[Blase2020](#Blase2020)],

$$ \begin{align}
\Psi_\text{excitation}^{(n)}(\mathbf{r}_e,\mathbf{r}_h) = \sum_{ia} X_{ia}^{(n)} \varphi_i(\mathbf{r}_h) \varphi_a(\mathbf{r}_e) + Y_{ia}^{(n)} \varphi_i(\mathbf{r}_e) \varphi_a(\mathbf{r}_h) \quad ,
\end{align}$$
i.e. $X_{ia}^{(n)}$ and $Y_{ia}^{(n)}$ describe the transition amplitude between occupied orbital $\varphi_i$ and empty orbital $\varphi_a$ of the $n$-th excitation.



If the $B$ matrix is small, the Tamm-Dancoff approximation (TDA) can be applied which neglects the $B$ matrix. 
In case $A$ is positive definite, and excitation energies $\Omega^{(n)}>0$, we have $\mathbf{Y}=0$ and $\mathbf{X}$ can be computed from
$\color{red}\text{Referenz?}$

$$ A \mathbf{X}^{(n)}_\mathbf{TDA} = \Omega^{(n)}_\mathbf{TDA} \mathbf{X}^{(n)}_\mathbf{TDA} \quad .$$


Diagonalizing $A$ in TDA, or the full block-matrix $ABBA$, takes approximately $(N_{occ} N_{empty})^3$ floating point operations.
$\color{red}\text{Anders worden - ABBA Diagonalisierung braucht mindestens 3 Diags, also anderen Vorfaktor}$
This translates to a computational scaling of $O(N^6)$ in the system size $N$, e.g. the number of electrons. 


## 2. BSE input

% JW comments: First, describe BSE input parameters
Standard:
         
```
          &GW
            &BSE
              BSE_APPROX BOTH      
            &END BSE
          &END GW
```

Zahl: Genauigkeit zu aims (MAE angeben für TZ), zu best estimates, average über Moleküle und Excitation energies







The parameters defining the BSE calculation have to be specified in the [BSE]-subsection when
running a [RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) `ENERGY` calculation. As highlighted above, ensure
a converged DFT and [GW] calculation before. The most important keywords are:

- [BSE_APPROX](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_APPROX): Option to
  switch on/off the TDA, i.e. $B=0$. Can either be the full BSE, TDA or both.
- [SPIN_CONFIG](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.SPIN_CONFIG): Specifies
  the desired spin configuration of the excitation. This determines the scalar factor
  $\alpha^\mathrm{(S/T)}$, which is 2 for Singlet excitations and 0 for Triplet excitations.


```{note}
The accuracy of the BSE relies heavily on well-converged settings in the prior DFT and GW steps (BSE@$G_0W_0$@DFT). 
For example, the chosen [XC_FUNCTIONAL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.SECTION_PARAMETERS) and the parameters for the analytic continutation in GW ([QUADRATURE_POINTS](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.QUADRATURE_POINTS), [NPARAM_PADE](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.NPARAM_PADE) and [OMEGA_MAX_FIT](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.OMEGA_MAX_FIT)) can have a profound influence on the excitation energies.
In particular, all MO's included in the BSE have to be corrected in GW by setting [CORR_MOS_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.CORR_MOS_OCC) and [CORR_MOS_VIRT](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.CORR_MOS_VIRT) to a sufficiently large number. By default, the invocation of the BSE section corrects all orbitals.
```

Relevant parameters of the BSE calculation which crucially determine the computational time and the numerical precision: 
- basis set, recommondation: aug-cc-DZVP should be good, check with aug-cc-TZVP, all-electron, retrieve from EMSL database
- ENERGY_CUTOFF 
  - [ENERGY_CUTOFF_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_OCC):
    Cutoff for the occupied MO's. If the energy of MO with index $i$ is more than `ENERGY_CUTOFF_OCC`
    away from the HOMO (highest occupied molecular orbital) energy, i.e.
    $\epsilon_{\mathrm{HOMO}}-\epsilon_i>\mathtt{ENERGY\_CUTOFF\_OCC}$, then all occupied MO's with
    index $\le i$ are not included in the matrices $A_{ia,jb}$ and $B_{ia,jb}$.
  - [ENERGY_CUTOFF_VIRT](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_VIRT):
    Cutoff for the unoccupied MO's. If the energy of MO with index $a$ is more than
    `ENERGY_CUTOFF_VIRT` away from the LUMO (lowest unoccupied molecular orbital) energy, i.e.
    $\epsilon_a-\epsilon_{\mathrm{LUMO}}>\mathtt{ENERGY\_CUTOFF\_VIRT}$, then all unoccupied MO's with
    index $\ge a$ are not included in the matrices $A_{ia,jb}$ and $B_{ia,jb}$.

```{note}
Usage of [ENERGY_CUTOFF_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_OCC) and [ENERGY_CUTOFF_VIRT](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_VIRT) requires careful checks of convergence!
```

            &BSE
              ENERGY_CUTOFF_OCC    ! in eV
              ENERGY_CUTOFF_VIR    ! in eV
              BSE_APPROX BOTH           ! In this case, full BSE and TDA are calculated
            &END BSE

- [BSE_APPROX](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_APPROX)

The memory consumption of the BSE algorithm is large, it is approximately $100 Nocc^2Nvir^2$ B. You can see Nocc, Nvir and the required memory from the output file (heißt...). 
The BSE implementation is well parallelized, i.e. you can use several nodes that can fit the memory. 

The relevant model parameters that influence the BSE results: Starting DFT xc functional and GW self-consistency: recommondation for molecules: G0W0@PBE0 or evGW0@PBE, Verweis auf https://pubs.acs.org/doi/10.1021/acs.jctc.8b00014



## 3. Minimal example for a BSE calculation

### 3.1 Input file

In this section, we provide a minimal example on a BSE calculation on H$_2$, for the calculation you need the input file BSE_H2.inp ($\color{red}\text{here}$) and the aug-cc-DZVP basis ($\color{red}\text{here}$). 

Please copy both files into your working directory and run CP2K by

```none
mpirun -n 1 cp2k.psmp BSE_H2.inp
```

which requires 4.5 GB RAM and takes roughly 90 seconds on 1 core. You can download the output file $\color{red}\text{here}$.
$\color{red}\text{Poisson solver?}$
```none
&GLOBAL
  PROJECT  H2
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS-aug               ! Custom Basis set file (aug-cc-pVDZ and aug-cc-pVDZ-RIFIT from EMSL database)
    POTENTIAL_FILE_NAME POTENTIAL     
    &QS
      METHOD GAPW                               ! All electron calculation
      EPS_DEFAULT 1.0E-16
      EPS_PGF_ORB 1.0E-16
    &END QS
    &POISSON
      PERIODIC NONE
      PSOLVER MT
    &END
    &SCF
      SCF_GUESS RESTART
      EPS_SCF 1e-7
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE                        ! Choice of functional has a profound influence on the results
      &END XC_FUNCTIONAL
      &WF_CORRELATION
        &RI_RPA                                 ! In the RI_RPA and the GW section, additional numerical parameters, e.g.
          &GW                                   ! QUADRATURE_POINTS or NPARAM_PADE, can be specified.
            &BSE
              BSE_APPROX BOTH                   ! In this case, full BSE and TDA are calculated
            &END BSE
          &END GW
        &END RI_RPA
      &END WF_CORRELATION
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 20 20 20
      PERIODIC NONE
    &END CELL
    &COORD
      H 0.0000 0.0000 0.0000                    ! H2 molecule geometry from GW100 Paper
      H 0.0000 0.0000 0.74144
    &END COORD
    &TOPOLOGY
      &CENTER_COORDINATES
      &END
    &END TOPOLOGY
    &KIND H
      BASIS_SET ORB    aug-cc-pVDZ              ! For production runs, the basis set should be checked for convergence.
      BASIS_SET RI_AUX aug-cc-pVDZ-RIFIT        ! In general, pVDZ should be a solid choice.
      POTENTIAL ALL
    &END KIND
    &KIND O
      BASIS_SET ORB    aug-cc-pVDZ
      BASIS_SET RI_AUX aug-cc-pVDZ-RIFIT
      POTENTIAL ALL
    &END KIND
  &END SUBSYS
&END FORCE_EVAL

```

The basis sets `aug-cc-pVDZ` and `aug-cc-pVDZ-RIFIT` in `BASIS-aug` can be obtained from the Basis
Set Exchange Library:
<a href="https://www.basissetexchange.org/basis/aug-cc-pvdz/format/cp2k/?version=1&elements=1" target="_blank">aug-cc-pVDZ</a>,
<a href="https://www.basissetexchange.org/basis/aug-cc-pvdz-rifit/format/cp2k/?version=1&elements=1" target="_blank">aug-cc-pVDZ-RIFIT</a>.
The geometry for $H_2$ was taken from the [GW100](https://doi.org/10.1021/acs.jctc.5b00453)-Paper.

### 3.2 Output

In the resulting output file (cf. $\color{red}\text{here}$),
the [BSE]-section itself starts with a banner after the [GW]-section. Therein, all lines a formatted
with a trailing `BSE|`.

At first, characteristics of the BSE-run are printed, as described in Sec. [2](#2-bse-input), which can be used to estimate memory consumption and computational cost.

Afterwards, depending on the chosen
[BSE_APPROX](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_APPROX), a banner
signalizes the start of the respective results section. Afterwards, the most important formulas and
quantities are summarized before the excitation energies and the single-particle transitions, i.e.
the eigenvector elements $X_{ia}^{(n)}$, are printed up to the given
[EPS_X](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.EPS_X). For the full solution (no
TDA applied), the first lines of the output for the energies of the requested singlet excitation
look like

```none
 BSE| Excitation energies from solving the BSE without the TDA:
 BSE|
 BSE|     Excitation n        Multiplet  TDA/full BSE   Excitation energy Ω (eV)
 BSE|                1    Singlet State        -full-                    11.4669
 BSE|                2    Singlet State        -full-                    12.4840
 BSE|                3    Singlet State        -full-                    15.2848
```

with the columns excitation index $n$, the requested multiplet, the
[BSE_APPROX](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_APPROX) and the energy
in eV. The single-particle transitions are displayed like that (for the first three excitations):

```none
 BSE| Excitations are built up by the following single-particle transitions,
 BSE| neglecting contributions where |X_ia^n| <  0.10 :
 BSE|         -- Quick reminder: HOMO i =    1 and LUMO a =    2 --
 BSE|
 BSE| Excitation n      i =>     a               TDA/full BSE           |X_ia^n|
 BSE|
 BSE|            1      1 =>     2                     -full-             0.6682
 BSE|            1      1 =>     4                     -full-             0.2459
 BSE|
 BSE|            2      1 =>     3                     -full-             0.7060
 BSE|
 BSE|            3      1 =>     5                     -full-             0.7077
```

Here, some reminders for [EPS_X](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.EPS_X)
and the indices of HOMO and LUMO are printed. The columns display the excitation index $n$, the
single-particle indices $i$ and $a$ contributing to the $n$-th excitation, the
[BSE_APPROX](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_APPROX) and the absolute
value of the eigenvector entry $|X_{ia}^{(n)}|$.

In the case of the $\mathrm{H}_2$, the first excitation is mainly built up by a transition from the first MO
to the second MO, i.e. the LUMO, but also contains a considerable contribution from the 1=>4
(HOMO=>LUMO+2) transition. The remaining contributions of the normalized $\mathbf{X}^{(n)}$ are
smaller than `0.10` and are therefore not printed.

## 4. References
<a id="Graml2024b">[Graml2024b]</a> M. Graml, J. Wilhelm, Manuscript in preparation

<a id="Blase2018">[Blase2018]</a>  X. Blase, I. Duchemin, D. Jacquemin, The Bethe–Salpeter equation in chemistry: relations with TD-DFT, applications and challenges. Chem. Soc. Rev., 47, 1022. doi: 10.1039/c7cs00049a (2018)

<a id="Blase2020">[Blase2020]</a> X. Blase, I. Duchemin, D. Jacquemin, P.-F. Loos, The Bethe−Salpeter Equation Formalism: From Physics to Chemistry.  J. Phys. Chem. Lett., 11, 7371−7382. doi: 10.1021/acs.jpclett.0c01875 (2020)

<a id="Bruneval2015">[Bruneval2015]</a> F. Bruneval, S. M. Hamed, J. B. Neaton, A systematic benchmark of the ab initio Bethe-Salpeter equation approach for low-lying optical excitations of small organic molecules. J. Chem. Phys. 142, 244101. doi: 10.1063/1.4922489 (2015)

<a id="Sander2015">[Sander2015]</a> T. Sander, E. Maggio and G. Kresse, Beyond the Tamm-Dancoff approximation for extended systems using exact diagonalization. Physical Review B 92, 045209. doi: 10.1103/PhysRevB.92.045209 (2015)

<a id="Golze2019">[Golze2019]</a> D. Golze, M. Dvorak and P. Rinke, The GW Compendium: A Practical Guide to Theoretical Photoemission Spectroscopy. Front. Chem. 7:377. doi: 10.3389/fchem.2019.00377 (2019)


[bse]: #CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE
[gw]: #CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW
