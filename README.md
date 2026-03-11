# PdFM-Abaqus-FEM


Finite element model used for the electromechanical nanoindentation simulations reported in:

**Cho et al., "Nanoscale frictional imaging of ferroelectric domains."**

This repository contains the Abaqus finite element model used to reproduce the electromechanical simulations reported in the paper. The code analyzes the **polarity-dependent indentation response of ferroelectric domains** in polarization-derived friction microscopy (PdFM).

---

## 🔬 Physical Background

Polarization-derived friction microscopy (PdFM) visualizes ferroelectric domains through friction contrast measured during AFM scanning.

When an AFM tip indents a ferroelectric surface, a highly localized **strain gradient** is generated beneath the tip.  
This strain gradient induces **flexoelectric polarization**, which interacts with the intrinsic **ferroelectric polarization** through electromechanical coupling.

As a result:

- In **down-polarized domains**, flexoelectric and ferroelectric polarization are aligned.
- In **up-polarized domains**, the two contributions partially oppose each other.

This difference leads to **polarity-dependent electromechanical stiffness**, producing measurable differences in indentation response such as:

- **Indentation force**
- **Contact area**
- **Effective stiffness**

These mechanical differences contribute to the **friction asymmetry observed in PdFM experiments**.

The polarity-dependent indentation response arises from the competition between **piezoelectric polarization** and **strain-gradient–induced flexoelectric polarization** beneath the AFM tip.

The present FEM simulations were performed to quantitatively analyze this mechanism.

---

## ⚙️ Model Description

The simulations model **axisymmetric nanoindentation of a ferroelectric half-space by a rigid spherical indenter**.

The continuum formulation includes the following coupled physical effects:

- **Linear elasticity**
- **Dielectric response**
- **Piezoelectric coupling**
- **Flexoelectric coupling**
- **Strain-gradient elasticity**

The governing equations are implemented in **Abaqus using a user-defined element (UEL)**.

Two domain polarizations are simulated:

| Domain | Subroutine |
|------|------|
| Up-polarized domain | `Subroutine_UEL_up.for` |
| Down-polarized domain | `Subroutine_UEL_down.for` |

The domain polarity is introduced by **reversing the sign of the piezoelectric tensor**, while keeping all other material parameters unchanged.

---

## 📂 Repository Structure


```text
PdFM-Abaqus-FEM
│
├── abaqus_model
│   ├── Abaqus_inp.inp
│   ├── Subroutine_UEL_up.for
│   └── Subroutine_UEL_down.for
│
├── LICENSE
└── README.md
```

---

## 📁 Files

### `Abaqus_inp.inp`

Abaqus input file defining the **axisymmetric indentation model**, including:

- mesh definition  
- boundary conditions  
- contact setup  
- material parameters  
- user element call  

### `Subroutine_UEL_up.for`

User-defined element (UEL) implementation for the **up-polarized ferroelectric domain**.

### `Subroutine_UEL_down.for`

User-defined element (UEL) implementation for the **down-polarized ferroelectric domain**.

The two subroutines differ only in the **sign of the piezoelectric coupling**, representing opposite ferroelectric polarization states.

---

## ⚙️ Simulation Setup

The indentation model corresponds to the configuration described in the paper.

**Key features:**

- Axisymmetric finite element model  
- Rigid spherical indenter  
- Nanoscale indentation depth  
- Hard contact in the normal direction  
- Frictionless tangential contact  

The simulations evaluate:

- **Indentation force–depth relationship**
- **Polarity-dependent electromechanical stiffness**
- **Differences in contact mechanics between ferroelectric domains**

---

## ▶️ Running the Simulation

The simulations are executed in **Abaqus with a user-defined element (UEL)**.

Example command:

```bash
abaqus job=example_job input=Abaqus_inp.inp user=Subroutine_UEL_up.for
