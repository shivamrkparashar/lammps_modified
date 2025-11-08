/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(mcemc,FixMCEMC);
// clang-format on
#else

#ifndef LMP_FIX_MCEMC_H
#define LMP_FIX_MCEMC_H

#include "fix.h"
#include "fix_gcmc.h"

namespace LAMMPS_NS {

class FixMCEMC : public FixGCMC {
 public:
  FixMCEMC(class LAMMPS *, int, char **);    // constructor

 protected:
  int &nmcemc_type;           // reference to base class type member
  int &mcemc_nmax;            // reference to base max atoms
  double gaugecell_volume;    // volume of the gauge cell
  int ntotal;                 // total number of atoms in system

  // protected methods
  void attempt_atomic_deletion() override;           // attempt to delete a single atom
  void attempt_atomic_insertion() override;          // attempt to insert a single atom
  void attempt_molecule_deletion() override;         // attempt to delete a molecule
  void attempt_molecule_insertion() override;        // attempt to insert a molecule
  void attempt_atomic_deletion_full() override;      // delete atom with full energy calculation
  void attempt_atomic_insertion_full() override;     // insert atom with full energy calculation
  void attempt_molecule_deletion_full() override;    // delete molecule with full energy calculation
  void
  attempt_molecule_insertion_full() override;    // insert molecule with full energy calculation
};

}    // namespace LAMMPS_NS

#endif
#endif
