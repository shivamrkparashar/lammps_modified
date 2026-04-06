// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Shivam Parashar (Rutgers University)
------------------------------------------------------------------------- */

#include "fix_mcemc.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "citeme.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "region.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

// large energy value used to signal overlap

static constexpr double MAXENERGYSIGNAL = 1.0e100;

// this must be lower than MAXENERGYSIGNAL
// by a large amount, so that it is still
// less than total energy when negative
// energy contributions are added to MAXENERGYSIGNAL

static constexpr double MAXENERGYTEST = 1.0e50;

enum { EXCHATOM, EXCHMOL };          // exchmode
enum { NONE, MOVEATOM, MOVEMOL };    // movemode

static const char cite_fix_mcemc[] =
  "fix mcemc command: doi:10.1016/j.jcis.2024.06.083\n\n"
  "@Article{Parashar,\n"
  "author = {Parashar, S. and Neimark, A. V.},\n"
  "title = {Understanding the Origins of Reversible and Hysteretic Pathways \n"
  "         of Adsorption Phase Transitions in Metal-Organic Frameworks},\n"
  "journal = {Journal of Colloid And Interface Science},\n"
  "year = {2024},\n"
  "doi = {10.1016/j.jcis.2024.06.083},\n"
  "}\n\n";


/* ---------------------------------------------------------------------- */
// Derive from no-op constructor using initializer list
FixMCEMC::FixMCEMC(LAMMPS *lmp, int narg, char **arg) : FixGCMC(lmp, narg, arg, NoOpTag()),
  nmcemc_type(ngcmc_type),   // bind reference to base member
  mcemc_nmax(gcmc_nmax)      // bind reference to base member

{
  // check minimum # of arguments
  if (narg < 12) utils::missing_cmd_args(FLERR, "fix mcemc", error);

  // add citation to log file
  if (lmp->citeme) lmp->citeme->add(cite_fix_mcemc);

  if (atom->molecular == Atom::TEMPLATE)
    error->all(FLERR, Error::NOPOINTER, "Fix mcemc does not (yet) work with atom_style template");

  dynamic_group_allow = 1;

  vector_flag = 1;
  size_vector = 8;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  ngroups = 0;
  ngrouptypes = 0;
  triclinic = domain->triclinic;

  // required args

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  nexchanges = utils::inumeric(FLERR, arg[4], false, lmp);
  nmcmoves = utils::inumeric(FLERR, arg[5], false, lmp);
  nmcemc_type = utils::expand_type_int(FLERR, arg[6], Atom::ATOM, lmp);
  seed = utils::inumeric(FLERR, arg[7], false, lmp);
  reservoir_temperature = utils::numeric(FLERR, arg[8], false, lmp);
  gaugecell_volume = utils::numeric(FLERR, arg[9], false, lmp);
  ntotal = utils::inumeric(FLERR, arg[10], false, lmp);
  displace = utils::numeric(FLERR, arg[11], false, lmp);

  if (nevery <= 0) error->all(FLERR, 3, "Fix mcemc nevery value must be > 0");
  if (nexchanges < 0) error->all(FLERR, 4, "Fix mcemc nexchanges value must be >= 0");
  if (nmcmoves < 0) error->all(FLERR, 5, "Fix mcemc nmcmoves value must be >= 0");
  if (seed <= 0) error->all(FLERR, 7, "Fix mcemc random seed must be > 0");
  if (reservoir_temperature < 0.0)
    error->all(FLERR, 8, "Fix mcemc gas reservoir temperature must be >= 0");
  if (gaugecell_volume < 0.0) error->all(FLERR, 9, "Fix mcemc gauge volume must be >= 0");
  if (ntotal < 0.0) error->all(FLERR, 10, "Fix mcemc ntotal value must be >= 0");
  if (displace < 0.0)
    error->all(FLERR, 11, "Fix mcemc translation displacement distance must be >= 0");

  // read options from end of input line

  options(narg - 12, &arg[12]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp, seed);

  // random number generator, not the same for all procs

  random_unequal = new RanPark(lmp, seed);

  // error checks on region and its extent being inside simulation box

  region_xlo = region_xhi = region_ylo = region_yhi = region_zlo = region_zhi = 0.0;
  if (region) {
    if (region->bboxflag == 0)
      error->all(FLERR, Error::NOPOINTER, "Fix mcemc region {} does not support a bounding box",
                 idregion);
    if (region->dynamic_check())
      error->all(FLERR, Error::NOPOINTER, "Fix mcemc region {} cannot be dynamic", idregion);

    region_xlo = region->extent_xlo;
    region_xhi = region->extent_xhi;
    region_ylo = region->extent_ylo;
    region_yhi = region->extent_yhi;
    region_zlo = region->extent_zlo;
    region_zhi = region->extent_zhi;

    // estimate region volume using MC trials

    double coord[3];
    int inside = 0;
    int attempts = 10000000;
    for (int i = 0; i < attempts; i++) {
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      if (region->match(coord[0],coord[1],coord[2]) != 0)
        inside++;
    }

    double max_region_volume = (region_xhi - region_xlo) *
      (region_yhi - region_ylo) * (region_zhi - region_zlo);

    region_volume = max_region_volume * static_cast<double>(inside) / static_cast<double>(attempts);
  }

  // error check and further setup for exchmode = EXCHMOL

  if (exchmode == EXCHMOL) {
    if (onemols[imol]->xflag == 0)
      error->all(FLERR,"Fix mcemc molecule must have coordinates");
    if (onemols[imol]->typeflag == 0)
      error->all(FLERR,"Fix mcemc molecule must have atom types");
    if (nmcemc_type != 0)
      error->all(FLERR,"Atom type must be zero in fix mcemc mol command");
    if (onemols[imol]->qflag == 1 && atom->q == nullptr)
      error->all(FLERR,"Fix mcemc molecule has charges, but atom style does not");

    if (atom->molecular == Atom::TEMPLATE && onemols != atom->avec->onemols)
      error->all(FLERR,"Fix mcemc molecule template ID must be same "
                 "as atom_style template ID");
    onemols[imol]->check_attributes();
  }

  if (charge_flag && atom->q == nullptr)
    error->all(FLERR,"Fix mcemc atom has charge, but atom style does not");

  if (rigidflag && exchmode == EXCHATOM)
    error->all(FLERR,"Cannot use fix mcemc rigid and not molecule");
  if (shakeflag && exchmode == EXCHATOM)
    error->all(FLERR,"Cannot use fix mcemc shake and not molecule");
  if (rigidflag && shakeflag)
    error->all(FLERR,"Cannot use fix mcemc rigid and shake");
  if (rigidflag && (nmcmoves > 0))
    error->all(FLERR,"Cannot use fix mcemc rigid with MC moves");
  if (shakeflag && (nmcmoves > 0))
    error->all(FLERR,"Cannot use fix mcemc shake with MC moves");

  // setup of array of coordinates for molecule insertion
  // also used by rotation moves for any molecule

  if (exchmode == EXCHATOM) natoms_per_molecule = 1;
  else natoms_per_molecule = onemols[imol]->natoms;
  nmaxmolatoms = natoms_per_molecule;
  grow_molecule_arrays(nmaxmolatoms);

  // compute the number of MC cycles that occur nevery timesteps

  ncycles = nexchanges + nmcmoves;

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters

  mc_active = 0;

  ntranslation_attempts = 0.0;
  ntranslation_successes = 0.0;
  nrotation_attempts = 0.0;
  nrotation_successes = 0.0;
  ndeletion_attempts = 0.0;
  ndeletion_successes = 0.0;
  ninsertion_attempts = 0.0;
  ninsertion_successes = 0.0;

  mcemc_nmax = 0;
  local_gas_list = nullptr;
}



/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixMCEMC::attempt_atomic_deletion()
{
  ndeletion_attempts += 1.0;

  if (ngas == 0 || ngas <= min_ngas) return;

  int i = pick_random_gas_atom();

  int success = 0;
  if (i >= 0) {
    double deletion_energy = energy(i,nmcemc_type,-1,atom->x[i]);
    if (random_unequal->uniform() <
        gaugecell_volume*ngas*exp(beta*deletion_energy)/((ntotal-ngas+1)*volume)) {
      atom->avec->copy(atom->nlocal-1,i,1);
      atom->nlocal--;
      success = 1;
    }
  }

  int success_all = 0;
  MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);

  if (success_all) {
    atom->natoms--;
    if (atom->tag_enable) {
      if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    }
    atom->nghost = 0;
    if (triclinic) domain->x2lamda(atom->nlocal);
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    update_gas_atoms_list();
    ndeletion_successes += 1.0;
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixMCEMC::attempt_atomic_insertion()
{
  double lamda[3];

  ninsertion_attempts += 1.0;

  if (ngas >= max_ngas) return;

  // pick coordinates for insertion point

  double coord[3];
  if (region) {
    int region_attempt = 0;
    coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
    coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
    coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
    while (region->match(coord[0],coord[1],coord[2]) == 0) {
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      region_attempt++;
      if (region_attempt >= max_region_attempts) return;
    }
    if (triclinic) domain->x2lamda(coord,lamda);
  } else {
    if (triclinic == 0) {
      coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
      coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
      coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
    } else {
      lamda[0] = random_equal->uniform();
      lamda[1] = random_equal->uniform();
      lamda[2] = random_equal->uniform();

      // wasteful, but necessary

      if (lamda[0] == 1.0) lamda[0] = 0.0;
      if (lamda[1] == 1.0) lamda[1] = 0.0;
      if (lamda[2] == 1.0) lamda[2] = 0.0;

      domain->lamda2x(lamda,coord);
    }
  }

  int proc_flag = 0;
  if (triclinic == 0) {
    domain->remap(coord);
    if (!domain->inside(coord))
      error->one(FLERR,"Fix mcemc put atom outside box");
    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) proc_flag = 1;
  } else {
    if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
        lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
        lamda[2] >= sublo[2] && lamda[2] < subhi[2]) proc_flag = 1;
  }

  int success = 0;
  if (proc_flag) {
    int ii = -1;
    if (charge_flag) {
      ii = atom->nlocal + atom->nghost;
      if (ii >= atom->nmax) atom->avec->grow(0);
      atom->q[ii] = charge;
    }
    double insertion_energy = energy(ii,nmcemc_type,-1,coord);

    if (insertion_energy < MAXENERGYTEST &&
        random_unequal->uniform() <
        (ntotal-ngas)*volume*exp(-beta*insertion_energy)/(ngas+1)/(gaugecell_volume)) {
      atom->avec->create_atom(nmcemc_type,coord);
      int m = atom->nlocal - 1;

      // add to groups
      // optionally add to type-based groups

      atom->mask[m] = groupbitall;
      for (int igroup = 0; igroup < ngrouptypes; igroup++) {
        if (nmcemc_type == grouptypes[igroup])
          atom->mask[m] |= grouptypebits[igroup];
      }

      atom->v[m][0] = random_unequal->gaussian()*sigma;
      atom->v[m][1] = random_unequal->gaussian()*sigma;
      atom->v[m][2] = random_unequal->gaussian()*sigma;
      modify->create_attribute(m);

      success = 1;
    }
  }

  int success_all = 0;
  MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);

  if (success_all) {
    atom->natoms++;
    if (atom->tag_enable) {
      atom->tag_extend();
      if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    }
    atom->nghost = 0;
    if (triclinic) domain->x2lamda(atom->nlocal);
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    update_gas_atoms_list();
    ninsertion_successes += 1.0;
  }
}


/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixMCEMC::attempt_molecule_deletion()
{
  ndeletion_attempts += 1.0;

  if (ngas == 0 || ngas <= min_ngas) return;

  // work-around to avoid n=0 problem with fix rigid/nvt/small

  if (rigidflag && ngas == natoms_per_molecule) return;

  tagint deletion_molecule = pick_random_gas_molecule();
  if (deletion_molecule == -1) return;

  double deletion_energy_sum = molecule_energy(deletion_molecule);

  if (random_equal->uniform() <
      gaugecell_volume*ngas*exp(beta*deletion_energy_sum)/((ntotal-ngas+1)*volume*natoms_per_molecule)) {
    int i = 0;
    while (i < atom->nlocal) {
      if (atom->molecule[i] == deletion_molecule) {
        atom->avec->copy(atom->nlocal-1,i,1);
        atom->nlocal--;
      } else i++;
    }
    atom->natoms -= natoms_per_molecule;
    if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    atom->nghost = 0;
    if (triclinic) domain->x2lamda(atom->nlocal);
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    update_gas_atoms_list();
    ndeletion_successes += 1.0;
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixMCEMC::attempt_molecule_insertion()
{
  double lamda[3];
  ninsertion_attempts += 1.0;

  if (ngas >= max_ngas) return;

  double com_coord[3];
  if (region) {
    int region_attempt = 0;
    com_coord[0] = region_xlo + random_equal->uniform() *
      (region_xhi-region_xlo);
    com_coord[1] = region_ylo + random_equal->uniform() *
      (region_yhi-region_ylo);
    com_coord[2] = region_zlo + random_equal->uniform() *
      (region_zhi-region_zlo);
    while (region->match(com_coord[0],com_coord[1],
                                           com_coord[2]) == 0) {
      com_coord[0] = region_xlo + random_equal->uniform() *
        (region_xhi-region_xlo);
      com_coord[1] = region_ylo + random_equal->uniform() *
        (region_yhi-region_ylo);
      com_coord[2] = region_zlo + random_equal->uniform() *
        (region_zhi-region_zlo);
      region_attempt++;
      if (region_attempt >= max_region_attempts) return;
    }
    if (triclinic) domain->x2lamda(com_coord,lamda);
  } else {
    if (triclinic == 0) {
      com_coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
      com_coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
      com_coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
    } else {
      lamda[0] = random_equal->uniform();
      lamda[1] = random_equal->uniform();
      lamda[2] = random_equal->uniform();

      // wasteful, but necessary

      if (lamda[0] == 1.0) lamda[0] = 0.0;
      if (lamda[1] == 1.0) lamda[1] = 0.0;
      if (lamda[2] == 1.0) lamda[2] = 0.0;

      domain->lamda2x(lamda,com_coord);
    }
  }

  // generate point in unit cube
  // then restrict to unit sphere

  double r[3],rotmat[3][3],quat[4];
  double rsq = 1.1;
  while (rsq > 1.0) {
    r[0] = 2.0*random_equal->uniform() - 1.0;
    r[1] = 2.0*random_equal->uniform() - 1.0;
    r[2] = 2.0*random_equal->uniform() - 1.0;
    rsq = MathExtra::dot3(r, r);
  }

  double theta = random_equal->uniform() * MY_2PI;
  MathExtra::norm3(r);
  MathExtra::axisangle_to_quat(r,theta,quat);
  MathExtra::quat_to_mat(quat,rotmat);

  double insertion_energy = 0.0;
  auto *procflag = new bool[natoms_per_molecule];

  for (int i = 0; i < natoms_per_molecule; i++) {
    MathExtra::matvec(rotmat,onemols[imol]->x[i],molcoords[i]);
    molcoords[i][0] += com_coord[0];
    molcoords[i][1] += com_coord[1];
    molcoords[i][2] += com_coord[2];

    // use temporary variable for remapped position
    // so unmapped position is preserved in molcoords

    double xtmp[3];
    xtmp[0] = molcoords[i][0];
    xtmp[1] = molcoords[i][1];
    xtmp[2] = molcoords[i][2];
    domain->remap(xtmp);
    if (!domain->inside(xtmp))
      error->one(FLERR,"Fix mcemc put atom outside box");

    procflag[i] = false;
    if (triclinic == 0) {
      if (xtmp[0] >= sublo[0] && xtmp[0] < subhi[0] &&
          xtmp[1] >= sublo[1] && xtmp[1] < subhi[1] &&
          xtmp[2] >= sublo[2] && xtmp[2] < subhi[2]) procflag[i] = true;
    } else {
      domain->x2lamda(xtmp,lamda);
      if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
          lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
          lamda[2] >= sublo[2] && lamda[2] < subhi[2]) procflag[i] = true;
    }

    if (procflag[i]) {
      int ii = -1;
      if (onemols[imol]->qflag == 1) {
        ii = atom->nlocal + atom->nghost;
        if (ii >= atom->nmax) atom->avec->grow(0);
        atom->q[ii] = onemols[imol]->q[i];
      }
      insertion_energy += energy(ii,onemols[imol]->type[i],-1,xtmp);
    }
  }

  double insertion_energy_sum = 0.0;
  MPI_Allreduce(&insertion_energy,&insertion_energy_sum,1,
                MPI_DOUBLE,MPI_SUM,world);

  if (insertion_energy_sum < MAXENERGYTEST &&
      random_equal->uniform() < (ntotal-ngas)*volume*natoms_per_molecule*
      exp(-beta*insertion_energy_sum)/(ngas + natoms_per_molecule)/(gaugecell_volume)) {

    tagint maxmol = 0;
    for (int i = 0; i < atom->nlocal; i++) maxmol = MAX(maxmol,atom->molecule[i]);
    tagint maxmol_all;
    MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
    maxmol_all++;
    if (maxmol_all >= MAXTAGINT)
      error->all(FLERR,"Fix mcemc ran out of available molecule IDs");

    tagint maxtag = 0;
    for (int i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag,atom->tag[i]);
    tagint maxtag_all;
    MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

    int nlocalprev = atom->nlocal;

    double vnew[3];
    vnew[0] = random_equal->gaussian()*sigma;
    vnew[1] = random_equal->gaussian()*sigma;
    vnew[2] = random_equal->gaussian()*sigma;

    for (int i = 0; i < natoms_per_molecule; i++) {
      if (procflag[i]) {
        atom->avec->create_atom(onemols[imol]->type[i],molcoords[i]);
        int m = atom->nlocal - 1;

        // add to groups
        // optionally add to type-based groups

        atom->mask[m] = groupbitall;
        for (int igroup = 0; igroup < ngrouptypes; igroup++) {
          if (nmcemc_type == grouptypes[igroup])
            atom->mask[m] |= grouptypebits[igroup];
        }

        atom->image[m] = imagezero;
        domain->remap(atom->x[m],atom->image[m]);
        atom->molecule[m] = maxmol_all;
        if (maxtag_all+i+1 >= MAXTAGINT)
          error->all(FLERR,"Fix mcemc ran out of available atom IDs");
        atom->tag[m] = maxtag_all + i + 1;
        atom->v[m][0] = vnew[0];
        atom->v[m][1] = vnew[1];
        atom->v[m][2] = vnew[2];

        atom->add_molecule_atom(onemols[imol],i,m,maxtag_all);
        modify->create_attribute(m);
      }
    }

    // FixRigidSmall::set_molecule stores rigid body attributes
    // FixShake::set_molecule stores shake info for molecule

    for (int submol = 0; submol < nmol; ++submol) {
      if (rigidflag)
        fixrigid->set_molecule(nlocalprev,maxtag_all,submol,com_coord,vnew,quat);
      else if (shakeflag)
        fixshake->set_molecule(nlocalprev,maxtag_all,submol,com_coord,vnew,quat);
    }
    atom->natoms += natoms_per_molecule;
    if (atom->natoms < 0)
      error->all(FLERR,"Too many total atoms");
    atom->nbonds += onemols[imol]->nbonds;
    atom->nangles += onemols[imol]->nangles;
    atom->ndihedrals += onemols[imol]->ndihedrals;
    atom->nimpropers += onemols[imol]->nimpropers;
    if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    atom->nghost = 0;
    if (triclinic) domain->x2lamda(atom->nlocal);
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    update_gas_atoms_list();
    ninsertion_successes += 1.0;
  }
  delete[] procflag;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixMCEMC::attempt_atomic_deletion_full()
{
  double q_tmp;
  const int q_flag = atom->q_flag;

  ndeletion_attempts += 1.0;

  if (ngas == 0 || ngas <= min_ngas) return;

  double energy_before = energy_stored;

  const int i = pick_random_gas_atom();

  int tmpmask;
  if (i >= 0) {
    tmpmask = atom->mask[i];
    atom->mask[i] = exclusion_group_bit;
    if (q_flag) {
      q_tmp = atom->q[i];
      atom->q[i] = 0.0;
    }
  }
  if (force->kspace) force->kspace->qsum_qsq();
  if (force->pair->tail_flag) force->pair->reinit();
  double energy_after = energy_full();

  if (random_equal->uniform() <
      gaugecell_volume*ngas*exp(beta*(energy_before - energy_after))/((ntotal-ngas+1)*volume)) {
    if (i >= 0) {
      atom->avec->copy(atom->nlocal-1,i,1);
      atom->nlocal--;
    }
    atom->natoms--;
    if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    ndeletion_successes += 1.0;
    energy_stored = energy_after;
  } else {
    if (i >= 0) {
      atom->mask[i] = tmpmask;
      if (q_flag) atom->q[i] = q_tmp;
    }
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
    energy_stored = energy_before;
  }
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixMCEMC::attempt_atomic_insertion_full()
{
  double lamda[3];
  ninsertion_attempts += 1.0;

  if (ngas >= max_ngas) return;

  double energy_before = energy_stored;

  double coord[3];
  if (region) {
    int region_attempt = 0;
    coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
    coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
    coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
    while (region->match(coord[0],coord[1],coord[2]) == 0) {
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      region_attempt++;
      if (region_attempt >= max_region_attempts) return;
    }
    if (triclinic) domain->x2lamda(coord,lamda);
  } else {
    if (triclinic == 0) {
      coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
      coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
      coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
    } else {
      lamda[0] = random_equal->uniform();
      lamda[1] = random_equal->uniform();
      lamda[2] = random_equal->uniform();

      // wasteful, but necessary

      if (lamda[0] == 1.0) lamda[0] = 0.0;
      if (lamda[1] == 1.0) lamda[1] = 0.0;
      if (lamda[2] == 1.0) lamda[2] = 0.0;

      domain->lamda2x(lamda,coord);
    }
  }

  int proc_flag = 0;
  if (triclinic == 0) {
    domain->remap(coord);
    if (!domain->inside(coord))
      error->one(FLERR,"Fix mcemc put atom outside box");
    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) proc_flag = 1;
  } else {
    if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
        lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
        lamda[2] >= sublo[2] && lamda[2] < subhi[2]) proc_flag = 1;
  }

  if (proc_flag) {
    atom->avec->create_atom(nmcemc_type,coord);
    int m = atom->nlocal - 1;

    // add to groups
    // optionally add to type-based groups

    atom->mask[m] = groupbitall;
    for (int igroup = 0; igroup < ngrouptypes; igroup++) {
      if (nmcemc_type == grouptypes[igroup])
        atom->mask[m] |= grouptypebits[igroup];
    }

    atom->v[m][0] = random_unequal->gaussian()*sigma;
    atom->v[m][1] = random_unequal->gaussian()*sigma;
    atom->v[m][2] = random_unequal->gaussian()*sigma;
    if (charge_flag) atom->q[m] = charge;
    modify->create_attribute(m);
  }

  atom->natoms++;
  if (atom->tag_enable) {
    atom->tag_extend();
    if (atom->map_style != Atom::MAP_NONE) atom->map_init();
  }
  atom->nghost = 0;
  if (triclinic) domain->x2lamda(atom->nlocal);
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  if (force->kspace) force->kspace->qsum_qsq();
  if (force->pair->tail_flag) force->pair->reinit();
  double energy_after = energy_full();

  if (energy_after < MAXENERGYTEST &&
      random_equal->uniform() <
      (ntotal-ngas)*volume*exp(beta*(energy_before - energy_after))/(ngas+1)/(gaugecell_volume)) {

    ninsertion_successes += 1.0;
    energy_stored = energy_after;
  } else {
    atom->natoms--;
    if (proc_flag) atom->nlocal--;
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
    energy_stored = energy_before;
  }
  update_gas_atoms_list();
}


/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixMCEMC::attempt_molecule_deletion_full()
{
  ndeletion_attempts += 1.0;

  if (ngas == 0 || ngas <= min_ngas) return;

  // work-around to avoid n=0 problem with fix rigid/nvt/small

  if (rigidflag && ngas == natoms_per_molecule) return;

  tagint deletion_molecule = pick_random_gas_molecule();
  if (deletion_molecule == -1) return;

  double energy_before = energy_stored;

  // check nmolq, grow arrays if necessary

  int nmolq = 0;
  for (int i = 0; i < atom->nlocal; i++)
    if (atom->molecule[i] == deletion_molecule)
      if (atom->q_flag) nmolq++;

  if (nmolq > nmaxmolatoms)
    grow_molecule_arrays(nmolq);

  int m = 0;
  int *tmpmask = new int[atom->nlocal];
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->molecule[i] == deletion_molecule) {
      tmpmask[i] = atom->mask[i];
      atom->mask[i] = exclusion_group_bit;
      toggle_intramolecular(i);
      if (atom->q_flag) {
        molq[m] = atom->q[i];
        m++;
        atom->q[i] = 0.0;
      }
    }
  }
  if (force->kspace) force->kspace->qsum_qsq();
  if (force->pair->tail_flag) force->pair->reinit();
  double energy_after = energy_full();

  // energy_before corrected by energy_intra

  double deltaphi = gaugecell_volume*ngas*exp(beta*((energy_before - energy_intra) - energy_after))/((ntotal-ngas+1)*volume*natoms_per_molecule);

  if (random_equal->uniform() < deltaphi) {
    int i = 0;
    while (i < atom->nlocal) {
      if (atom->molecule[i] == deletion_molecule) {
        atom->avec->copy(atom->nlocal-1,i,1);
        atom->nlocal--;
      } else i++;
    }
    atom->natoms -= natoms_per_molecule;
    if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    ndeletion_successes += 1.0;
    energy_stored = energy_after;
  } else {
    energy_stored = energy_before;
    int m = 0;
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->molecule[i] == deletion_molecule) {
        atom->mask[i] = tmpmask[i];
        toggle_intramolecular(i);
        if (atom->q_flag) {
          atom->q[i] = molq[m];
          m++;
        }
      }
    }
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
  }
  update_gas_atoms_list();
  delete[] tmpmask;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixMCEMC::attempt_molecule_insertion_full()
{
  double lamda[3];
  ninsertion_attempts += 1.0;

  if (ngas >= max_ngas) return;

  double energy_before = energy_stored;

  tagint maxmol = 0;
  for (int i = 0; i < atom->nlocal; i++) maxmol = MAX(maxmol,atom->molecule[i]);
  tagint maxmol_all;
  MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
  maxmol_all++;
  if (maxmol_all >= MAXTAGINT)
    error->all(FLERR,"Fix mcemc ran out of available molecule IDs");
  int insertion_molecule = maxmol_all;

  tagint maxtag = 0;
  for (int i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag,atom->tag[i]);
  tagint maxtag_all;
  MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

  int nlocalprev = atom->nlocal;

  double com_coord[3];
  if (region) {
    int region_attempt = 0;
    com_coord[0] = region_xlo + random_equal->uniform() *
      (region_xhi-region_xlo);
    com_coord[1] = region_ylo + random_equal->uniform() *
      (region_yhi-region_ylo);
    com_coord[2] = region_zlo + random_equal->uniform() *
      (region_zhi-region_zlo);
    while (region->match(com_coord[0],com_coord[1],
                                           com_coord[2]) == 0) {
      com_coord[0] = region_xlo + random_equal->uniform() *
        (region_xhi-region_xlo);
      com_coord[1] = region_ylo + random_equal->uniform() *
        (region_yhi-region_ylo);
      com_coord[2] = region_zlo + random_equal->uniform() *
        (region_zhi-region_zlo);
      region_attempt++;
      if (region_attempt >= max_region_attempts) return;
    }
    if (triclinic) domain->x2lamda(com_coord,lamda);
  } else {
    if (triclinic == 0) {
      com_coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
      com_coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
      com_coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
    } else {
      lamda[0] = random_equal->uniform();
      lamda[1] = random_equal->uniform();
      lamda[2] = random_equal->uniform();

      // wasteful, but necessary

      if (lamda[0] == 1.0) lamda[0] = 0.0;
      if (lamda[1] == 1.0) lamda[1] = 0.0;
      if (lamda[2] == 1.0) lamda[2] = 0.0;

      domain->lamda2x(lamda,com_coord);
    }

  }

  // generate point in unit cube
  // then restrict to unit sphere

  double r[3],rotmat[3][3],quat[4];
  double rsq = 1.1;
  while (rsq > 1.0) {
    r[0] = 2.0*random_equal->uniform() - 1.0;
    r[1] = 2.0*random_equal->uniform() - 1.0;
    r[2] = 2.0*random_equal->uniform() - 1.0;
    rsq = MathExtra::dot3(r, r);
  }

  double theta = random_equal->uniform() * MY_2PI;
  MathExtra::norm3(r);
  MathExtra::axisangle_to_quat(r,theta,quat);
  MathExtra::quat_to_mat(quat,rotmat);

  double vnew[3];
  vnew[0] = random_equal->gaussian()*sigma;
  vnew[1] = random_equal->gaussian()*sigma;
  vnew[2] = random_equal->gaussian()*sigma;

  for (int i = 0; i < natoms_per_molecule; i++) {
    double xtmp[3];
    MathExtra::matvec(rotmat,onemols[imol]->x[i],xtmp);
    xtmp[0] += com_coord[0];
    xtmp[1] += com_coord[1];
    xtmp[2] += com_coord[2];

    // need to adjust image flags in remap()

    imageint imagetmp = imagezero;
    domain->remap(xtmp,imagetmp);
    if (!domain->inside(xtmp))
      error->one(FLERR,"Fix mcemc put atom outside box");

    int proc_flag = 0;
    if (triclinic == 0) {
      if (xtmp[0] >= sublo[0] && xtmp[0] < subhi[0] &&
          xtmp[1] >= sublo[1] && xtmp[1] < subhi[1] &&
          xtmp[2] >= sublo[2] && xtmp[2] < subhi[2]) proc_flag = 1;
    } else {
      domain->x2lamda(xtmp,lamda);
      if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
          lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
          lamda[2] >= sublo[2] && lamda[2] < subhi[2]) proc_flag = 1;
    }

    if (proc_flag) {
      atom->avec->create_atom(onemols[imol]->type[i],xtmp);
      int m = atom->nlocal - 1;

      // add to groups
      // optionally add to type-based groups

      atom->mask[m] = groupbitall;
      for (int igroup = 0; igroup < ngrouptypes; igroup++) {
        if (nmcemc_type == grouptypes[igroup])
          atom->mask[m] |= grouptypebits[igroup];
      }

      atom->image[m] = imagetmp;
      atom->molecule[m] = insertion_molecule;
      if (maxtag_all+i+1 >= MAXTAGINT)
        error->all(FLERR,"Fix mcemc ran out of available atom IDs");
      atom->tag[m] = maxtag_all + i + 1;
      atom->v[m][0] = vnew[0];
      atom->v[m][1] = vnew[1];
      atom->v[m][2] = vnew[2];

      atom->add_molecule_atom(onemols[imol],i,m,maxtag_all);
      modify->create_attribute(m);
    }
  }

  // FixRigidSmall::set_molecule stores rigid body attributes
  // FixShake::set_molecule stores shake info for molecule

  for (int submol = 0; submol < nmol; ++submol) {
    if (rigidflag)
      fixrigid->set_molecule(nlocalprev,maxtag_all,submol,com_coord,vnew,quat);
    else if (shakeflag)
      fixshake->set_molecule(nlocalprev,maxtag_all,submol,com_coord,vnew,quat);
  }
  atom->natoms += natoms_per_molecule;
  if (atom->natoms < 0)
    error->all(FLERR,"Too many total atoms");
  atom->nbonds += onemols[imol]->nbonds;
  atom->nangles += onemols[imol]->nangles;
  atom->ndihedrals += onemols[imol]->ndihedrals;
  atom->nimpropers += onemols[imol]->nimpropers;
  if (atom->map_style != Atom::MAP_NONE) atom->map_init();
  atom->nghost = 0;
  if (triclinic) domain->x2lamda(atom->nlocal);
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  if (force->kspace) force->kspace->qsum_qsq();
  if (force->pair->tail_flag) force->pair->reinit();
  double energy_after = energy_full();

  // energy_after corrected by energy_intra

  double deltaphi = (ntotal-ngas)*volume*natoms_per_molecule*
    exp(beta*(energy_before - (energy_after - energy_intra)))/(ngas + natoms_per_molecule)/(gaugecell_volume);

  if (energy_after < MAXENERGYTEST &&
      random_equal->uniform() < deltaphi) {

    ninsertion_successes += 1.0;
    energy_stored = energy_after;

  } else {

    atom->nbonds -= onemols[imol]->nbonds;
    atom->nangles -= onemols[imol]->nangles;
    atom->ndihedrals -= onemols[imol]->ndihedrals;
    atom->nimpropers -= onemols[imol]->nimpropers;
    atom->natoms -= natoms_per_molecule;

    energy_stored = energy_before;
    int i = 0;
    while (i < atom->nlocal) {
      if (atom->molecule[i] == insertion_molecule) {
        atom->avec->copy(atom->nlocal-1,i,1);
        atom->nlocal--;
      } else i++;
    }
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
  }
  update_gas_atoms_list();
}

