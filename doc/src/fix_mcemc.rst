.. index:: fix mcemc

fix mcemc command
=================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID mcemc N X M type seed T Vgauge Ntotal displace keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* mcemc = style name of this fix command
* N = invoke this fix every N steps
* X = average number of MCEMC exchanges to attempt every N steps
* M = average number of MC moves to attempt every N steps
* type = atom type (1-Ntypes or type label) for inserted atoms (must be 0 if mol keyword used)
* seed = random # seed (positive integer)
* T = temperature of the ideal gas reservoir (temperature units)
* Vgauge = gauge volume for MCEMC exchanges (volume units)
* Ntotal = total number of atoms in the system (positive integer)
* displace = maximum Monte Carlo translation distance (length units)
* zero or more keyword/value pairs may be appended to args

  .. parsed-literal::

     keyword = *mol*, *region*, *maxangle*, *full_energy*, *charge*, *group*, *grouptype*, *intra_energy*, *tfac_insert*, or *overlap_cutoff*
       *mol* value = template-ID
         template-ID = ID of molecule template specified in a separate :doc:`molecule <molecule>` command
       *mcmoves* values = Patomtrans Pmoltrans Pmolrotate
         Patomtrans = proportion of atom translation MC moves
         Pmoltrans = proportion of molecule translation MC moves
         Pmolrotate = proportion of molecule rotation MC moves
       *rigid* value = fix-ID
         fix-ID = ID of :doc:`fix rigid/small <fix_rigid>` command
       *shake* value = fix-ID
         fix-ID = ID of :doc:`fix shake <fix_shake>` command
       *region* value = region-ID
         region-ID = ID of region where MCEMC exchanges and MC moves are allowed
       *maxangle* value = maximum molecular rotation angle (degrees)
       *full_energy* = compute the entire system energy when performing MCEMC exchanges and MC moves
       *charge* value = charge of inserted atoms (charge units)
       *group* value = group-ID
         group-ID = group-ID for inserted atoms (string)
       *grouptype* values = type group-ID
         type = atom type (1-Ntypes or type label)
         group-ID = group-ID for inserted atoms (string)
       *intra_energy* value = intramolecular energy (energy units)
       *tfac_insert* value = scale up/down temperature of inserted atoms (unitless)
       *overlap_cutoff* value = maximum pair distance for overlap rejection (distance units)
       *max* value = Maximum number of atoms allowed in the fix group (and region)
       *min* value = Minimum number of atoms allowed in the fix group (and region)

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 gas mcemc 10 1000 1000 2 29494 298.0 1200 100 0.01
   fix 3 water mcemc 10 100 100 0 3456543 3.0 1200 100 mol my_one_water maxangle 180 full_energy
   fix 4 my_gas mcemc 1 10 10 1 123456543 300.0 1200 100 1.0 region disk

Description
"""""""""""

.. versionadded:: TBD

This fix performs mesocanonical Monte Carlo (MCEMC) also known as gauge
cell simulation by exchanging particles with a finite volume ideal gas
reservoir (gauge cell) at the same temperature as the system as
discussed in :ref:`(Parashar) <Parashar>`.  It also attempts Monte Carlo
moves (translations and rotations) of particles within the simulation
cell.  Specific uses of this fix are to compute the adsorption isotherm
in porous materials or to compute the vapor-liquid equilibrium of
fluids.  This fix is complementary to the :doc:`fix gcmc <fix_gcmc>`
command, which performs grand canonical Monte Carlo (GCMC) by exchanging
particles with an infinite chemical potential reservoir.  MCEMC and GCMC
give identical adsorption isotherms for microporous materials.  But for
large pores (> 2 nm), GCMC gives a hysteretic adsorption/desorption
isotherm, while MCEMC gives a reversible S- shaped van der Waals type
isotherm, as discussed in :ref:`(Parashar) <Parashar>`.  The MCEMC
isotherm spans the stable and meta-stable states, while GCMC samples
only the stable states.  The MCEMC method is a middle ground between the
grand canonical ensemble which permits unlimited fluctuations, and the
canonical ensemble, which considers a closed system.  MCEMC simulations
generate adsorption isotherms equivalent to the canonical ensemble
isotherms with an accuracy of one molecule.  MCEMC is equivalent to GCMC
when the gauge cell volume is infinite and is equivalent to the
canonical ensemble when gauge cell volume is zero.

Every *N* timesteps the fix attempts both MCEMC exchanges (insertions or
deletions) and MC moves of gas atoms or molecules.  On those timesteps,
the average number of attempted MCEMC exchanges is *X*, while the average
number of attempted MC moves is *M*.  For MCEMC exchanges of either
molecular or atomic gases, these exchanges can be either deletions or
insertions, with equal probability.

The possible choices for MC moves are translation of an atom,
translation of a molecule, and rotation of a molecule.  The relative
amounts of each are determined by the optional *mcmoves* keyword (see
below).  The default behavior is as follows.  If the *mol* keyword is
used, only molecule translations and molecule rotations are performed
with equal probability.  Conversely, if the *mol* keyword is not used,
only atom translations are performed.  *M* should typically be chosen to
be approximately equal to the expected number of gas atoms or molecules
of the given type within the simulation cell or region, which will
result in roughly one MC move per atom or molecule per MC cycle.

All inserted particles are always added to two groups: the default group
"all" and the fix group specified in the fix command.  In addition,
particles are also added to any groups specified by the *group* and
*grouptype* keywords.  If inserted particles are individual atoms, they
are assigned the atom type given by the type argument.  If they are
molecules, the type argument has no effect and must be set to zero.
Instead, the type of each atom in the inserted molecule is specified in
the file read by the :doc:`molecule <molecule>` command.

.. note::

   Care should be taken to apply *fix mcemc* to a group that contains
   *only* those atoms and molecules that you wish to manipulate using
   Monte Carlo.  Hence it is generally not a good idea to specify the
   default group "all" in the fix command, although it is allowed.

This command may optionally use the *region* keyword to define an
exchange and move volume.  The specified region must have been
previously defined with a :doc:`region <region>` command.  It must be
defined with side = *in*\ .  Insertion attempts occur only within the
specified region.  For non-rectangular regions, random trial points are
generated within the rectangular bounding box until a point is found
that lies inside the region.  If no valid point is generated after 1000
trials, no insertion is performed, but it is counted as an attempted
insertion.  Move and deletion attempt candidates are selected from gas
atoms or molecules within the region.  If there are no candidates, no
move or deletion is performed, but it is counted as an attempt move or
deletion.  If an attempted move places the atom or molecule
center-of-mass outside the specified region, a new attempted move is
generated.  This process is repeated until the atom or molecule
center-of-mass is inside the specified region.

Note that neighbor lists are re-built every timestep that this fix is
invoked, so you should not set *N* to be too small.  However, periodic
rebuilds are necessary in order to avoid dangerous rebuilds and missed
interactions.  Specifically, avoid performing so many MC translations
per timestep that atoms can move beyond the neighbor list skin distance.
See the :doc:`neighbor <neighbor>` command for details.

When an atom or molecule is to be inserted, its coordinates are chosen
at a random position within the current simulation cell or region, and
new atom velocities are randomly chosen from the specified temperature
distribution given by *T*.  The effective temperature for new atom
velocities can be increased or decreased using the optional keyword
*tfac_insert* (see below).  Relative coordinates for atoms in a molecule
are taken from the template molecule provided by the user.  The center
of mass of the molecule is placed at the insertion point.  The
orientation of the molecule is chosen at random by rotating about this
point.

Individual atoms are inserted, unless the *mol* keyword is used.  It
specifies a *template-ID* previously defined using the :doc:`molecule
<molecule>` command, which reads a file that defines the molecule.  The
coordinates, atom types, charges, etc., as well as any bonding and
special neighbor information for the molecule can be specified in the
molecule file.  See the :doc:`molecule <molecule>` command for details.
The only settings required to be in this file are the coordinates and
types of atoms in the molecule.

When not using the *mol* keyword, you should ensure you do not delete
atoms that are bonded to other atoms, or LAMMPS will soon generate an
error when it tries to find bonded neighbors.  LAMMPS will warn you if
any of the atoms eligible for deletion have a non-zero molecule ID, but
does not check for this at the time of deletion.

If you wish to insert molecules using the *mol* keyword that will be
treated as rigid bodies, use the *rigid* keyword, specifying as its
value the ID of a separate :doc:`fix rigid/small <fix_rigid>` command
which also appears in your input script.

.. note::

   If you wish the new rigid molecules (and other rigid molecules) to be
   thermostatted correctly via :doc:`fix rigid/small/nvt <fix_rigid>` or
   :doc:`fix rigid/small/npt <fix_rigid>`, then you need to use the
   :doc:`fix_modify dynamic/dof yes <fix_modify>` command for the rigid
   fix.  This is to inform that fix that the molecule count will vary
   dynamically.

If you wish to insert molecules via the *mol* keyword, that will have
their bonds or angles constrained via SHAKE, use the *shake* keyword,
specifying as its value the ID of a separate :doc:`fix shake
<fix_shake>` command which also appears in your input script.

Optionally, users may specify the relative amounts of different MC moves
using the *mcmoves* keyword. The values *Patomtrans*, *Pmoltrans*,
*Pmolrotate* specify the average proportion of atom translations,
molecule translations, and molecule rotations, respectively. The values
must be non-negative integers or real numbers, with at least one
non-zero value. For example, (10,30,0) would result in 25% of the MC
moves being atomic translations, 75% molecular translations, and no
molecular rotations.

Optionally, users may specify the maximum rotation angle for molecular
rotations using the *maxangle* keyword and specifying the angle in
degrees. Rotations are performed by generating a random point on the
unit sphere and a random rotation angle on the range [0,maxangle).  The
molecule is then rotated by that angle about an axis passing through the
molecule center of mass. The axis is parallel to the unit vector defined
by the point on the unit sphere.  The same procedure is used for
randomly rotating molecules when they are inserted, except that the
maximum angle is 360 degrees.

Note that *fix mcemc* does not use configurational bias MC or any other
kind of sampling of intramolecular degrees of freedom.  Inserted
molecules can have different orientations, but they will all have the
same intramolecular configuration, which was specified in the molecule
command input.

For atomic gasses, inserted atoms have the specified atom type, but
deleted atoms are any atoms that have been inserted or that already
belong to the fix group.  For molecular gasses, exchanged molecules use
the same atom types as in the template molecule supplied by the user.
In both cases, exchanged atoms/molecules are assigned to two groups: the
default group "all" and the fix group (which can also be "all").

During the MCEMC exchange move, the particles are exchanged between the
system and a finite ideal gas reservoir (gauge cell) such that the total
number of particles in the combined system (system + gauge cell) remains
constant.

.. math::

   N_{total} = N_{system} + N_{gauge}

where :math:`N_{total}` is the total number of particles in the combined
system, :math:`N_{system}` is the number of particles in the simulation
system, and :math:`N_{gauge}` is the number of particles in the ideal
gas reservoir (gauge cell). The gauge cell has a fixed volume
(:math:`V_{gauge}`) and is maintained at the same temperature (*T*) as
the simulation system. The combined system is in thermal and chemical
equilibrium, hence the chemical potential of the system is equal to that
of the gauge cell.  The chemical potential of the gauge cell (and hence
the system) is given by:

.. math::
    \mu^{id} = k_{B}T ln(\frac{N_{gauge}\Lambda^{3}}{V_{gauge}})

where :math:`k_{B}` is the Boltzmann constant, :math:`\Lambda` is the
thermal de Broglie wavelength of the ideal gas particles at temperature
*T*, :math:`V_{gauge}` is the volume of gauge cell, and
:math:`N_{gauge}` is the average number of particles in the gauge cell.
The constant :math:`\Lambda` is required for dimensional consistency.
For all unit styles except *lj* it is defined as the thermal de Broglie
wavelength.

.. math::

   \Lambda = \sqrt{ \frac{h^2}{2 \pi m k_B T}}

where *h* is Planck's constant, and *m* is the mass of the exchanged
atom or molecule.  For unit style *lj*, :math:`\Lambda` is simply set to
unity.

During an MCEMC insertion move, a particle is randomly selected from the
gauge cell and inserted into the simulation system. During an MCEMC
deletion move, a particle is randomly selected from the simulation
system and moved to the gauge cell. The acceptance probability for
particle addition to the system is given by

.. math::

   acc(N \rightarrow N+1) = \min\left(1, \frac{V N_{gauge}}{(N+1) V_{gauge}} \exp(-\beta[E(N+1)-E(N)])\right)

The acceptance probability for particle deletion from system is given by:

.. math::

   acc(N \rightarrow N-1) = \min\left(1, \frac{V_{gauge} N}{(N_{gauge}+1) V} \exp(-\beta[E(N)-E(N-1)])\right)

In GCMC, the chemical potential of the infinite reservoir is imposed on
the system using the exchange move, while in MCEMC, the gauge cell
(which is assumed to be an ideal gas) is used to measure the chemical
potential of the system. In GCMC, the only input is the chemical
potential of the infinite reservoir and one can obtain the number of
particles in the system as the average number of particles observed in
the system. In MCEMC, there are two inputs: :math:`N_{total}` and
:math:`V_{gauge}`.  From these two inputs, the chemical potential (and
hence fugacity) of the system is calculated based on ideal gas chemical
potential in the gauge cell. The number of particles adsorbed
is simply the average number of particles observed in the system. Note
that the choice of :math:`N_{total}` and :math:`V_{gauge}` is not
arbitrary. The gauge cell should be sufficiently smaller that it can
stabilize the fluid configuration within the system but should be sufficiently
large for accurate measurement of chemical potential. The recommendation is to 
choose :math:`V_{gauge}` such that the gauge cell contain roughly
70-80 particles. Generating a GCMC isotherm beforehand can help you
choose an appropriate value of :math:`N_{total}`.

The *full_energy* option means that the fix calculates the total
potential energy of the entire simulated system, instead of just the
energy of the part that is changed.  The total system energy before and
after the proposed MCEMC exchange or MC move is then used in the
Metropolis criterion to determine whether or not to accept the proposed
change.  By default, this option is off, in which case only partial
energies are computed to determine the energy difference due to the
proposed change.

The *full_energy* option is needed for systems with complicated
potential energy calculations, including the following:

* long-range electrostatics (kspace)
* many-body pair styles
* hybrid pair styles
* eam pair styles
* tail corrections
* need to include potential energy contributions from other fixes

In these cases, LAMMPS will automatically apply the *full_energy*
keyword and issue a warning message.

When the *mol* keyword is used, the *full_energy* option also includes
the intramolecular energy of inserted and deleted molecules, whereas
this energy is not included when *full_energy* is not used.  If this is
not desired, the *intra_energy* keyword can be used to define an amount
of energy that is subtracted from the final energy when a molecule is
inserted, and subtracted from the initial energy when a molecule is
deleted.  For molecules that have a non-zero intramolecular energy, this
will ensure roughly the same behavior whether or not the *full_energy*
option is used.

Inserted atoms and molecules are assigned random velocities based on the
specified temperature *T*.  Because the relative velocity of all atoms
in the molecule is zero, this may result in inserted molecules that are
systematically too cold.  In addition, the intramolecular potential
energy of the inserted molecule may cause the kinetic energy of the
molecule to quickly increase or decrease after insertion.  The
*tfac_insert* keyword allows the user to counteract these effects by
changing the temperature used to assign velocities to inserted atoms and
molecules by a constant factor.  For a particular application, some
experimentation may be required to find a value of *tfac_insert* that
results in inserted molecules that equilibrate quickly to the correct
temperature.

Some fixes have an associated potential energy. Examples of such fixes
include: :doc:`efield <fix_efield>`, :doc:`gravity <fix_gravity>`,
:doc:`addforce <fix_addforce>`, :doc:`langevin <fix_langevin>`,
:doc:`restrain <fix_restrain>`, :doc:`temp/berendsen
<fix_temp_berendsen>`, :doc:`temp/rescale <fix_temp_rescale>`, and
:doc:`wall fixes <fix_wall>`.  For that energy to be included in the
total potential energy of the system (the quantity used when performing
MCEMC exchange and MC moves), you MUST enable the :doc:`fix_modify
<fix_modify>` *energy* option for that fix.  The doc pages for
individual :doc:`fix <fix>` commands specify if this should be done.

Use the *charge* option to insert atoms with a user-specified point
charge. Note that doing so will cause the system to become non-neutral.
LAMMPS issues a warning when using long-range electrostatics (kspace)
with non-neutral systems. See the :doc:`compute group/group
<compute_group_group>` documentation for more details about simulating
non-neutral systems with kspace on.

Use of this fix typically will cause the number of atoms to fluctuate,
therefore, you will want to use the :doc:`compute_modify dynamic/dof
<compute_modify>` command to ensure that the current number of atoms is
used as a normalizing factor each time temperature is computed. A simple
example of this is:

.. code-block:: LAMMPS

   compute_modify thermo_temp dynamic/dof yes

A more complicated example is listed earlier on this page
in the context of NVT dynamics.

.. note::

   If the density of the cell is initially very small or zero, and
   increases to a much larger density after a period of equilibration,
   then certain quantities that are only calculated once at the start
   (kspace parameters) may no longer be accurate.  The solution is to
   start a new simulation after the equilibrium density has been
   reached.

With some pair_styles, such as :doc:`Buckingham <pair_buck>`,
:doc:`Born-Mayer-Huggins <pair_born>` and :doc:`ReaxFF <pair_reaxff>`,
two atoms placed close to each other may have an arbitrary large,
negative potential energy due to the functional form of the potential.
While these unphysical configurations are inaccessible to typical
dynamical trajectories, they can be generated by Monte Carlo moves. The
*overlap_cutoff* keyword suppresses these moves by effectively assigning
an infinite positive energy to all new configurations that place any
pair of atoms closer than the specified overlap cutoff distance.

The *max* and *min* keywords allow for the restriction of the number of
atoms in the fix group (and region in case the *region* keyword is
used).  They automatically reject all insertion or deletion moves that
would take the system beyond the set boundaries.  Should the system
already be beyond the boundary, only moves that bring the system closer
to the bounds may be accepted.

The *group* keyword adds all inserted atoms to the :doc:`group <group>`
of the group-ID value. The *grouptype* keyword adds all inserted atoms
of the specified type to the :doc:`group <group>` of the group-ID value.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the fix to :doc:`binary restart files
<restart>`.  This includes information about the random number
generator seed, the next timestep for MC exchanges, the number of MC
step attempts and successes etc.  See the :doc:`read_restart
<read_restart>` command for info on how to re-specify a fix in an
input script that reads a restart file, so that the operation of the
fix continues in an uninterrupted fashion.

.. note::

   For this to work correctly, the timestep must **not** be changed
   after reading the restart with :doc:`reset_timestep <reset_timestep>`.
   The fix will try to detect it and stop with an error.

None of the :doc:`fix_modify <fix_modify>` options are relevant to
this fix.

This fix computes a global vector of length 8, which can be accessed
by various :doc:`output commands <Howto_output>`.  The vector values are
the following global cumulative quantities:

  #. translation attempts
  #. translation successes
  #. insertion attempts
  #. insertion successes
  #. deletion attempts
  #. deletion successes
  #. rotation attempts
  #. rotation successes

The vector values calculated by this fix are "intensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the MC package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

This fix style requires an :doc:`atom style <atom_style>` with per atom
type masses.

Do not set "neigh_modify once yes" or else this fix will never be
called.  Reneighboring is **required**.

This fix style is only usable for 3D simulations.

This fix can be run in parallel, but aspects of the MCEMC part will not
scale well in parallel.  Currently, molecule translations and rotations
are not supported with more than one MPI process.  It is still possible
to do parallel molecule exchange without translation and rotation moves
by setting MC moves to zero and/or by using the *mcmoves* keyword with
*Pmoltrans* = *Pmolrotate* = 0 .

When using *fix mcemc* in combination with :doc:`fix shake <fix_shake>`
or :doc:`fix rigid <fix_rigid>`, only MCEMC exchange moves are
supported, so the argument *M* must be zero.

When using *fix mcemc* in combination with :doc:`fix rigid <fix_rigid>`,
deletion of the last remaining molecule is not allowed for technical
reasons, and so the molecule count will never drop below 1, regardless
of the specified chemical potential.

Note that very lengthy simulations involving insertions/deletions of
billions of gas molecules may run out of atom or molecule IDs and
trigger an error, so it is better to run multiple shorter-duration
simulations.  The :doc:`reset_atoms <reset_atoms>` command can be used
to "compress" the atom and molecule IDs between runs.  Likewise, very
large molecules have not been tested and may turn out to be problematic.

Use of multiple *fix mcemc* commands in the same input script can be
problematic if using a template molecule.  The issue is that the
user-referenced template molecule in the second *fix mcemc* command may
no longer exist since it might have been deleted by the first *fix
mcemc* command.  An existing template molecule will need to be referenced
by the user for each subsequent *fix mcemc* command.

Related commands
""""""""""""""""

:doc:`fix gcmc <fix_gcmc>`,
:doc:`fix atom/swap <fix_atom_swap>`,
:doc:`fix nvt <fix_nh>`, :doc:`neighbor <neighbor>`,
:doc:`fix deposit <fix_deposit>`, :doc:`fix evaporate <fix_evaporate>`,
:doc:`delete_atoms <delete_atoms>`

Defaults
""""""""

The option defaults are mol = no, maxangle = 10, overlap_cutoff = 0.0,
intra_energy = 0.0, tfac_insert = 1.0.
(Patomtrans, Pmoltrans, Pmolrotate) = (1, 0, 0) for mol = no and
(0, 1, 1) for mol = yes. full_energy = no,
except for the situations where full_energy is required, as
listed above.

----------

.. _Parashar:

**(Parashar)** Parashar, S., Neimark, A.V., Journal of Colloid And Interface Science, 2024. DOI: 10.1016/j.jcis.2024.06.083

.. _Neimark:

**(Neimark)** Neimark, A.V., Vishnyakov, A., Physical Review E, 2000. DOI: 10.1103/PhysRevE.62.4611
