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
* keyword/value pairs = See :doc:`fix gcmc <fix_gcmc>` for descriptions of keywords such as
   *mol*, *region*, *maxangle*, *full_energy*, *charge*, *group*, *grouptype*,
   *intra_energy*, *tfac_insert*, *overlap_cutoff*, *max*, *min*, *mcmoves*,
   *rigid*, and *shake* (optional). *pressure* and *fugacity_coeff* keywords are not
   valid for this fix.

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 gas mcemc 10 1000 1000 2 29494 298.0 1200 100 0.01
   fix 3 water mcemc 10 100 100 0 3456543 3.0 1200 100 mol my_one_water maxangle 180 full_energy
   fix 4 my_gas mcemc 1 10 10 1 123456543 300.0 1200 100 1.0 region disk

Description
"""""""""""

.. versionadded:: TBD   

This fix performs mesocanonical Monte Carlo (MCEMC), also known as gauge cell
simulation, by exchanging particles with a finite-volume ideal gas
reservoir (gauge cell) at the same temperature as the system, as
discussed in :ref:`(Parashar) <Parashar>`. It also attempts Monte Carlo
moves (translations and rotations) of particles within the simulation
cell, similar to :doc:`fix gcmc <fix_gcmc>`. Specific uses of this fix
include computing adsorption isotherms in porous materials and
computing vapor–liquid equilibrium of fluids. MCEMC is complementary to
:doc:`fix gcmc <fix_gcmc>`, which performs grand canonical Monte Carlo
(GCMC) by exchanging particles with an infinite chemical-potential
reservoir. MCEMC and GCMC give identical adsorption isotherms for
microporous materials. For larger pores (> 2 nm), GCMC often shows a
hysteretic adsorption/desorption isotherm, while MCEMC yields a
reversible S-shaped van der Waals type isotherm (see :ref:`(Parashar)
<Parashar>`). The MCEMC isotherm spans stable and meta-stable states,
whereas GCMC samples only stable states. MCEMC lies between the grand
canonical ensemble (unrestricted fluctuations) and the canonical
ensemble (closed system). MCEMC reduces to GCMC when the gauge-cell volume
is infinite and to the canonical ensemble when the gauge-cell volume is
zero.

The syntax and operation are mostly identical to :doc:`fix gcmc <fix_gcmc>`,
except that fix mcemc requires *Vgauge* and *Ntotal* in place of the *mu*
(chemical potential) argument used by *fix gcmc*. For operational details
that are common to both commands (particle insertion, MC moves, region
handling, keyword usage), see the :doc:`fix gcmc <fix_gcmc>` documentation.

MCEMC Theory
"""""""""""""

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
    \mu^{id} = k_{B}T \ln(\frac{N_{gauge}\Lambda^{3}}{V_{gauge}})

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

Where :math:`V` is the volume of the simulation system, :math:`N` is the
current number of particles in the simulation system. 
In GCMC, the chemical potential of the infinite ideal gas reservoir is imposed on
the system using the exchange move, while in MCEMC, the gauge cell
(also assumed to be an ideal gas) is used to measure the chemical
potential of the system and it is not explicitly simulated. 
For calculating adsorption isotherm using GCMC, the input is the chemical
potential of the infinite reservoir and the output is the average number of particles observed in
the system. In MCEMC, there are two inputs: :math:`N_{total}` and
:math:`V_{gauge}`.  From these two inputs, the chemical potential (and
hence fugacity) of the system is calculated based on ideal gas chemical
potential in the gauge cell. The number of particles adsorbed 
is simply the average number of particles observed in the system. Note
that the choice of :math:`N_{total}` and :math:`V_{gauge}` is not
arbitrary. 

The gauge cell should be small enough to stabilize the fluid configuration
within the system yet large enough for accurate measurement of chemical 
potential. The recommendation is to calculate :math:`V_{gauge}` using the
ideal gas equation, such that the gauge cell contain roughly 70-80 particles
during a simulation. Generating a GCMC isotherm beforehand can help you choose an 
appropriate value of :math:`N_{total}`.


Restrictions
""""""""""""

Same restrictions as for :doc:`fix gcmc <fix_gcmc>` apply to this fix.


Related commands
""""""""""""""""

:doc:`fix gcmc <fix_gcmc>`,

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
