"""Bulk Na+ in SPC/E water with the fitted Na+ / water pair terms (OpenMM).

Water-water and long-range electrostatics come from SPC/E with PME. The fitted
Na+ / water interaction enters as a short-range correction: inside the cutoff a
CustomNonbondedForce adds (fitted terms) minus (SPC/E Coulomb with the Na+ charge
already handled by PME), and the Na+ Lennard-Jones term is switched off. The
correction is pair-additive, so many-body induction is only included on average.
"""

import time

import numpy as np

from .atomtypes import TYPE_INDEX
from .terms import K_COUL

SPCE_CHARGE = {"O": -0.8476, "H": 0.4238}
SPCE_VISCOSITY = 7.29e-4  # Pa s, SPC/E at 298 K (Gonzalez and Abascal 2010)
KB = 1.380649e-23
YEH_HUMMER = 2.837297

EXPRESSION = (
    "elst - coul + exch + ind + disp;"
    "elst = K*q*(1-exp(-beta*rA)*(1+beta*rA))/rA;"
    "coul = K*qs/rA;"
    "exch = A*exp(-B*rA);"
    "ind = -0.5*K*alpha/rA^4*(1-exp(-xi)*(1+xi+xi^2/2+xi^3/6+xi^4/24));"
    "disp = -C6/rA^6*(1-exp(-xd)*(1+xd+xd^2/2+xd^3/6+xd^4/24+xd^5/120+xd^6/720));"
    "xi = c*rA; xd = d*rA; rA = 10*r;"
    "q = q1*q2; qs = qs1*qs2; beta = beta1*beta2; A = A1*A2; B = B1*B2;"
    "alpha = alpha1*alpha2; c = c1*c2; C6 = C61*C62; d = d1*d2"
)
PER_PARTICLE = ("q", "qs", "beta", "A", "B", "alpha", "c", "C6", "d")


def correction_force(system, topology, params, na_indices, cutoff_nm=1.0, switch_nm=0.9):
    """Add the Na+ / water correction. Na+ carries 1 in every slot, so products pick the water value."""
    import openmm as mm

    f = mm.CustomNonbondedForce(EXPRESSION)
    f.addGlobalParameter("K", K_COUL)
    for name in PER_PARTICLE:
        f.addPerParticleParameter(name)
    na = set(na_indices)
    water = []
    for atom in topology.atoms():
        i = atom.index
        if i in na:
            f.addParticle([1.0] * len(PER_PARTICLE))
            continue
        t = TYPE_INDEX["OW" if atom.element.symbol == "O" else "HW"]
        v = params.values
        f.addParticle([v["q"][t], SPCE_CHARGE[atom.element.symbol], v["beta"][t], np.exp(v["logA"][t]),
                       v["B"][t], v["alpha"][t], v["c"][t], v["C6"][t], v["d"][t]])
        water.append(i)
    nb = next(x for x in system.getForces() if isinstance(x, mm.NonbondedForce))
    for k in range(nb.getNumExceptions()):  # the CPU platform needs identical exclusions in all nonbonded forces
        i, j = nb.getExceptionParameters(k)[:2]
        f.addExclusion(i, j)
    f.addInteractionGroup(sorted(na), water)
    f.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    f.setCutoffDistance(cutoff_nm)
    f.setUseSwitchingFunction(True)
    f.setSwitchingDistance(switch_nm)
    system.addForce(f)
    return f


def build_system(params, box_nm=3.0, n_na=4, seed=1):
    """SPC/E box with n_na waters replaced by Na+. Returns (system, topology, positions, na_indices)."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u

    ff = app.ForceField("amber14/spce.xml")
    m = app.Modeller(app.Topology(), [] * u.nanometer)
    m.addSolvent(ff, model="spce", boxSize=mm.Vec3(box_nm, box_nm, box_nm) * u.nanometer, neutralize=False)
    waters = [r for r in m.topology.residues() if r.name == "HOH"]
    rng = np.random.default_rng(seed)
    chosen = [waters[i] for i in rng.choice(len(waters), n_na, replace=False)]
    pos = m.positions
    sites = [pos[next(a.index for a in r.atoms() if a.element.symbol == "O")] for r in chosen]
    m.delete(chosen)
    top = app.Topology()
    chain = top.addChain()
    for _ in sites:
        res = top.addResidue("NA", chain)
        top.addAtom("NA", app.element.sodium, res)
    m.add(top, sites)
    system = ff.createSystem(m.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometer,
                             constraints=app.HBonds, rigidWater=True)
    na_idx = [a.index for a in m.topology.atoms() if a.residue.name == "NA"]
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    for i in na_idx:
        q, _, _ = nb.getParticleParameters(i)
        nb.setParticleParameters(i, q, 0.3 * u.nanometer, 0.0 * u.kilojoule_per_mole)
    nb.setUseDispersionCorrection(False)
    correction_force(system, m.topology, params, na_idx)
    return system, m.topology, m.positions, na_idx


def run(params, out_prefix, box_nm=3.0, n_na=4, temperature=298.15, equil_ps=200.0, prod_ps=1000.0,
        stride_ps=0.5, platform=None, seed=1):
    """NPT equilibration, then NVE production. Saves Na+ trajectory and Na-O distances to <prefix>.npz."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u

    system, top, pos, na_idx = build_system(params, box_nm, n_na, seed)
    water_o = np.array([a.index for a in top.atoms() if a.residue.name == "HOH" and a.element.symbol == "O"])
    barostat = mm.MonteCarloBarostat(1.0 * u.bar, temperature * u.kelvin, 25)
    system.addForce(barostat)
    integ = mm.LangevinMiddleIntegrator(temperature * u.kelvin, 1.0 / u.picosecond, 0.002 * u.picoseconds)
    integ.setRandomNumberSeed(seed)
    plat = mm.Platform.getPlatformByName(platform) if platform else None
    sim = app.Simulation(top, system, integ, plat) if plat else app.Simulation(top, system, integ)
    sim.context.setPositions(pos)
    sim.minimizeEnergy(maxIterations=500)
    sim.context.setVelocitiesToTemperature(temperature * u.kelvin, seed)
    sim.step(int(equil_ps / 0.002))

    # production without thermostat or barostat so that dynamics are not perturbed
    state = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    box = state.getPeriodicBoxVectors()
    system.removeForce(system.getNumForces() - 1)  # barostat was added last
    verlet = mm.VerletIntegrator(0.002 * u.picoseconds)
    sim2 = app.Simulation(top, system, verlet, plat) if plat else app.Simulation(top, system, verlet)
    sim2.context.setPeriodicBoxVectors(*box)
    sim2.context.setPositions(state.getPositions())
    sim2.context.setVelocities(state.getVelocities())
    n_frames = int(prod_ps / stride_ps)
    per = int(stride_ps / 0.002)
    traj, dist = [], []
    L = box[0][0].value_in_unit(u.nanometer)
    t0 = time.time()
    for i in range(n_frames):
        sim2.step(per)
        x = sim2.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(u.nanometer)
        traj.append(x[na_idx])
        d = x[water_o][None, :, :] - x[na_idx][:, None, :]
        d -= L * np.round(d / L)
        dist.append(np.linalg.norm(d, axis=-1).astype(np.float32))
        if (i + 1) % 100 == 0:
            np.savez_compressed(str(out_prefix) + "_partial", na=np.array(traj), dist=np.array(dist), box=L,
                                stride_ps=stride_ps)
            print(f"production frame {i + 1}/{n_frames}, {time.time() - t0:.0f} s", flush=True)
    np.savez_compressed(out_prefix, na=np.array(traj), dist=np.array(dist), box=L, stride_ps=stride_ps)
    return str(out_prefix) + ".npz"


def unwrap(traj, box):
    """Unwrap a periodic (frames, ions, 3) trajectory using minimum-image steps."""
    step = np.diff(traj, axis=0)
    step -= box * np.round(step / box)
    return np.concatenate([traj[:1], traj[:1] + np.cumsum(step, axis=0)])


def msd_fft(x):
    """Mean squared displacement for every lag, averaged over time origins and ions.

    x has shape (frames, ions, 3) and must be unwrapped. FFT autocorrelation
    keeps the cost at O(N log N) per ion.
    """
    n, ions = x.shape[:2]
    sq = np.sum(x**2, axis=-1)
    pad = np.concatenate([sq, np.zeros((1, ions))])
    f = np.fft.rfft(x, n=2 * n, axis=0)
    acf = np.fft.irfft(f * np.conj(f), n=2 * n, axis=0)[:n].sum(axis=-1)
    cross = acf / (n - np.arange(n))[:, None]
    q = 2.0 * pad.sum(axis=0)
    first = np.empty((n, ions))
    for m in range(n):
        q = q - pad[m - 1] - pad[n - m]
        first[m] = q / (n - m)
    return (first - 2.0 * cross).mean(axis=1)


def diffusion(traj, box, stride_ps, fit_from_ps=5.0, fit_to_frac=0.2, fit_to_max_ps=100.0, n_blocks=4,
              temperature=298.15):
    """Diffusion coefficient (m^2/s) from the Einstein relation with Yeh-Hummer correction.

    The MSD slope is fitted between fit_from_ps (after the ballistic and caging regime) and the
    smaller of fit_to_frac of the run length and fit_to_max_ps, since long lags are statistically poor.
    Returns (D_raw, D_corrected, standard error from n_blocks independent stretches of the run).
    """
    x = unwrap(traj, box)

    def d_of(seg):
        n = len(seg)
        msd = msd_fft(seg)[: n // 2]
        t = np.arange(n // 2) * stride_ps
        t_hi = min(fit_to_max_ps, fit_to_frac * n * stride_ps)
        t_lo = min(fit_from_ps, 0.25 * t_hi)
        sel = (t >= t_lo) & (t <= t_hi)
        slope = np.polyfit(t[sel], msd[sel], 1)[0]  # nm^2/ps
        return slope / 6.0 * 1e-6  # nm^2/ps to m^2/s

    d_raw = d_of(x)
    d_blocks = np.array([d_of(b) for b in np.array_split(x, n_blocks)])
    corr = YEH_HUMMER * KB * temperature / (6 * np.pi * SPCE_VISCOSITY * box * 1e-9)
    return d_raw, d_raw + corr, float(d_blocks.std(ddof=1) / np.sqrt(n_blocks))


def rdf_and_cn(dist, box, r_max=0.8, n_bins=160):
    """Na-O g(r) in nm and the running coordination number, from stored distances (frames, ions, waters)."""
    edges = np.linspace(0.0, r_max, n_bins + 1)
    counts = np.histogram(dist.ravel(), bins=edges)[0] / (dist.shape[0] * dist.shape[1])
    r = 0.5 * (edges[1:] + edges[:-1])
    rho = dist.shape[2] / box**3
    shell = 4.0 * np.pi * r**2 * (edges[1] - edges[0])
    g = counts / (rho * shell)
    cn = np.cumsum(counts)
    return r, g, cn


def first_shell(r, g, cn):
    """Position of the first g(r) maximum, first minimum after it, and CN at that minimum."""
    i_max = int(np.argmax(g[: len(g) // 2]))
    tail = g[i_max:]
    i_min = i_max + int(np.argmin(tail[: len(tail) // 2]))
    return r[i_max], r[i_min], cn[i_min]
