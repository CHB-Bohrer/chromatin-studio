"""Read-only physical binding estimates, separate from simulated kinetics.

The simulator's ``on_rate`` is a finite-radius, per-eligible-pair hazard in
inverse reduced time. It is not a molar association constant. This calculator
instead assumes dilute, independent 3D spherical reactants and an explicitly
supplied length/time calibration. It does not set or infer a thermodynamic Kd
for the interacting polymer simulation.

For unlike species, the Smoluchowski pair kernel is also the complex formation
coefficient. For identical species, unordered collisions give half that kernel
as the coefficient of [A]^2 in d[A2]/dt. Kd is defined as [A]^2/[A2], not a
site-labeled or monomer-depletion convention.
"""

import math


AVOGADRO = 6.02214076e23  # Exact SI definition, mol^-1.
LITERS_PER_UM3 = 1e-15

SOURCES = (
    {
        "title": "Donev, Yang & Kim (2018), Efficient reactive Brownian dynamics",
        "url": "https://doi.org/10.1063/1.5009464",
        "supports": "Finite-radius Doi hazards, diffusion-limited association, identical-particle counting, and detailed-balance requirements.",
    },
    {
        "title": "Vijaykumar, Bolhuis & ten Wolde (2016), The intrinsic rate constants in diffusion-influenced reactions",
        "url": "https://doi.org/10.1039/C6FD00104A",
        "supports": "Relative diffusion Da+Db and the distinction between intrinsic contact kinetics and effective association/dissociation rates.",
    },
    {
        "title": "Lawley & Keener (2016), Including rebinding reactions in well-mixed models of distributive biochemical reactions",
        "url": "https://doi.org/10.1016/j.bpj.2016.10.008",
        "supports": "Rapid rebinding can alter apparent kinetics and requires more than directly replacing a rate constant.",
    },
)

LIMITATIONS = (
    "Estimate only: no configuration, trajectory, reaction rate, or physical calibration is changed.",
    "The length and time scales are user-supplied assumptions; genomic bp per bead does not determine a physical length or time scale.",
    "4π(Da+Db)R is the dilute, unconfined, independent-sphere 3D encounter limit without interaction potentials. It is not a universal upper bound in attractive fields.",
    "Efficiency is a phenomenological multiplier of this encounter limit, not the simulator's contact probability or on_rate.",
    "The existing on_rate is an inverse-time hazard for each eligible nearby pair, not M^-1 s^-1. Its value does not enter this diffusion-limited estimate.",
    "The nominal koff enforces the requested Kd only in an ideal two-state, well-mixed association model. It is not a calibrated replacement for the current microscopic bond-breaking hazard.",
    "Binding, polymer tethering, crowding, repeated rebinding, and multivalent assemblies change measured diffusion and residence times. A slowed bound/cluster diffusion coefficient does not define an intrinsic Kd; apparent multivalent avidity is not a single-pair Kd.",
)


def _number(value, name, *, positive=True, maximum=None):
    """Validate JSON numbers without accepting booleans, strings or infinities."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number.")
    try:
        value = float(value)
    except (OverflowError, ValueError):
        raise ValueError(f"{name} must be a finite number.") from None
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number.")
    if value <= 0 if positive else value < 0:
        relation = "greater than zero" if positive else "nonnegative"
        raise ValueError(f"{name} must be {relation}.")
    if maximum is not None and value > maximum:
        raise ValueError(f"{name} must be at most {maximum:g}.")
    return value


def _result(value, name, *, positive=True):
    """Reject overflow/underflow rather than emitting non-JSON numbers."""
    if not math.isfinite(value) or (value <= 0 if positive else value < 0):
        raise ValueError(
            f"The supplied calibration gives an unrepresentable {name}; "
            "use less extreme physical length, time, Kd, or efficiency values."
        )
    return value


def estimate_binding(config, request):
    """Estimate ideal association and nominal dissociation without mutation.

    ``config`` is the validated, applied model configuration. ``request`` must
    provide species IDs ``a`` and ``b``, ``length_nm`` per reduced length unit,
    ``time_s`` per reduced time unit, and a positive ``target_kd_um`` in μM.
    Optional ``efficiency`` is in (0, 1], defaulting to the absorbing limit 1.

    ``b='chromatin'`` selects a hypothetical unmodified chromatin-bead target,
    using the exact engine capture cutoff protein.radius + chromatin.bond_length.
    This is an isolated-sphere proxy using bare bead mobility, not a prediction
    of tethered DNA-site kinetics. A species actually named 'chromatin' makes
    that selector ambiguous and is rejected explicitly.

    All invalid inputs raise ValueError. No simulator, RNG, or mutable global
    state is accessed. Output is JSON serializable with finite numerical values.
    """
    if not isinstance(config, dict):
        raise ValueError("Applied configuration must be an object.")
    if not isinstance(request, dict):
        raise ValueError("Kinetics request must be an object.")
    a, b = request.get("a"), request.get("b")
    if not isinstance(a, str) or not a or not isinstance(b, str) or not b:
        raise ValueError("a and b must be nonempty species identifiers.")
    inputs = {
        "length_nm": _number(request.get("length_nm"), "length_nm"),
        "time_s": _number(request.get("time_s"), "time_s"),
        "target_kd_um": _number(request.get("target_kd_um"), "target_kd_um"),
        "efficiency": _number(request.get("efficiency", 1.0), "efficiency", maximum=1.0),
    }
    proteins = {p["id"]: p for p in config.get("proteins", [])}
    if a not in proteins:
        raise ValueError(f"Unknown applied protein species a: {a}.")
    protein_a = proteins[a]
    diffusion_a = _number(protein_a.get("diffusion"), f"{a}.diffusion", positive=False)
    limitations = list(LIMITATIONS)
    model = {"diffusion_a_reduced": diffusion_a,
             "diffusion_a_basis": "free molecule center-of-mass diffusion"}
    if b == "chromatin":
        if b in proteins:
            raise ValueError("The species ID 'chromatin' conflicts with the chromatin-target selector; rename that species to use this calculator.")
        chromatin = config.get("chromatin", {})
        diffusion_b = _number(chromatin.get("diffusion"), "chromatin.diffusion", positive=False)
        radius = (_number(protein_a.get("radius"), f"{a}.radius")
                  + _number(chromatin.get("bond_length"), "chromatin.bond_length"))
        on_rate = _number(protein_a.get("on_rate"), f"{a}.on_rate", positive=False)
        off_rate = _number(protein_a.get("off_rate"), f"{a}.off_rate", positive=False)
        baseline = _number(protein_a.get("chromatin_affinity", {}).get("unmodified", 0.0),
                           f"{a}.chromatin_affinity.unmodified", positive=False)
        model.update(radius_source="protein.radius + chromatin.bond_length",
                     diffusion_b_basis="bare chromatin bead diffusion, not polymer effective diffusion",
                     unmodified_affinity_multiplier=baseline)
        kind, b_name, homodimer = "protein_chromatin_proxy", "Unmodified chromatin bead (sphere proxy)", False
        limitations.append(
            "Chromatin target: an independent spherical-bead proxy using the exact engine encounter cutoff and bare bead diffusion. A tethered polymer site is not an independently diffusing sphere; this does not infer a DNA-site Kd. Histone-mark multipliers are not converted into physical affinities."
        )
        if baseline == 0 or on_rate == 0:
            limitations.append("The applied model currently disables unmodified chromatin association for this species; the physical estimate remains hypothetical.")
        involved_proteins = [protein_a]
    else:
        if b not in proteins:
            raise ValueError(f"Unknown applied protein species b: {b}.")
        protein_b = proteins[b]
        rule = next((rule for rule in config.get("interactions", [])
                     if sorted((rule["a"], rule["b"])) == sorted((a, b))), None)
        if rule is None:
            raise ValueError(f"No applied protein interaction rule exists for {a} and {b}.")
        diffusion_b = _number(protein_b.get("diffusion"), f"{b}.diffusion", positive=False)
        radius = _number(rule.get("radius"), "interaction.radius")
        on_rate = _number(rule.get("on_rate"), "interaction.on_rate", positive=False)
        off_rate = _number(rule.get("off_rate"), "interaction.off_rate", positive=False)
        model.update(radius_source="interaction.radius",
                     diffusion_b_basis="free molecule center-of-mass diffusion")
        kind, b_name, homodimer = "protein_protein", protein_b.get("name", b), a == b
        involved_proteins = [protein_a, protein_b]
        if on_rate == 0:
            limitations.append("The applied interaction has zero association hazard; the physical estimate remains hypothetical.")
    if any(p.get("representation") == "idr" for p in involved_proteins):
        limitations.append(
            "IDR approximation: diffusion uses each whole molecule's configured free center-of-mass coefficient, while the capture radius acts between sticker/domain beads in the simulator. An extended multivalent chain is not a single reactive sphere; no bead-count multiplier is applied to the molecular diffusion."
        )
    if homodimer:
        limitations.append(
            "Homodimer convention: d[A2]/dt = (k_pair/2)[A]^2 - koff[A2], with Kd=[A]^2/[A2]. The factor 1/2 counts identical unordered pairs once; monomer depletion and site-labeled conventions use different coefficients."
        )
    diffusion_sum = _result(diffusion_a + diffusion_b, "positive relative diffusion")
    model.update(diffusion_b_reduced=diffusion_b, relative_diffusion_reduced=diffusion_sum,
                 radius_reduced=radius, on_rate_reduced=on_rate, off_rate_reduced=off_rate)

    length_um = _result(inputs["length_nm"] / 1000.0, "length scale")
    time_s = inputs["time_s"]
    diffusion_scale = _result(length_um * length_um / time_s, "diffusion scale")
    da = _result(diffusion_a * diffusion_scale, "diffusion of a", positive=diffusion_a > 0)
    db = _result(diffusion_b * diffusion_scale, "diffusion of b", positive=diffusion_b > 0)
    relative_d = _result(da + db, "relative diffusion")
    physical_radius = _result(radius * length_um, "capture radius")
    volume_kernel = _result(4.0 * math.pi * relative_d * physical_radius, "encounter kernel")
    molar_kernel_limit = _result(volume_kernel * (AVOGADRO * LITERS_PER_UM3), "molar encounter kernel")
    molar_kernel = _result(inputs["efficiency"] * molar_kernel_limit, "efficiency-adjusted encounter kernel")
    factor = 0.5 if homodimer else 1.0
    formation_rate = _result(factor * molar_kernel, "complex formation coefficient")
    kd_m = _result(inputs["target_kd_um"] * 1e-6, "target Kd")
    koff = _result(kd_m * formation_rate, "nominal dissociation rate")
    koff_reduced = _result(koff * time_s, "nominal reduced dissociation rate")
    pair = {"a": a, "b": b, "a_name": protein_a.get("name", a), "b_name": b_name,
            "kind": kind, "homodimer": homodimer,
            "kd_convention": "[A]^2/[A2]" if homodimer else "[A][B]/[AB]",
            "formation_equation": "d[A2]/dt = macroscopic_kon*[A]^2 - nominal_koff*[A2]" if homodimer else "d[AB]/dt = macroscopic_kon*[A]*[B] - nominal_koff*[AB]"}
    return {
        "status": "estimate_only", "pair": pair, "inputs": inputs, "model": model,
        "length_um_per_reduced_unit": length_um,
        "physical_D_a_um2_s": da, "physical_D_b_um2_s": db,
        "physical_relative_D_um2_s": relative_d, "physical_radius_um": physical_radius,
        "kon_limit_um3_s": volume_kernel,
        "kon_limit_M_inv_s": molar_kernel_limit, "kon_M_inv_s": molar_kernel,
        "association_stoichiometric_factor": factor,
        "macroscopic_kon_M_inv_s": formation_rate,
        "nominal_koff_s_inv": koff, "nominal_koff_reduced": koff_reduced,
        "limitations": limitations, "sources": [dict(source) for source in SOURCES],
    }
