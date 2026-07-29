#!/usr/bin/env python3
"""
Generate a MICM mechanism configuration (v0 "camp-data" JSON) and the
reaction map companion file for a compiled CAM chemistry mechanism, for use
with the -micm CAM configure option (mo_micm.F90).

Every compiled reaction is emitted as a rate-injected slot (USER_DEFINED, or
EMISSION when no reactant is a solution species): CAM computes all rate
constants (setrxt/usrrxt/photolysis/adjrxt) exactly as for the built-in
solvers and mo_micm injects them each timestep, so MICM performs only the
implicit integration. Reactants and products keep only solution species --
invariant (Fixed) reactants and products are dropped, because adjrxt folds
invariant concentrations into the injected rate constant and invariants are
not part of the solved system.

Additionally, a FIRST_ORDER_LOSS slot is emitted for every solution species
(heterogeneous washout rates, selected at runtime via gas_wetdep_list) and
an EMISSION slot for every external forcing species (extfrc). Unused slots
are injected as zero at runtime.

The reaction map companion file gives, for each compiled reaction in global
index order, the number of solution-species reactant molecules (the unit
conversion exponent used by mo_micm) and the MICM rate parameter name.

Usage: gen_micm_config.py <pp_mechanism_dir> [output_dir]
       (default output_dir: <pp_mechanism_dir>/micm)
"""

import json
import os
import re
import sys


def parse_chem_mech(path):
    """Parse a chem_proc mechanism input file (chem_mech.in)."""
    solution = []   # solution species names, in order (= solsym)
    fixed = []      # invariant species names (= inv_lst)
    photolysis = [] # (tag, reactants, products) in order
    reactions = []  # (tag_or_None, reactants, products) in order
    ext_forcing = []# external forcing species names, in order (= extfrc_lst)

    section = None
    with open(path) as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()
        if not line or line.startswith('*'):
            continue
        upper = line.upper()
        # section transitions (End before Begin, since "End X" contains "X")
        if upper.startswith('END ') or upper.startswith('END'):
            token = upper.replace(' ', '')
            if token in ('ENDSOLUTION', 'ENDFIXED', 'ENDPHOTOLYSIS',
                         'ENDREACTIONS', 'ENDEXTFORCING'):
                section = None
                continue
        if upper == 'SOLUTION':
            section = 'solution'
            continue
        if upper == 'FIXED':
            section = 'fixed'
            continue
        if upper == 'PHOTOLYSIS':
            section = 'photolysis'
            continue
        if upper == 'REACTIONS':
            section = 'reactions'
            continue
        if upper == 'EXT FORCING':
            section = 'ext'
            continue

        if section == 'solution' or section == 'fixed':
            for entry in line.split(','):
                entry = entry.strip()
                if not entry:
                    continue
                name = entry.split('->')[0].strip()
                if name:
                    (solution if section == 'solution' else fixed).append(name)
        elif section in ('photolysis', 'reactions'):
            target = photolysis if section == 'photolysis' else reactions
            if '->' not in line and not line.startswith('['):
                # continuation of the previous reaction's product list
                if not target:
                    raise ValueError(f'unexpected line in {section} section: {line}')
                target[-1] += ' ' + line
            else:
                target.append(line)
        elif section == 'ext':
            for entry in line.split(','):
                name = entry.split('<-')[0].strip()
                if name:
                    ext_forcing.append(name)

    photolysis = [parse_reaction(ln) for ln in photolysis]
    reactions = [parse_reaction(ln) for ln in reactions]
    return solution, fixed, photolysis, reactions, ext_forcing


def parse_reaction(line):
    """Parse one reaction line into (tag_or_None, reactants, products)."""
    tag = None
    if line.startswith('['):
        close = line.index(']')
        tag = line[1:close]
        # strip alias/annotation syntax: [jsoa_a1->,.0004*jno2],
        # [O1D_N2,cph=189.81], [jo2_a=userdefined,]
        tag = tag.split('->')[0].split(',')[0].split('=')[0].strip()
        line = line[close + 1:].strip()
    if '->' not in line:
        raise ValueError(f'not a reaction line: {line}')
    # rate coefficients after ';' are not needed (rates are injected)
    eqn = line.split(';')[0]
    lhs, _, rhs = eqn.partition('->')
    return tag, parse_side(lhs), parse_side(rhs)


def parse_side(side):
    """Parse one side of a reaction equation into [(coefficient, species)]."""
    out = []
    for term in side.split('+'):
        term = term.strip()
        if not term or term == 'hv':
            continue
        # coefficient and species separated by '*' or (older mechanisms)
        # just whitespace: "0.5*SO2", ".25 GLYOXAL"
        m = re.match(r'^(?:([0-9.]+(?:[eE][+-]?\d+)?)\s*(?:\*\s*)?)?(\{?[A-Za-z][A-Za-z0-9_]*\}?)$', term)
        if m is None:
            raise ValueError(f'cannot parse reaction term: {term!r}')
        coeff = float(m.group(1)) if m.group(1) else 1.0
        name = m.group(2).strip('{}')  # {X} marks non-transported bookkeeping
        out.append((coeff, name))
    return out


def cross_validate(pp_dir, solution, all_reactions, ext_forcing):
    """Check the parse against the mechanism's generated Fortran: counts from
    chem_mods.F90, and reaction tags at their global indices from
    mo_sim_dat.F90's rxt_tag_lst/rxt_tag_map (complete for fully-tagged
    mechanisms, spot checks otherwise)."""
    # chem_mods declares its parameters as one continued multi-name list
    # (name = value, & per line)
    chem_mods = open(os.path.join(pp_dir, 'chem_mods.F90')).read()
    counts = dict(re.findall(r'(\w+)\s*=\s*(\d+)\s*,?\s*&?\s*(?:!|$)',
                             chem_mods, re.M))
    for key, have in (('rxntot', len(all_reactions)),
                      ('gas_pcnst', len(solution)),
                      ('extcnt', len(ext_forcing))):
        want = int(counts[key])
        if have != want:
            raise ValueError(
                f'parsed {key} = {have} but chem_mods.F90 has {want}')

    sim_dat = open(os.path.join(pp_dir, 'mo_sim_dat.F90')).read()
    tags = []
    for chunk in re.findall(r'rxt_tag_lst\([^)]*\)\s*=\s*\(/(.*?)/\)',
                            sim_dat, re.S):
        tags += [t.strip() for t in re.findall(r"'([^']*)'", chunk)]
    tag_map = []
    for chunk in re.findall(r'rxt_tag_map\([^)]*\)\s*=\s*\(/(.*?)/\)',
                            sim_dat, re.S):
        tag_map += [int(n) for n in re.findall(r'\d+', chunk)]
    if len(tags) != len(tag_map):
        raise ValueError('cannot parse rxt_tag_lst/rxt_tag_map from mo_sim_dat.F90')
    for tag, gidx in zip(tags, tag_map):
        parsed_tag = all_reactions[gidx - 1][0]
        if parsed_tag != tag:
            raise ValueError(
                f'reaction {gidx}: parsed tag {parsed_tag!r} does not match '
                f'mo_sim_dat rxt_tag_lst {tag!r}')


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    pp_dir = sys.argv[1].rstrip('/')
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(pp_dir, 'micm')
    mech_name = os.path.basename(pp_dir)

    solution, fixed, photolysis, reactions, ext_forcing = \
        parse_chem_mech(os.path.join(pp_dir, 'chem_mech.in'))
    solution_set = set(solution)

    # global reaction index order: photolysis block first, then reactions
    all_reactions = [(tag, r, p) for (tag, r, p) in photolysis] + reactions
    cross_validate(pp_dir, solution, all_reactions, ext_forcing)

    species_json = {'camp-data': []}
    for name in solution:
        species_json['camp-data'].append({
            'name': name,
            'type': 'CHEM_SPEC',
            'absolute tolerance': 1e-12,
        })

    mech_reactions = []
    # map entries: (index, n_solution_reactants, micm_rate_parameter_name,
    # yield). A reaction with no solution-species reactants becomes one
    # EMISSION slot per solution product (mo_micm injects rate*yield into
    # each); a reaction with no solution reactants AND no solution products
    # is invisible to the solved system and maps to the no-op name NONE.
    rxt_map = []
    for idx, (tag, reactants, products) in enumerate(all_reactions, start=1):
        label = tag if tag else f'rxn{idx}'
        if re.search(r'[,\s]', label):
            # mo_micm's list-directed map read cannot represent these
            raise ValueError(f'reaction label contains comma/whitespace: {label!r}')
        sol_reactants = [(c, s) for (c, s) in reactants if s in solution_set]
        sol_products = [(c, s) for (c, s) in products if s in solution_set]
        n_sol = int(round(sum(c for (c, s) in sol_reactants)))
        if abs(n_sol - sum(c for (c, s) in sol_reactants)) > 1e-9:
            raise ValueError(f'non-integer reactant count in reaction {label}')

        if n_sol == 0:
            if not sol_products:
                rxt_map.append((idx, 0, 'NONE', 0.0))
                continue
            for ip, (coeff, name) in enumerate(sol_products, start=1):
                sub_label = label if len(sol_products) == 1 else f'{label}_p{ip}'
                mech_reactions.append({
                    'type': 'EMISSION',
                    'MUSICA name': sub_label,
                    'species': name,
                })
                rxt_map.append((idx, n_sol, f'EMIS.{sub_label}', coeff))
        else:
            entry = {
                'type': 'USER_DEFINED',
                'MUSICA name': label,
                'reactants': {},
                'products': {},
            }
            rcounts = {}
            for coeff, name in sol_reactants:
                rcounts[name] = rcounts.get(name, 0) + int(round(coeff))
            for name, qty in rcounts.items():
                entry['reactants'][name] = {'qty': qty} if qty != 1 else {}
            pcounts = {}
            for coeff, name in sol_products:
                pcounts[name] = pcounts.get(name, 0.0) + coeff
            for name, yld in pcounts.items():
                entry['products'][name] = {'yield': yld} if yld != 1.0 else {}
            mech_reactions.append(entry)
            rxt_map.append((idx, n_sol, f'USER.{label}', 1.0))

    # first-order loss slot per solution species (washout / sethet)
    for name in solution:
        mech_reactions.append({
            'type': 'FIRST_ORDER_LOSS',
            'MUSICA name': name,
            'species': name,
        })
    # emission slot per external forcing species (setext)
    for name in ext_forcing:
        if name not in solution_set:
            raise ValueError(f'external forcing species {name} is not a solution species')
        mech_reactions.append({
            'type': 'EMISSION',
            'MUSICA name': f'ext_{name}',
            'species': name,
        })

    reactions_json = {'camp-data': [{
        'name': mech_name,
        'type': 'MECHANISM',
        'reactions': mech_reactions,
    }]}

    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'species.json'), 'w') as f:
        json.dump(species_json, f, indent=2)
        f.write('\n')
    with open(os.path.join(out_dir, 'reactions.json'), 'w') as f:
        json.dump(reactions_json, f, indent=2)
        f.write('\n')
    with open(os.path.join(out_dir, 'config.json'), 'w') as f:
        json.dump({'camp-files': ['species.json', 'reactions.json']}, f, indent=2)
        f.write('\n')
    with open(os.path.join(out_dir, 'rxt_map.txt'), 'w') as f:
        f.write(f'# MICM reaction map for {mech_name} -- generated by gen_micm_config.py\n')
        f.write('# rxntot gas_pcnst extcnt n_entries\n')
        f.write(f'{len(all_reactions)} {len(solution)} {len(ext_forcing)} {len(rxt_map)}\n')
        f.write('# reaction_index n_solution_reactants micm_rate_parameter_name yield\n')
        for idx, n_sol, name, yld in rxt_map:
            f.write(f'{idx} {n_sol} {name} {yld}\n')

    n_rxt_slots = len(mech_reactions) - len(solution) - len(ext_forcing)
    print(f'{mech_name}: {len(all_reactions)} reactions '
          f'({len(photolysis)} photolysis), {len(solution)} solution species, '
          f'{len(fixed)} invariants, {len(ext_forcing)} external forcings')
    print(f'rate parameters: {n_rxt_slots} reaction + {len(solution)} loss '
          f'+ {len(ext_forcing)} emission = {len(mech_reactions)}')
    print(f'wrote {out_dir}/{{config,species,reactions}}.json and rxt_map.txt')


if __name__ == '__main__':
    main()
