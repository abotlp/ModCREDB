#!/usr/bin/env python3

import statistics
import sys
from pathlib import Path

import gemmi


def main() -> int:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} MODEL.cif", file=sys.stderr)
        return 2

    path = Path(sys.argv[1])
    if not path.is_file():
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        return 1

    structure = gemmi.read_structure(str(path))
    structure.setup_entities()

    if len(structure) == 0:
        print("ERROR: structure contains no models", file=sys.stderr)
        return 1

    model = structure[0]

    protein_chain_indices = set()
    nucleic_chain_indices = set()

    print("FILE:", path)
    print("CHAINS:")

    for chain_index, chain in enumerate(model):
        polymer = chain.get_polymer()

        if len(polymer) == 0:
            print(
                f"  {chain.name}: non-polymer/empty "
                f"residues={len(chain)} atoms={chain.count_atom_sites()}"
            )
            continue

        polymer_type = str(polymer.check_polymer_type())

        if "Peptide" in polymer_type:
            protein_chain_indices.add(chain_index)

        if "Dna" in polymer_type or "Rna" in polymer_type:
            nucleic_chain_indices.add(chain_index)

        print(
            f"  {chain.name}: {polymer_type} "
            f"polymer_residues={len(polymer)} atoms={chain.count_atom_sites()}"
        )

    protein_names = [model[i].name for i in sorted(protein_chain_indices)]
    nucleic_names = [model[i].name for i in sorted(nucleic_chain_indices)]

    print("PROTEIN CHAINS:", ",".join(protein_names) or "NONE")
    print("DNA/RNA CHAINS:", ",".join(nucleic_names) or "NONE")

    if not protein_chain_indices:
        print("RESULT: FAIL_NO_PROTEIN")
        return 1

    if not nucleic_chain_indices:
        print("RESULT: FAIL_NO_DNA")
        return 0

    cutoff = 4.5
    neighbor_search = gemmi.NeighborSearch(
        model,
        structure.cell,
        cutoff,
    ).populate(include_h=False)

    atom_contacts = set()
    protein_contact_residues = set()
    nucleic_contact_residues = set()
    ca_plddt = []

    for protein_chain_index in sorted(protein_chain_indices):
        chain_residues = set()
    nucleic_contact_residues = set()
    ca_plddt = []

    for protein_chain_index in sorted(protein_chain_indices):
        chain = model[protein_chain_index]

        for protein_residue_index, residue in enumerate(chain):
            for protein_atom_index, atom in enumerate(residue):
                if atom.is_hydrogen():
                    continue

                if atom.name.strip() == "CA":
                    ca_plddt.append(float(atom.b_iso))

                marks = neighbor_search.find_neighbors(
                    atom,
                    min_dist=0.1,
                    max_dist=cutoff,
                )

                for mark in marks:
                    if mark.image_idx != 0:
                        continue

                    if mark.chain_idx not in nucleic_chain_indices:
                        continue

                    dna_cra = mark.to_cra(model)

                    atom_contacts.add(
                        (
                            protein_chain_index,
                            protein_residue_index,
                            protein_atom_index,
                            mark.chain_idx,
                            mark.residue_idx,
                            mark.atom_idx,
                        )
                    )

                    protein_contact_residues.add(
                        (
                            chain.name,
                            str(residue.seqid),
                            residue.name,
                        )
                    )

                    nucleic_contact_residues.add(
                        (
                            dna_cra.chain.name,
                            str(dna_cra.residue.seqid),
                            dna_cra.residue.name,
                        )
                    )

    print(f"CUTOFF_A: {cutoff}")
    print("ATOM_CONTACTS:", len(atom_contacts))
    print("PROTEIN_INTERFACE_RESIDUES:", len(protein_contact_residues))
    print("DNA_INTERFACE_RESIDUES:", len(nucleic_contact_residues))

    if ca_plddt:
        print(
            "PROTEIN_CA_PLDDT:",
            f"mean={statistics.mean(ca_plddt):.2f}",
            f"min={min(ca_plddt):.2f}",
            f"max={max(ca_plddt):.2f}",
            f"n={len(ca_plddt)}",
        )

    if atom_contacts:
        print("RESULT: PASS_NONEMPTY_PROTEIN_DNA_INTERFACE")
    else:
        print("RESULT: FAIL_EMPTY_PROTEIN_DNA_INTERFACE")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
