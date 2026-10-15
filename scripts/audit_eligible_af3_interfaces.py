#!/usr/bin/env python3

import csv
import statistics
from collections import Counter
from pathlib import Path

import gemmi


INPUT = Path("external/baldo_model_inventory/af3_candidate_db_eligibility.tsv")
OUTPUT = Path("external/baldo_model_inventory/af3_interface_audit.tsv")
CUTOFF = 4.5

FIELDS = [
    "tf_id",
    "directory_name",
    "variant",
    "model_cif",
    "result",
    "protein_chains",
    "nucleic_chains",
    "atom_contacts",
    "protein_interface_residues",
    "nucleic_interface_residues",
    "protein_ca_plddt_mean",
    "protein_ca_plddt_min",
    "protein_ca_plddt_max",
    "protein_ca_count",
    "error",
]




def audit_model(tf_id: str, path: Path) -> dict:
    directory_name = path.parent.name

    if directory_name.lower().endswith("_dimer"):
        variant = "dimer"
    elif directory_name.lower().endswith("_monomer"):
        variant = "monomer"
    else:
        variant = "base"

    record = {
        "tf_id": tf_id,
        "directory_name": directory_name,
        "variant": variant,
        "model_cif": str(path),
        "result": "",
        "protein_chains": "",
        "nucleic_chains": "",
        "atom_contacts": 0,
        "protein_interface_residues": 0,
        "nucleic_interface_residues": 0,
        "protein_ca_plddt_mean": "",
        "protein_ca_plddt_min": "",
        "protein_ca_plddt_max": "",
        "protein_ca_count": 0,
        "error": "",
    }

    if not path.is_file():
        record["result"] = "ERROR_MISSING_FILE"
        record["error"] = "File does not exist"
        return record

    try:
        structure = gemmi.read_structure(str(path))
        structure.setup_entities()

        if len(structure) == 0:
            record["result"] = "ERROR_NO_MODEL"
            record["error"] = "Structure contains no models"
            return record

        model = structure[0]
        protein_chain_indices = set()
        nucleic_chain_indices = set()

        for chain_index, chain in enumerate(model):
            polymer = chain.get_polymer()
            if len(polymer) == 0:
                continue

            polymer_type = str(polymer.check_polymer_type())

            if "Peptide" in polymer_type:
                protein_chain_indices.add(chain_index)

            if "Dna" in polymer_type or "Rna" in polymer_type:
                nucleic_chain_indices.add(chain_index)

        protein_names = [
            model[index].name for index in sorted(protein_chain_indices)
        ]
        nucleic_names = [
            model[index].name for index in sorted(nucleic_chain_indices)
        ]

        record["protein_chains"] = ",".join(protein_names)
        record["nucleic_chains"] = ",".join(nucleic_names)

        if not protein_chain_indices:
            record["result"] = "FAIL_NO_PROTEIN"
            return record

        if not nucleic_chain_indices:
            record["result"] = "FAIL_NO_DNA"
            return record

        neighbor_search = gemmi.NeighborSearch(
            model,
            structure.cell,
            CUTOFF,
        ).populate(include_h=False)

        atom_contacts = set()
        protein_contact_residues = set()
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

                    for mark in neighbor_search.find_neighbors(
                        atom,
                        min_dist=0.1,
                        max_dist=CUTOFF,
                    ):
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

        record["atom_contacts"] = len(atom_contacts)
        record["protein_interface_residues"] = len(protein_contact_residues)
        record["nucleic_interface_residues"] = len(nucleic_contact_residues)

        if ca_plddt:
            record["protein_ca_plddt_mean"] = f"{statistics.mean(ca_plddt):.3f}"
            record["protein_ca_plddt_min"] = f"{min(ca_plddt):.3f}"
            record["protein_ca_plddt_max"] = f"{max(ca_plddt):.3f}"
            record["protein_ca_count"] = len(ca_plddt)

        if atom_contacts:
            record["result"] = "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"
        else:
            record["result"] = "FAIL_EMPTY_PROTEIN_DNA_INTERFACE"

    except Exception as exc:
        record["result"] = "ERROR"
        record["error"] = f"{type(exc).__name__}: {exc}"

    return record


def main() -> None:
    if not INPUT.is_file():
        raise SystemExit(f"Missing input: {INPUT}")

    candidates = []

    with INPUT.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["eligible_for_af3_import"].lower() != "yes":
                continue

            for value in row["selected_cifs"].split(";"):
                value = value.strip()
                if value:
                    candidates.append((row["tf_id"], Path(value)))

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()

        for index, (tf_id, path) in enumerate(candidates, start=1):
            record = audit_model(tf_id, path)
            writer.writerow(record)
            handle.flush()

            counts[record["result"]] += 1

            if index % 25 == 0 or index == len(candidates):
                print(f"Processed {index}/{len(candidates)}")

    print("\nModels audited:", len(candidates))
    print("Unique eligible TFs:", len({tf_id for tf_id, _ in candidates}))

    for result, count in sorted(counts.items()):
        print(f"{result}: {count}")

    print("\nWritten:", OUTPUT)


if __name__ == "__main__":
    main()
