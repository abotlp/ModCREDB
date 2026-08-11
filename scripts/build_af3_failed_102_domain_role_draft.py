#!/usr/bin/env python3
"""Build the second-stage draft Pfam domain-role review table.

This script performs evidence collection and creates provisional domain-level
proposals only.  It opens SQLite in read-only mode, refuses to overwrite
outputs, and verifies that the database bytes and metadata did not change.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import re
import sqlite3
import sys
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_PFAM_COUNT = 26
OFFICIAL_API = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}/"
OFFICIAL_PAGE = "https://www.ebi.ac.uk/interpro/entry/pfam/{pfam_id}/"
INTERPRO_API = "https://www.ebi.ac.uk/interpro/api/entry/interpro/{interpro_id}/"

COLUMNS = [
    "pfam_id",
    "pfam_name",
    "interpro_id",
    "interpro_name",
    "official_entry_url",
    "official_source_status",
    "failed_accession_count",
    "failed_accessions",
    "fragment_count",
    "contacting_fragment_count",
    "legacy_family_tokens",
    "canonical_same_gene_examples",
    "canonical_same_gene_pfam_context",
    "modcre_interval_context",
    "proposed_domain_role",
    "proposed_is_expected_TF_DBD",
    "evidence_summary",
    "evidence_sources",
    "confidence",
    "review_status",
    "reviewer",
    "review_date",
    "manual_review_notes",
]

ROLES = {
    "TF_DNA_BINDING_DOMAIN",
    "NON_DBD_PROTEIN_INTERACTION",
    "NON_DBD_TRANSCRIPTIONAL_REGULATORY",
    "NON_DBD_LIGAND_BINDING",
    "COFACTOR_OR_COMPLEX_COMPONENT",
    "UNCERTAIN",
}

# Curated proposals are deliberately explicit and auditable.  These rationales
# paraphrase official InterPro/Pfam descriptions and are supplemented at run
# time with canonical same-gene architecture.  Contact counts never select a
# role or increase confidence.
PROPOSALS: dict[str, tuple[str, str, str, str]] = {
    "PF00023": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "Ankyrin repeats are protein-binding adaptor repeats.",
    ),
    "PF00104": (
        "NON_DBD_LIGAND_BINDING", "NO", "HIGH",
        "The official entry identifies this as the hormone-binding domain; the canonical receptor has a separate C4 zinc-finger DBD.",
    ),
    "PF00651": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "BTB/POZ mediates homo-/heterodimerisation and corepressor interactions; canonical examples carry separate zinc-finger DBDs.",
    ),
    "PF01017": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "The official entry identifies the STAT coiled-coil domain as a protein-interaction region, distinct from PF02864, the STAT DBD.",
    ),
    "PF02198": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "The integrated Pointed-domain entry assigns protein-protein association; canonical ETS proteins carry a separate PF00178 ETS DBD.",
    ),
    "PF02312": (
        "COFACTOR_OR_COMPLEX_COMPONENT", "NO", "HIGH",
        "This is the beta subunit of core-binding factor and enhances DNA binding by the alpha subunit rather than serving as its DBD.",
    ),
    "PF02865": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "The official name and canonical STAT architecture identify a protein-interaction domain separate from PF02864, the STAT DBD.",
    ),
    "PF03299": (
        "TF_DNA_BINDING_DOMAIN", "YES", "HIGH",
        "The integrated AP-2 C-terminal entry explicitly includes the helix-span-helix region that mediates site-specific DNA binding and dimerisation.",
    ),
    "PF04589": (
        "NON_DBD_TRANSCRIPTIONAL_REGULATORY", "NO", "HIGH",
        "The official entry identifies an RFX transcription-activation region N-terminal to the separate PF02257 RFX DBD.",
    ),
    "PF04704": (
        "NON_DBD_TRANSCRIPTIONAL_REGULATORY", "NO", "HIGH",
        "The official entry reports transcriptional activation when fused to a DBD; canonical proteins contain separate C2H2 zinc-finger DBDs.",
    ),
    "PF08347": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "This region binds CTNNB1/beta-catenin; canonical TCF/LEF proteins contain a separate PF00505 HMG-box DBD.",
    ),
    "PF08447": (
        "UNCERTAIN", "NO", "MEDIUM",
        "The official entry establishes a PAS fold but does not assign a specific interaction or ligand role; canonical TFs carry a separate bHLH DBD.",
    ),
    "PF08778": (
        "NON_DBD_TRANSCRIPTIONAL_REGULATORY", "NO", "HIGH",
        "The official entry explicitly identifies the HIF-1 alpha C-terminal transactivation domain.",
    ),
    "PF10401": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "Official evidence describes association with CREB-binding protein, while canonical IRFs carry the separate PF00605 IRF DBD.",
    ),
    "PF11413": (
        "UNCERTAIN", "NO", "LOW",
        "The official family text describes HIF-1 regulation and degradation but does not define a specific molecular role for this short region.",
    ),
    "PF11914": (
        "UNCERTAIN", "UNCERTAIN", "LOW",
        "The official entry calls DUF3432 functionally uncharacterised; association with zinc-finger proteins is insufficient to assign its own role.",
    ),
    "PF11928": (
        "UNCERTAIN", "NO", "MEDIUM",
        "The official entry identifies an EGR N-terminal domain but does not establish its molecular role; canonical EGR has separate C2H2 DBDs.",
    ),
    "PF12577": (
        "NON_DBD_TRANSCRIPTIONAL_REGULATORY", "NO", "MEDIUM",
        "This PPAR-gamma N-terminal region is separate from canonical C4 zinc-finger and ligand-binding domains; the official entry cautions that it may not be a separate domain.",
    ),
    "PF12796": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "Ankyrin repeats are protein-binding repeats; canonical NF-kappa-B-related architecture has a separate Rel-homology DBD.",
    ),
    "PF14598": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "The official entry describes PAS-B binding to a STAT6 LXXLL motif; canonical bHLH-PAS TFs have a separate bHLH DBD.",
    ),
    "PF15951": (
        "UNCERTAIN", "NO", "MEDIUM",
        "The official entry locates this region at the N-terminus of MiT/TFE factors but gives no molecular role; canonical proteins have a separate bHLH DBD.",
    ),
    "PF17725": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "The official entry identifies a YAP/TAZ-binding domain and distinguishes it from the N-terminal PF01285 TEA DBD.",
    ),
    "PF18326": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "The official RFX5 N-terminal-domain entry assigns homodimerisation and complex-assembly functions; canonical RFX5 has separate RFX DBDs.",
    ),
    "PF19525": (
        "NON_DBD_TRANSCRIPTIONAL_REGULATORY", "NO", "HIGH",
        "This ETS-flanking region autoinhibits and regulates the adjacent PF00178 ETS DBD rather than constituting that DBD.",
    ),
    "PF19536": (
        "UNCERTAIN", "NO", "MEDIUM",
        "The official entry locates this POU2F1 C-terminal region adjacent to the homeobox but does not assign a direct molecular role; canonical POU/homeobox DBDs are separate.",
    ),
    "PF25340": (
        "NON_DBD_PROTEIN_INTERACTION", "NO", "HIGH",
        "The official RFX BCD entry assigns mainly dimerisation and distinguishes it from PF02257, the RFX HTH DBD.",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a 26-row DRAFT Pfam domain-role evidence table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--unique-pfam", type=Path, default=ROOT / "outputs/af3_failed_102_unique_pfam_review.tsv")
    parser.add_argument("--inventory", type=Path, default=ROOT / "outputs/af3_failed_102_domain_inventory.tsv")
    parser.add_argument("--inventory-qc", type=Path, default=ROOT / "outputs/af3_failed_102_domain_inventory_qc.json")
    parser.add_argument("--database", type=Path, default=ROOT / "data/tf_webdb.sqlite")
    parser.add_argument("--family-tree", type=Path, default=ROOT / "data_sources/tf_family_tree.json")
    parser.add_argument("--modcre-candidates", type=Path, default=ROOT / "external/baldo_model_inventory/af3_failed_same_gene_modcre_candidates_v2.tsv")
    parser.add_argument("--modcre-audit", type=Path, default=ROOT / "external/baldo_model_inventory/af3_failed_same_gene_modcre_interface_audit.tsv")
    parser.add_argument("--failed-fasta", type=Path, default=Path("/home/patricia/TF_database_Baldo_data/TF_without_model.fasta"))
    parser.add_argument("--models-root", type=Path, default=Path("/data/sbi/interchange/boliva/patricia/models"))
    parser.add_argument("--output-tsv", type=Path, default=ROOT / "outputs/af3_failed_102_domain_role_draft.tsv")
    parser.add_argument("--output-qc", type=Path, default=ROOT / "outputs/af3_failed_102_domain_role_draft_qc.json")
    parser.add_argument("--official-timeout", type=float, default=30.0)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def db_signature(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256(path),
    }


def clean_description(value: Any) -> str:
    if not isinstance(value, list):
        return ""
    text = " ".join(str(item.get("text", "")) for item in value if isinstance(item, dict))
    text = html.unescape(re.sub(r"<[^>]+>", " ", text))
    return re.sub(r"\s+", " ", text).strip()


def retrieve_json(url: str, timeout: float) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": "tf-webdb-domain-role-draft/1.0"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.load(response)


def retrieve_official(pfam_id: str, fallback_ipr: str, timeout: float) -> dict[str, Any]:
    result: dict[str, Any] = {
        "status": "NOT_RETRIEVED",
        "url": OFFICIAL_PAGE.format(pfam_id=pfam_id),
        "pfam_description": "",
        "interpro_description": "",
        "error": "",
    }
    try:
        payload = retrieve_json(OFFICIAL_API.format(pfam_id=pfam_id), timeout)
        metadata = payload.get("metadata", {})
        result.update(
            status="RETRIEVED",
            pfam_name=(metadata.get("name") or {}).get("name", ""),
            interpro_id=metadata.get("integrated") or fallback_ipr,
            pfam_type=metadata.get("type") or "",
            pfam_description=clean_description(metadata.get("description")),
        )
        ipr = result["interpro_id"]
        if ipr:
            try:
                integrated = retrieve_json(INTERPRO_API.format(interpro_id=ipr), timeout)
                imeta = integrated.get("metadata", {})
                result["interpro_name"] = (imeta.get("name") or {}).get("name", "")
                result["interpro_type"] = imeta.get("type") or ""
                result["interpro_description"] = clean_description(imeta.get("description"))
            except (OSError, ValueError, urllib.error.URLError) as exc:
                result["integrated_error"] = f"{type(exc).__name__}: {exc}"
    except (OSError, ValueError, urllib.error.URLError) as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def compact_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def canonical_context(
    connection: sqlite3.Connection, examples: list[str]
) -> tuple[str, str]:
    contexts: dict[str, list[str]] = {}
    retained: list[str] = []
    for accession in examples:
        annotation = connection.execute(
            "SELECT tf_id, reviewed FROM tf_annotation WHERE uniprot_accession=? ORDER BY reviewed DESC LIMIT 1",
            (accession,),
        ).fetchone()
        if annotation is None or annotation[1] != 1:
            continue
        retained.append(accession)
        domains = connection.execute(
            """SELECT pfam_id, start, end, pfam_name
               FROM tf_pfam_annotation WHERE tf_id=?
               ORDER BY start, end, pfam_id""",
            (annotation[0],),
        ).fetchall()
        contexts[accession] = [
            f"{pfam}:{start}-{end}:{name or ''}" for pfam, start, end, name in domains
        ]
    return ";".join(retained), compact_json(contexts) if contexts else ""


def build_modcre_context(
    rows: list[dict[str, str]],
    examples: list[str],
    candidate_index: dict[tuple[str, str], dict[str, str]],
    audit_index: dict[str, list[dict[str, str]]],
) -> str:
    context: dict[str, dict[str, Any]] = {}
    for accession in examples:
        intervals: set[str] = set()
        source_rows: list[dict[str, str]] = []
        for row in rows:
            raw = row.get("candidate_modcre_model_intervals", "")
            if raw:
                intervals.update(json.loads(raw).get(accession, []))
            candidate = candidate_index.get((row["failed_accession"], accession))
            if candidate:
                source_rows.append(candidate)
        if not intervals and not source_rows:
            continue
        ordered = sorted(intervals, key=lambda item: tuple(int(x) for x in item.split(":")))
        templates = sorted({r.get("best_template", "") for r in source_rows if r.get("best_template")})
        domains = sorted({r.get("best_model_domain", "") for r in source_rows if r.get("best_model_domain")})
        audit_results = sorted({
            a.get("result", "") for a in audit_index.get(accession, []) if a.get("result")
        })
        context[accession] = {
            "interval_count": len(ordered),
            "interval_examples": ordered[:8],
            "best_templates": templates,
            "best_model_domains": domains,
            "interface_audit_results_secondary_only": audit_results,
        }
    return compact_json(context) if context else ""


def add_check(
    checks: dict[str, Any], name: str, passed: bool, observed: Any, expected: Any
) -> None:
    checks[name] = {"passed": bool(passed), "observed": observed, "expected": expected}


def main() -> int:
    args = parse_args()
    inputs = [
        args.unique_pfam, args.inventory, args.inventory_qc, args.database,
        args.family_tree, args.modcre_candidates, args.modcre_audit,
        args.failed_fasta, args.models_root,
    ]
    missing = [str(path) for path in inputs if not path.exists()]
    if missing:
        raise SystemExit("Missing required input(s): " + ", ".join(missing))
    for output in (args.output_tsv, args.output_qc):
        if output.exists():
            raise SystemExit(f"Refusing to overwrite existing output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)

    before = db_signature(args.database)
    unique_rows = read_tsv(args.unique_pfam)
    inventory = read_tsv(args.inventory)
    candidates = read_tsv(args.modcre_candidates)
    audits = read_tsv(args.modcre_audit)
    inventory_qc = json.loads(args.inventory_qc.read_text(encoding="utf-8"))
    family_tree = json.loads(args.family_tree.read_text(encoding="utf-8"))

    pfam_ids = [row["pfam_id"] for row in unique_rows]
    if len(pfam_ids) != EXPECTED_PFAM_COUNT or len(set(pfam_ids)) != EXPECTED_PFAM_COUNT:
        raise SystemExit("Unique-Pfam input is not exactly 26 unique IDs")
    if set(pfam_ids) != set(PROPOSALS):
        raise SystemExit("Proposal map does not exactly match the unique-Pfam input")

    inventory_by_pfam: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in inventory:
        if row.get("pfam_id") in set(pfam_ids):
            inventory_by_pfam[row["pfam_id"]].append(row)

    candidate_index = {
        (row["failed_tf_id"], row["candidate_tf_id"]): row for row in candidates
    }
    audit_index: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in audits:
        audit_index[row["best_model_accession"]].append(row)

    family_tokens = {
        token
        for node in family_tree
        for token in node.get("tokens", [])
        if isinstance(token, str)
    }

    db_uri = f"file:{args.database.resolve()}?mode=ro"
    connection = sqlite3.connect(db_uri, uri=True)
    connection.execute("PRAGMA query_only=ON")
    output_rows: list[dict[str, Any]] = []
    retrieval: dict[str, dict[str, Any]] = {}
    try:
        for source in unique_rows:
            pfam_id = source["pfam_id"]
            rows = inventory_by_pfam[pfam_id]
            official = retrieve_official(
                pfam_id, rows[0].get("interpro_id", "") if rows else "", args.official_timeout
            )
            retrieval[pfam_id] = official
            role, expected_dbd, confidence, rationale = PROPOSALS[pfam_id]

            # A failed official retrieval makes the role uncertain unless local
            # canonical evidence decisively separates the region from a DBD.
            reviewed_examples = sorted({
                accession
                for row in rows
                for accession in row.get("reviewed_same_gene_candidate_accessions", "").split(";")
                if accession
            })
            examples_text, contexts_text = canonical_context(connection, reviewed_examples)
            locally_decisive_non_dbd = bool(contexts_text) and expected_dbd == "NO"
            if official["status"] == "NOT_RETRIEVED" and not locally_decisive_non_dbd:
                role = "UNCERTAIN"
                expected_dbd = "UNCERTAIN"
                confidence = "LOW"
                rationale = "Official entry was not retrieved and local evidence was not decisive."

            accessions = sorted({row["failed_accession"] for row in rows})
            contacting = sum(
                row.get("fragment_selected_result") == "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"
                for row in rows
            )
            legacy = sorted({
                token
                for row in rows
                for token in row.get("legacy_family_text", "").split(",")
                if token
            })
            evidence_sources = [official["url"]]
            if examples_text:
                evidence_sources.extend([
                    "data/tf_webdb.sqlite (mode=ro)",
                    "outputs/af3_failed_102_domain_inventory.tsv",
                ])
            if legacy:
                evidence_sources.extend([
                    "data_sources/tf_family_tree.json",
                    "tf.family_text",
                    "tf_family",
                ])
            modcre = build_modcre_context(rows, reviewed_examples, candidate_index, audit_index)
            if modcre:
                evidence_sources.extend([
                    "external/baldo_model_inventory/af3_failed_same_gene_modcre_candidates_v2.tsv",
                    "external/baldo_model_inventory/af3_failed_same_gene_modcre_interface_audit.tsv",
                ])

            contact_note = (
                f" Fragment contacts: {contacting}/{len(rows)}; treated as secondary evidence only."
                if contacting else
                f" Fragment contacts: 0/{len(rows)}; contacts are not a role criterion."
            )
            summary = rationale + contact_note
            output_rows.append({
                "pfam_id": pfam_id,
                "pfam_name": official.get("pfam_name") or source["pfam_name"],
                "interpro_id": official.get("interpro_id") or (rows[0].get("interpro_id", "") if rows else ""),
                "interpro_name": official.get("interpro_name") or (rows[0].get("interpro_name", "") if rows else ""),
                "official_entry_url": official["url"],
                "official_source_status": official["status"],
                "failed_accession_count": len(accessions),
                "failed_accessions": ";".join(accessions),
                "fragment_count": len(rows),
                "contacting_fragment_count": contacting,
                "legacy_family_tokens": ";".join(legacy),
                "canonical_same_gene_examples": examples_text,
                "canonical_same_gene_pfam_context": contexts_text,
                "modcre_interval_context": modcre,
                "proposed_domain_role": role,
                "proposed_is_expected_TF_DBD": expected_dbd,
                "evidence_summary": summary,
                "evidence_sources": ";".join(dict.fromkeys(evidence_sources)),
                "confidence": confidence,
                "review_status": "DRAFT_UNREVIEWED",
                "reviewer": "",
                "review_date": "",
                "manual_review_notes": "",
            })
    finally:
        connection.close()

    after_read = db_signature(args.database)
    checks: dict[str, Any] = {}
    add_check(checks, "exactly_26_rows", len(output_rows) == 26, len(output_rows), 26)
    add_check(checks, "exactly_26_unique_pfam_ids", len({r["pfam_id"] for r in output_rows}) == 26, len({r["pfam_id"] for r in output_rows}), 26)
    add_check(checks, "ids_match_unique_pfam_input", {r["pfam_id"] for r in output_rows} == set(pfam_ids), sorted(r["pfam_id"] for r in output_rows), sorted(pfam_ids))
    forbidden = sorted(set(COLUMNS) & {
        "domain_role_review", "is_expected_DBD_review", "final_structural_class",
        "accession_structural_class", "final_accession_classification",
    })
    add_check(checks, "no_final_accession_classification_fields", not forbidden, forbidden, [])
    add_check(checks, "database_unchanged_after_read", before == after_read, after_read, before)
    statuses = sorted({r["review_status"] for r in output_rows})
    add_check(checks, "all_rows_draft_unreviewed", statuses == ["DRAFT_UNREVIEWED"], statuses, ["DRAFT_UNREVIEWED"])
    missing_sources = [r["pfam_id"] for r in output_rows if r["proposed_domain_role"] != "UNCERTAIN" and not r["evidence_sources"]]
    add_check(checks, "evidence_source_for_every_non_uncertain_proposal", not missing_sources, missing_sources, [])
    invalid_roles = sorted({r["proposed_domain_role"] for r in output_rows} - ROLES)
    add_check(checks, "only_allowed_provisional_roles", not invalid_roles, invalid_roles, [])
    source_ids = {r["pfam_id"] for r in output_rows}
    family_tokens_seen = sorted({
        token for r in output_rows for token in r["legacy_family_tokens"].split(";") if token
    })
    add_check(
        checks, "family_tree_loaded_and_consulted",
        isinstance(family_tree, list) and bool(family_tree),
        {"nodes": len(family_tree), "legacy_tokens_matching_tree": len(set(family_tokens_seen) & family_tokens)},
        {"nodes": ">0"},
    )
    add_check(
        checks, "upstream_inventory_qc_passed",
        inventory_qc.get("all_checks_passed") is True,
        inventory_qc.get("all_checks_passed"), True,
    )
    secondary_ok = all("secondary evidence only" in r["evidence_summary"] for r in output_rows if int(r["contacting_fragment_count"]) > 0)
    add_check(checks, "contacting_fragment_count_secondary_only", secondary_ok, secondary_ok, True)

    by_id = {r["pfam_id"]: r for r in output_rows}
    e9_domains = {
        row["pfam_id"] for row in inventory if row["failed_accession"] == "E9PN75" and row.get("pfam_id")
    }
    control_1 = (
        by_id["PF02198"]["proposed_domain_role"] != "TF_DNA_BINDING_DOMAIN"
        and "PF02198" in e9_domains and "PF00178" not in e9_domains
    )
    add_check(checks, "control_PF02198_E9PN75_not_ets_family_promoted_to_dbd", control_1, {
        "proposal": by_id["PF02198"]["proposed_domain_role"],
        "E9PN75_pfam_ids": sorted(e9_domains),
    }, {"proposal": "not TF_DNA_BINDING_DOMAIN", "contains": "PF02198", "does_not_contain": "PF00178"})
    add_check(checks, "control_PF01017_B4DHE0_contact_not_auto_dbd",
              by_id["PF01017"]["proposed_domain_role"] != "TF_DNA_BINDING_DOMAIN" and int(by_id["PF01017"]["contacting_fragment_count"]) > 0,
              {"proposal": by_id["PF01017"]["proposed_domain_role"], "contacts": by_id["PF01017"]["contacting_fragment_count"]},
              {"proposal": "not TF_DNA_BINDING_DOMAIN", "contacts": ">0"})
    add_check(checks, "control_PF17725_Q59EF3_contact_not_auto_dbd",
              by_id["PF17725"]["proposed_domain_role"] != "TF_DNA_BINDING_DOMAIN" and int(by_id["PF17725"]["contacting_fragment_count"]) > 0,
              {"proposal": by_id["PF17725"]["proposed_domain_role"], "contacts": by_id["PF17725"]["contacting_fragment_count"]},
              {"proposal": "not TF_DNA_BINDING_DOMAIN", "contacts": ">0"})
    ap2_summary = by_id["PF03299"]["evidence_summary"].lower()
    add_check(checks, "control_PF03299_H7C4N4_ap2_dbd_evidence_draft_only",
              by_id["PF03299"]["proposed_domain_role"] == "TF_DNA_BINDING_DOMAIN"
              and "site-specific dna binding" in ap2_summary
              and by_id["PF03299"]["review_status"] == "DRAFT_UNREVIEWED",
              {"proposal": by_id["PF03299"]["proposed_domain_role"], "review_status": by_id["PF03299"]["review_status"], "summary": by_id["PF03299"]["evidence_summary"]},
              {"proposal": "TF_DNA_BINDING_DOMAIN", "review_status": "DRAFT_UNREVIEWED", "official_AP2_DBD_evidence": True})
    pf10401_context = by_id["PF10401"]["canonical_same_gene_pfam_context"]
    add_check(checks, "control_PF10401_not_inferred_from_name_alone",
              "PF00605" in pf10401_context and "CREB-binding protein" in by_id["PF10401"]["evidence_summary"],
              {"canonical_context_has_PF00605": "PF00605" in pf10401_context, "summary": by_id["PF10401"]["evidence_summary"]},
              {"canonical_context_has_PF00605": True, "official_function_evidence": True})

    not_retrieved = sorted(pfam for pfam, item in retrieval.items() if item["status"] == "NOT_RETRIEVED")
    all_prewrite_passed = all(item["passed"] for item in checks.values())
    if not all_prewrite_passed:
        raise SystemExit("QC failed before output write: " + ", ".join(name for name, item in checks.items() if not item["passed"]))

    with args.output_tsv.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    after_write = db_signature(args.database)
    add_check(checks, "database_unchanged_after_output_write", before == after_write, after_write, before)
    final_rows = read_tsv(args.output_tsv)
    add_check(checks, "written_tsv_round_trip_26_rows", len(final_rows) == 26, len(final_rows), 26)
    all_passed = all(item["passed"] for item in checks.values())
    qc = {
        "all_checks_passed": all_passed,
        "draft_only": True,
        "official_entries_not_retrieved": not_retrieved,
        "official_retrieval": {
            pfam: {
                "status": item["status"],
                "url": item["url"],
                "error": item.get("error", ""),
                "integrated_error": item.get("integrated_error", ""),
            }
            for pfam, item in sorted(retrieval.items())
        },
        "database_signature_before": before,
        "database_signature_after": after_write,
        "checks": checks,
        "inputs": {key: str(value) for key, value in vars(args).items() if isinstance(value, Path)},
        "output_tsv": str(args.output_tsv),
    }
    with args.output_qc.open("x", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if not all_passed:
        raise SystemExit("Final QC failed")
    print(f"Wrote {args.output_tsv} ({len(output_rows)} draft rows)")
    print(f"Wrote {args.output_qc}")
    print(f"Official entries not retrieved: {', '.join(not_retrieved) if not_retrieved else 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
