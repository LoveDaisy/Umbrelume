"""Schema gate for the prior-source survey products (task-prior-source-survey).

Light-weight contract checks only: the pinned-source SHA format, the
unavailable-row search-trail requirement, the coverage-matrix verdict
vocabulary, and the dimension-value alignment against the committed
inversion product. These are data-contract gates, not physics tests; the
survey changed no default and this file must stay removable without
touching the data products.
"""

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SURVEY = ROOT / "validation" / "results" / "prior_source_survey"
INVERSION = ROOT / "validation" / "results" / "atmosphere_inversion" / "inversion.json"

SHA_RE = re.compile(r"^[0-9a-f]{64}$")
VERDICTS = {"covered", "partial", "unavailable"}


def _load(path):
    with open(path) as fh:
        return json.load(fh)


def test_pinned_sources_carry_wellformed_sha_and_location():
    sources = _load(SURVEY / "sources.json")
    assert sources["scope"] and sources["source_convention"], "header contract fields must be present"
    checked = 0
    for entry in sources["sources"]:
        if entry["status"] in ("available", "partial"):
            if entry.get("source") is not None:
                assert SHA_RE.match(entry["source"]["sha256"]), entry["name"]
                assert entry["source"]["location"], entry["name"]
                assert entry["source"]["bytes"] > 0, entry["name"]
                checked += 1
            for f in entry.get("files") or []:
                assert SHA_RE.match(f["sha256"]), f"{entry['name']}/{f['file_name']}"
                assert f["location"], f"{entry['name']}/{f['file_name']}"
                assert f["bytes"] > 0, f"{entry['name']}/{f['file_name']}"
                checked += 1
    assert checked >= 20, "the survey must pin at least the OPAC paper + component files + site export + climatology papers"
    # every numeric value in the humidity growth table carries a pinned source
    growth = _load(SURVEY / "humidity_growth.json")
    for comp in growth["components"]:
        for rh, row in comp["per_rh"].items():
            assert SHA_RE.match(row["source"]["sha256"]), f"{comp['component']}/{rh}"
            assert row["source"]["location"], f"{comp['component']}/{rh}"


def test_unavailable_rows_carry_search_trail_pointers():
    sources = _load(SURVEY / "sources.json")
    for entry in sources["sources"]:
        if entry["status"] == "unavailable":
            notes = entry.get("search_notes") or ""
            assert "search_log.md" in notes, entry["name"]
    matrix = _load(SURVEY / "coverage_matrix.json")
    for dim in matrix["matrix"]:
        for cell in dim["cells"]:
            if cell["verdict"] == "unavailable":
                blob = (cell.get("basis") or "") + (cell.get("location") or "")
                assert "search_log.md" in blob or "W1." in blob or "W2." in blob or "W3." in blob, (
                    f"{dim['json_path']}/{cell['source']}"
                )


def test_coverage_matrix_verdict_vocabulary_and_dimension_alignment():
    matrix = _load(SURVEY / "coverage_matrix.json")
    inversion = _load(INVERSION)
    # verdict vocabulary gate: no fourth verdict, no justification language
    for dim in matrix["matrix"]:
        for cell in dim["cells"]:
            assert cell["verdict"] in VERDICTS, f"{dim['json_path']}/{cell['source']}"
            blob = (cell.get("basis") or "") + json.dumps(dim.get("value", ""))
            assert "justified" not in blob.lower(), f"justification language in {dim['json_path']}"
    # dimension values must equal the committed inversion product (json-path resolved)
    for dim in matrix["matrix"]:
        path = dim["json_path"]
        if not path.startswith("retrieval_spec."):
            continue
        node = inversion
        for part in path.split("."):
            node = node[part]
        assert dim["value"] == node, path
    # branch-rule status mirrors the committed flag
    assert matrix["branch_rule_status"]["value"] == inversion["retrieval_spec"]["branch_rule_fired"]


def test_survey_products_present_and_default_untouched_claims_consistent():
    for name in ("sources.json", "search_log.md", "humidity_growth.json", "coverage_matrix.json", "three_paths.json"):
        assert (SURVEY / name).exists(), name
    lanes = _load(SURVEY / "three_paths.json")
    assert {lane["lane"] for lane in lanes["lanes"]} == {"fixed_higher_prior", "documented_range", "keep_current"}
    for lane in lanes["lanes"]:
        assert lane["existing_evidence"] and lane["gaps"] and lane["what_this_lane_still_needs"], lane["lane"]
        # no recommendation column: lane objects must not carry verdict/recommendation keys
        assert not any(k.startswith("recommend") for k in lane), lane["lane"]


def test_coverage_matrix_complete_over_frozen_sources():
    matrix = _load(SURVEY / "coverage_matrix.json")
    sources = _load(SURVEY / "sources.json")
    frozen = {e["name"] for e in sources["sources"]}
    dims = matrix["matrix"]
    assert len(dims) == 6, "five retrieval-spec dimensions plus the humidity-variant dimension"
    paths = [d["json_path"] for d in dims]
    assert sum(p.startswith("retrieval_spec.") for p in paths) == 5
    assert sum(p.startswith("humidity_variants") for p in paths) == 1
    for d in dims:
        got = [c["source"] for c in d["cells"]]
        assert len(got) == len(set(got)), f"duplicate source cell in {d['json_path']}"
        assert set(got) == frozen, f"incomplete source set in {d['json_path']}: missing {frozen - set(got)}"


def test_worked_derived_examples_recompute():
    """The delivered conversions and geometric derivations must recompute from
    their own stated formula and inputs (numeric contract, not just SHA format)."""
    import math

    matrix = _load(SURVEY / "coverage_matrix.json")
    examples = {}
    for entry in matrix["derived_comparisons"]:
        for w in entry.get("worked_examples", []):
            examples[w["id"]] = w
    assert {"lanai_low", "lanai_high", "gsfc_median", "h63_target_layer"} <= set(examples)

    for eid in ("lanai_low", "lanai_high", "gsfc_median"):
        ex = examples[eid]
        tau = ex["inputs"]["tau"]
        wl = ex["inputs"]["wavelength_nm"]
        a0, a1 = ex["inputs"]["alpha_range"]
        recomputed = sorted(tau * (550.0 / wl) ** (-a) for a in (a0, a1))
        stored = ex["outputs"]["tau550_range"]
        for got, want in zip(recomputed, stored):
            assert abs(got - want) <= 5e-4 * max(1.0, want), f"{eid}: {got} vs {want}"

    ex = examples["h63_target_layer"]
    H, top = ex["inputs"]["H_m"], ex["inputs"]["top_m"]
    h63 = -H * math.log(1 - 0.63 * (1 - math.exp(-top / H)))
    assert abs(ex["outputs"]["h63_m"] - h63) <= 1.0, f"h63: {ex['outputs']['h63_m']} vs {h63}"
