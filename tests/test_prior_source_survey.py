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
