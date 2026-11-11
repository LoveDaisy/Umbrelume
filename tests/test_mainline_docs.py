"""MAINLINE evidence addresses: text/JSON only, no observational solver imports."""
import json
from pathlib import Path

import pytest

from docs_references import check_local_links, resolve_product_fields

ROOT = Path(__file__).resolve().parents[1]
MAINLINE = ROOT / "docs/MAINLINE.zh.md"


def test_mainline_product_fields():
    assert resolve_product_fields(MAINLINE.read_text(), ROOT / "validation/results")


def test_mainline_local_links():
    assert check_local_links(MAINLINE.read_text(), MAINLINE.parent)


def test_field_resolver_nested_arrays_and_all_wildcard_branches(tmp_path):
    (tmp_path / "sample.json").write_text(json.dumps({
        "events": {"a": {"rows": [{"value": 1}]}, "b": {"rows": [{"value": 2}]}}
    }))
    text = "`sample.json:events.*.rows[0].value`"
    assert resolve_product_fields(text, tmp_path) == [text.strip("`")]
    # A good first branch must not hide a bad later branch.
    with pytest.raises(AssertionError, match="missing field"):
        resolve_product_fields("`sample.json:events.*.rows[1].value`", tmp_path)
    (tmp_path / "sample.json").write_text(json.dumps({
        "events": {"a": {"value": 1}, "b": {}}
    }))
    with pytest.raises(AssertionError, match="missing field value"):
        resolve_product_fields("`sample.json:events.*.value`", tmp_path)


@pytest.mark.parametrize("content,pointer", [
    ({"rows": []}, "rows[0]"), ({"rows": {}}, "rows.*"),
    ({"rows": 1}, "rows.value"), ({"rows": {}}, "rows.no_such_field"),
])
def test_field_resolver_rejects_bad_fields(tmp_path, content, pointer):
    (tmp_path / "sample.json").write_text(json.dumps(content))
    with pytest.raises(AssertionError):
        resolve_product_fields(f"`sample.json:{pointer}`", tmp_path)


def test_field_resolver_rejects_missing_and_ambiguous_files(tmp_path):
    reference = "`sample.json:value`"
    with pytest.raises(AssertionError, match="exactly one"):
        resolve_product_fields(reference, tmp_path)
    for folder in ("first", "second"):
        directory = tmp_path / folder
        directory.mkdir()
        (directory / "sample.json").write_text('{"value": 1}')
    with pytest.raises(AssertionError, match="exactly one"):
        resolve_product_fields(reference, tmp_path)


def test_links_inline_images_references_and_url_components(tmp_path):
    (tmp_path / "with space.md").write_text("target")
    text = '''[inline](with%20space.md?view=raw#section)
![image](<with space.md>)
[full][target]
[target][]
[target]
[target]: <with space.md> "A title"
[page](#section) [web](https://example.invalid/missing.md)
[mail](mailto:someone@example.invalid)
```bash
[example](missing.md)
```
'''
    assert check_local_links(text, tmp_path) == ["with space.md"] * 5


@pytest.mark.parametrize("text", [
    "[bad](missing.md)", "![bad](missing.png)",
    "[bad][ref]\n[ref]: missing.md", "[bad][undefined]",
])
def test_links_reject_missing_files_and_undefined_references(tmp_path, text):
    with pytest.raises(AssertionError):
        check_local_links(text, tmp_path)
