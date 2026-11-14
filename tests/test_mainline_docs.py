"""MAINLINE evidence addresses: text/JSON only, no observational solver imports."""
import json
from pathlib import Path
import re

import pytest

from docs_references import check_local_links, resolve_product_fields

ROOT = Path(__file__).resolve().parents[1]
MAINLINE = ROOT / "docs/MAINLINE.zh.md"


def test_mainline_product_fields():
    assert resolve_product_fields(MAINLINE.read_text(), ROOT / "validation/results")


def test_mainline_local_links():
    assert check_local_links(MAINLINE.read_text(), MAINLINE.parent)


EXPECTED_IMAGES = {
    "../figures/mainline/l1_molecular.png",
    "../figures/mainline/l2_free_troposphere.png",
    "../figures/mainline/l25_stratosphere.png",
    "../figures/mainline/transmission_ladder.png",
    "../figures/mainline/shadow_profile.png",
    "../figures/multiple_events/comparison.png",
    "../figures/filter_response_2011/comparison.png",
    "../figures/dec2011/comparison.png",
}


def check_mainline_images(text):
    # Reuse the local-link resolver (including missing-image failures), but
    # require actual embeds too: ordinary links or zero images cannot pass.
    check_local_links(text, MAINLINE.parent)
    targets = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text)
    assert set(targets) == EXPECTED_IMAGES
    assert len(targets) == len(EXPECTED_IMAGES)
    assert re.findall(r"\*\*图 (\d+) ·", text) == [str(i) for i in range(1, 9)]


def test_mainline_embeds_expected_figures_and_numbered_captions():
    check_mainline_images(MAINLINE.read_text())
    assert "equal-Y 色条已逐半径移除亮度" in MAINLINE.read_text()
    assert "历史 v0.2" in MAINLINE.read_text()
    assert "假设响应包络" in MAINLINE.read_text()
    assert "不是置信区间" in MAINLINE.read_text()


def test_mainline_images_missing_embed_and_bad_path_fail():
    text = MAINLINE.read_text()
    for broken in (text.replace("![", "["), text.replace(
            "../figures/mainline/l1_molecular.png", "../figures/mainline/missing.png")):
        with pytest.raises(AssertionError):
            check_mainline_images(broken)


def test_mainline_figure_readme_links():
    readme = ROOT / "figures/mainline/README.md"
    assert check_local_links(readme.read_text(), readme.parent)


BILINGUAL_PAIRS = [
    ((ROOT / "README.md"), (ROOT / "README.zh.md")),
    ((ROOT / "docs/ROADMAP.zh.md"), (ROOT / "docs/ROADMAP.en.md")),
    ((ROOT / "docs/EXTENDING.zh.md"), (ROOT / "docs/EXTENDING.en.md")),
    ((ROOT / "docs/MAINLINE.zh.md"), (ROOT / "docs/MAINLINE.en.md")),
]


def test_bilingual_pairs_exist_and_cross_link():
    """双语集合闸：成员存在、与孪生实链互指、全部本地链接可解析。

    文件名后缀只标语言，权威版本逐对声明（docs/README.md）；内容同步靠定期
    人工检查，此闸只拦截「忘了译本存在/互链断裂」这类失效。
    """
    for a, b in BILINGUAL_PAIRS:
        for source, twin in ((a, b), (b, a)):
            assert source.exists(), f"missing bilingual file: {source}"
            text = source.read_text()
            assert f"]({twin.name})" in text, \
                f"{source.name} does not link its twin {twin.name}"
            assert check_local_links(text, source.parent)


@pytest.mark.parametrize("path", sorted(
    p for p in (ROOT / "docs").rglob("*.md") if p.name != "MAINLINE.zh.md"
))
def test_docs_tree_local_links_resolve(path):
    """docs 树（含 records/、physics/）全部本地链接可解析。

    MAINLINE.zh.md 由更强的图片/字段闸单独覆盖；根目录双语 README 由
    pairing 闸覆盖。figures/ 与 validation/ 的 README 索引暂不在本闸射程。
    """
    # 无本地链接的文件返回空列表也是合法的——这里只消费闸的抛错路径。
    check_local_links(path.read_text(), path.parent)


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
