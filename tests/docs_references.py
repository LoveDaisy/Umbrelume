"""Lightweight, offline readers for the repository's documented evidence links."""
import json
import re
from urllib.parse import unquote, urlsplit


FIELD_REFERENCE = re.compile(r"`([a-z_0-9]+\.json:[^`]+)`")


def resolve_product_fields(text, results_root):
    """Resolve every cited basename:dot.path, rejecting missing/ambiguous files.

    Brackets are path separators; '*' expands every branch of a nonempty dict.
    Return the references so callers can reject accidentally empty coverage.
    """
    references = FIELD_REFERENCE.findall(text)
    for reference in references:
        file_name, _, pointer = reference.partition(":")
        hits = sorted(results_root.rglob(file_name))
        assert len(hits) == 1, f"{reference}: expected exactly one {file_name}, found {hits}"
        normalized = re.sub(r"\[([^\]]+)\]", r".\1", pointer)
        keys = [key for key in normalized.split(".") if key]
        assert keys, reference
        branches = [json.loads(hits[0].read_text())]
        for key in keys:
            expanded = []
            for branch in branches:
                if key == "*":
                    assert isinstance(branch, dict) and branch, reference
                    expanded.extend(branch.values())
                else:
                    try:
                        expanded.append(branch[int(key)] if key.isdigit() else branch[key])
                    except (KeyError, IndexError, TypeError) as error:
                        raise AssertionError(f"{reference}: missing field {key}") from error
            branches = expanded
        assert branches, f"{reference}: resolved to nothing"
    return references


def local_link_targets(text):
    """Read inline/image and reference links, not examples in fenced code.

    This is the small Markdown subset used by MAINLINE, not a Markdown linter.
    External URLs and page-local anchors do not denote local files.
    """
    text = re.sub(r"^```[^\n]*\n.*?^```[^\n]*$", "", text, flags=re.M | re.S)
    definition = re.compile(r"^\s{0,3}\[([^\]\n]+)\]:\s*(<[^>\n]+>|\S+).*$", re.M)

    def label(value):
        return " ".join(value.casefold().split())

    definitions = {label(key): target.strip("<>")
                   for key, target in definition.findall(text)}
    text = definition.sub("", text)
    targets = [target.strip("<>") for target in
               re.findall(r"!?\[[^\]\n]*\]\(\s*(<[^>\n]+>|[^\s)]+)", text)]
    reference = re.compile(r"!?\[([^\]\n]+)\]\[([^\]\n]*)\]")
    for title, key in reference.findall(text):
        key = label(key or title)
        assert key in definitions, f"undefined reference link: {key}"
        targets.append(definitions[key])
    # Shortcut references are links only when their definition exists.
    text = reference.sub("", text)
    for key in re.findall(r"!?\[([^\]\n]+)\](?![\[(])", text):
        if label(key) in definitions:
            targets.append(definitions[label(key)])
    local = []
    for target in targets:
        parsed = urlsplit(target)
        if not parsed.scheme and not parsed.netloc and parsed.path:
            local.append(unquote(parsed.path))
    return local


def check_local_links(text, document_dir):
    targets = local_link_targets(text)
    for target in targets:
        assert (document_dir / target).exists(), f"missing local link: {target}"
    return targets
