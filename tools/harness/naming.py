"""C1: filename / folder parsing into setting labels, plus a template tool for writing new patterns.

Patterns live in config.json under "naming". Every pattern is searched (case-insensitive) against
the file's path relative to the library root. Fields from earlier patterns win, so a folder pattern
(e.g. "Short RT/") and a filename pattern (e.g. "RT1.2_LF0.9") combine.
A pattern may carry "scale": {field: factor} for implied-decimal fields (0.1 turns "26" into 2.6).
"""
import re
from collections import Counter

FIELDS = ("program", "rt", "lf", "hf", "delay", "amp", "lp", "span", "early", "late")
NUMERIC = ("rt", "lf", "hf", "delay", "amp", "lp", "span", "early", "late")


def parse_name(relpath, patterns):
    out = {f"label_{k}": None for k in FIELDS}
    matched = []
    for p in patterns:
        m = re.search(p["regex"], relpath, flags=re.IGNORECASE)
        if not m:
            continue
        matched.append(p["name"])
        for k, v in p.get("set", {}).items():  # constants implied by the pattern itself
            if out.get(f"label_{k}") is None:
                out[f"label_{k}"] = v
        for k, v in m.groupdict().items():
            if v is None or out.get(f"label_{k}") is not None:
                continue
            v = p.get("map", {}).get(k, {}).get(v.lower(), v)
            if k in NUMERIC:
                try:
                    v = float(v) * p.get("scale", {}).get(k, 1)  # e.g. 0.1 for implied-decimal fields
                except ValueError:
                    continue
                v = round(v, 6)
            elif k == "program":
                v = v.upper()
            out[f"label_{k}"] = v
    out["label_patterns"] = "+".join(matched)
    return out


def name_templates(relpaths, top=40):
    """Collapse digit runs to '#' and count, to show the library's naming scheme at a glance."""
    tmpl = Counter()
    example = {}
    for r in relpaths:
        t = re.sub(r"\d+(?:[.,]\d+)?", "#", r)
        tmpl[t] += 1
        example.setdefault(t, r)
    return [(t, n, example[t]) for t, n in tmpl.most_common(top)]
