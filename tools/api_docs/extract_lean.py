#!/usr/bin/env python3
"""Literate Lean chapters: Hugo pages written from commented Lean source.

A repository whose specification lives in commented Lean source lists modules
in ``<repo>/docs/lean/config.toml`` (schema in README.md, §6). Each chapter is a
Hugo section: a landing page plus one page per module, written under the
chapter's ``out`` directory (inside ``docs/content/generated/``, never
committed). A module page reads top to bottom like the source file:

* module docs (``/-! … -/``) become prose; Lean doc comments are Markdown;
* each declaration with a doc comment (``/-- … -/``) becomes a block: its kind,
  its name (an anchor), a link to the source at the building commit, the
  signature as a ``lean`` code block (a structure or inductive type in full),
  and the doc text;
* a theorem is marked ``proved`` when neither its own proof nor any file of the
  Lean package contains ``sorry`` or ``admit`` and the package declares no
  ``axiom``.

Only the source text is read: no Lean toolchain, Mathlib or build is needed.

    extract_lean.py --repo-root .              # every chapter in docs/lean/config.toml
    extract_lean.py --repo-root . --sha abc123 # commit for source links (default: git HEAD)

Python 3.11 or newer (tomllib; `pip install tomli` on older versions), standard
library only.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib  # type: ignore

CONFIG = Path("docs") / "lean" / "config.toml"

# Declaration commands, longest spellings first so `class inductive` beats `class`.
DECL_KEYWORDS = (
    "class inductive", "theorem", "lemma", "def", "abbrev", "instance", "structure",
    "class", "inductive", "example", "opaque", "axiom", "macro_rules", "macro",
    "syntax", "notation", "infixl", "infixr", "infix", "prefix", "postfix", "elab",
    "declare_syntax_cat",
)
MODIFIERS = ("private", "protected", "noncomputable", "partial", "unsafe", "nonrec",
             "scoped", "local")
# Shown whole: the body is the definition (fields, constructors, the rule).
WHOLE = {"structure", "class", "inductive", "class inductive", "notation", "infixl",
         "infixr", "infix", "prefix", "postfix", "syntax", "declare_syntax_cat",
         "macro_rules"}
THEOREMS = {"theorem", "lemma"}
ANONYMOUS = {"example", "macro_rules", "notation", "infixl", "infixr", "infix", "prefix",
             "postfix", "syntax", "macro", "elab"}
OPEN = "([{⟨⦃"
CLOSE = ")]}⟩⦄"
SENTINEL = "\x01"  # marks where a comment or string starts in the masked text


class LeanError(Exception):
    pass


def _is_ident(c: str) -> bool:
    return c.isalnum() or c in "_'!?"


# ---- lexing ---------------------------------------------------------------


@dataclass
class Span:
    kind: str  # "moddoc", "doc", "comment", "string"
    start: int
    end: int  # exclusive


def scan(text: str) -> tuple[list[Span], str]:
    """Comments and string literals of Lean source, and a masked copy of it.

    The masked copy has the same length and line breaks; every comment and
    string is blanked to spaces except its first character, which becomes
    SENTINEL. Searching it for brackets, `:=` or `sorry` therefore sees code
    only, and a line whose first column is not blank starts a new command.
    """
    spans: list[Span] = []
    n = len(text)
    i = 0
    char_lit = re.compile(r"'(?:\\(?:x[0-9a-fA-F]{2}|u\{[0-9a-fA-F]+\}|.)|[^\\'\n])'")
    raw_lit = re.compile(r'r(#*)"')
    while i < n:
        c = text[i]
        prev = text[i - 1] if i else ""
        if text.startswith("--", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            spans.append(Span("comment", i, j))
            i = j
        elif text.startswith("/-", i):
            kind = "comment"
            if text.startswith("/-!", i):
                kind = "moddoc"
            elif text.startswith("/--", i) and not text.startswith("/--/", i):
                kind = "doc"
            depth, j = 1, i + 2
            while j < n and depth:
                if text.startswith("/-", j):
                    depth += 1
                    j += 2
                elif text.startswith("-/", j):
                    depth -= 1
                    j += 2
                else:
                    j += 1
            if depth:
                raise LeanError(f"unterminated block comment at offset {i}")
            spans.append(Span(kind, i, j))
            i = j
        elif c == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 2 if text[j] == "\\" else 1
            spans.append(Span("string", i, min(j + 1, n)))
            i = j + 1
        elif c == "r" and not _is_ident(prev) and (m := raw_lit.match(text, i)):
            close = '"' + m.group(1)
            j = text.find(close, m.end())
            j = n if j < 0 else j + len(close)
            spans.append(Span("string", i, j))
            i = j
        elif c == "'" and not _is_ident(prev) and (m := char_lit.match(text, i)):
            spans.append(Span("string", i, m.end()))
            i = m.end()
        else:
            i += 1
    masked = list(text)
    for s in spans:
        for k in range(s.start, s.end):
            if masked[k] != "\n":
                masked[k] = " "
        masked[s.start] = SENTINEL
    return spans, "".join(masked)


def comment_body(raw: str) -> str:
    """The Markdown inside `/-! … -/` or `/-- … -/`, dedented."""
    inner = raw[3:-2]
    first, _, rest = inner.partition("\n")
    body = first.strip() + ("\n" + textwrap.dedent(rest) if rest else "")
    return body.strip("\n").rstrip()


# ---- modules ---------------------------------------------------------------


@dataclass
class Decl:
    kind: str
    name: str  # as written ("" when anonymous)
    full: str  # namespace-qualified ("" when anonymous)
    anchor: str
    signature: str
    doc: str
    line: int
    end_line: int
    escape: bool = False  # the proof text contains sorry/admit


@dataclass
class Module:
    name: str
    path: str  # relative to the repo root
    title: str = ""
    summary: str = ""
    imports: list[str] = field(default_factory=list)
    items: list = field(default_factory=list)  # str (prose) or Decl, in source order
    namespaces: list[str] = field(default_factory=list)  # namespace at each item
    undocumented_theorems: int = 0
    slug: str = ""

    @property
    def all_namespaces(self) -> tuple[str, ...]:
        """Every namespace the module's items sit in, first seen first."""
        return tuple(dict.fromkeys(n for n in self.namespaces if n))

    @property
    def decls(self) -> list[Decl]:
        return [i for i in self.items if isinstance(i, Decl)]


_CMD = re.compile(
    r"(?P<attrs>(?:@\[[^\]\n]*\]\s*)*)"
    r"(?P<mods>(?:(?:" + "|".join(MODIFIERS) + r")\s+)*)"
    r"(?P<kw>" + "|".join(k.replace(" ", r"\s+") for k in DECL_KEYWORDS) + r")(?![\w'])"
)
_NAME = re.compile(r"\s+((?:«[^»]*»|[^\s:({\[⦃⟨]+))")
_SCOPE = re.compile(
    r"^(?:noncomputable\s+)?(?P<kw>namespace|section|mutual|end)\b[ \t]*(?P<name>[^\s-]*)", re.M
)
_THEOREM_LINE = re.compile(
    r"^(?:@\[[^\]\n]*\]\s*)*(?:(?:" + "|".join(MODIFIERS) + r")\s+)*(?:theorem|lemma)\s", re.M
)
_ESCAPE = re.compile(r"(?<![\w.'])(?:sorry|admit)(?![\w'!?])")
_AXIOM = re.compile(r"^(?:(?:" + "|".join(MODIFIERS) + r")\s+)*axiom\s", re.M)


def summary(para: str) -> str:
    """The first sentence of a paragraph, for a table of contents."""
    para = " ".join(para.split())
    m = re.search(r"[.:;](?=\s+[A-Z`(])", para)
    first = para[: m.start() + 1] if m else para
    return first[:-1] + "." if first.endswith((":", ";")) else first


def _line_of(text: str, pos: int) -> int:
    return text.count("\n", 0, pos) + 1


def _command_end(masked: str, start: int) -> int:
    """Where the command starting at `start` ends: the next line that starts a
    command (first column not blank, and not a `|` alternative or a
    `deriving` clause), or the end of the text."""
    pos = masked.find("\n", start)
    while pos >= 0:
        line_start = pos + 1
        if line_start >= len(masked):
            return len(masked)
        c = masked[line_start]
        if not c.isspace():
            rest = masked[line_start:line_start + 40]
            continuation = c == "|" or (
                re.match(r"deriving\s", rest) and not re.match(r"deriving\s+instance\b", rest)
            ) or re.match(r"(termination_by|decreasing_by)\b", rest)
            if not continuation:
                return line_start
        pos = masked.find("\n", line_start)
    return len(masked)


def _signature_end(masked: str, kw: str, start: int, end: int) -> tuple[int, int | None]:
    """End of the shown signature, and where the body (proof) starts, if any."""
    if kw in WHOLE:
        return end, None
    depth = 0
    i = start
    while i < end:
        c = masked[i]
        if c in OPEN:
            depth += 1
        elif c in CLOSE:
            depth = max(0, depth - 1)
        elif depth == 0:
            if masked.startswith(":=", i):
                return i, i
            if kw in ("macro", "elab") and masked.startswith("=>", i):
                return i, i
            if masked.startswith("where", i) and not _is_ident(masked[i - 1]) and (
                i + 5 >= end or not _is_ident(masked[i + 5])
            ):
                return i, i
            if c == "\n":
                m = re.match(r"\n[ \t]*\|", masked[i:i + 80])
                if m:
                    return i, i
        i += 1
    return end, None


def _anchor(text: str) -> str:
    return re.sub(r"\s+", "_", text.replace("«", "").replace("»", ""))


def parse_module(text: str, name: str, path: str) -> Module:
    """Read one Lean source file into prose blocks and documented declarations."""
    spans, masked = scan(text)
    mod = Module(name=name, path=path)
    mod.imports = re.findall(r"^import\s+(\S+)", masked, re.M)

    events: list[tuple[int, str, object]] = []
    for s in spans:
        if s.kind in ("moddoc", "doc"):
            events.append((s.start, s.kind, s))
    for m in _SCOPE.finditer(masked):
        events.append((m.start(), m.group("kw"), m.group("name")))
    events.sort(key=lambda e: e[0])

    stack: list[tuple[str, str]] = []  # (kind, name)
    consumed = 0
    used_anchors: set[str] = set()
    documented_theorems = 0
    for pos, kind, payload in events:
        if pos < consumed:
            continue
        if kind == "namespace":
            stack.append(("namespace", payload))
            continue
        if kind in ("section", "mutual"):
            stack.append((kind, payload))
            continue
        if kind == "end":
            if stack:
                stack.pop()
            continue
        ns = ".".join(n for k, n in stack if k == "namespace" and n)
        span: Span = payload
        body = comment_body(text[span.start:span.end])
        if kind == "moddoc":
            if not mod.title and not mod.items:
                m = re.match(r"#\s+(.+)\n?", body)
                if m:
                    mod.title = m.group(1).strip()
                    body = body[m.end():].lstrip("\n")
                    para = re.split(r"\n\s*\n", body, maxsplit=1)[0].strip()
                    if para and not para.startswith(("#", "```", "|", "*", "-")):
                        mod.summary = summary(para)
            if body:
                mod.items.append(body)
                mod.namespaces.append(ns)
            continue
        # A declaration doc comment: is a command right after it?
        after = span.end
        while after < len(masked) and masked[after].isspace():
            after += 1
        m = _CMD.match(masked, after)
        if not m:
            continue  # a field or constructor doc, shown inside its declaration
        kw = re.sub(r"\s+", " ", m.group("kw"))
        cmd_end = _command_end(masked, after)
        sig_end, body_start = _signature_end(masked, kw, m.end(), cmd_end)
        decl_name = ""
        if kw not in ANONYMOUS:
            nm = _NAME.match(masked, m.end())
            decl_name = nm.group(1) if nm else ""
        full = ""
        if decl_name:
            full = decl_name[len("_root_."):] if decl_name.startswith("_root_.") else (
                f"{ns}.{decl_name}" if ns else decl_name)
        line = _line_of(text, after)
        base = _anchor(full) if full else f"{name}.L{line}"
        anchor, k = base, 2
        while anchor in used_anchors:
            anchor, k = f"{base}-{k}", k + 1
        used_anchors.add(anchor)
        signature = text[after:sig_end].rstrip()
        escape = False
        if kw in THEOREMS:
            documented_theorems += 1
            proof = masked[body_start:cmd_end] if body_start is not None else ""
            escape = bool(_ESCAPE.search(proof))
        end_line = _line_of(text, len(text[:cmd_end].rstrip()))
        mod.items.append(Decl(kind=kw, name=decl_name, full=full, anchor=anchor,
                              signature=signature, doc=body, line=line,
                              end_line=max(line, end_line), escape=escape))
        mod.namespaces.append(ns)
        consumed = cmd_end
    mod.undocumented_theorems = len(_THEOREM_LINE.findall(masked)) - documented_theorems
    if not mod.title:
        mod.title = name
    return mod


def package_is_clean(root: Path) -> bool:
    """No file of the Lean package contains `sorry`/`admit` or declares an axiom."""
    for p in sorted(root.rglob("*.lean")):
        if ".lake" in p.parts:
            continue
        _, masked = scan(p.read_text(encoding="utf-8"))
        if _ESCAPE.search(masked) or _AXIOM.search(masked):
            return False
    return True


# ---- module lists -----------------------------------------------------------


def module_name(root: Path, path: Path) -> str:
    return ".".join(path.relative_to(root).with_suffix("").parts)


def expand_modules(root: Path, patterns: list[str]) -> tuple[list[str], list[str]]:
    """Module names in configuration order, and warnings.

    `A.B` is one module (it must exist), `A.B.*` its direct submodules and
    `A.B.**` every module below it, each set sorted; a module already listed is
    not repeated, so an explicit list can put a few modules first and a glob
    collect the rest."""
    seen: list[str] = []
    warnings: list[str] = []
    for pat in patterns:
        if pat.endswith((".**", ".*")):
            deep = pat.endswith(".**")
            base = root.joinpath(*pat.rsplit(".", 1)[0].split("."))
            found = sorted(base.rglob("*.lean") if deep else base.glob("*.lean"))
            names = [module_name(root, p) for p in found if ".lake" not in p.parts]
            if not names:
                warnings.append(f"module pattern {pat!r} matches no module under {root}")
        else:
            if not root.joinpath(*pat.split(".")).with_suffix(".lean").is_file():
                raise LeanError(f"module {pat} not found ({root.joinpath(*pat.split('.'))}.lean)")
            names = [pat]
        for nm in names:
            if nm not in seen:
                seen.append(nm)
    return seen, warnings


# ---- rendering --------------------------------------------------------------


_FENCE = re.compile(r"^\s{0,3}(`{3,}|~{3,})")
_CODE_SPAN = re.compile(r"(?<![`\[])`([^`\n]+)`(?!`)(?!\]\()")
_IDENT = re.compile(r"^[^\W\d][\w.'!?]*$")


class Linker:
    """Turns code spans naming a chapter declaration or module into links."""

    def __init__(self, modules: list[Module], url_of_page):
        self.decls: dict[str, str] = {}
        self.modules: dict[str, str] = {}
        for m in modules:
            page = url_of_page(m)
            self.modules[m.name] = page
            for d in m.decls:
                if d.full:
                    self.decls.setdefault(d.full, f"{page}#{d.anchor}")

    def resolve(self, name: str, ns: str, page: str, also: tuple[str, ...] = ()) -> str | None:
        """The URL of the module or declaration `name` means at namespace `ns`
        (then at the namespaces in `also`), or None."""
        if name in self.modules:
            return self.modules[name]
        for space in (ns, *also):
            parts = space.split(".") if space else []
            for k in range(len(parts), -1, -1):
                url = self.decls.get(".".join(parts[:k] + [name]))
                if url:
                    return url[len(page):] if url.startswith(page + "#") else url
        return None

    def link(self, md: str, ns: str, page: str, also: tuple[str, ...] = ()) -> str:
        out: list[str] = []
        fence: str | None = None
        for line in md.splitlines(keepends=True):
            m = _FENCE.match(line)
            if fence:
                out.append(line)
                if m and m.group(1)[0] == fence[0] and len(m.group(1)) >= len(fence):
                    fence = None
                continue
            if m:
                fence = m.group(1)
                out.append(line)
                continue

            def repl(mm: re.Match) -> str:
                name = mm.group(1)
                if not _IDENT.match(name):
                    return mm.group(0)
                url = self.resolve(name, ns, page, also)
                return f"[`{name}`]({url})" if url else mm.group(0)

            out.append(_CODE_SPAN.sub(repl, line))
        return "".join(out)


def _yaml_str(s: str) -> str:
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _plain(s: str) -> str:
    return s.replace("`", "")


@dataclass
class Chapter:
    title: str
    out: Path
    url: str  # project-root-relative URL of the chapter ("/generated/…")
    modules: list[Module]
    description: str = ""
    intro: str = ""


class Renderer:
    def __init__(self, repo_url: str, sha: str, clean: bool):
        self.repo_url = repo_url.rstrip("/")
        self.sha = sha
        self.clean = clean

    def source_url(self, path: str, line: int | None = None, end: int | None = None) -> str:
        url = f"{self.repo_url}/blob/{self.sha}/{path}"
        if line:
            url += f"#L{line}"
            if end and end != line:
                url += f"-L{end}"
        return url

    def decl_block(self, mod: Module, d: Decl, ns: str, linker: Linker, page: str) -> str:
        cls = "lean-decl lean-theorem" if d.kind in THEOREMS else "lean-decl"
        shown = d.name or ""
        head = [f'<span class="lean-kind">{d.kind}</span>']
        if shown:
            head.append(f'<a class="lean-name" href="#{d.anchor}" title="{_html(d.full)}">'
                        f"<code>{_html(shown)}</code></a>")
        else:
            head.append(f'<a class="lean-name" href="#{d.anchor}">¶</a>')
        if d.kind in THEOREMS:
            if d.escape:
                head.append('<span class="lean-badge lean-badge-sorry" title="the proof text '
                            'contains sorry or admit">sorry</span>')
            elif self.clean:
                head.append('<span class="lean-badge lean-badge-proved" title="no sorry or '
                            'admit in the proof or anywhere in the Lean package, and no axiom">'
                            "proved</span>")
        head.append(f'<a class="lean-src" href="{self.source_url(mod.path, d.line, d.end_line)}">'
                    "source</a>")
        fence = "```"
        while fence in d.signature:
            fence += "`"
        parts = [
            f'<div class="{cls}" id="{_html(d.anchor)}">',
            f'<div class="lean-decl-head">{" ".join(head)}</div>',
            "",
            f"{fence}lean\n{d.signature}\n{fence}",
            "",
        ]
        if d.doc:
            parts += [linker.link(d.doc, ns, page, mod.all_namespaces), ""]
        parts.append("</div>")
        return "\n".join(parts)

    def module_page(self, chapter: Chapter, mod: Module, weight: int, linker: Linker) -> str:
        page = f"{chapter.url}/{mod.slug}"
        names = {m.name: m for m in chapter.modules}
        imports = []
        for imp in mod.imports:
            if imp in names:
                imports.append(f'<a href="{chapter.url}/{names[imp].slug}"><code>{imp}</code></a>')
            else:
                imports.append(f"<code>{imp}</code>")
        meta = (f'<p class="lean-meta">Lean module <code>{mod.name}</code>, '
                f'<a href="{self.source_url(mod.path)}"><code>{mod.path}</code></a> '
                f"at <code>{self.sha[:12]}</code>")
        if imports:
            meta += f'<br>Imports {", ".join(imports)}'
        meta += "</p>"
        out = [
            "---",
            f"title: {_yaml_str(_plain(mod.title))}",
            "type: lean",
            f"lean_module: {_yaml_str(mod.name)}",
            f"weight: {weight}",
            "---",
            "",
            meta,
            "",
        ]
        for item, ns in zip(mod.items, mod.namespaces):
            if isinstance(item, Decl):
                out += [self.decl_block(mod, item, ns, linker, page), ""]
            else:
                out += [linker.link(item, ns, page, mod.all_namespaces), ""]
        if mod.undocumented_theorems > 0:
            n = mod.undocumented_theorems
            out += [(f'<p class="lean-meta">This module also proves {n} '
                     f'{"lemma" if n == 1 else "lemmas"} without doc comments; they are in '
                     f'<a href="{self.source_url(mod.path)}">the source</a>.</p>'), ""]
        return "\n".join(out).rstrip() + "\n"

    def landing_page(self, chapter: Chapter) -> str:
        theorems = sum(1 for m in chapter.modules for d in m.decls if d.kind in THEOREMS)
        out = [
            "---",
            f"title: {_yaml_str(chapter.title)}",
            "type: lean",
        ]
        if chapter.description:
            out.append(f"description: {_yaml_str(chapter.description)}")
        out += ["---", ""]
        if chapter.intro:
            out += [chapter.intro.strip(), ""]
        badge = ("Every theorem marked **proved** has no `sorry` or `admit` in its proof, and "
                 "no file of the Lean package contains either or declares an axiom. "
                 if self.clean else "")
        out += [
            (f"Generated from the Lean source at commit [`{self.sha[:12]}`]"
             f"({self.repo_url}/tree/{self.sha}), without running Lean. Module comments are the "
             "prose; each documented declaration shows its statement, its doc comment and a "
             f"link to its source. {len(chapter.modules)} modules, {theorems} documented "
             f"theorems. {badge}").rstrip(),
            "",
        ]
        group = None
        for m in chapter.modules:
            parent = m.name.rsplit(".", 1)[0]
            if parent != group:
                group = parent
                out += [f"## `{parent}`", ""]
            entry = f"- [{_plain(m.title)}]({chapter.url}/{m.slug}) (`{m.name}`)"
            if m.summary:
                entry += f": {m.summary}"
            out.append(entry)
            nxt = chapter.modules.index(m) + 1
            if nxt == len(chapter.modules) or chapter.modules[nxt].name.rsplit(".", 1)[0] != group:
                out.append("")
        return "\n".join(out).rstrip() + "\n"

    def write(self, chapter: Chapter) -> int:
        linker = Linker(chapter.modules, lambda m: f"{chapter.url}/{m.slug}")
        if chapter.out.exists():
            shutil.rmtree(chapter.out)
        chapter.out.mkdir(parents=True)
        (chapter.out / "_index.md").write_text(self.landing_page(chapter), encoding="utf-8")
        for w, m in enumerate(chapter.modules, start=1):
            (chapter.out / f"{m.slug}.md").write_text(
                self.module_page(chapter, m, w, linker), encoding="utf-8")
        return len(chapter.modules)


def _html(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def slug(module: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", module.lower()).strip("-")


# ---- driver -------------------------------------------------------------------


def load_config(repo_root: Path, path: Path | None = None) -> dict:
    path = path or repo_root / CONFIG
    if not path.is_file():
        raise LeanError(f"{path}: not found (see tools/api_docs/README.md, §6)")
    with path.open("rb") as f:
        cfg = tomllib.load(f)
    if not cfg.get("root"):
        raise LeanError(f"{path}: `root` (the Lean package directory) is required")
    chapters = cfg.get("chapter") or []
    if not chapters:
        raise LeanError(f"{path}: no [[chapter]] tables")
    for i, ch in enumerate(chapters):
        for key in ("title", "out", "modules"):
            if not ch.get(key):
                raise LeanError(f"{path}: chapter[{i}] needs `{key}`")
    return cfg


def site_info(repo_root: Path) -> dict:
    p = repo_root / "docs" / "site.toml"
    if not p.is_file():
        return {}
    with p.open("rb") as f:
        return tomllib.load(f)


def git_sha(repo_root: Path) -> str:
    try:
        return subprocess.run(["git", "-C", str(repo_root), "rev-parse", "HEAD"], check=True,
                              capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return os.environ.get("GITHUB_SHA", "") or "main"


def build(repo_root: Path, sha: str | None = None, repo_url: str | None = None,
          config: Path | None = None) -> tuple[list[Chapter], list[str]]:
    repo_root = repo_root.resolve()
    cfg = load_config(repo_root, config)
    site = site_info(repo_root)
    if not repo_url:
        repo_url = cfg.get("repo_url") or (
            f"https://github.com/{site['repo']}" if site.get("repo") else "")
    if not repo_url:
        raise LeanError("no repository URL: set `repo` in docs/site.toml, `repo_url` in "
                        "docs/lean/config.toml, or pass --repo-url")
    sha = sha or git_sha(repo_root)
    root = repo_root / cfg["root"]
    content = (repo_root / site.get("content", "docs/content")).resolve()
    clean = package_is_clean(root)
    renderer = Renderer(repo_url, sha, clean)
    warnings: list[str] = []
    chapters: list[Chapter] = []
    for ch in cfg["chapter"]:
        out = (repo_root / ch["out"]).resolve()
        try:
            rel = out.relative_to(content)
        except ValueError:
            raise LeanError(f"chapter out {ch['out']!r} is not inside the content "
                            f"directory {content}") from None
        names, warn = expand_modules(root, ch["modules"])
        warnings += warn
        modules = []
        for nm in names:
            p = root.joinpath(*nm.split(".")).with_suffix(".lean")
            mod = parse_module(p.read_text(encoding="utf-8"), nm,
                               p.relative_to(repo_root).as_posix())
            mod.slug = slug(nm)
            modules.append(mod)
        intro = ""
        if ch.get("intro"):
            ip = repo_root / ch["intro"]
            if not ip.is_file():
                raise LeanError(f"chapter intro {ch['intro']} not found")
            intro = ip.read_text(encoding="utf-8")
        chapter = Chapter(title=ch["title"], out=out, url="/" + rel.as_posix(), modules=modules,
                          description=ch.get("description", ""), intro=intro)
        renderer.write(chapter)
        chapters.append(chapter)
    return chapters, warnings


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--config", type=Path, default=None,
                    help="chapter configuration (default: <repo>/docs/lean/config.toml)")
    ap.add_argument("--sha", default=None, help="commit for source links (default: git HEAD)")
    ap.add_argument("--repo-url", default=None,
                    help="default: https://github.com/<repo in site.toml>")
    ap.add_argument("--strict", action="store_true", help="warnings fail the run")
    args = ap.parse_args(argv)
    try:
        chapters, warnings = build(args.repo_root, args.sha, args.repo_url, args.config)
    except LeanError as e:
        print(f"extract_lean.py: error: {e}", file=sys.stderr)
        return 1
    for w in warnings:
        print(f"extract_lean.py: warning: {w}", file=sys.stderr)
    for ch in chapters:
        n = sum(len(m.decls) for m in ch.modules)
        print(f"extract_lean.py: {ch.title}: {len(ch.modules)} modules, {n} declarations "
              f"-> {ch.out}")
    return 1 if (warnings and args.strict) else 0


if __name__ == "__main__":
    sys.exit(main())
