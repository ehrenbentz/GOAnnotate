"""
Name and symbol filtering.

The purpose of this module is to stop designations that belong to another
genome or another annotation pipeline from being carried into the genome
being annotated, and to keep names that carry no information from being
applied as products. Text in a database hit falls into three classes:

  A. Genome-specific identifiers (LOC numbers, locus tags, RIKEN clone IDs,
     C#orf# names, KIAA/FLJ/MGC clone IDs). These are removed from names as
     tokens; the name is discarded only when nothing informative remains.
     As gene symbols they are discarded, but only after the symbol has had
     the chance to resolve to a reference-species ortholog symbol.
  B. Genome-specific qualifiers ("-like", "isoform X1", "LOW QUALITY
     PROTEIN:", "(Fragment)", "putative"). These are stripped and the rest
     of the name is used.
  C. Family-level descriptors ("Olfactory receptor", "Zinc finger protein").
     These are correct and transferable, only less specific; they are kept
     and used as the product when no more specific name is available.

Pattern types in the pattern file (all case-insensitive, one per line):

  qualifier_regex:   class B, removed from names first (alias: strip_regex:)
  identifier_regex:  class A, removed from names as tokens
  placeholder:       whole names that carry no information once A and B
                     are removed ("uncharacterized protein", "protein")
  exact:             discard names equal to this text
  contains:          discard names containing this text
  regex:             discard names matching this regex
  generic_regex:     class C, family-level names (kept, deprioritized)
  normalize_regex:   removed only when names are compared, never from output
  symbol_regex:      discard gene symbols matching this regex

Plain lines without a prefix are accepted for compatibility with 1.x files:
6 characters or fewer are treated as exact:, longer lines as contains:.
"""

import logging
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def default_bad_names_path() -> Path:
    """Path of the bad_names.txt file shipped inside the package."""
    return Path(__file__).resolve().parent / "data" / "bad_names.txt"


# Connector words left dangling at either end of a name after a token was removed
_ORPHAN_LEAD = re.compile(r"^(?:like|of|and|or|the|with|for|to)\s+", re.I)
_ORPHAN_TAIL = re.compile(r"\s+(?:like|of|and|or|the|with|for|to)$", re.I)

# Structural test for family-level names (class C) that complements generic_regex:
# a short name ending in a protein-class word, with no designator token
# (digits, a short uppercase code, a Greek letter name, a Roman numeral).
_CLASS_WORDS = {
    "receptor", "kinase", "protein", "transporter", "channel", "hydrolase", "transferase",
    "ligase", "phosphatase", "peptidase", "protease", "proteinase", "reductase", "oxidase",
    "synthase", "synthetase", "dehydrogenase", "factor", "enzyme", "esterase", "lipase",
    "isomerase", "lyase", "helicase", "nuclease", "polymerase", "methyltransferase",
    "acetyltransferase", "glycosyltransferase", "deaminase", "decarboxylase", "carboxylase",
    "transposase", "integrase", "exchanger", "symporter", "antiporter", "regulator",
    "activator", "inhibitor", "repressor", "cyclase", "chaperone", "pump", "carrier",
    "mutase", "epimerase", "aminotransferase", "peroxidase", "catalase", "oxygenase",
    "hydroxylase", "monooxygenase", "dioxygenase", "glucosidase", "galactosidase",
    "mannosidase", "glucuronidase", "sulfotransferase", "sulfatase", "lectin", "ubiquitin-protein ligase",
}
_DESIGNATOR = re.compile(
    r"\d|^[A-Z]{1,3}$|^[A-Z]+\d|^(?:alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|"
    r"mu|nu|xi|omicron|pi|rho|sigma|tau|upsilon|phi|chi|psi|omega)(?:-\d+)?$|^[IVX]{1,4}$")


class NameAudit:
    """Counters of what the filter did to the names and symbols of one run."""

    def __init__(self):
        self.hits_seen = 0
        self.hits_kept = 0
        self.qualifier_strips = Counter()        # pattern -> names changed
        self.identifier_removals = Counter()     # pattern -> names changed
        self.identifier_examples = defaultdict(Counter)   # pattern -> raw name -> n
        self.discards = Counter()                # reason label -> hits discarded
        self.discard_examples = defaultdict(Counter)      # reason label -> name -> n
        self.kept_names = Counter()              # cleaned name -> hits
        self.symbols_raw_bad = 0                 # GN= symbols matching an identifier pattern
        self.symbols_rescued = 0                 # ... of which resolved to a reference symbol
        self.symbol_discards = Counter()         # pattern -> symbols blanked after resolution
        self.symbol_discard_examples = defaultdict(Counter)
        self.kept_symbols = Counter()            # final symbol -> hits


class BadNameFilter:
    """
    Cleans protein names of genome-specific identifiers and qualifiers,
    discards names that carry no information, recognizes family-level
    names, and discards genome-specific gene symbols.

    Method names are kept from 1.x: clean_name(), is_bad_name(),
    is_bad_symbol(). New: is_generic_name(), normalize_for_comparison(),
    explain_discard().
    """

    def __init__(self, patterns_file: str, extra_files: list = None):
        self.qualifier_patterns = []     # (compiled, source text)
        self.identifier_patterns = []
        self.placeholders = set()
        self.exact_matches = set()
        self.substrings = []
        self.regex_patterns = []
        self.generic_patterns = []
        self.normalize_patterns = []
        self.symbol_regex_patterns = []
        self.product_patterns = []       # (compiled, replacement, source text), applied to output products
        self._legacy_plain = 0

        files = [patterns_file] + list(extra_files or [])
        for path in files:
            if not Path(path).exists():
                raise FileNotFoundError(f"Bad names pattern file not found: {path}")
            with open(path) as f:
                self._parse_patterns(f.read(), path)

        if self._legacy_plain:
            logger.info(f"BadNameFilter: {self._legacy_plain} pattern lines without a prefix were read "
                        f"as exact:/contains: (1.x compatibility); prefixes are recommended")
        logger.info(
            f"BadNameFilter: {len(self.qualifier_patterns)} qualifier, {len(self.identifier_patterns)} identifier, "
            f"{len(self.placeholders)} placeholder, {len(self.exact_matches)} exact, {len(self.substrings)} contains, "
            f"{len(self.regex_patterns)} regex, {len(self.generic_patterns)} generic, "
            f"{len(self.normalize_patterns)} normalize, {len(self.symbol_regex_patterns)} symbol, "
            f"{len(self.product_patterns)} product patterns from {', '.join(repr(str(p)) for p in files)}")

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------

    def _compile(self, pattern: str, kind: str, source: str):
        try:
            return re.compile(pattern, re.IGNORECASE)
        except re.error as e:
            logger.warning(f"Invalid {kind} pattern '{pattern}' in {source}: {e}")
            return None

    def _parse_patterns(self, text: str, source: str = ""):
        prefixes = (
            ("qualifier_regex:", "qualifier"), ("strip_regex:", "qualifier"),
            ("identifier_regex:", "identifier"), ("placeholder:", "placeholder"),
            ("exact:", "exact"), ("contains:", "contains"), ("regex:", "regex"),
            ("generic_regex:", "generic"), ("normalize_regex:", "normalize"),
            ("symbol_regex:", "symbol"), ("product_regex:", "product"),
        )
        for line in text.split("\n"):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            kind = None
            body = line
            for prefix, k in prefixes:
                if line.lower().startswith(prefix):
                    kind, body = k, line[len(prefix):].strip()
                    break
            if kind is None:
                # 1.x compatibility: unprefixed lines
                self._legacy_plain += 1
                kind = "exact" if len(line) <= 6 else "contains"
                body = line
            if not body:
                continue

            if kind == "qualifier":
                pat = self._compile(body, "qualifier_regex", source)
                if pat:
                    self.qualifier_patterns.append((pat, body))
            elif kind == "identifier":
                pat = self._compile(body, "identifier_regex", source)
                if pat:
                    self.identifier_patterns.append((pat, body))
            elif kind == "placeholder":
                self.placeholders.add(self._placeholder_key(body))
            elif kind == "exact":
                self.exact_matches.add(body.lower())
            elif kind == "contains":
                self.substrings.append(body.lower())
            elif kind == "regex":
                pat = self._compile(body, "regex", source)
                if pat:
                    self.regex_patterns.append(pat)
            elif kind == "generic":
                pat = self._compile(body, "generic_regex", source)
                if pat:
                    self.generic_patterns.append(pat)
            elif kind == "normalize":
                pat = self._compile(body, "normalize_regex", source)
                if pat:
                    self.normalize_patterns.append(pat)
            elif kind == "symbol":
                pat = self._compile(body, "symbol_regex", source)
                if pat:
                    self.symbol_regex_patterns.append(pat)
            elif kind == "product":
                # "PATTERN => REPLACEMENT"; without "=>" the match is removed
                pattern, sep, repl = body.partition("=>")
                pattern = pattern.strip() if sep else body
                repl = repl.strip() if sep else ""
                pat = self._compile(pattern, "product_regex", source)
                if pat:
                    self.product_patterns.append((pat, repl, body))

    @staticmethod
    def _placeholder_key(name: str) -> str:
        s = name.lower()
        s = re.sub(r"[^a-z0-9]+", " ", s)
        return " ".join(s.split())

    # ------------------------------------------------------------------
    # Names
    # ------------------------------------------------------------------

    @staticmethod
    def _tidy(s: str) -> str:
        """Repair punctuation and whitespace left behind by token removal."""
        s = re.sub(r"\(\s*\)", "", s)
        s = re.sub(r"\[\s*\]", "", s)
        s = re.sub(r"\s*,\s*,", ",", s)
        s = re.sub(r"\s+,", ",", s)
        s = re.sub(r"\s+", " ", s).strip()
        s = re.sub(r"^[\s,;:/\-]+|[\s,;:/\-]+$", "", s)
        if _ORPHAN_LEAD.match(s + " "):
            return ""
        for _ in range(3):
            before = s
            s = _ORPHAN_LEAD.sub("", s)
            s = _ORPHAN_TAIL.sub("", s)
            s = re.sub(r"^[\s,;:/\-]+|[\s,;:/\-]+$", "", s)
            if s == before:
                break
        return s

    def clean_name(self, name: str, audit: NameAudit = None) -> Optional[str]:
        """
        Remove genome-specific qualifiers (class B) and identifiers (class A)
        from a protein name. Returns the cleaned name, or None when nothing
        is left. Runs before is_bad_name() so that informative names
        obscured by identifiers or qualifiers are preserved.
        """
        if not name:
            return None
        cleaned = name
        for pat, text in self.qualifier_patterns:
            new = pat.sub("", cleaned)
            if audit is not None and new != cleaned:
                audit.qualifier_strips[text] += 1
            cleaned = new
        for pat, text in self.identifier_patterns:
            new = pat.sub(" ", cleaned)
            if new != cleaned:
                if audit is not None:
                    audit.identifier_removals[text] += 1
                    if len(audit.identifier_examples[text]) < 50:
                        audit.identifier_examples[text][name[:70]] += 1
                cleaned = new
        cleaned = self._tidy(cleaned)
        return cleaned or None

    def is_placeholder(self, name: str) -> bool:
        return self._placeholder_key(name) in self.placeholders

    def is_bad_name(self, name: str) -> bool:
        """True when a (cleaned) protein name carries no information or matches a discard pattern."""
        return self.explain_discard(name) is not None

    def explain_discard(self, name: str) -> Optional[str]:
        """The reason a (cleaned) name is discarded, or None when it is kept."""
        if not name:
            return "(empty)"
        if self.is_placeholder(name):
            return f"placeholder: {self._placeholder_key(name)}"
        name_lower = name.lower().strip()
        if name_lower in self.exact_matches:
            return f"exact: {name_lower}"
        for substring in self.substrings:
            if substring in name_lower:
                return f"contains: {substring}"
        for pattern in self.regex_patterns:
            if pattern.search(name):
                return f"regex: {pattern.pattern}"
        return None

    def process_name(self, name: str, audit: NameAudit = None) -> Optional[str]:
        """clean_name() followed by the discard test; returns the kept name or None. Records to the audit."""
        if audit is not None:
            audit.hits_seen += 1
        cleaned = self.clean_name(name, audit)
        if cleaned is None:
            if audit is not None:
                audit.discards["(nothing left after removing identifiers and qualifiers)"] += 1
                ex = audit.discard_examples["(nothing left after removing identifiers and qualifiers)"]
                if len(ex) < 50:
                    ex[name[:70]] += 1
            return None
        reason = self.explain_discard(cleaned)
        if reason is not None:
            if audit is not None:
                audit.discards[reason] += 1
                if len(audit.discard_examples[reason]) < 50:
                    audit.discard_examples[reason][cleaned[:70]] += 1
            return None
        if audit is not None:
            audit.hits_kept += 1
            audit.kept_names[cleaned] += 1
        return cleaned

    def is_generic_name(self, name: str) -> bool:
        """
        True for family-level names (class C): matched by a generic_regex
        pattern, or a short name ending in a protein-class word with no
        designator token (a number, a code such as 4K36 or PLK1, a Greek
        letter name), and at most two words.
        """
        if not name:
            return False
        for pattern in self.generic_patterns:
            if pattern.search(name):
                return True
        tokens = re.findall(r"[A-Za-z0-9/\-]+", name)
        if not tokens or len(tokens) > 2:
            return False
        if tokens[-1].lower() not in _CLASS_WORDS:
            return False
        return not any(_DESIGNATOR.search(t) for t in tokens)

    def normalize_for_comparison(self, name: str) -> str:
        """Lowercased name with normalize_regex patterns removed and whitespace collapsed."""
        if not name:
            return ""
        s = name
        for pat in self.normalize_patterns:
            s = pat.sub(" ", s)
        return " ".join(s.lower().split())

    def clean_product(self, name: str, audit: dict = None) -> str:
        """
        Apply the product_regex rules (product naming conventions) to an
        output product name: EC numbers, molecular weights, species names,
        "<SYMBOL> protein" word order, doubled words. Records the rules that
        fired in audit (rule text -> count) when given.
        """
        if not name:
            return ""
        s = name
        for pat, repl, text in self.product_patterns:
            new = pat.sub(repl, s)
            if new != s:
                if audit is not None:
                    audit[text] = audit.get(text, 0) + 1
                s = new
        s = self._tidy(s)
        return s or name

    # ------------------------------------------------------------------
    # Symbols
    # ------------------------------------------------------------------

    def is_bad_symbol(self, symbol: str) -> bool:
        """True when a gene symbol is a genome-specific identifier (symbol_regex patterns)."""
        return self.explain_bad_symbol(symbol) is not None

    def explain_bad_symbol(self, symbol: str) -> Optional[str]:
        if not symbol:
            return "(empty)"
        for pattern in self.symbol_regex_patterns:
            if pattern.search(symbol):
                return pattern.pattern
        return None

    def filter_names(self, names: list) -> list:
        """Return only names that are NOT bad."""
        return [n for n in names if not self.is_bad_name(n)]


# ----------------------------------------------------------------------
# Audit report
# ----------------------------------------------------------------------

_SUSPECT_SYMBOL = re.compile(r"\d{4,}|_|^[A-Z]{1,2}\d{5,}$|^(?:GM|LOC|SI:|ZGC:|WU:|IM:)\d*", re.I)


def write_name_audit(audit: NameAudit, bad_name_filter: BadNameFilter, path, annotations: dict = None,
                     top: int = 40) -> None:
    """Write the per-run name audit: what was removed, what was discarded, what survived."""
    def section(fh, title):
        fh.write(f"\n{title}\n{'-' * len(title)}\n")

    with open(path, "w") as fh:
        fh.write("GOAnnotate name audit\n=====================\n")
        fh.write(f"Hits examined: {audit.hits_seen:,}; kept {audit.hits_kept:,} "
                 f"({100 * audit.hits_kept / max(1, audit.hits_seen):.1f}%); "
                 f"discarded {audit.hits_seen - audit.hits_kept:,}\n")

        section(fh, "Discards by reason")
        for label, n in audit.discards.most_common():
            ex = ", ".join(f"{k} ({v})" for k, v in audit.discard_examples[label].most_common(3))
            fh.write(f"  {n:>9,}  {label:<70} e.g. {ex}\n")

        section(fh, "Identifiers removed from names (class A), by pattern")
        for text, n in audit.identifier_removals.most_common():
            ex = ", ".join(f"{k}" for k, _ in audit.identifier_examples[text].most_common(2))
            fh.write(f"  {n:>9,}  {text:<60} e.g. {ex}\n")
        for _pat, text in bad_name_filter.identifier_patterns:
            if text not in audit.identifier_removals:
                fh.write(f"  {0:>9,}  {text}\n")

        section(fh, "Qualifiers stripped from names (class B), by pattern")
        for text, n in audit.qualifier_strips.most_common():
            fh.write(f"  {n:>9,}  {text}\n")
        for _pat, text in bad_name_filter.qualifier_patterns:
            if text not in audit.qualifier_strips:
                fh.write(f"  {0:>9,}  {text}\n")

        section(fh, f"Most frequent surviving names (top {top}; G = family-level name)")
        for name, n in audit.kept_names.most_common(top):
            flag = "G" if bad_name_filter.is_generic_name(name) else " "
            fh.write(f"  {n:>9,}  {flag} {name}\n")

        section(fh, f"Surviving names that look like identifiers or placeholders (top {top})")
        suspect = re.compile(
            r"\b(cDNA|clone|EST|contig|scaffold|unnamed|unknown|similar to|hypothetical|predicted|"
            r"open reading frame|readthrough|fragment)\b|\bprotein [A-Z]\d{1,2}[A-Z]?\d*$|^\S+ protein$|"
            r"\b(?:UPF|DUF)\d+\b|\b[A-Z]{1,4}\d{4,}\b", re.I)
        found = [(n, name) for name, n in audit.kept_names.items() if suspect.search(name)]
        found.sort(reverse=True)
        fh.write(f"  {len(found):,} distinct names, {sum(n for n, _ in found):,} hits\n")
        for n, name in found[:top]:
            fh.write(f"  {n:>9,}  {name}\n")

        section(fh, "Symbols")
        fh.write(f"  GN= symbols matching an identifier pattern: {audit.symbols_raw_bad:,}; "
                 f"of these {audit.symbols_rescued:,} resolved to a reference-species ortholog symbol\n")
        fh.write("  Symbols blanked after resolution, by pattern:\n")
        for text, n in audit.symbol_discards.most_common():
            ex = ", ".join(k for k, _ in audit.symbol_discard_examples[text].most_common(4))
            fh.write(f"  {n:>9,}  {text:<50} e.g. {ex}\n")

        section(fh, f"Surviving symbols that look like identifiers (top {top})")
        sus = [(n, s) for s, n in audit.kept_symbols.items() if _SUSPECT_SYMBOL.search(s)]
        sus.sort(reverse=True)
        fh.write(f"  {len(sus):,} distinct symbols, {sum(n for n, _ in sus):,} hits\n")
        for n, s in sus[:top]:
            fh.write(f"  {n:>9,}  {s}\n")

        if annotations:
            section(fh, "Products assigned to genes")
            products = Counter(a.consensus_name for a in annotations.values() if a.consensus_name)
            n_generic = sum(n for name, n in products.items() if bad_name_filter.is_generic_name(name))
            n_total = sum(products.values())
            fh.write(f"  genes with a product: {n_total:,}; with a family-level product: {n_generic:,} "
                     f"({100 * n_generic / max(1, n_total):.1f}%)\n")
            fh.write(f"  Most frequent family-level products (top {top}):\n")
            gen = [(n, name) for name, n in products.items() if bad_name_filter.is_generic_name(name)]
            gen.sort(reverse=True)
            for n, name in gen[:top]:
                fh.write(f"  {n:>9,}  {name}\n")
    logger.info(f"Name audit written to {path}")
