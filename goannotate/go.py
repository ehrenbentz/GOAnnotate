"""GO hierarchy (OBO) parsing and GO OBO retrieval."""

import logging
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Optional

from goannotate import GOAnnotateError, __version__

logger = logging.getLogger(__name__)

GO_OBO_URL = "https://purl.obolibrary.org/obo/go.obo"
GO_OBO_MAX_AGE_DAYS = 30

# The three ontology roots; assigned by ND ("no biological data") annotations only
GO_ROOT_TERMS = frozenset({"GO:0008150", "GO:0003674", "GO:0005575"})

# GO evidence codes by class (http://geneontology.org/docs/guide-go-evidence-codes/)
EVIDENCE_CLASSES = {
    "experimental": frozenset({"EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                               "HTP", "HDA", "HMP", "HGI", "HEP"}),
    "phylogenetic": frozenset({"IBA", "IBD", "IKR", "IRD"}),
    "computational": frozenset({"ISS", "ISO", "ISA", "ISM", "IGC", "RCA"}),
    "author": frozenset({"TAS", "NAS"}),
    "curator": frozenset({"IC"}),
    "no_data": frozenset({"ND"}),
    "electronic": frozenset({"IEA"}),
}
ELECTRONIC_CODES = EVIDENCE_CLASSES["electronic"]
NO_DATA_CODES = EVIDENCE_CLASSES["no_data"]

# --go-evidence policies: the codes a term needs at least one of to be used.
# ND is never used (its terms are the roots). Terms without recorded codes
# (databases built by earlier versions) are used under every policy.
EVIDENCE_POLICIES = {
    "all": frozenset().union(*(v for k, v in EVIDENCE_CLASSES.items() if k != "no_data")),
    "non-electronic": frozenset().union(*(v for k, v in EVIDENCE_CLASSES.items()
                                          if k not in ("no_data", "electronic"))),
    "experimental": EVIDENCE_CLASSES["experimental"],
}


def evidence_class(code: str) -> str:
    for name, codes in EVIDENCE_CLASSES.items():
        if code in codes:
            return name
    return "other"


def best_evidence_class(codes) -> str:
    """The strongest class among a set of codes: experimental > phylogenetic > computational > author > curator > electronic."""
    order = ("experimental", "phylogenetic", "computational", "author", "curator", "electronic", "no_data")
    classes = {evidence_class(c) for c in codes}
    for name in order:
        if name in classes:
            return name
    return "unknown" if not codes else "other"


class GOHierarchy:
    """
    Parse GO OBO file and track parent-child relationships.
    Used to filter GO terms to the most specific (leaf-like) terms, and to
    propagate taxon constraints to subtypes and parts.
    """

    def __init__(self, obo_file: str):
        self.names = {}              # GO ID -> name
        self.namespace = {}          # GO ID -> namespace (BP/MF/CC)
        self.parents = defaultdict(set)   # GO ID -> set of parent IDs (is_a)
        self.children = defaultdict(set)  # GO ID -> set of child IDs (is_a)
        self.part_of_children = defaultdict(set)  # whole GO ID -> set of part GO IDs (part_of)
        self.alt_ids = {}            # alt_id -> primary GO ID

        self._parse_obo(obo_file)
        logger.info(
            f"Loaded GO hierarchy: {len(self.names)} terms "
            f"({sum(1 for ns in self.namespace.values() if ns == 'biological_process')} BP, "
            f"{sum(1 for ns in self.namespace.values() if ns == 'molecular_function')} MF, "
            f"{sum(1 for ns in self.namespace.values() if ns == 'cellular_component')} CC)"
        )

    def _parse_obo(self, obo_file: str):
        """Parse the OBO file."""
        current_id = None
        current_name = None
        current_ns = None
        current_alt_ids = []
        is_obsolete = False

        def save():
            if current_id and current_name and current_ns and not is_obsolete:
                self.names[current_id] = current_name
                self.namespace[current_id] = current_ns
                for alt_id in current_alt_ids:
                    self.alt_ids[alt_id] = current_id

        with open(obo_file) as f:
            for line in f:
                line = line.strip()

                if line == "[Term]":
                    save()
                    current_id = None
                    current_name = None
                    current_ns = None
                    current_alt_ids = []
                    is_obsolete = False
                    continue

                if line.startswith("[") and line != "[Term]":
                    save()
                    current_id = None
                    continue

                if line.startswith("id: GO:"):
                    current_id = line.split()[1]
                elif line.startswith("alt_id: GO:"):
                    current_alt_ids.append(line.split()[1])
                elif line.startswith("name:"):
                    current_name = line[5:].strip()
                elif line.startswith("namespace:"):
                    current_ns = line.split()[1]
                elif line.startswith("is_a: GO:"):
                    parent_id = line.split()[1]
                    if current_id:
                        self.parents[current_id].add(parent_id)
                        self.children[parent_id].add(current_id)
                elif line.startswith("relationship: part_of GO:"):
                    whole_id = line.split()[2]
                    if current_id:
                        self.part_of_children[whole_id].add(current_id)
                elif line == "is_obsolete: true":
                    is_obsolete = True

            save()

    def normalize_id(self, go_id: str) -> Optional[str]:
        """Normalize a GO ID, resolving alt_ids to primary IDs."""
        if go_id in self.names:
            return go_id
        if go_id in self.alt_ids:
            return self.alt_ids[go_id]
        return None

    def get_ancestors(self, go_id: str) -> set:
        """Get all ancestors of a GO term."""
        ancestors = set()
        stack = list(self.parents.get(go_id, []))
        while stack:
            parent = stack.pop()
            if parent not in ancestors:
                ancestors.add(parent)
                stack.extend(self.parents.get(parent, []))
        return ancestors

    def get_descendants(self, go_id: str, include_part_of: bool = False) -> set:
        """
        Get all descendants of a GO term through is_a, and through part_of
        as well when include_part_of is set (the relations along which
        taxon constraints propagate).
        """
        descendants = set()
        stack = list(self.children.get(go_id, []))
        if include_part_of:
            stack.extend(self.part_of_children.get(go_id, []))
        while stack:
            child = stack.pop()
            if child not in descendants:
                descendants.add(child)
                stack.extend(self.children.get(child, []))
                if include_part_of:
                    stack.extend(self.part_of_children.get(child, []))
        return descendants

    def filter_to_specific(self, go_terms: set) -> set:
        """
        Filter a set of GO terms to retain only the most specific terms.
        Removes any term that is an ancestor of another term in the set.
        """
        if not go_terms:
            return set()

        normalized = set()
        for term in go_terms:
            norm = self.normalize_id(term)
            if norm:
                normalized.add(norm)

        if not normalized:
            return set()

        specific = set()
        for term in normalized:
            descendants = self.get_descendants(term)
            if not descendants.intersection(normalized):
                specific.add(term)

        return specific

    def filter_by_namespace(self, go_terms: set, namespaces: list) -> set:
        """Filter GO terms to only include those in specified namespaces."""
        ns_map = {
            'BP': 'biological_process',
            'MF': 'molecular_function',
            'CC': 'cellular_component',
            'biological_process': 'biological_process',
            'molecular_function': 'molecular_function',
            'cellular_component': 'cellular_component'
        }
        target_ns = {ns_map.get(ns, ns) for ns in namespaces}
        return {t for t in go_terms if self.namespace.get(t) in target_ns}


def download_go_obo(output_path: str):
    """Download the current GO OBO hierarchy file from the Gene Ontology."""
    import urllib.request
    import urllib.error

    logger.info(f"Downloading GO OBO file from {GO_OBO_URL} ...")

    try:
        req = urllib.request.Request(GO_OBO_URL,
                                     headers={"User-Agent": f"GOAnnotate/{__version__}"})
        with urllib.request.urlopen(req) as response, open(output_path, 'wb') as f:
            f.write(response.read())
    except urllib.error.URLError as e:
        raise GOAnnotateError(f"Failed to download GO OBO file: {e}") from e

    logger.info(f"GO OBO file saved to: {output_path}")


def resolve_go_obo(go_obo_path: Optional[str]) -> str:
    """
    Ensure a valid, up-to-date GO OBO file is available.

    If a path is provided and the file exists and is less than 30 days old,
    it is used as-is. If the file is older than 30 days, a fresh copy is
    re-downloaded to the same path. If no path is provided, defaults to
    ./go.obo in the current working directory.
    """
    if not go_obo_path:
        go_obo_path = "go.obo"

    obo = Path(go_obo_path)

    if obo.exists():
        age_days = (datetime.now() - datetime.fromtimestamp(obo.stat().st_mtime)).days
        if age_days <= GO_OBO_MAX_AGE_DAYS:
            logger.info(f"Using existing GO OBO file: {go_obo_path} ({age_days} days old)")
            return go_obo_path
        logger.warning(f"GO OBO file is {age_days} days old (>{GO_OBO_MAX_AGE_DAYS}); "
                       f"re-downloading to {go_obo_path}")
    else:
        logger.info(f"GO OBO file not found at {go_obo_path}; downloading")

    download_go_obo(str(obo))
    return str(obo)
