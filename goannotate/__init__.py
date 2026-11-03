"""GOAnnotate: functional annotation of gene models with Diamond searches against UniProt."""

__version__ = "1.0.0"


class GOAnnotateError(Exception):
    """Raised for user-facing errors (missing files, failed external commands)."""
