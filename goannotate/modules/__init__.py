"""Optional annotation modules wrapping external tools.

Modules are selected with --modules (names below, or auto for every module
whose tool is found). Each adds columns to the annotation tables,
attributes to the annotated GFF3 and a report. A module with
runs_before_vote (orthology) runs before the Diamond vote and may steer it
(prior()) and change symbols afterwards; the others run after the core's
symbol assignment and add to it.
"""

import logging
import time

from goannotate.modules.base import Module, ModuleContext, ModuleError, ModuleSkipped, ToolStatus

logger = logging.getLogger(__name__)


def registry() -> dict:
    """{module name: class}, in the order modules run."""
    from goannotate.modules.interproscan import InterProScanModule
    from goannotate.modules.orthology import OrthologyModule
    from goannotate.modules.signalp import SignalPModule
    # Orthology runs before InterProScan, so that genes named from relatives' orthologs are not
    # given a domain-based name first
    return {
        SignalPModule.name: SignalPModule,
        OrthologyModule.name: OrthologyModule,
        InterProScanModule.name: InterProScanModule,
    }


def parse_selection(value: str) -> list:
    """
    Module names from the --modules value: comma- or space-separated names,
    'auto' (every registered module, run when available), or none/'' for
    no modules. Unknown names raise ValueError.
    """
    if not value or value.strip().lower() in ("none", "off"):
        return []
    known = registry()
    if value.strip().lower() == "auto":
        return ["auto"]
    names = [n.strip().lower() for n in value.replace(",", " ").split() if n.strip()]
    unknown = [n for n in names if n not in known]
    if unknown:
        raise ValueError(f"unknown module(s): {', '.join(unknown)}; available: {', '.join(known)}")
    seen = []
    for n in names:
        if n not in seen:
            seen.append(n)
    return seen


def select_modules(args) -> list:
    """Module instances for the --modules value (auto: all registered)."""
    names = parse_selection(getattr(args, "modules", None))
    if not names:
        return []
    known = registry()
    if names == ["auto"]:
        return [cls(args) for cls in known.values()]
    return [known[n](args) for n in names]


def detect_all(args) -> list:
    """ToolStatus of every registered module, for --check-tools."""
    statuses = []
    for cls in registry().values():
        statuses.append(cls(args).detect())
    return statuses


def _start(module: Module, ctx: ModuleContext, auto: bool) -> bool:
    """Detect and run one module; False when it was skipped or failed (reason in skipped_reason)."""
    status = module.detect()
    if module.needs_proteins and ctx.search_mode != "blastp":
        module.skipped_reason = "needs protein queries (nucleotide input)"
        logger.info(f"Module {module.name} skipped: {module.skipped_reason}")
        return False
    if not status.available:
        module.skipped_reason = status.reason
        if auto:
            logger.info(f"Module {module.name} not run: {status.reason}")
        else:
            logger.warning(f"Module {module.name} requested but not available: {status.reason}")
        return False
    logger.info(f"Module {module.name}: {status.describe()}")
    module._t0 = time.time()
    try:
        module.run(ctx)
    except ModuleSkipped as e:
        module.skipped_reason = str(e)
        logger.warning(f"Module {module.name} not run: {e}")
        return False
    except ModuleError as e:
        module.skipped_reason = str(e)
        logger.error(f"Module {module.name} failed: {e}")
        return False
    return True


def _finish(module: Module, ctx: ModuleContext) -> bool:
    try:
        module.integrate(ctx)
    except ModuleError as e:
        module.skipped_reason = str(e)
        logger.error(f"Module {module.name} failed: {e}")
        return False
    logger.info(f"Module {module.name} done in {(time.time() - module._t0) / 60:.1f} min")
    return True


def run_modules(modules: list, ctx: ModuleContext, auto: bool = False) -> list:
    """
    Detect, run and integrate each module in turn. A module whose tool is
    missing is skipped: silently under auto, with a warning otherwise. A
    module that fails is reported and skipped; the pipeline continues.
    Returns the modules that ran.
    """
    return [m for m in modules if _start(m, ctx, auto) and _finish(m, ctx)]


def start_modules(modules: list, ctx: ModuleContext, auto: bool = False) -> list:
    """Detect and run (not integrate) the modules that run before the vote; returns those that ran."""
    return [m for m in modules if _start(m, ctx, auto)]


def integrate_modules(modules: list, ctx: ModuleContext) -> list:
    """Integrate modules started with start_modules once the annotations exist; returns those that completed."""
    return [m for m in modules if _finish(m, ctx)]


__all__ = ["Module", "ModuleContext", "ModuleError", "ModuleSkipped", "ToolStatus", "registry", "parse_selection",
           "select_modules", "detect_all", "run_modules", "start_modules", "integrate_modules"]
