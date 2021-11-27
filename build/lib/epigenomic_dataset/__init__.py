"""Module offering methods to retrieve data and tasks for CRR predictions."""
from .build import build
from .load_epigenomes import load_epigenomes
from .load_tasks import (
    active_enhancers_vs_active_promoters,
    active_enhancers_vs_inactive_enhancers,
    inactive_enhancers_vs_inactive_promoters,
    active_promoters_vs_inactive_promoters,
    load_all_tasks
)
from .logging import logger

__all__ = [
    "build", "load_epigenomes", "logger",
    "active_enhancers_vs_active_promoters",
    "active_enhancers_vs_inactive_enhancers",
    "inactive_enhancers_vs_inactive_promoters",
    "active_promoters_vs_inactive_promoters",
    "load_all_tasks"
]
