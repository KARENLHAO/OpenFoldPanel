"""Shared helpers for report-facing aggregation and sorting."""

from __future__ import annotations

from openfoldpanel.models import ContactEntry, ContactHit


def scope_sort_key(scope_name: str) -> tuple[int, int]:
    """Return a stable ordering key for global and cluster scopes."""

    if scope_name == "global":
        return (0, 0)
    if scope_name.startswith("cluster_"):
        suffix = scope_name.removeprefix("cluster_")
        if suffix.isdigit():
            return (1, int(suffix))
    return (2, 10**9)


def contact_entry_has_hit(entry: ContactEntry) -> bool:
    """Return whether a contact entry contains one or more resolved hits."""

    if entry.all_contacts:
        return True
    return (
        entry.partner_type is not None
        and entry.min_distance is not None
        and entry.symbol is not None
        and entry.strength_category is not None
    )


def contact_hits_for_entry(entry: ContactEntry) -> list[ContactHit]:
    """Return all concrete contact hits for one residue position."""

    if entry.all_contacts:
        return list(entry.all_contacts)
    if not contact_entry_has_hit(entry):
        return []
    return [
        ContactHit(
            partner_type=entry.partner_type,
            partner_chain=entry.partner_chain,
            partner_resname=entry.partner_resname or "",
            partner_resid=entry.partner_resid or "",
            min_distance=entry.min_distance,
            symbol=entry.symbol,
            strength_category=entry.strength_category,
        )
    ]
