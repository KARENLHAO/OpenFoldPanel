"""Write track data as JSON."""

from __future__ import annotations

import json
from pathlib import Path

from openfoldpanel.models import JobReportData, dataclass_to_dict


def write_tracks_json(report_data: JobReportData, path: Path) -> None:
    """Serialize panel data to JSON."""

    payload = dataclass_to_dict(report_data)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
