"""Subprocess wrapper with readable error reporting."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


class MissingExecutableError(RuntimeError):
    """Raised when an expected executable is unavailable."""


class ExternalToolError(RuntimeError):
    """Raised when an external command exits unsuccessfully."""


@dataclass(slots=True)
class CommandResult:
    command: list[str]
    returncode: int
    stdout: str
    stderr: str


def which(executable: str) -> str | None:
    """Return the resolved executable path if present."""

    return shutil.which(executable)


def run_command(
    command: Sequence[str],
    *,
    cwd: Path | None = None,
    check: bool = True,
    input_text: str | None = None,
) -> CommandResult:
    """Run a subprocess and capture text output."""

    if not command:
        raise ValueError("command must not be empty")
    if which(command[0]) is None:
        raise MissingExecutableError(f"Executable not found: {command[0]}")

    completed = subprocess.run(
        list(command),
        cwd=str(cwd) if cwd else None,
        input=input_text,
        text=True,
        capture_output=True,
        check=False,
    )
    result = CommandResult(
        command=list(command),
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
    )
    if check and completed.returncode != 0:
        raise ExternalToolError(
            f"Command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return result
