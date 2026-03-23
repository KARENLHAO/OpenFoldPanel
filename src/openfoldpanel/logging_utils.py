"""Logging helpers for console and per-job log files."""

from __future__ import annotations

import logging
from pathlib import Path


def configure_root_logger(verbose: bool = False) -> logging.Logger:
    """Configure the root logger once."""

    logger = logging.getLogger("openfoldpanel")
    if logger.handlers:
        logger.setLevel(logging.DEBUG if verbose else logging.INFO)
        return logger

    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.addHandler(handler)
    logger.propagate = False
    return logger


def attach_file_logger(logger: logging.Logger, log_path: Path) -> logging.Handler:
    """Attach a file handler to the shared project logger."""

    log_path.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s | %(levelname)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)
    return file_handler


def detach_handler(logger: logging.Logger, handler: logging.Handler) -> None:
    """Detach and close a previously attached handler."""

    logger.removeHandler(handler)
    handler.close()
