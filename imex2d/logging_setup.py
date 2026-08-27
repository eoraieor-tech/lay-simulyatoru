"""Mərkəzi loglama — C4.

Səbəb: uzun işə salmalarda nə baş verdiyini bilmək lazımdır. `print`
konsola gedir və itir, `QMessageBox` isə skript rejimində ümumiyyətlə
işləmir. Bu modul hər ikisini əvəz edir.

İki hədəf:
    · konsol  — inkişaf zamanı
    · fayl    — sessiya tarixçəsi (istəyə görə, rotasiya ilə)

UI qatı `QtLogHandler` vasitəsilə eyni axını jurnal tabında göstərir,
yəni istifadəçinin gördüyü mətn ilə fayldakı mətn eynidir.
"""

from __future__ import annotations

import logging
import logging.handlers
import os
import sys
from typing import Optional

LOGGER_NAME = "imex2d"
DEFAULT_FORMAT = "%(asctime)s  %(levelname)-8s %(name)-28s %(message)s"
DATE_FORMAT = "%H:%M:%S"

_configured = False


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Modul üçün logger. `get_logger(__name__)` şəklində işlədilir."""
    if name is None or name == LOGGER_NAME:
        return logging.getLogger(LOGGER_NAME)
    if name.startswith("imex2d."):
        return logging.getLogger(name)
    return logging.getLogger(f"{LOGGER_NAME}.{name}")


def configure(level: int = logging.INFO,
              log_file: Optional[str] = None,
              console: bool = True,
              max_bytes: int = 2_000_000,
              backup_count: int = 3,
              force: bool = False) -> logging.Logger:
    """Loglamanı bir dəfə qurur. Təkrar çağırışlar `force` olmadan keçilir."""
    global _configured
    logger = logging.getLogger(LOGGER_NAME)
    if _configured and not force:
        return logger

    logger.setLevel(level)
    logger.propagate = False
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    formatter = logging.Formatter(DEFAULT_FORMAT, DATE_FORMAT)

    if console:
        stream = logging.StreamHandler(sys.stderr)
        stream.setFormatter(formatter)
        logger.addHandler(stream)

    if log_file:
        directory = os.path.dirname(os.path.abspath(log_file))
        if directory:
            os.makedirs(directory, exist_ok=True)
        rotating = logging.handlers.RotatingFileHandler(
            log_file, maxBytes=max_bytes, backupCount=backup_count,
            encoding="utf-8")
        rotating.setFormatter(formatter)
        logger.addHandler(rotating)

    _configured = True
    return logger


def add_handler(handler: logging.Handler) -> None:
    """UI kimi xarici dinləyicilər üçün."""
    handler.setFormatter(logging.Formatter(DEFAULT_FORMAT, DATE_FORMAT))
    logging.getLogger(LOGGER_NAME).addHandler(handler)


def remove_handler(handler: logging.Handler) -> None:
    logging.getLogger(LOGGER_NAME).removeHandler(handler)


def reset() -> None:
    """Testlər üçün — vəziyyəti təmizləyir."""
    global _configured
    logger = logging.getLogger(LOGGER_NAME)
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()
    _configured = False
