from pathlib import Path

from .PRJEB75711 import PRJEB75711
from .PRJNA256013 import PRJNA256013
from .PRJNA382632 import PRJNA382632
from .PRJNA637636 import PRJNA637636

ROOT = Path(__file__).parent

by_assembly = {
    "CHM13v2": [
        PRJNA256013,
        PRJNA382632,
        PRJNA637636,
        PRJEB75711
    ],
}

__all__ = [
    "ROOT", "by_assembly", "PRJNA382632", "PRJNA256013", "PRJNA637636", "PRJEB75711"
]
