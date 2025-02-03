from pathlib import Path

from .B256178 import B256178  # IAV infection of MEFs (batch 1)
from .B261790 import B261790  # HSV-1 infection of MEFs (batch 1)
from .B319096 import B319096  # HSV-1 infection of MEFs (batch 2)
from .B831009 import B831009  # HSV-1/IAV RIP from HT-29s

ROOT = Path(__file__).parent

by_assembly = {
    "GRCm39": [
        B256178,  # IAV infection of MEFs (batch 1, Z22 RIP-seq)
        B261790,  # HSV-1 infection of MEFs (batch 1, Z22 RIP-seq)
        B319096,  # HSV-1 infection of MEFs (batch 2, FLAG RIP-seq)
    ],
    "CHM13v2": [
        B831009,  # HSV-1/IAV RIP from HT-29s (FLAG + Z22 RIP-seq)
    ]
}

__all__ = ["ROOT", "by_assembly", "B319096", "B256178", "B261790", "B831009"]
