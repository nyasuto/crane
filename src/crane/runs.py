"""data/runs/<timestamp>_<tag>/ 出力ディレクトリ管理。"""

from datetime import datetime
from pathlib import Path


def new_run_dir(tag: str, base: Path | str = "data/runs") -> Path:
    """タイムスタンプ付き run ディレクトリを作って返す。上書きしない。"""
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = Path(base) / f"{stamp}_{tag}"
    path.mkdir(parents=True, exist_ok=False)
    return path
