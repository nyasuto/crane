"""data/runs/<timestamp>_<tag>/ 出力ディレクトリ管理。"""

from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

# run ディレクトリ名は手元の日付で並べたいので JST 固定にする。
JST = ZoneInfo("Asia/Tokyo")


def new_run_dir(tag: str, base: Path | str = "data/runs") -> Path:
    """タイムスタンプ付き run ディレクトリを作って返す。上書きしない。"""
    stamp = datetime.now(JST).strftime("%Y%m%d_%H%M%S")
    path = Path(base) / f"{stamp}_{tag}"
    path.mkdir(parents=True, exist_ok=False)
    return path
