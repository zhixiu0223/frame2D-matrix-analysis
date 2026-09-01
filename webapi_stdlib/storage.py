"""
存檔/讀檔的最小抽象層。現在只有 LocalFileStorage(存在本機磁碟的
saved_models/ 資料夾),以後部署到 GCP Cloud Run 時(容器檔案系統是
暫時性的,不能拿來長期存檔),只要照同樣的介面(list/save/load/delete)
另外寫一個 GCSStorage(存 Cloud Storage bucket),換掉 main.py/server.py
裡 instantiate storage 的那一行就好,路由邏輯完全不用動。

檔名做了基本的白名單檢查(只允許英數字/底線/連字號/中文),避免路徑
穿越(../, /)之類的問題。
"""
import json
import re
from pathlib import Path

_NAME_RE = re.compile(r'^[\w\-\u4e00-\u9fff]{1,64}$')


class InvalidNameError(ValueError):
    pass


class NotFoundError(FileNotFoundError):
    pass


class LocalFileStorage:
    def __init__(self, base_dir):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def _path(self, name: str) -> Path:
        if not _NAME_RE.match(name):
            raise InvalidNameError(
                f"檔名 '{name}' 不合法: 只能用英數字、底線、連字號、中文,長度 1-64")
        return self.base_dir / f"{name}.json"

    def list(self):
        return sorted(p.stem for p in self.base_dir.glob("*.json"))

    def save(self, name: str, data: dict):
        path = self._path(name)
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")

    def load(self, name: str) -> dict:
        path = self._path(name)
        if not path.exists():
            raise NotFoundError(f"找不到存檔 '{name}'")
        return json.loads(path.read_text(encoding="utf-8"))

    def delete(self, name: str):
        path = self._path(name)
        if not path.exists():
            raise NotFoundError(f"找不到存檔 '{name}'")
        path.unlink()
