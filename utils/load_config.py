from pathlib import Path
import yaml


class Config:
    def __init__(self, data, base_dir="."):
        self._data = self._convert(data, Path(base_dir).resolve())

    def _convert(self, obj, base):
        if isinstance(obj, dict):
            return {k: self._convert(v, base) for k, v in obj.items()}

        elif isinstance(obj, list):
            return [self._convert(v, base) for v in obj]

        elif isinstance(obj, str):
            p = Path(obj).expanduser()

            # convert likely paths only
            if "/" in obj or obj.startswith("."):
                if not p.is_absolute():
                    p = base / p
                return p.resolve()

            return obj

        return obj

    def __getattr__(self, name):
        value = self._data[name]

        if isinstance(value, dict):
            return Config(value)

        return value

    def __getitem__(self, key):
        return self._data[key]

    def as_dict(self):
        return self._data


def load_config(file):
    with open(file) as f:
        data = yaml.safe_load(f)

    return Config(data, base_dir=Path(file).parent)
