"""Fetch only public numeric inputs, enforcing the checksums used in P03."""
import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    "koa_k2_rv_activity.csv": (
        "https://koa.ipac.caltech.edu/data/Contributed/K2ExoDiversity/I2_spectra/K2_Key_Project_Observation_name_and_RV__activity_Compilation.csv",
        "4cee16fbe7d1a5f8fabc21cbab57fd915fb30b42797800a69fb25a44be0e6d5b"),
    "koa_k2_rv_activity.json": (
        "https://koa.ipac.caltech.edu/Datasets/K2ExoDiversity/K2_RV_Activity.json",
        "8ab68e741e51ac711bce37fac288e711fab891993736a2ea8e8150e600834015"),
    "rvs_toi141.dat": (
        "https://raw.githubusercontent.com/nespinoza/juliet/763f49c4/docs/tutorials/rvs_toi141.dat",
        "a0e0e9d0d218ba971fed61bcabe01c3c368180c2528ed8c2824b6f72b5291506"),
    "KOA_README.txt": (
        "https://koa.ipac.caltech.edu/data/Contributed/K2ExoDiversity/I2_spectra/README",
        "d458bd903f59ad72fcfae070b756618725a503fd38aaba425df9f648418d6201"),
}


def main():
    out = ROOT / "data/raw"
    out.mkdir(parents=True, exist_ok=True)
    for name, (url, expected) in SOURCES.items():
        path = out / name
        if path.exists():
            content = path.read_bytes()
        else:
            with urllib.request.urlopen(url, timeout=90) as response:
                content = response.read()
        found = hashlib.sha256(content).hexdigest()
        if found != expected:
            raise RuntimeError(f"Input changed: {url}; expected {expected}, got {found}. Audit before using.")
        if not path.exists():
            path.write_bytes(content)
        print(json.dumps({"path": str(path.relative_to(ROOT)), "url": url, "sha256": found}))


if __name__ == "__main__":
    main()
