"""Fetch pinned ExoMol data separately from the MIT software source."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen


def main():
    folder=Path(__file__).resolve().parents[1]/'data/exomol/CO_Li2015'
    manifest=json.loads((folder/'manifest.json').read_text())
    for item in manifest['files']:
        path=folder/item['name']
        if path.exists():payload=path.read_bytes()
        else:
            with urlopen(item['url'],timeout=60) as response:payload=response.read()
        if hashlib.sha256(payload).hexdigest()!=item['sha256']:
            raise RuntimeError(f"Checksum mismatch: {item['name']}; do not silently replace pinned data")
        if not path.exists():path.write_bytes(payload)
        print('VERIFIED',item['name'])


if __name__=='__main__':main()
