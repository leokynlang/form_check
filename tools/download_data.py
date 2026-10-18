"""Download and verify the repository's published data bundles (Python stdlib)."""
from __future__ import annotations
import argparse
import hashlib
import json
import shutil
import stat
import zipfile
from pathlib import Path
from urllib.request import Request, urlopen
ROOT=Path(__file__).resolve().parents[1]
def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):h.update(block)
    return h.hexdigest()
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle",action="append",help="Bundle name; repeat to select several. Default: all.")
    parser.add_argument("--list",action="store_true")
    parser.add_argument("--verify-only",action="store_true")
    parser.add_argument("--overwrite",action="store_true",help="Explicitly allow replacement of different local data.")
    parser.add_argument("--manifest",default="data_release/manifest.json",
                        help="Release manifest to restore from; use data_release/manifest_r3.json for the "
                             "R3 record-level bundles.")
    args=parser.parse_args()
    manifest_path=Path(args.manifest)
    if not manifest_path.is_absolute():manifest_path=ROOT/manifest_path
    manifest=json.loads(manifest_path.read_text(encoding="utf-8"))
    available={a["name"]:a for a in manifest["assets"]}
    names=args.bundle or list(available)
    if set(names)-set(available):parser.error("Unknown bundle: "+", ".join(set(names)-set(available)))
    for name in names:
        asset=available[name]
        print(f"{name}: {asset['file_count']} files, {asset['size']/1024/1024:.1f} MiB compressed")
        if args.list:continue
        failures=[]
        if args.verify_only:
            for row in asset["files"]:
                p=ROOT/row["path"]
                if not p.is_file() or p.stat().st_size!=row["size"] or digest(p)!=row["sha256"]:failures.append(row["path"])
            if failures:raise SystemExit("Missing/different files: "+"; ".join(failures[:20]))
            print("  all recorded files match SHA-256")
            continue
        cache=ROOT/"data_release/.downloads";cache.mkdir(parents=True,exist_ok=True)
        target=cache/asset["filename"]
        if not target.exists() or digest(target)!=asset["sha256"]:
            partial=target.with_suffix(target.suffix+".partial")
            with urlopen(Request(asset["url"],headers={"User-Agent":"form-check-data-download"}),timeout=180) as response,partial.open("wb") as out:
                shutil.copyfileobj(response,out,1024*1024)
            if partial.stat().st_size!=asset["size"] or digest(partial)!=asset["sha256"]:raise SystemExit("Archive checksum failed: "+name)
            partial.replace(target)
        expected={r["path"]:r for r in asset["files"]}
        with zipfile.ZipFile(target) as archive:
            members={entry.filename:entry for entry in archive.infolist() if entry.filename!="DATA_BUNDLE_INFO.json"}
            if set(members)!=set(expected):raise SystemExit("Unexpected ZIP member list")
            for rel,entry in members.items():
                p=(ROOT/rel).resolve()
                if not p.is_relative_to(ROOT) or stat.S_ISLNK(entry.external_attr>>16):raise SystemExit("Unsafe ZIP path")
                row=expected[rel]
                if p.exists():
                    if p.is_file() and p.stat().st_size==row["size"] and digest(p)==row["sha256"]:continue
                    if not args.overwrite:raise SystemExit(f"Different local file: {rel}; move it aside or explicitly use --overwrite")
                p.parent.mkdir(parents=True,exist_ok=True)
                staged=p.with_name(p.name+".download-part")
                h=hashlib.sha256();count=0
                with archive.open(entry) as source,staged.open("wb") as out:
                    for chunk in iter(lambda:source.read(1024*1024),b""):
                        h.update(chunk);count+=len(chunk);out.write(chunk)
                if count!=row["size"] or h.hexdigest()!=row["sha256"]:raise SystemExit("Member checksum failed: "+rel)
                staged.replace(p)
        print("  restored and verified")
if __name__=="__main__":main()
