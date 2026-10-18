"""Package R3 record-level research data for a public GitHub release (stdlib only).

Companion to tools/build_data_release.py, which packaged the R1 snapshot
(tag data-2026-10-18). This tool packages the R3 addition: the P0-1 v2 and P0-2
record-level parquet files. .gitignore keeps those out of the repository
(runs/**/*.parquet) because they are large; they are not otherwise published.

Archive layout and integrity metadata are identical to the R1 snapshot
(DATA_BUNDLE_INFO.json plus a per-member SHA-256 record), so tools/download_data.py
restores either snapshot with the same verification path.

Only files that git reports as *untracked and ignored* are packaged, so this can
never duplicate content already committed to the repository.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = "r3-2026-10-18"
OWNER_REPO = "leokynlang/form_check"
SKIP_PARTS = {"__pycache__", ".git", ".pytest_cache", ".ruff_cache", "tmp"}
BUNDLES = {
    "p01-v2-records": ["runs/p01_v2"],
    "p02-fixed-domain-records": ["runs/p02_fixed_domain"],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def release_files(paths: list[str]) -> list[Path]:
    """Return the untracked-and-ignored files under `paths` (the ones not in the repo)."""
    listed = subprocess.run(
        ["git", "ls-files", "--others", "--ignored", "--exclude-standard", "--", *paths],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout.splitlines()
    files = []
    for rel in listed:
        p = (ROOT / rel).resolve()
        if not p.is_file() or p.is_symlink() or not p.is_relative_to(ROOT):
            continue
        if any(part in SKIP_PARTS for part in p.relative_to(ROOT).parts):
            continue
        files.append(p)
    return sorted(files)


def build(name, paths, out, created_at, store, reuse_existing):
    dest = out / (name + ".zip")
    if dest.exists():
        if not reuse_existing:
            raise RuntimeError(f"Refusing to overwrite existing archive: {dest}")
        with zipfile.ZipFile(dest) as archive:
            records = json.loads(archive.read("DATA_BUNDLE_INFO.json"))["files"]
            actual = {p.relative_to(ROOT).as_posix(): p for p in paths}
            if set(actual) != {row["path"] for row in records}:
                raise RuntimeError("Archive source set changed: " + name)
            for row in records:
                p = actual[row["path"]]
                if p.stat().st_size != row["size"] or sha256(p) != row["sha256"]:
                    raise RuntimeError("Archive source changed: " + row["path"])
            if archive.testzip():
                raise RuntimeError("Archive CRC failed: " + name)
        print(f"{name}: verified and reused", flush=True)
    else:
        records = []
        members = []
        for p in paths:
            rel = p.relative_to(ROOT).as_posix()
            info = zipfile.ZipInfo.from_file(p, arcname=rel)
            info.compress_type = zipfile.ZIP_STORED if store else zipfile.ZIP_DEFLATED
            members.append((p, rel, info))
        mode = zipfile.ZIP_STORED if store else zipfile.ZIP_DEFLATED
        with zipfile.ZipFile(dest, "w", compression=mode, compresslevel=None if store else 6,
                             allowZip64=True) as archive:
            for index, (p, rel, info) in enumerate(members):
                before = p.stat()
                h = hashlib.sha256()
                count = 0
                with p.open("rb") as source, archive.open(info, "w", force_zip64=True) as target:
                    for chunk in iter(lambda: source.read(1024 * 1024), b""):
                        h.update(chunk)
                        count += len(chunk)
                        target.write(chunk)
                after = p.stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    raise RuntimeError(f"File changed while packaging: {rel}; rerun into a new directory")
                records.append({"path": rel, "size": count, "sha256": h.hexdigest()})
                if index % 500 == 0:
                    print(f"{name}: {index + 1}/{len(paths)} files", flush=True)
            archive.writestr("DATA_BUNDLE_INFO.json", json.dumps(
                {"bundle": name, "created_at": created_at, "tag": TAG, "files": records}, indent=2))
        with zipfile.ZipFile(dest) as archive:
            bad = archive.testzip()
            if bad:
                raise RuntimeError("ZIP CRC failed: " + bad)
        print(f"{name}: wrote {dest.stat().st_size / 1024 / 1024:.1f} MiB", flush=True)
    return {
        "name": name,
        "filename": dest.name,
        "size": dest.stat().st_size,
        "sha256": sha256(dest),
        "url": f"https://github.com/{OWNER_REPO}/releases/download/{TAG}/{dest.name}",
        "file_count": len(records),
        "uncompressed_size": sum(row["size"] for row in records),
        "files": records,
    }


def verify(name, paths, out):
    """Independently re-read the archive and recompute every member hash."""
    dest = out / (name + ".zip")
    if not dest.is_file():
        raise SystemExit("Archive missing: " + str(dest))
    actual = {p.relative_to(ROOT).as_posix(): p for p in paths}
    with zipfile.ZipFile(dest) as archive:
        info = json.loads(archive.read("DATA_BUNDLE_INFO.json"))
        if info["bundle"] != name or info["tag"] != TAG:
            raise SystemExit(f"{name}: DATA_BUNDLE_INFO.json disagrees on bundle/tag")
        records = info["files"]
        expected = {row["path"]: row for row in records}
        if len(expected) != len(records):
            raise SystemExit(f"{name}: duplicate member path")
        members = {e.filename for e in archive.infolist() if e.filename != "DATA_BUNDLE_INFO.json"}
        if members != set(expected):
            raise SystemExit(f"{name}: ZIP member set differs from the recorded file list")
        if set(actual) != set(expected):
            raise SystemExit(f"{name}: archive does not match the current untracked source set")
        for index, row in enumerate(records):
            h = hashlib.sha256()
            size = 0
            with archive.open(row["path"]) as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    h.update(block)
                    size += len(block)
            if size != row["size"] or h.hexdigest() != row["sha256"]:
                raise SystemExit(f"{name}: member hash mismatch: {row['path']}")
            source = actual[row["path"]]
            if source.stat().st_size != row["size"] or sha256(source) != row["sha256"]:
                raise SystemExit(f"{name}: source file differs from the archived member: {row['path']}")
            if index % 1000 == 0:
                print(f"{name}: verified {index + 1}/{len(records)} members", flush=True)
        bad = archive.testzip()
        if bad:
            raise SystemExit(f"{name}: ZIP CRC failed: {bad}")
    print(f"{name}: OK - {len(records)} members, {dest.stat().st_size / 1024 / 1024:.1f} MiB, "
          f"sha256 {sha256(dest)}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true",
                        help="Re-verify existing archives and sources without rebuilding.")
    parser.add_argument("--bundle", action="append", choices=sorted(BUNDLES),
                        help="Bundle to build; repeat to select several. Default: all.")
    parser.add_argument("--store", action="store_true",
                        help="Store members uncompressed (parquet is already compressed).")
    parser.add_argument("--reuse-existing", action="store_true",
                        help="Reuse an intact archive only when every current source hash still matches.")
    args = parser.parse_args()

    out = args.output.resolve()
    if not out.is_relative_to(ROOT / "tmp"):
        raise SystemExit("Output must be inside this repository's ignored tmp directory")
    out.mkdir(parents=True, exist_ok=True)

    names = args.bundle or sorted(BUNDLES)
    if args.verify_only:
        for name in names:
            paths = release_files(BUNDLES[name])
            if not paths:
                raise SystemExit(f"No untracked files found for bundle {name}; cannot verify")
            verify(name, paths, out)
        return
    created_at = datetime.now(timezone.utc).isoformat()
    assets = []
    for name in names:
        paths = release_files(BUNDLES[name])
        if not paths:
            raise SystemExit(f"No untracked files found for bundle {name}; refusing to publish an empty archive")
        raw = sum(p.stat().st_size for p in paths) / 1024 / 1024
        print(f"{name}: {len(paths)} files, {raw:.1f} MiB raw", flush=True)
        assets.append(build(name, paths, out, created_at, args.store, args.reuse_existing))

    manifest = {
        "schema_version": 1,
        "repository": OWNER_REPO,
        "release_tag": TAG,
        "created_at": created_at,
        "assets": assets,
        "notes": [
            "R3 record-level data, kept out of the repository by the runs/**/*.parquet ignore rule.",
            "P0-1 v2: independently calibrated joint release; 58,112 fresh datasets over 20 calibration rounds.",
            "P0-2: fixed-domain full grid; 246,240 datasets (46,080 reference, 92,160 calibration, 108,000 test).",
            "Largest member is runs/p02_fixed_domain/stage2/decisions.parquet (about 540 MiB).",
            "Current scientific interpretation is paper/research_R3/manuscript.md.",
            "Integrity: DATA_BUNDLE_INFO.json inside each archive plus SHA256SUMS_r3.txt here.",
        ],
    }
    content = json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
    for target in [ROOT / "data_release/manifest_r3.json", out / "manifest_r3.json"]:
        target.write_text(content, encoding="utf-8")
    sums = "".join(a["sha256"] + "  " + a["filename"] + "\n" for a in assets)
    for target in [ROOT / "data_release/SHA256SUMS_r3.txt", out / "SHA256SUMS_r3.txt"]:
        target.write_text(sums, encoding="utf-8")
    print(json.dumps({"created_at": created_at,
                      "assets": [{k: v for k, v in a.items() if k != "files"} for a in assets]},
                     indent=2), flush=True)


if __name__ == "__main__":
    main()