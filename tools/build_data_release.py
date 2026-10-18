"""Package stable local research data for a public GitHub release (stdlib only)."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import re
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
TAG = "data-2026-10-18"
OWNER_REPO = "leokynlang/form_check"
SKIP_PARTS = {"__pycache__", ".git", ".pytest_cache", ".ruff_cache"}
INCLUDED_SOURCES = {"uci", "nist_strd", "usgs", "usgs_eq", "open_meteo", "epa_fueleconomy"}
SOURCE_TERMS = {
    "uci": ("CC BY 4.0", "https://archive.ics.uci.edu/"),
    "nist_strd": ("NIST non-SRD data terms", "https://www.nist.gov/open/license"),
    "usgs": ("USGS public-domain data", "https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits"),
    "usgs_eq": ("USGS public-domain data", "https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits"),
    "open_meteo": ("CC BY 4.0; Open-Meteo / Copernicus ERA5 attribution", "https://open-meteo.com/en/licence"),
    "epa_fueleconomy": ("DOE/EPA FuelEconomy.gov public data", "https://www.fueleconomy.gov/feg/download.shtml"),
    "webbook": ("not mirrored: SRD 69 copyright; obtain from source", "https://webbook.nist.gov/chemistry/"),
    "rdatasets": ("not mirrored: original data licensing not definitive", "https://vincentarelbundock.github.io/Rdatasets/"),
    "openml": ("not mirrored: metadata says Public, without an unambiguous standard license", "https://www.openml.org/"),
}

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

def eligible_files(folder):
    return sorted(p for p in folder.rglob("*") if p.is_file() and
                  not any(part in SKIP_PARTS for part in p.parts) and
                  not p.is_symlink() and p.resolve().is_relative_to(ROOT))

def uci_license_check(uid):
    url = f"https://archive.ics.uci.edu/dataset/{uid}"
    try:
        with urlopen(Request(url, headers={"User-Agent":"form-check-research-publication"}), timeout=35) as response:
            page = response.read().decode("utf-8", errors="replace")
        found = bool(re.search(r"CC(?: |&#32;|&nbsp;|\\u0020|-)*BY(?: |-|&nbsp;)*4\.0|creativecommons\.org/licenses/by/4\.0", page, re.I))
        return str(uid), {"url":url, "cc_by_4_0_found":found}
    except Exception as exc:
        return str(uid), {"url":url,"cc_by_4_0_found":False,"error":type(exc).__name__}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reuse-existing", action="store_true",
                        help="Reuse an intact archive only when all current source hashes still match.")
    parser.add_argument("--reuse-license-checks", action="store_true",
                        help="Reuse previously successful checks and documented manual source reviews.")
    args = parser.parse_args()
    out = args.output.resolve()
    if not out.is_relative_to(ROOT / "tmp"):
        raise SystemExit("Output must be inside this repository's ignored tmp directory")
    out.mkdir(parents=True, exist_ok=True)
    dataset_manifest = json.loads((ROOT/"real_data/clean/manifest.json").read_text(encoding="utf-8"))
    meta_files = sorted((ROOT/"real_data/raw/uci").glob("*_meta.json"))
    uci_ids = {int(p.stem.split("_")[0]) for p in meta_files} | {291}
    license_checks = {}
    check_path=ROOT/"data_release/license_checks.json"
    if args.reuse_license_checks and check_path.exists():
        license_checks={k:v for k,v in json.loads(check_path.read_text(encoding="utf-8")).items()
                        if v.get("cc_by_4_0_found")}
        manual_path=ROOT/"data_release/source_license_reviews.json"
        if manual_path.exists():
            license_checks.update(json.loads(manual_path.read_text(encoding="utf-8")))
    remaining=[uid for uid in sorted(uci_ids) if str(uid) not in license_checks]
    with ThreadPoolExecutor(max_workers=6) as executor:
        license_checks.update(dict(executor.map(uci_license_check, remaining)))
    failed_uci = [uid for uid, entry in license_checks.items() if not entry["cc_by_4_0_found"]]
    print("UCI license pages checked; unresolved IDs:", failed_uci, flush=True)
    (ROOT/"data_release/license_checks.json").write_text(json.dumps(license_checks,indent=2),encoding="utf-8")
    source_rows = []
    real_paths = {ROOT/"real_data/clean/manifest.json"}
    for name, info in sorted(dataset_manifest.items()):
        source = info["source"]
        mirrored = source in INCLUDED_SOURCES
        if source == "uci":
            uid = re.search(r"/dataset/(\d+)", info.get("source_url",""))
            mirrored = bool(uid and license_checks.get(uid.group(1),{}).get("cc_by_4_0_found"))
        source_rows.append({"dataset":name,"source":source,"source_url":info.get("source_url"),
                            "redistributed":mirrored,"terms":SOURCE_TERMS[source][0],
                            "terms_url":SOURCE_TERMS[source][1],
                            "changes":"See clean manifest notes and real_data/prepare_data.py",
                            "raw_files":info.get("raw_files",[])})
        if mirrored:
            clean = ROOT/"real_data/clean"/(name+".csv")
            if clean.exists(): real_paths.add(clean)
            for item in info.get("raw_files",[]):
                raw = ROOT/"real_data"/item["path"].replace("\\","/")
                if raw.exists(): real_paths.add(raw)
    # Include original raw snapshots for accepted source families; UCI is per-ID checked.
    for source in sorted(INCLUDED_SOURCES):
        folder = ROOT/"real_data/raw"/source
        if not folder.exists(): continue
        for p in eligible_files(folder):
            if source == "uci":
                rel = p.relative_to(folder).parts[0]
                match = re.match(r"(\d+)", rel)
                if not match or not license_checks.get(match.group(1),{}).get("cc_by_4_0_found"): continue
            real_paths.add(p)
    (ROOT/"data_release/dataset_sources.json").write_text(json.dumps(source_rows,indent=2,ensure_ascii=False),encoding="utf-8")
    real_paths.update([ROOT/"data_release/dataset_sources.json",ROOT/"data_release/license_checks.json",
                       ROOT/"data_release/THIRD_PARTY_DATA.md",ROOT/"real_data/README.md",
                       ROOT/"real_data/prepare_data.py"])
    core = eligible_files(ROOT/"runs")
    traces = eligible_files(ROOT/"outputs")
    scratch = ROOT/"_tmp_pysr"
    if scratch.exists():
        traces += [p for p in eligible_files(scratch) if p.suffix in {".py",".md",".json",".txt"}]
    bundles = [("experiment-records",core),("engineering-data",sorted(real_paths)),
               ("pysr-search",sorted(traces))]
    created_at = datetime.now(timezone.utc).isoformat()
    def build(bundle):
        name, paths = bundle
        dest = out / (name+".zip")
        if dest.exists():
            if not args.reuse_existing:
                raise RuntimeError(f"Refusing to overwrite existing archive: {dest}")
            with zipfile.ZipFile(dest) as archive:
                records=json.loads(archive.read("DATA_BUNDLE_INFO.json"))["files"]
                actual={p.relative_to(ROOT).as_posix():p for p in paths}
                if set(actual)!={r["path"] for r in records}:
                    raise RuntimeError("Archive source set changed: "+name)
                for row in records:
                    p=actual[row["path"]]
                    if p.stat().st_size!=row["size"] or sha256(p)!=row["sha256"]:
                        raise RuntimeError("Archive source changed: "+row["path"])
                if archive.testzip():raise RuntimeError("Archive CRC failed: "+name)
            print(name+": verified and reused",flush=True)
            return {"name":name,"filename":dest.name,"size":dest.stat().st_size,"sha256":sha256(dest),
                    "url":f"https://github.com/{OWNER_REPO}/releases/download/{TAG}/{dest.name}",
                    "file_count":len(records),"uncompressed_size":sum(r["size"] for r in records),
                    "files":records}
        records = []
        with zipfile.ZipFile(dest,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as archive:
            for index,p in enumerate(paths):
                if p.is_symlink() or not p.resolve().is_relative_to(ROOT): raise RuntimeError("Unsafe input path")
                before=p.stat()
                h=hashlib.sha256()
                count=0
                rel=p.relative_to(ROOT).as_posix()
                zinfo=zipfile.ZipInfo.from_file(p,arcname=rel)
                zinfo.compress_type=zipfile.ZIP_DEFLATED
                with p.open("rb") as source, archive.open(zinfo,"w",force_zip64=True) as target:
                    for chunk in iter(lambda:source.read(1024*1024),b""):
                        h.update(chunk);count+=len(chunk);target.write(chunk)
                after=p.stat()
                if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
                    raise RuntimeError(f"File changed while packaging: {rel}; rerun into a new output directory")
                records.append({"path":rel,"size":count,"sha256":h.hexdigest()})
                if index%500==0: print(f"{name}: {index+1}/{len(paths)} files",flush=True)
            archive.writestr("DATA_BUNDLE_INFO.json",json.dumps({"bundle":name,"created_at":created_at,
                  "tag":TAG,"files":records},indent=2))
        with zipfile.ZipFile(dest) as archive:
            bad=archive.testzip()
            if bad: raise RuntimeError("ZIP CRC failed: "+bad)
        return {"name":name,"filename":dest.name,"size":dest.stat().st_size,"sha256":sha256(dest),
                "url":f"https://github.com/{OWNER_REPO}/releases/download/{TAG}/{dest.name}",
                "file_count":len(records),"uncompressed_size":sum(r["size"] for r in records),
                "files":records}
    with ThreadPoolExecutor(max_workers=3) as executor:
        assets=list(executor.map(build,bundles))
    manifest={"schema_version":1,"repository":OWNER_REPO,"release_tag":TAG,"created_at":created_at,
              "assets":assets,"engineering_datasets":{"total":len(source_rows),
              "redistributed":sum(row["redistributed"] for row in source_rows),
              "source_download_only":[row["dataset"] for row in source_rows if not row["redistributed"]]},
              "notes":["Recorded experiment outputs are snapshots, including historical and pilot runs.",
                       "Current scientific interpretation is paper/research_R2/manuscript.md.",
                       "Third-party data retain source-specific terms; consult dataset_sources.json."]}
    content=json.dumps(manifest,indent=2,ensure_ascii=False)+"\n"
    (ROOT/"data_release/manifest.json").write_text(content,encoding="utf-8")
    (out/"manifest.json").write_text(content,encoding="utf-8")
    sums="".join(a["sha256"]+"  "+a["filename"]+"\n" for a in assets)
    (ROOT/"data_release/SHA256SUMS.txt").write_text(sums,encoding="utf-8")
    (out/"SHA256SUMS.txt").write_text(sums,encoding="utf-8")
    print(json.dumps({"created_at":created_at,"assets":[{k:v for k,v in a.items() if k!="files"} for a in assets],
                      "engineering_datasets":manifest["engineering_datasets"]},indent=2),flush=True)
if __name__=="__main__":main()
