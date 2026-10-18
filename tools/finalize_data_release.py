"""Finalize integrity metadata from already-built, frozen ZIP snapshots."""
import argparse
import hashlib
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG="data-2026-10-18"
def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):h.update(block)
    return h.hexdigest()
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    out=args.output.resolve()
    if not out.is_relative_to(ROOT/"tmp"):raise SystemExit("Expected ignored tmp publication directory")
    assets=[]
    for name in ["experiment-records","engineering-data","pysr-search"]:
        path=out/(name+".zip")
        with zipfile.ZipFile(path) as archive:
            info=json.loads(archive.read("DATA_BUNDLE_INFO.json"))
            records=info["files"]
            expected={r["path"]:r for r in records}
            if len(expected)!=len(records):raise SystemExit("Duplicate member")
            if set(archive.namelist())!=set(expected)|{"DATA_BUNDLE_INFO.json"}:raise SystemExit("Member set differs")
            for index,row in enumerate(records):
                h=hashlib.sha256();size=0
                with archive.open(row["path"]) as stream:
                    for block in iter(lambda:stream.read(1024*1024),b""):
                        h.update(block);size+=len(block)
                if size!=row["size"] or h.hexdigest()!=row["sha256"]:raise SystemExit("Member SHA-256 mismatch: "+row["path"])
                if index%1000==0:print(f"{name}: verified {index+1}/{len(records)} members",flush=True)
            assets.append({"name":name,"filename":path.name,"size":path.stat().st_size,
              "sha256":digest(path),"archive_created_at":info["created_at"],
              "url":f"https://github.com/leokynlang/form_check/releases/download/{TAG}/{path.name}",
              "file_count":len(records),"uncompressed_size":sum(r["size"] for r in records),"files":records})
        print(name+": full member and archive hashes verified",flush=True)
    source_rows=json.loads((ROOT/"data_release/dataset_sources.json").read_text(encoding="utf-8"))
    manifest={"schema_version":1,"repository":"leokynlang/form_check","release_tag":TAG,
      "created_at":datetime.now(timezone.utc).isoformat(),"assets":assets,
      "engineering_datasets":{"total":len(source_rows),"redistributed":sum(r["redistributed"] for r in source_rows),
      "source_download_only":[r["dataset"] for r in source_rows if not r["redistributed"]]},
      "notes":["Frozen recorded snapshots; each archive preserves its own creation time.",
               "Historical and pilot records are not all part of the main manuscript's analysis pool.",
               "Third-party redistribution follows data_release/THIRD_PARTY_DATA.md."]}
    value=json.dumps(manifest,indent=2,ensure_ascii=False)+"\n"
    for target in [ROOT/"data_release/manifest.json",out/"manifest.json"]:target.write_text(value,encoding="utf-8")
    sums="".join(a["sha256"]+"  "+a["filename"]+"\n" for a in assets)
    for target in [ROOT/"data_release/SHA256SUMS.txt",out/"SHA256SUMS.txt"]:target.write_text(sums,encoding="utf-8")
    print(json.dumps({k:v for k,v in manifest.items() if k!="assets"},indent=2),flush=True)
if __name__=="__main__":main()
