#!/usr/bin/env python3
"""./tam.py RUN_FILE validate|run
settings.json: base_url, api_key.  RUN_FILE: type, jobName, input (task json), settings (extra hyperparameters).
Everything about the design itself is read from the input task json. Responses are printed as Tamarind sends them."""
import json, sys, urllib.error, urllib.request
from pathlib import Path

if len(sys.argv) != 3 or sys.argv[2] not in ("validate", "run"):
    sys.exit(__doc__)
here = Path(__file__).resolve().parent
cfg = json.loads((here / "settings.json").read_text())
run_path = Path(sys.argv[1]).resolve()
run = json.loads(run_path.read_text())
src = (run_path.parent / run["input"]).resolve()
task = next(iter(json.loads(src.read_text()).values()))  # first entry of the input json
cif = src.parent / Path(task["input"]).name              # structure named by the input json, sitting next to it


class Follow308(urllib.request.HTTPRedirectHandler):  # urllib will not re-send a PUT body on 307/308 (curl -L does)
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return urllib.request.Request(newurl, data=req.data, headers=req.headers, method=req.get_method())


def call(method, route, body, ctype):
    req = urllib.request.Request(cfg["base_url"] + route, data=body, method=method,
                                 headers={"x-api-key": cfg["api_key"], "Content-Type": ctype})
    try:
        with urllib.request.build_opener(Follow308).open(req) as r:
            return r.status, r.read().decode()
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode()


if run["type"] == "rfdiffusion3":  # Tamarind's "json" task takes this file as is
    uploads, derived = [src, cif], {"jsonFile": src.name, "inputFiles": [cif.name]}
elif run["type"] == "rfd4-proteina":
    lo, hi = task["length"].split("-")
    chains = {}
    for res in task["unindex"].split(","):
        chains.setdefault(res[0], []).append(res[1:])
    uploads, derived = [cif], {
        "pdbFile": cif.name, "ligandCode": task["ligand"],
        "minScaffoldLength": int(lo), "maxScaffoldLength": int(hi),
        "motifResidues": {c: ",".join(n) for c, n in chains.items()},
        "catalyticAtomPicks": " ".join(f"{k}:{v}" for k, v in task["select_fixed_atoms"].items())}
else:
    sys.exit(f"unsupported type {run['type']}")

for f in uploads:
    print("upload", f.name, *call("PUT", f"/upload/{f.name}", f.read_bytes(), "application/octet-stream"))
payload = json.dumps({"jobName": run["jobName"], "type": run["type"],
                      "settings": {**run["settings"], **derived}}).encode()
code, text = call("POST", "/validate-job", payload, "application/json")
print("validate", code, text)
if sys.argv[2] == "run":
    if code != 200 or not json.loads(text).get("valid"):
        sys.exit("NOT SUBMITTED")
    print("submit", *call("POST", "/submit-job", payload, "application/json"))
