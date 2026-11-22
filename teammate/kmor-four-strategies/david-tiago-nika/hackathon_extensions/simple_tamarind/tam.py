#!/usr/bin/env python3
"""Minimal Tamarind Bio client. Standard library only.

This script contains NO job values. Everything that reaches Tamarind comes from ONE
user-defined JSON file (the "params file"):

{
  "base_url": "https://app.tamarind.bio/api",   # Tamarind API root
  "key_file": ".tamarind_key",                  # where the API key is stored
  "type":     "rfd4-proteina",                  # tool name
  "jobName":  "my_job",                         # name of the job
  "files":    {"pdbFile": "path/to/input.cif"}, # setting name -> local file to upload
  "settings": { ... }                           # sent to Tamarind exactly as written
}

Relative paths (key_file and the paths in "files") are resolved against the folder
that contains the params file.

For every entry in "files", the local file is uploaded and settings[<name>] is set
to the uploaded file name. That is the only thing the script adds to "settings".

Commands:
  setup-key  --params P   save your Tamarind API key into key_file (hidden prompt)
  schema     --params P   print the real settings schema of the tool in "type"
  validate   --params P   upload files, then POST /validate-job (no cost, no job)
  run        --params P   validate; ONLY if valid, POST /submit-job with the
                          identical payload (this one costs credits)
"""
import argparse
import getpass
import hashlib
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

REQUIRED_KEYS = ("base_url", "key_file", "type", "jobName", "settings")


# ----------------------------------------------------------------------------- params
def load_params(params_path):
    path = Path(params_path).resolve()
    with open(path, "r", encoding="utf-8") as f:
        params = json.load(f)
    missing = [k for k in REQUIRED_KEYS if k not in params]
    if missing:
        sys.exit(f"ERROR: params file {path} is missing required key(s): {missing}")
    params.setdefault("files", {})
    params["_dir"] = path.parent
    return params


def resolve(params, rel):
    p = Path(rel)
    return p if p.is_absolute() else params["_dir"] / p


def read_key(params):
    kf = resolve(params, params["key_file"])
    if not kf.is_file():
        sys.exit(f"ERROR: no API key at {kf}. Run: setup-key --params <file>")
    key = kf.read_text(encoding="utf-8").strip()
    if not key:
        sys.exit(f"ERROR: {kf} is empty. Run: setup-key --params <file>")
    return key


# ------------------------------------------------------------------------------- http
class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *a, **k):
        return None


_opener = urllib.request.build_opener(_NoRedirect)


def http(method, url, key, body=None, content_type=None, timeout=60):
    """Return (status, parsed_body, headers). Never follows redirects."""
    headers = {"x-api-key": key}
    data = None
    if body is not None:
        data = body
        headers["Content-Type"] = content_type
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with _opener.open(req, timeout=timeout) as r:
            status, raw, hdrs = r.status, r.read(), r.headers
    except urllib.error.HTTPError as e:
        status, raw, hdrs = e.code, e.read(), e.headers
    text = raw.decode("utf-8", errors="replace")
    try:
        return status, json.loads(text), hdrs
    except ValueError:
        return status, text, hdrs


def api(params):
    return params["base_url"].rstrip("/")


# ----------------------------------------------------------------------------- upload
def upload_files(params, key):
    """Upload each local file in params['files']; return {setting_name: uploaded_name}."""
    uploaded = {}
    for setting_name, local in params["files"].items():
        lp = resolve(params, local)
        if not lp.is_file():
            sys.exit(f"ERROR: file for '{setting_name}' not found: {lp}")
        name = lp.name
        status, body, hdrs = http(
            "PUT", f"{api(params)}/upload/{name}", key,
            body=lp.read_bytes(), content_type="application/octet-stream",
        )
        if status in (200, 201):
            print(f"[upload] {name}: uploaded")
        elif status in (301, 302, 307, 308):
            # Observed behaviour: when the name already exists in the account, Tamarind
            # answers with a redirect to its CDN. The existing file is NOT replaced.
            s2, listing, _ = http("GET", f"{api(params)}/files", key)
            names = listing if isinstance(listing, list) else (listing.get("files", []) if isinstance(listing, dict) else [])
            names = [str(n.get("name", n)) if isinstance(n, dict) else str(n) for n in names]
            if name in names:
                print(f"[upload] {name}: ALREADY in your account (not replaced). "
                      f"If the local content changed, rename the file or delete the remote one first.")
            else:
                sys.exit(f"ERROR: upload of {name} was redirected and the file is not in the account.")
        else:
            sys.exit(f"ERROR: upload of {name} failed: HTTP {status}: {body}")
        uploaded[setting_name] = name
    return uploaded


# ---------------------------------------------------------------------------- payload
def build_payload(params, uploaded):
    settings = dict(params["settings"])
    for setting_name, name in uploaded.items():
        settings[setting_name] = name
    return {"jobName": params["jobName"], "type": params["type"], "settings": settings}


def fingerprint(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def validate(params, key, payload):
    body = {"type": payload["type"], "settings": payload["settings"], "jobName": payload["jobName"]}
    status, resp, _ = http("POST", f"{api(params)}/validate-job", key,
                           body=json.dumps(body).encode("utf-8"), content_type="application/json")
    if status != 200 or not isinstance(resp, dict):
        sys.exit(f"ERROR: validate-job HTTP {status}: {resp}")
    return resp


def print_verdict(resp):
    print(f"[validate] valid = {resp.get('valid')}")
    if not resp.get("valid"):
        print(f"[validate] error: {resp.get('error')}")
        for e in resp.get("errors", []) or []:
            print(f"  - {e.get('field')}: {e.get('message')}")
        if resp.get("missing_fields"):
            print(f"  missing: {resp['missing_fields']}")
    else:
        print("[validate] normalized payload (your settings + defaults Tamarind fills in):")
        sent = resp.get("normalized") or {}
        print(json.dumps(sent, indent=2))


# --------------------------------------------------------------------------- commands
def cmd_setup_key(a):
    params = load_params(a.params)
    key = getpass.getpass("Tamarind API key (hidden): ") if sys.stdin.isatty() else sys.stdin.readline().strip()
    if not key:
        sys.exit("ERROR: empty key")
    status, body, _ = http("GET", f"{api(params)}/files", key)
    if status != 200:
        sys.exit(f"ERROR: Tamarind rejected this key (HTTP {status}). Nothing saved.")
    kf = resolve(params, params["key_file"])
    kf.write_text(key + "\n", encoding="utf-8")
    try:
        os.chmod(kf, 0o600)
    except OSError:
        pass
    print(f"Key verified against Tamarind and saved to {kf}")


def cmd_schema(a):
    params = load_params(a.params)
    key = read_key(params)
    status, tools, _ = http("GET", f"{api(params)}/tools", key)
    if status != 200:
        sys.exit(f"ERROR: GET /tools HTTP {status}")
    tool = next((t for t in tools if t.get("name") == params["type"]), None)
    if tool is None:
        sys.exit(f"ERROR: tool '{params['type']}' not found. Available: {[t['name'] for t in tools]}")
    print(json.dumps(tool, indent=2))


def cmd_validate(a):
    params = load_params(a.params)
    key = read_key(params)
    payload = build_payload(params, upload_files(params, key))
    print(f"[payload] fingerprint {fingerprint(payload)}")
    print(json.dumps(payload, indent=2))
    resp = validate(params, key, payload)
    print_verdict(resp)
    sys.exit(0 if resp.get("valid") else 1)


def cmd_run(a):
    params = load_params(a.params)
    key = read_key(params)
    payload = build_payload(params, upload_files(params, key))
    fp = fingerprint(payload)
    print(f"[payload] fingerprint {fp}")
    resp = validate(params, key, payload)
    print_verdict(resp)
    if not resp.get("valid"):
        sys.exit("NOT SUBMITTED: validation failed.")
    if fingerprint(payload) != fp:  # paranoia: payload must be byte-identical to the validated one
        sys.exit("NOT SUBMITTED: payload changed after validation.")
    status, body, _ = http("POST", f"{api(params)}/submit-job", key,
                           body=json.dumps(payload).encode("utf-8"), content_type="application/json")
    print(f"[submit] HTTP {status}")
    print(body if isinstance(body, str) else json.dumps(body, indent=2))
    sys.exit(0 if status == 200 else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("setup-key", cmd_setup_key), ("schema", cmd_schema),
                     ("validate", cmd_validate), ("run", cmd_run)):
        p = sub.add_parser(name)
        p.add_argument("--params", required=True, help="path to the user-defined params JSON")
        p.set_defaults(fn=fn)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
