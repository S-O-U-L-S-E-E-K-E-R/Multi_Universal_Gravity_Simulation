"""Look up sky coordinates for SPARC galaxies with the CDS Sesame name resolver.
Usage: python3 analysis/resolve_coords.py <pergal.json> <out_coords.json>"""
import json, re, subprocess, sys, time, urllib.parse
names = [x["name"] for x in json.load(open(sys.argv[1]))]
def variants(n):
    v = [n]
    m = re.match(r"^([A-Za-z]+)0*(\d+.*)$", n)
    if m: v.append(f"{m.group(1)} {m.group(2)}")
    m = re.match(r"^ESO(\d+)-G(\d+)$", n)
    if m: v += [f"ESO {int(m.group(1))}-{int(m.group(2))}", f"ESO {m.group(1)}-G{m.group(2)}"]
    m = re.match(r"^F(\d+)-(\w+)$", n)
    if m: v += [f"LSBC F{m.group(1)}-{m.group(2).zfill(2)}", f"[BKM96] F{m.group(1)}-{m.group(2)}"]
    m = re.match(r"^UGCA(\d+)$", n)
    if m: v.append(f"UGCA {int(m.group(1))}")
    m = re.match(r"^KK98-(\d+)$", n)
    if m: v.append(f"[KK98] {m.group(1)}")
    return v
out = {}
for n in names:
    got = None
    for v in variants(n):
        url = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/SNV?" + urllib.parse.quote(v)
        try:
            txt = subprocess.run(["curl", "-sS", "--max-time", "20", url], capture_output=True, text=True).stdout
        except Exception:
            txt = ""
        m = re.search(r"^%J ([\d.]+) ([+-][\d.]+)", txt, re.M)
        if m: got = (float(m.group(1)), float(m.group(2)), v); break
        time.sleep(0.2)
    out[n] = got
    print(n, got, flush=True)
json.dump(out, open(sys.argv[2], "w"))
print("resolved", sum(v is not None for v in out.values()), "of", len(out))
