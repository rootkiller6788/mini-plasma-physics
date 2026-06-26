
import os
BASE = r"F:
ano-everything\mini-pure-physics\8. mini-plasma-physics\mini-industrial-plasma"

def w(relpath, content):
    path = os.path.join(BASE, relpath)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    lines = content.count(chr(10)) + 1
    print(f"  OK {relpath}: {lines} lines, {len(content)} bytes")

print("Building remaining files...")
