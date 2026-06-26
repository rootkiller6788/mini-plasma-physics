import os

base = r"F:\nano-everything\mini-pure-physics\8. mini-plasma-physics\mini-laser-plasma"

def write_file(rel_path, content):
    full = os.path.join(base, rel_path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  Created {rel_path} ({len(content)} bytes)")

print("Building mini-laser-plasma module...")
print(f"Base: {base}")
