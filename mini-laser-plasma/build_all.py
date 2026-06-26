import os
base = r"F:
ano-everything\mini-pure-physics\8. mini-plasma-physics\mini-laser-plasma"
def W(rel, content):
    full = os.path.join(base, rel)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  {rel}: {content.count(chr(10))} lines")
print("Building mini-laser-plasma...")

