#!/usr/bin/env python3
import os, math
BASE = os.path.dirname(os.path.abspath(__file__))
total_lines = 0
def w(rel, content):
    global total_lines
    path = os.path.join(BASE, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    n = content.count(chr(10))
    total_lines += n
    print(f"  {rel}: {n} lines")
    return n
