
import os
base = r"F:\nano-everything\mini-pure-physics\8. mini-plasma-physics\mini-industrial-plasma"
path = os.path.join(base, 'src', 'plasma_chemistry.c')
append_content = open(os.path.join(base, '_append.txt'), 'r').read()
with open(path, 'a') as f:
    f.write(append_content)
print('Appended', len(append_content), 'bytes')
