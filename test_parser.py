import re

def parse_vina_log(content):
    # Patrón robusto para tabla de Vina 1.1.x y 1.2.x
    # Busca una línea que empiece con espacios, luego '1', luego espacios, luego el número flotante
    # Ejemplo: "   1       -8.903          0          0"
    match = re.search(r'^\s*1\s+([-+]?\d+\.\d+)', content, re.MULTILINE)
    if match:
        return float(match.group(1))
    
    # Fallback: buscar "Estimated Free Energy of Binding" (algunas versiones)
    match = re.search(r'Estimated Free Energy of Binding\s*:\s*([-+]?\d+\.\d+)', content)
    if match:
        return float(match.group(1))

    return None

log_text = """
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -8.903          0          0
   2       -8.819      2.647      9.288
"""

print(f"Resultado: {parse_vina_log(log_text)}")
