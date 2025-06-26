import subprocess
import sys
import os
from Bio.Seq import Seq

def calcular_energia(seq, base, nearest_neighbor_energy = "nearest_neighbor_energy.txt", dangling_3_energy = "dangling_3_energy.txt"):
    # Caminho absoluto para o perl.exe portátil
    perl_exe = r"C:/strawberry-perl/perl/bin/perl.exe"
    script_dir = os.path.dirname(__file__)
    script_path = os.path.join(script_dir, "run_thermo.pl")

    result = subprocess.run(
        [perl_exe, "-I.", script_path, seq, base, nearest_neighbor_energy, dangling_3_energy],
        cwd=script_dir,  # <- garante que o script rode no mesmo diretório
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print("Erro ao executar o script Perl:")
        print(result.stderr)
        sys.exit(1)
    energia = result.stdout.strip()
    return float(energia)


def gerar_complementar_reversa(seq):
    seq_obj = Seq(seq)
    return str(seq_obj.reverse_complement())


def calcular_diferenca_energia(seq_sense):
    """Calcula a diferença de energia livre entre sense e antisense para uma sequência dada."""
    if len(seq_sense) < 6:
        raise ValueError("Sequência muito curta para extrair base complementar.")

    base = seq_sense[5]  # Ajuste conforme seu critério para a base dangling
    energia_sense = calcular_energia(seq_sense, base)
    if energia_sense is None:
        raise RuntimeError("Erro ao calcular energia da sequência sense.")

    seq_antisense = gerar_complementar_reversa(seq_sense)
    energia_antisense = calcular_energia(seq_antisense, base)
    if energia_antisense is None:
        raise RuntimeError("Erro ao calcular energia da sequência antisense.")

    diff = round((energia_sense - energia_antisense), 2)
    return energia_antisense, diff

if __name__ == "__main__":
    seq_exemplo = "AGCTGAC"
    try:
        energia_antisense, diff = calcular_diferenca_energia(seq_exemplo)
        print(f"Diferença de energia livre (sense - antisense): {diff}")
    except Exception as e:
        print(f"Erro: {e}")