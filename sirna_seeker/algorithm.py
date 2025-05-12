# Importando
#--------------------#
## Tratamento dos dados
#--------------------#
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import json
import time

## Filtros
#--------------------#
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
import pandas as pd
from ViennaRNA import RNA

## Blast
#--------------------#
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML

from concurrent.futures import ThreadPoolExecutor, TimeoutError

import re

#sequence = 'seqence.fasta'

# Algoritmo
#-------------------------------------------------------------------#

#-------------------------------------------------------------------#
#                        Meta data                                  #
#-------------------------------------------------------------------#
def extract_gene_symbol_from_header(header):
    # Expressão regular para capturar siglas de genes entre parênteses.
    match = re.search(r'\(([^)]+)\)', header)
    if match:
        return match.group(1)
    return None

def meta_data(sequence, tamanho, autor,sequence_tag):
  sequencias = list(SeqIO.parse(sequence, 'fasta'))
  # Gene Symbol
  header = sequencias[0].description
  gene_symbol = extract_gene_symbol_from_header(header)
  # RefSeq acenssion number
  refseq_numbers = sequencias[0].name
  # Original gene Size
  gene_size = len(sequencias[0].seq)
  # dic
  meta = {
          #'Form_id': form_id, #enter data
          'gene_symbols': gene_symbol, #run funciton
          'refseq_numbers': refseq_numbers, #run function
          'autor':autor, #enter data
          'gene_size': gene_size, #run function
          'sirna_size':tamanho, #enter data
          'query_title':sequence_tag
  }
  return meta

# ----------------------------------------------------------------- #
#                      Trascribing and selecting                    #
# ----------------------------------------------------------------- #

# Trascreve a sequencia dada para rna
# ----------------------------------------------------------------- #
def transcrever(arquivo):
    sequencias = list(SeqIO.parse(arquivo, 'fasta'))
    for sequencia in sequencias:
        name = sequencia.name
        # Converte as letras T para U
        sequencia_convertida = Bio.Seq.transcribe(sequencia.seq)
    return sequencia_convertida

# Achando todas as possiveis sequencias de siRNA dado == sequencias trascritas
# ----------------------------------------------------------------- #
def possiveis_siRNA(dado, tamanho=21): 
    # definindo variaveis
    possiveis_siRNA = []
    tuplas = []
    # iterando
    for index, _ in enumerate(dado):
        f = index + tamanho
        sequence_p = Seq(dado[index:f])  # convertendo para objeto Seq
        possiveis_siRNA.append(str(sequence_p))
        tuplas.append((str(sequence_p), f'{index}:{index+tamanho}'))

    siRNA = possiveis_siRNA[:-tamanho]

    return siRNA, tuplas

# ----------------------------------------------------------------- #
#                            Free Energy                            #
# ----------------------------------------------------------------- #
from Bio.Seq import Seq

# Parâmetros de energia do modelo nearest-neighbor (valores refinados do usuário)
nearest_neighbor_params = {
    'GC': (-16.52, -42.13),  # GC/CG
    'CG': (-9.61, -23.46),   # CG/GC
    'GG': (-13.94, -34.41),  # CC/GG
    'CC': (-13.94, -34.41),  # CC/GG
    'GA': (-13.75, -36.53),  # GA/CU
    'UC': (-13.75, -36.53),  # GA/CU (reverso complementar)
    'AC': (-11.98, -31.37),  # AC/UG
    'UG': (-11.98, -31.37),  # AC/UG (reverso complementar)
    'CA': (-10.47, -27.08),  # CA/GU
    'GU': (-10.47, -27.08),  # CA/GU (reverso complementar)
    'AG': (-9.34, -23.66),   # AG/UC
    'CU': (-9.34, -23.66),   # AG/UC (reverso complementar)
    'UA': (-9.16, -25.40),   # UA/AU
    'AU': (-8.91, -25.22),   # AU/UA
    'AA': (-7.44, -20.98),   # AA/UU
    'UU': (-7.44, -20.98)    # AA/UU (reverso complementar)
}

# Correções adicionais
INITIATION_ENTHALPY = 4.66     # kcal/mol
INITIATION_ENTROPY = 1.78      # cal/mol·K
SYMMETRY_CORRECTION = -1.38    # cal/mol·K

# Função para verificar se a sequência é auto-complementar
def is_self_complementary(seq):
    complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    rev_comp = ''.join(complement.get(base, base) for base in reversed(seq))
    return seq == rev_comp

# Função principal que recebe fita senso (DNA) e retorna a energia livre da complementar
def free_energy(seq):
    
    def calc_dG(rna_seq, temp_C=37.0):
        seq = rna_seq.upper()
        total_dh = INITIATION_ENTHALPY
        total_ds = INITIATION_ENTROPY

        for i in range(len(seq) - 1):
            pair = seq[i:i+2]
            if pair in nearest_neighbor_params:
                dh, ds = nearest_neighbor_params[pair]
                total_dh += dh
                total_ds += ds
            else:
                raise ValueError(f"Par não encontrado nos parâmetros: {pair}")

        if is_self_complementary(seq):
            total_ds += SYMMETRY_CORRECTION

        temp_K = temp_C + 273.15
        delta_G = total_dh - (temp_K * total_ds / 1000)
        return round(delta_G, 2)


    sense_rna = Seq(seq)

    # Gera a fita antissenso como complementar reversa da fita senso (transcrita para RNA)
    antisense_rna = sense_rna.reverse_complement().transcribe()

    # Calcula energia livre de cada fita
    dG_sense = calc_dG(str(sense_rna[:7]))
    dG_antisense = calc_dG(str(antisense_rna[:7]))

    # Calcula diferença ΔG_senso - ΔG_antissenso e força resultado a ser negativo
    dG_diff = round(dG_sense - dG_antisense, 2)
    dG_diff = -abs(dG_diff)

    # Retorna ΔG da fita complementar e diferença
    return dG_antisense, dG_diff

# ----------------------------------------------------------------- #
#                             Classifing
# ----------------------------------------------------------------- #

# Reynolds 
# ----------------------------------------------------------------- #
def reynolds (sequence):
  score = 0
  falha = []
  # Instabilidade na position 15 - 19
  #----------------------------------#
  check_estabilidade = 0
  for letra in sequence[14:19]:
     if letra == "A" or letra == "G":
        check_estabilidade += 1
  if check_estabilidade >= 1:
     score += 2
  else:
     falha.append(str("internal stability"))

  # Caracteristicas da Fita
  # position 13
  #----------------------------------#
  if sequence[12] != "G":
     score += 1
  else:
     falha.append(str("position 13"))
  # position 19
  #----------------------------------#
  if sequence [18] != "G" and sequence [18] != "C":
     score += 1
  else:
     falha.append(str("position 19"))
  if sequence [18] == "A":
     score += 1
  else:
     falha.append(str("position 19!= A"))
  # position 10
  #----------------------------------#
  if sequence [9] == "U":
     score += 1
  else:
     falha.append(str("position 10"))

  # position 3
  #----------------------------------#
  if sequence [2] == "A":
     score += 1
  else:
     falha.append(str("position 3"))

  return score, falha

# Ui-Tei
# ----------------------------------------------------------------- #
def Ui_Tei (sequence, extremidade=7):
  # Variaveis
  #----------------------------------#
  score = 0
  falha = []
  antisenso = Seq(sequence).complement_rna()

  #Senso
  #----------------------------------#
  check_senso = 0

  for letra in sequence[:extremidade]:
      if letra == "G" or letra == "C":
         check_senso += 1
  if check_senso > 0:
     score += 1
  else:
     falha.append(str("no C/G at the 5' sense end"))

  #Antisenso
  #----------------------------------#
  check_antisenso = 0

  for letra in antisenso[:extremidade]:
      if letra == "A" or letra == "U":
         check_antisenso += 1
  if check_antisenso > 0:
     score += 1
  else:
     falha.append(str("no A/U at the 5' antisense end"))

  # Presença de no minimo 5 A/U nas posições [:7]
  #----------------------------------#
  check_min_5 = 0

  for letra in antisenso[-extremidade:]:
      if letra == "A" or letra == "U":
         check_min_5 += 1
  if check_min_5 >= 5:
     score += 2
  else:
     falha.append(str("Less than 5pb A/U in antisense 5'"))

  # Mais de 9 repetições
  #----------------------------------#
  segment_size = 9
  check = True

  for index in range(len(sequence) - segment_size + 1):
        segmento = sequence[index:index + segment_size]
        count = 0
        for letra in segmento:
            if letra == 'C' or letra == 'G':
                count += 1
        if count >= segment_size:
            check = False

  if check == True:
    score += 1
  else:
    falha.append(str("Repetition of more than 9 C/G in a row"))
    score -= 2

  return score, falha

# Amarzguioui
# ----------------------------------------------------------------- #
def Amarzguioui(sequence):
  score = 0
  falha = []
  check_assimetria_5 = 0
  check_assimetria_3 = 0
  antisenso = Seq(sequence).complement_rna()

  # position 1
  #----------------------------------#
  if sequence [0] != "U":
     score += 1
  else:
     falha.append(str("position 1"))
  
  if sequence [0] == "C" and sequence [0] == "G":
     score += 1
  else:
     falha.append(str("position 1 has no C/G"))

  # position 6
  #----------------------------------#
  if sequence [5] == "A":
     score += 1
  else:
     falha.append(str("position 6"))
  # position 19
  #----------------------------------#
  if sequence [18] != "G":
     score += 1
  else:
     falha.append(str("position 19 is G"))
  
  if sequence [18] == "A" or sequence [18] == "U":
     score += 1
  else:
     falha.append(str("position 19 is C/G"))
  # Assimetria
  #----------------------------------#
  for letra in sequence[:3]:
      if letra == "A" or letra == "U": # baixo
         check_assimetria_5 += 1
  for letra in antisenso[-3:]:
      if letra == "A" or letra == "U": # alto
         check_assimetria_3 += 1
  if check_assimetria_3 > check_assimetria_5:
     score += 1
  else:
     falha.append(str("asymmetry"))

  return score, falha

# Testando a funcionalidade de um siRNA
# ----------------------------------------------------------------- #
def siRNA_score (sequence, tuplas,
                 autor=['reynolds', 'ui-tei', 'amarzguioui'],
                 tm=True, tmmax = 21.5):
        """
        Avalia o score e a conformidade do siRNA com base em diversos parâmetros.

        :param sequence: Sequência do siRNA para avaliação.
        :param autor: Lista de autores para o método de scoring.
        :param tm: Flag para calcular a temperatura de fusão (melt temperature).
        :param tmmax: Valor máximo para a temperatura de fusão.
        :param tuplas: Lista de tuplas contendo siRNA e suas posições.
        :return: Score, lista de falhas, conteúdo GC, energia livre e temperatura de fusão.
        """
        # Definindo variaveis
        #--------------------------------------------------------------#
        score = 0
        falha = []
        posicao = 'Não encontrado'

        # Conteudo Baixo GC
        #--------------------------------------------------------------#
        conteudo_gc = round(gc_fraction(sequence)*100, 2)
        if 30 <= conteudo_gc <= 52:
            score += 1

        # tempmelt
        #--------------------------------------------------------------#
        tm_score = round(mt.Tm_GC(sequence[1:9]), 2)
        if tm == True:
          #if autor != 'reynolds':
              if 0 < tm_score <= tmmax:
                  score += 1
              else:
                  falha.append(str("tm"))

        # G°
        #--------------------------------------------------------------#
        energia_livre = free_energy(sequence)
        dG_antisense, dG_diff = energia_livre
        #if autor != 'ui-tei':
        if -13 < dG_antisense < -7:
                score += 2
        else:
                falha.append(str("free energy"))

        # Autores
        # Score reynolds total = 9
        #--------------------------------------------------------------#
        if autor == 'reynolds':
            r = reynolds(sequence)
            score += r[0]
            falha.extend(r[1] if isinstance(r[1], list) else [r[1]])

        # Score ui-tei total = 7
        #--------------------------------------------------------------#
        if autor == 'ui-tei':
            u = Ui_Tei(sequence)
            score += u[0]
            falha.extend(u[1] if isinstance(u[1], list) else [u[1]])

        # Score amarzguioui total = 4
        #--------------------------------------------------------------#
        if autor == 'amarzguioui':
            a = Amarzguioui(sequence)
            score += a[0]
            falha.extend(a[1] if isinstance(a[1], list) else [a[1]])

        # Encontrando a posição do siRNA nas tuplas
        #--------------------------------------------------------------#
        if tuplas:
            for seq, pos in tuplas:
                if sequence == seq:
                    posicao = pos
                    break
        #--------------------------------------------------------------#
        return score, falha, conteudo_gc, dG_diff, tm_score, posicao
# ----------------------------------------------------------------- #
# Selecionando os siRNA funcionais
# ----------------------------------------------------------------- #
def filtro_siRNA(sequences, tuplas, threshold=0.6,
                    autor=['reynolds', 'ui-tei', 'amarzguioui'],
                    tm=True, tmmax=21.5):
        """
        Filtra e avalia siRNAs com base em parâmetros de conformidade e scoring.

        :param sequences: Lista de sequências de siRNA para avaliação.
        :param conformidade: Nível de conformidade para filtrar siRNAs.
        :param autor: Lista de autores para métodos de scoring.
        :param tm: Flag para calcular a temperatura de fusão (melt temperature).
        :param tmmax: Valor máximo para a temperatura de fusão.
        :param tuplas: Lista de tuplas contendo siRNA e suas posições.
        :return: Dados JSON e lista de siRNAs verificados.
        """
        # Definindo Variaveis
        # --------------------------------------------------------------#
        siRNA_verificados = []
        posicao = []
        score = []
        falha = []
        conteudo_gc = []
        energia_livre = []
        TM_score = []
        
        if tm == True:
          total_reynolds = 11 * threshold
          total_uitei = 9 * threshold
          total_ama = 11 * threshold
        else:
          total_reynolds = 10 * threshold
          total_uitei = 8 * threshold
          total_ama = 10 * threshold


        # Iterador
        # --------------------------------------------------------------#
        for index, sequence in enumerate(sequences):
            # Veficando a qualidade
            # ----------------------------------#
            resultado = siRNA_score(sequence=sequence, tuplas=tuplas, autor=autor,
                                    tm=tm, tmmax=tmmax)

            # Exclui RNAs indesejadas
            # ----------------------------------#
            incluir_siRNA = True

            if 30 <= resultado[2] <= 52:
                incluir_siRNA = False

            if autor == "reynolds":
                if resultado[0] <= total_reynolds: #or not (0 < resultado[4] <= tmmax):
                    incluir_siRNA = False

            if autor == 'ui-tei':
              if resultado[0] <= total_uitei: #or not (-13 < resultado[3] <-7):
                    incluir_siRNA = False

            if autor == 'amarzguioui':
                if resultado[0] <= total_ama:
                    incluir_siRNA = False

            if incluir_siRNA:
                siRNA_verificados.append(sequence)
                score.append(resultado[0])
                falha.append(resultado[1])
                conteudo_gc.append(resultado[2])
                energia_livre.append(resultado[3])
                TM_score.append(resultado[4])
                posicao.append(resultado[5])

        # Organizando os resultados por pontuação
        resultados = list(zip(siRNA_verificados, score, TM_score, conteudo_gc, energia_livre, falha, posicao))
        resultados_ordenados = sorted(resultados, key=lambda x: x[1], reverse=True)[:50]
        siRNA_verificados = [resultado[0] for resultado in resultados_ordenados]
        score = [resultado[1] for resultado in resultados_ordenados]
        TM_score = [resultado[2] for resultado in resultados_ordenados]
        conteudo_gc = [resultado[3] for resultado in resultados_ordenados]
        energia_livre = [resultado[4] for resultado in resultados_ordenados]
        falha = [resultado[5] for resultado in resultados_ordenados]
        posicao = [resultado[6] for resultado in resultados_ordenados]

        #"score, falha, conteudo_gc, energia_livre, tm_score, posicao"

        # Preparando dados para JSON
        dados_json = [
            {
                "Sequencia'": resultado[0],
                "Pontuacao": resultado[1],
                "TM": resultado[2],
                "CG": resultado[3],
                "AG": resultado[4],
                "Falhas": resultado[5],
                "posicao":resultado[6],
            }
            for resultado in resultados_ordenados
        ]

        # Salvando em um arquivo JSON
        with open('resultados_siRNA.json', 'w') as arquivo_json:
            json.dump(dados_json, arquivo_json, indent=4)

        return dados_json, siRNA_verificados

# ----------------------------------------------------------------- #
#                              Blast                                #
# ----------------------------------------------------------------- #

# Fasta com as sequencias funcionais
# ---------------------------------------------- #    
def guardando_sequence (data):
        records = []
        for index, sequence in enumerate(data):
            # Criando o cabeçalho enumerado
            header = f"sequencia_{index+1}"

            # Criando um objeto SeqRecord com a sequência e o cabeçalho
            record = SeqRecord(Seq(sequence), id=header, description="")

            # Adicionando o registro à lista
            records.append(record)

        # Salvando os registros no arquivo FASTA
        fasta_filename = "minha_sequencia.fasta"
        with open(fasta_filename, "w") as output_file:
            SeqIO.write(records, output_file, "fasta")
        return fasta_filename

# Rodando o Blast
# ---------------------------------------------- #    
def blast_siRNA(sequence, sequence_tag, db=None, organismo="txid9606[ORGN]",
                df=None, identidad=0.78):
    """
    Executa o BLAST de uma sequência de siRNA contra um banco de dados especificado.

    Parâmetros:
    -----------
    sequence : str
        A sequência de siRNA a ser analisada.
    sequence_tag : str
        Identificador para a sequência.
    db : str, opcional
        Banco de dados BLAST a ser usado. Padrão é None.
    organismo : str, opcional
        Filtro taxonômico para o BLAST. Padrão é "txid9606[ORGN]" (Homo sapiens).
    df : str, opcional
        Programa BLAST a ser usado, como "blastn". Padrão é None.
    identidad : float, opcional
        Identidade mínima para considerar um alinhamento. Padrão é 0.78 (78%).

    Retorna:
    --------
    list of dict
        Uma lista de dicionários contendo os resultados do BLAST com os campos:
        - "description": Identificador do hit (nome do alinhamento).
        - "identity": Porcentagem de identidade do alinhamento.
        - "coverage": Cobertura da sequência alinhada.

    Exceções:
    ---------
    ValueError:
        Lançada se a identidade ou a sequência forem inválidas.

    """
    
    resultados = []  # Lista para armazenar os resultados

    try:
        print("Iniciando a requisição BLAST...")
        result_handle = NCBIWWW.qblast(program=df, database=db,
                                       sequence=sequence,
                                       entrez_query=organismo,
                                       perc_ident=identidad)
        print("Analisando os resultados...")
        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                hit = alignment.hsps[0]

                # Calcular porcentagem de identidade e cobertura
                p_identity = hit.identities / hit.align_length
                coverage = hit.align_length / blast_record.query_length

                # Adicionar o resultado como dicionário
                resultados.append({
                    "description": alignment.hit_id[18:],  # Nome do hit
                    "identity": p_identity,  # Porcentagem de identidade
                    "coverage": coverage  # Cobertura da sequência
                })
                print(resultados)
        
                time.sleep(20)

    except Exception as e:
        print(f"Erro durante o processamento do BLAST: {e}")
        raise

    return resultados

def identidade_siRNA(fasta_file, sequence_tag, db = "refseq_rna", organismo = "txid9606[ORGN]",
                 df=None, identidade=0.78):
  # Ler o arquivo fasta
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

  # Realizar o BLAST para cada sequência
    for i, sequence in enumerate(sequences, start=1):
        print(f"Processando sequência {i}...")

        # Blast
        resultados = blast_siRNA (sequence , sequence_tag, 
                                               db, organismo,
                                                df, identidade)

    return resultados

