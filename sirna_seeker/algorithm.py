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
def free_energy(seq1):
    """
    Calculates the Minimum Free Energy (MFE) and structure of an RNA duplex
    using the complementary sequence of the input RNA sequence.

    Parameters:
        seq1 (str): The input RNA sequence.

    Returns:
        tuple: A tuple containing the MFE (float) and the dot-bracket structure (str).
    """
    # Create a Seq object
    seq1_obj = Seq(seq1)

    # Generate the complementary sequence and convert it to a string
    complementary_seq_str = str(seq1_obj.complement())

    # Perform RNA duplex folding
    result = RNA.duplexfold(seq1, complementary_seq_str)

    # Extract the free energy and structure
    mfe = round(result.energy,3)  # Free energy in kcal/mol

    return mfe

# ----------------------------------------------------------------- #
#                             Classifing
# ----------------------------------------------------------------- #

# Reynolds 
# ----------------------------------------------------------------- #
def reynolds (sequence):
    score = 0
    falha = []
    # Instabilidade na posicao 15 - 19
    #----------------------------------#
    check_estabilidade = 0
    for letra in sequence[14:19]:
        if letra == "A" or letra == "G":
            check_estabilidade += 1
    if check_estabilidade >= 1:
        score += 1
    else:
        falha.append(str("estabilidade interna"))
    # posicao 13
    #----------------------------------#
    if sequence[12] != "G":
        score += 1
    else:
        falha.append(str("posicao 13"))
    # posicao 19
    #----------------------------------#
    if sequence [18] != "G" and sequence [18] != "C":
        score += 1
    else:
        falha.append(str("posicao 19"))
    if sequence [18] == "A":
        score += 1
    else:
        falha.append(str("posicao 19!= A"))
    # posicao 10
    #----------------------------------#
    if sequence [9] == "U":
        score += 1
    else:
        falha.append(str("posicao 10"))

    # posicao 3
    #----------------------------------#
    if sequence [2] == "A":
        score += 1
    else:
        falha.append(str("posicao 3"))

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
    if check_senso == 5:
        score += 1
    else:
        falha.append(str("Extremidade senso"))

    #Antisenso
    #----------------------------------#
    check_antisenso = 0

    for letra in antisenso[-extremidade:]:
        if letra == "A" or letra == "U":
            check_antisenso += 1
    if check_antisenso == extremidade:
        score += 1
    else:
        falha.append(str("Extremidade antisenso"))

    # Presença de no minimo 5 A/U nas posições [:7]
    #----------------------------------#
    check_min_5 = 0

    for letra in antisenso[-7:]:
        if letra == "A" or letra == "U":
            check_min_5 += 1
    if check_min_5 >= 5:
        score += 1
    else:
        falha.append(str("posições [:7]"))

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
        falha.append(str("Repetição de mais de 9 C/G seguidos"))

    return score, falha

# Amarzguioui
# ----------------------------------------------------------------- #
def Amarzguioui(sequence):
    score = 0
    falha = []
    check_assimetria_5 = 0
    check_assimetria_3 = 0
    antisenso = Seq(sequence).complement_rna()

    # posicao 1
    #----------------------------------#
    if sequence [0] != "U" and sequence [0] != "A":
        score += 1
    else:
        falha.append(str("posicao 1"))
    # posicao 6
    #----------------------------------#
    if sequence [5] == "A":
        score += 1
    else:
        falha.append(str("posicao 6"))
    # posicao 19
    #----------------------------------#
    if sequence [18] != "G" and sequence [18] != "C":
        score += 1
    else:
        falha.append(str("posicao 19"))
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
        falha.append(str("assimetria"))

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
        else:
            falha.append(str("Conteudo CG"))

        # tempmelt
        #--------------------------------------------------------------#
        if tm == True:
            tm_score = round(mt.Tm_GC(sequence[1:8]), 2)
            if tm_score <= tmmax:
                score += 1
            else:
                falha.append(str("tm"))

        # G°
        #--------------------------------------------------------------#
        energia_livre = free_energy(sequence)
        if -13 < energia_livre < -7:
            score += 2
        else:
            falha.append(str("Energia livre"))

        # Autores
        # Score reynolds total = 6
        #--------------------------------------------------------------#
        if autor == 'reynolds':
            r = reynolds(sequence)
            score += r[0]
            falha.extend(r[1] if isinstance(r[1], list) else [r[1]])

        # Score ui-tei total = 4
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
        return score, falha, conteudo_gc, energia_livre, tm_score, posicao

# ----------------------------------------------------------------- #
# Selecionando os siRNA funcionais
# ----------------------------------------------------------------- #
def filtro_siRNA(sequences, tuplas, conformidade=0.6,
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

        total_reynolds = round(10 * conformidade)
        total_uitei_ama = round(8 * conformidade)

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
            tm_value = resultado[4]

            if autor == "reynolds":
                if resultado[0] <= total_reynolds or not (0 < tm_value <= tmmax):
                    incluir_siRNA = False

            if autor == 'ui-tei' or autor == 'amarzguioui':
                if resultado[0] <= total_uitei_ama or not (0 < tm_value <= tmmax):
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




