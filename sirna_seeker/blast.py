## Importando
#--------------------#
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import requests
import xml.etree.ElementTree as ET
import time


# ----------------------------------------------------------------- #
# Blast
# ----------------------------------------------------------------- #


# Criando fasta com siRNA_verificados
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
import requests
import xml.etree.ElementTree as ET
import time
from Bio import SeqIO

# Função para submeter uma busca no BLAST
def submit_blast(sequence, organism, program='blastn', database='nt'):
    url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
    params = {
        'CMD': 'Put',
        'PROGRAM': program,
        'DATABASE': database,
        'QUERY': sequence,
        'FORMAT_TYPE': 'XML',
        'ENTREZ_QUERY': organism
    }
    response = requests.post(url, data=params)
    if response.status_code == 200:
        return response.text
    else:
        response.raise_for_status()

# Função para verificar o status de uma busca no BLAST
def check_status(rid):
    url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
    params = {
        'CMD': 'Get',
        'RID': rid,
        'FORMAT_OBJECT': 'SearchInfo'
    }
    response = requests.get(url, params=params)
    response.raise_for_status()  # Raise an exception for HTTP errors
    
    # Extrair o status da resposta
    status = response.text.split('Status=')[1].split('\n')[0].strip() 
    return status

# Função para recuperar os resultados do BLAST
def retrieve_results(rid):
    status = check_status(rid)
    while status == 'WAITING':
        print("A busca no BLAST ainda está em andamento. Aguardando...")
        time.sleep(60)  # Esperar 10 segundos antes de verificar novamente
        status = check_status(rid)

    if status == 'READY':
        url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
        params = {
            'CMD': 'Get',
            'RID': rid,
            'FORMAT_TYPE': 'XML'
        }
        response = requests.get(url, params=params)
        response.raise_for_status()  # Levantar uma exceção para erros HTTP
        return response.text
    else:
        print(f"Erro: A busca no BLAST falhou com status {status}")
        return None

# Função para analisar os resultados do BLAST
def parse_blast_results(xml_data, sequence, identidade, query_cover):
    root = ET.fromstring(xml_data)
    results = []
    for hit in root.findall('.//Hit'):
        accession = hit.find('Hit_accession').text  # Obtendo a Accession ao invés de Hit_def
        hit_len = hit.find('Hit_len').text
        hsp = hit.find('.//Hsp')
        identity = hsp.find('Hsp_identity').text
        query_from = int(hsp.find('Hsp_query-from').text)
        query_to = int(hsp.find('Hsp_query-to').text)
        query_len = int(root.find('BlastOutput_query-len').text)
        query_cover = (query_to - query_from + 1) / query_len * 100 

        if float(identity) >= identidade and query_cover >= query_cover:
            results.append({
                'Accession': accession,
                'Hit_len': hit_len,
                'Identity': identity,
                'Query_cover': query_cover
            })
    return results

# Função para realizar a busca de siRNA
def identidade_siRNA(sequence, organism, identidade, query_cover, db='refseq_rna'):
    response_text = submit_blast(sequence, organism, program='blastn', database=db)
    
    rid_start = response_text.find('RID = ') + len('RID = ')
    rid_end = response_text.find('\n', rid_start) 
    rid = response_text[rid_start:rid_end].strip()

    results = retrieve_results(rid)
    if results:
        parsed_results = parse_blast_results(results, sequence, identidade, query_cover,)
        return parsed_results
    else:
        return []

# Função principal para processar um arquivo FASTA e armazenar os resultados
def process_fasta_file(file_path, organism, identidade, query_cover, db='refseq_rna'):
    results = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        sequence_id = record.id
        blast_results = identidade_siRNA(sequence, organism, db, identidade, query_cover)
        results[sequence_id] = blast_results
    return results

"""# Exemplo de uso:
process_fasta_file(fasta_file_path, organism=organismo, db=db, identidade=identidade, query_cover=query_cover)"""