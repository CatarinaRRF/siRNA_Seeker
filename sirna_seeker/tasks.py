from celery import shared_task
from celery_progress.backend import ProgressRecorder
from django.core.cache import cache
from django.contrib.auth.models import User
from .algorithm import *
from .blast import *

@shared_task(bind=True)
def selection(self, user, sequence, sequence_tag, autor, size, include_tm,
              max_tm, run_blast, organism, database, identity, query_cover):
    
    print(f'include_tm:{include_tm} and run_blast:{run_blast}')
    user_id = User.objects.get(username=user)
    # Cache settings 
    cache.set(self.request.id, user_id, 300)
    print(f"**** {cache} ****")
    print(f"**** {user_id} ****")
    print(f"**** {self.request.id} ****")

    # Meta data 
    meta = meta_data(sequence, size, autor)
    
    # Starting the selectio 
    progress_recorder = ProgressRecorder(self)
    try:
 
        # Transcribing sequences
        progress_recorder.set_progress(10, 100, description='Transcribing sequences')
        trascribe = transcrever(sequence)
        candidates, tuplas = possiveis_siRNA(trascribe, tamanho=size)

        # Filtering siRNA candidates 
        progress_recorder.set_progress(30, 100, description='Filtering siRNA candidates')
        table, sirna_verified = filtro_siRNA(candidates, tuplas, conformidade=0.6,
                                             autor=autor,
                                             tm=include_tm, tmmax=max_tm)
        
        
        # Runing Blast
        if run_blast == True:
            progress_recorder.set_progress(50, 100, description='Running Blast')
            sirna_verified_fasta = guardando_sequence(sirna_verified) 
 
            
            sequences = list(SeqIO.parse(sirna_verified_fasta, "fasta"))

            # Blast for Every sequence
            for i, sequence in enumerate(sequences, start=1):
                    print(f"Processando sequência {i}...")

                    # Blast
                    tuplas_blast, identidade = blast_siRNA (sequence , sequence_tag, 
                                                                    database, organism,
                                                                        "blastn", identity)
                    
            print("siRNAs selected successfully")

            # Returning results
            progress_recorder.set_progress(100, 100, description='siRNAs selected successfully')
            return [{
                        'table': table,
                        'sirna_verified': sirna_verified
                    },
                    meta, 
                    {
                        'tuplas_blast': tuplas_blast,
                        'identidade': identidade
                    },
                   ]
        
        else: 
            # Returning results
            progress_recorder.set_progress(100, 100, description='siRNAs selected successfully without Blast')
            print("não rodei o blast")

            return [
                    {
                        'table': table,
                        'sirna_verified': sirna_verified 
                    },  
                    meta,
                    {}
                    ]
            
    
    # Configuration of error descripitions
    except Exception as e:
        progress_recorder.set_progress(100, 100, description=f'Error: {str(e)}')
        self.update_state(state='FAILURE', meta={
            'exc_type': type(e).__name__,
            'exc_message': str(e)
        })
        raise e