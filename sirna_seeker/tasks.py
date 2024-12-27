from celery import shared_task, current_task
from celery_progress.backend import ProgressRecorder
from django.core.cache import cache
from django.contrib.auth.models import User
from .algorithm import *
from .blast import *
from django_celery_results.models import TaskResult

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
    meta = meta_data(sequence, size, autor, sequence_tag)
    
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
            perc_identidy = identity/100
            result_dict = {}
            # Blast for Every sequence
            # Blast for the first 10 sequences
            for i, sequence in enumerate(sequences, start=1):
                if i > 10:  # Após 10 sequências, interrompe o loop
                    break
                print(f"Processando sequência {i}...")

                # Blast
                resultados = blast_siRNA (sequence , sequence_tag, database, organism, "blastn", perc_identidy)
                result_dict[str(sequence.seq)] = resultados

            print("siRNAs selected successfully")
        
           
            # Returning results
            progress_recorder.set_progress(100, 100, description='siRNAs selected successfully')
            return [{
                        'table': table,
                        'sirna_verified': sirna_verified
                    },
                    meta, 
                    result_dict,
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
    
