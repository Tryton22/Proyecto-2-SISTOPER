import sys
import threading
from utilities import *
from bio_seq import *
import time


# Semaforo para sincronizar la salida
output_lock = threading.Lock()


def create_object(fasta_dict, record):

    id_seq = record
    sequence = bio_seq(seq=fasta_dict[id_seq], label=id_seq)
    
    # Imprimir la informacion dentro del bloque protegido por el semaforo
    with output_lock:
        print(colored(sequence.get_seq_info()))
        print("Frecuencia de nucleotidos: ", sequence.nucleotide_frequency())
        print("Transcrito: ", colored(sequence.transcription()))
        print("Hebra complementaria: ", colored(sequence.reverse_complement()))
        print("Contenido GC:", sequence.gc_content())
        print("Traduccion: ", sequence.translate_seq())
        print("Frecuencia codon (L):", sequence.codon_usage('L'))
        for rf in sequence.gen_reading_frames():
            print("ORF: ",rf)
        print("Proteina: ", sequence.all_proteins_from_orfs())
        print("\n\n")

        # Hacer una pausa para simular el tiempo que tarda cada hebra en procesar la secuencia
        time.sleep(1)
    


def main(file_path):
    # Leer el archivo FASTA
    fasta_dict = read_FASTA(file_path)
    print(fasta_dict)

    for record in fasta_dict:
        # Crea hilo para cada secuencia
        t = threading.Thread(target=create_object, args=( fasta_dict, record, ))
        t.start()
    
    t.join()


# Verificar si se proporciona una ruta de archivo como argumento de linea de comandos
if len(sys.argv) > 1:
    # Obtener la ruta de archivo del argumento de linea de comandos
    file_path = sys.argv[1]
    # Llamar a la funcion principal con la ruta de archivo
    main(file_path)
else:
    print("No se proporciono una ruta de archivo.")


