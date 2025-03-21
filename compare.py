# Funzione per confrontare i valori con un errore delta
def confronta_valori(file1, file2, delta=0.01):
    # Carica i valori dai file, separando per virgola
    with open(file1, 'r') as f1:
        # Legge la prima riga e la separa con la virgola
        valori_file1 = [float(val.strip()) for val in f1.readline().split(', ')]

    with open(file2, 'r') as f2:
        # Legge la prima riga e la separa con la virgola
        valori_file2 = [float(val.strip()) for val in f2.readline().split(', ')]

    # Verifica che ogni valore di file1 sia presente in file2 entro un margine di errore delta
    for val1 in valori_file1:
        trovato = False
        for val2 in valori_file2:
            if abs(val1 - val2) <= delta:
                trovato = True
                break
        if not trovato:
            print(f"Valore {val1} non trovato in file2 con delta {delta}.")
            return False

    print("Tutti i valori di file1 sono presenti in file2 con errore delta accettabile.")
    return True

# Esegui la funzione con i nomi dei file e un delta specificato
file1 = "./py_vina.txt"
file2 = "./vina_out.txt"
delta = 0.0003

confronta_valori(file1, file2, delta)

