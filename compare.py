def confronta_file_con_ordine(file1, file2, delta=0.01):
    """
    Confronta due file che contengono valori float separati da virgole.
    Verifica che i valori corrispondano entro un errore massimo delta e che l'ordine sia mantenuto.
    """
    # Legge i file e separa i valori
    with open(file1, 'r') as f1:
        valori_file1 = [float(val.strip()) for val in f1.readline().split(',')]

    with open(file2, 'r') as f2:
        valori_file2 = [float(val.strip()) for val in f2.readline().split(',')]

    # Controlla che abbiano la stessa lunghezza
    if len(valori_file1) != len(valori_file2):
        print(f"Errore: il numero di valori nei file è diverso! ({len(valori_file1)} vs {len(valori_file2)})")
        return False

    # Confronta ogni valore nella stessa posizione
    for i, (val1, val2) in enumerate(zip(valori_file1, valori_file2)):
        if abs(val1 - val2) > delta:
            print(f"Errore alla posizione {i}: {val1} ≠ {val2} (differenza {abs(val1 - val2)}, fuori dal delta {delta})")
            return False

    print("Tutti i valori corrispondono e l'ordine è mantenuto!")
    return True

# Esegui la funzione con i file di output
file1 = "py_vina_out.txt"
file2 = "vina_out.txt"
delta = 0.0002

confronta_file_con_ordine(file1, file2, delta)
