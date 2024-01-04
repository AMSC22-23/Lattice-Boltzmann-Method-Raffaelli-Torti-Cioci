import numpy as np
import matplotlib.pyplot as plt

def profilo_alare_NACA(chord, altezza, spessore, camber, num_punti=100):
    theta = np.linspace(0, 2 * np.pi, num_punti)
    
    # Equazione di un profilo alare NACA
    x = 0.5 * (1 - np.cos(theta)) * chord
    yt = 5 * spessore * (0.2969 * np.sqrt(x/chord) -
                         0.126 * (x/chord) - 0.3516 * (x/chord)**2 +
                         0.2843 * (x/chord)**3 - 0.1015 * (x/chord)**4)

    yc = camber * (1 - np.cos(theta))
    
    # Calcola la coordinata y del profilo alare sopra e sotto la linea del profilo
    y_upper = yc + yt
    y_lower = yc - yt
    
    return x, y_upper, y_lower

def genera_griglia(chord, altezza, spessore, camber, num_punti=100, griglia_size=200):
    griglia = np.zeros((griglia_size, griglia_size), dtype=int)  # Inizializza con 0
    
    x, y_upper, y_lower = profilo_alare_NACA(chord, altezza, spessore, camber, num_punti)
    
    # Normalizza le coordinate alare tra -0.2 e 0.2 per centrare l'ala nella griglia
    x_normalized = (x - chord/2) / chord * 0.4
    y_upper_normalized = (y_upper - altezza/2) / altezza * 0.4
    y_lower_normalized = (y_lower - altezza/2) / altezza * 0.4
    
    # Trova le coordinate nella griglia
    i_upper = np.round((y_upper_normalized + 0.5) * (griglia_size - 1)).astype(int)
    i_lower = np.round((y_lower_normalized + 0.5) * (griglia_size - 1)).astype(int)
    j = np.round((x_normalized + 0.5) * (griglia_size - 1)).astype(int)
    
    # Imposta i valori della griglia sopra la linea del profilo a 1
    griglia[j, i_upper] = 1
    
    # Imposta i valori della griglia sotto la linea del profilo a 1
    griglia[j, i_lower] = 1
    
    return np.flipud(griglia)  # Inverti l'ordine delle colonne

def salva_su_file_e_plot(griglia, nome_file):
    plt.imshow(griglia, cmap='gray', origin='lower')
    plt.axis('off')
    plt.savefig(nome_file)
    plt.close()

    with open(nome_file.replace('.png', '.txt'), 'w') as file:
        for colonna in np.flipud(griglia.T):  # Inverti l'ordine delle colonne
            file.write(' '.join(map(lambda x: '1' if x == 0 else '0', colonna)) + '\n')

def replace_1_with_0(line):
    for i in range(len(line)):
        if line[i] == '1':
            if '0' in line[:i] and '0' in line[i+1:]:
                line = line[:i] + '0' + line[i+1:]
    return line

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    modified_lines = [replace_1_with_0(line) for line in lines]

    with open(output_file, 'w') as outfile:
        outfile.writelines(modified_lines)

def plot_da_file(nome_file):
    griglia_modificata = np.loadtxt(nome_file, dtype=int)
    plt.imshow(griglia_modificata, cmap='gray', origin='lower')
    plt.axis('off')

    # Salva l'immagine in un file PNG
    nome_file_immagine = nome_file.replace('.txt', '_plot.png')
    plt.savefig(nome_file_immagine)

    print(f"Immagine salvata su '{nome_file_immagine}'.")

if __name__ == "__main__":
    chord = 0.7        # Lunghezza della corda dell'ala
    altezza = 0.06     # Altezza dell'ala
    spessore = 0.01    # Spessore dell'ala
    camber = 0.01      # Camber dell'ala
    
    num_punti = 500    # Numero di punti sulla corda dell'ala
    griglia_size = 300 # Dimensione della griglia
    
    # Genera e salva la griglia originale
    griglia = genera_griglia(chord, altezza, spessore, camber, num_punti, griglia_size)
    nome_file = 'profilo_alare.png'
    salva_su_file_e_plot(griglia, nome_file)
    print(f"Risultato salvato su '{nome_file}'.")

    # Modifica il file di testo
    input_file_name = "profilo_alare.txt"
    output_file_name = "profilo_alare_new.txt"
    process_file(input_file_name, output_file_name)

    # Plotta l'immagine dal file modificato
    plot_da_file(output_file_name)
