# bioinformatica

El archivo .txt tiene los resultados impresos y ordenados por el número de apariciones en la cadena original,
siendo la cadena con el mayor número de "hits" (20): "GCGCACACAC". 
Además imprime las matrices de probabilidad de transición de un aminoácido a otro, y la matriz de observacion (que son las probabilidades de apración de los aminoácidos).

El archivo validation_kmers.py, al ejecutarse imprime en consola los resultados de las cadenas generadas, y su respectiva cantidad e "hits" con la cadena original.
Básicamente el programa produce cadenas "aleatorias" de tamaño n (10), a las cuales se le calcula la distancia Hamming con las ventanas (10 elementos) de la cadena original,
aquellas cadenas que tengan distancia Hamming de máximo 2 errores serán guardadas en un diccionario {cadena: "hits"}

Para producir las cadenas aleatorias, se genera un "aminoácido" raiz (con base a la probabilidad de aparicion de los aminoácidos en la cadena original),
y los aminoácidos subsecuentes obedecen la propiedad de Markov, ya que se generan con probabilidad que depende de la matriz de transicin generada a partir de la cadena original.

