1a) Escriba un código que genere secuencias de ADN random de un largo preestablecido por el usuario.

generate_random_DNA_sequence(length)
Parámetro: el largo de la secuencia.

1a) Modifique el código ahora para poder controlar además el contenido % de GC.

generate_random_DNA_weighted()
Sin parámetros, la función va a preguntar la longitud y la proporción de GC.

1b) Escriba un código que genere secuencias de proteínas random de un largo preestablecido por el usuario.

generate_random_protein_sequence(length)
Parámetro: el largo de la secuencia.

1b) modifique el mismo tal que sea posible indicar la frecuencia/probabilidad de ocurrencia de cada uno de los 20 aminoácidos.

generate_protein_sequence_weighted()
Sin parámetros, la función va a preguntar qué largo se quiere. 
Las proporciones de los aminoácidos se definen justo arriba de la función. 

1c) Escriba un código que dada una secuencia determine

i) la frecuencia de cada una de las bases de ADN. (graficoFrecBases(sequence,name))

graficoFrecBases(sequence,name=None)
Parámetros: la secuencia, nombre (opcional). 
Esta función a su vez llama a la función frecBases(sequence) que devuelve una lista con las frecuencias de las bases. 

ii) la frecuencia de cada uno de los 64 codones posibles. graficoFredCod(sequence,name=None)

graficoFredCod(sequence,name=None)
Parámetros: secuencia, nombre (opcional). 
Esta función llama a su vez a la función frecCod(sequence)

iii) la frecuencia de cada uno de los 20 aminoácidos (para una secuencia de proteínas) y para una secuencia de ADN a partir de los codones correspondientes.  Que el código sea capaz de mostrar los resultados en un gráfico de barras.

graficoFrecAa(aa_sequence,name=None)
Parámetros: secuencia, nombre (opcional). 
Si es una secuencia de nucleótidos, entonces se traduce con la función translate(sequence). 

una función que analiza todos estos aspectos de la secuencia: 
process_sequence(sequence,name=None):
Parámetros: secuencia, nombre (opcional). 

1d) Escriba un código que dada una secuencia de ADN determine la distribución de largo de los ORFs contenidos en la misma. (Nota: utilice SOLO 1 hebra y solo 1 marco de lectura x ahora)

count_orf_lengths(secuencia,name=None)
Parámetros: secuencia, nombre (opcional, aparece en el título del gráfico)
Output: una lista con las longitudes de los marcos de lectura abiertos. 

graficoOrfs(secuencia,name=None)
parámetros: secuencia, nombre (opcional, aparece en el título del gráfico). 
Esta función llama a su vez a la función count_orf_lengths(secuencia,name=None) para obtener una lista de longitudes
Output: histograma.

graficoOrfs6ML(seq,name=None)
parámetros: secuencia, nombre (opcional, aparece en el título del gráfico). 
Esta función encuentra las longitudes de los orfs para los 6 marcos de lectura, llamando a la función graficoOrf3ML(seq,name). 

Ejercicio 2) Aprovechando el manejo de secuencias de biopython, vuelva sobre su(s) secuencias de ADN al azar y ahora:
2a) Obtenga la reversa complementaria y compare las estadísticas de la misma con la original.

comparacionEjercicio2a(length) 
Parámetro: el largo de la secuencia. 
Output: 3 gráficos con la frecuencia de las bases, codones y aminoácidos. 

2b) Traduzca la secuencia y obtenga la distribución de largo de los ORF (compare con los resultados obtenidos anteriormente.

graficoOrfs_comparación(length)
Parámetro: el largo de la secuencia. 
Esta función llama a las funciones graficoFrecBases_comparacion(seq1,seq2), graficoFredCod_comparacion(seq1,seq2) y graficoFrecAa_comparacion(seq1,seq2)
Output: un histograma con los largos de los orfs. 

2b) Elija un registros de geneBank bajelo en formato fasta y gb, levantelos conBiopyhton y explore sus atributos (record, record.id, record.name, record.description, etc.) es capaz de obtener “solo” la secuencia y convertirla en un “string” para manipularla?

atributos(archivo_gb)
Parámetro: un url con una secuencia (archivo_gb ya está definida en el documento pero puede probar con otra)
Output: ID, nombre, descripción y los 100 primeros caracteres de la secuencia. 

2c) Combine lo aprendido en la clase previas (for / while loops) y/o la información disponible en biopython para cargar un conjunto de secuencias de ADN y/o proteínas. (Puede bajarlos manualmente a disco y realizar código para cargarlos y/o hacer el código que los baje directamente de internet). Obtenga solo las secuencias como cadenas de caracteres y guárdela en una lista tal que cada elemento de la lista sea una de las secuencias.

abrir_archivo(url) 
Parámetro: el url 
Output: una lista en la cual cada elemento es una secuencia. 

3a) Combine los códigos de los objetivos 1 y 2 del presente módulo y compare gráficamente las siguientes distribuciones

i) Distribución de bases en secuencia de ADN al azar vs secuencia de ADN de regiones codificantes (genes)

histogramasOrf(secuencia,name=None)
Parámetros: secuencia, nombre (opcional, aparece en el título del histograma)

histogramasJuntos(sec_cRNA,sec_RNA_azar)
Parámetros: la secuencia al real y la secuencia al azar. 
Definir la secuencia real: sec_cRNA = str(abrir_archivo(url))
Generar la secuencia al azar: sec_RNA_azar = generate_random_DNA_sequence(len(sec_cRNA))
Output: un histograma con las dos series juntas. 

ii) Distribución de aminoácidos de secuencia de proteínas al azar vs secuencias de proteínas reales

FrecAaJuntos(sec_aa_real,sec_aa_azar)
Parámetros: la secuencia de aminoácidos real y la secuencia al azar. 
Definir la secuencia real: sec_aa_real = str(abrir_archivo(url))
Generar la secuencia al azar: sec_aa_azar = generate_random_protein_sequence(len(sec_aa_real))
Output: un histograma con las dos series juntas. 

iii) Distribución de aminoácidos de secuencias de ADN al azar traducidas vs secuencias de proteínas reales

FrecAaJuntos(sec_aa_real,sec_aa_azar)
Parámetros: la secuencia de aminoácidos real y la secuencia al azar. 
Definir la secuencia real: sec_aa_real = str(abrir_archivo(url))
Generar la secuencia al azar traducida: 
sec_azar = translate(generate_random_DNA_sequence(0.3*len(sec_cRNA)))
Output: un histograma con las dos series juntas.




