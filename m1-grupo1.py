import random
import requests
from collections import Counter
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#from scipy.stats import chi2_contingency
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
from Bio import ExPASy
from Bio import SwissProt
import inspect
from weakref import ref


nucleotides = ['A', 'C', 'T', 'G']
aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
codones = {codon: 0 for codon in ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']}
genetic_code = {
  "TTT": "F", "TTC": "F",  # Fenilalanina
  "TTA": "L", "TTG": "L",  # Leucina
  "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",  # Leucina
  "ATT": "I", "ATC": "I", "ATA": "I",  # Isoleucina
  "ATG": "M",  # Metionina (START)
  "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",  # Valina
  "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",  # Serina
  "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",  # Prolina
  "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",  # Treonina
  "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",  # Alanina
  "TAT": "Y", "TAC": "Y",  # Tirosina
  "TAA":  "", "TAG": "",  # STOP
  "CAT": "H", "CAC": "H",  # Histidina
  "CAA": "Q", "CAG": "Q",  # Glutamina
  "AAT": "N", "AAC": "N",  # Asparagina
  "AAA": "K", "AAG": "K",  # Lisina
  "GAT": "D", "GAC": "D",  # Ácido aspártico
  "GAA": "E", "GAG": "E",  # Ácido glutámico
  "TGT": "C", "TGC": "C",  # Cisteína
  "TGA":  "",  # STOP
  "TGG": "W",  # Triptófano
  "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",  # Arginina
  "AGT": "S", "AGC": "S",  # Serina
  "AGA": "R", "AGG": "R",  # Arginina
  "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"  # Glicina
          }
aminoacidos = {aa: 0 for aa in ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']}

# 1a)
def generate_random_DNA_sequence(length): # Genera una secuencia de ntd al azar con iguales probabilidades para cada ntd
  sequence = ''
  for i in range(length):
    sequence += random.choice(nucleotides)
  return sequence

# 1a)
def generate_random_DNA_weighted(): # Genera una secuencia  (string) en la que el contenido GC lo define el usuario.
  # input= La longitud de la secuencia y la frecuencia de GC
  # output= string de la secuencia
  length=int(input('Seleccione el largo de la secuencia de nt: '))
  CG=float(input("Seleccione la frecuencia GC (entre 0 y 1): ")) #frecuencia de CG (entre 0 y 1)
  AT=1-CG #frecuencia de AT (se calcula sola)
  sequence=""
  frec=[AT/2,CG/2,AT/2,CG/2] # frecuencias de cada par de ntds
  for i in range(length):
    sequence += np.random.choice(nucleotides,p=frec) # agrega nucleotidos con la probabilidad definida
#  print("la secuencia ponderada es:",sequence)
  return(sequence)

# 1b)
def generate_random_protein_sequence(length):
  sequence_protein = ''
  for i in range(length):
    sequence_protein += random.choice(aminoacids)
  return sequence_protein

# 1b)
probabilidades=[0.00,0.1,0.00,0.1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]

def generate_protein_sequence_weighted(): # Genera un str de secuencia aminoacídica donde todas los aá tienen != probabilidad. Las frecuencias se definen arriba
  lista_num_azar=[]
  numeros=list(range(0,20))
  sec_prot=""
  len_sec_prot=int(input("Seleccione el largo de la secuencia de aminoácidos: "))
  for i in range(len_sec_prot):
    lista_num_azar.append(np.random.choice(numeros,p=probabilidades))
  for i in lista_num_azar:
    sec_prot+=aminoacids[int(i)-1]
  return(sec_prot)

# 1c) i)
def frecBases(sequence):
  lista_frec=[]
  cod_gen_unicos=set(sequence)
  lista_frec.append(sequence.count("A")*100/len(sequence))
  lista_frec.append(sequence.count("C")*100/len(sequence))
  lista_frec.append(sequence.count("T")*100/len(sequence))
  lista_frec.append(sequence.count("G")*100/len(sequence))
  return(lista_frec)

def graficoFrecBases(sequence,name=None): # Gráfico de frecuencias de cada base
  frec_bases=frecBases(sequence)
  if name is None:
    frame_info = inspect.currentframe().f_back # para obtener el nombre de la variable y pasarla al titulo del histograma
    name = [name for name, value in frame_info.f_locals.items() if value is sequence][0]
  x_values = nucleotides # Valores para el eje X
  y_values = frec_bases # Valores para el eje Y
  fig = go.Figure(data=go.Bar(x=x_values, y=y_values,marker_color='#576EF4'))
  fig.update_layout(
    title_text=f'Frecuencias de las bases en la secuencia: {name}',
    xaxis_title='Bases',
    yaxis_title='Porcentaje (%)',
    yaxis=dict(range=[0, 100]),  # Para que el eje y vaya de 0 a 100%
    template='plotly')
  fig.show()

# 1c) ii)
def frecCod(sequence,name=None): # Separa los codones y los cuenta. Devuelve
  for cd in range(0,len(sequence)-2,3):
    if sequence[cd:cd+3] in codones :
      codones[sequence[cd:cd+3]] +=1 #Cuenta las veces que aparece cada codon.
    codonesTot =sum(list(codones.values()))
    porcentajesCod = [(codon*100)/codonesTot for codon in list(codones.values())] #Valores para el eje Y.
    y_values_cod = porcentajesCod
  return(y_values_cod)

def graficoFredCod(sequence,name=None): # Gráfico de frecuencias de cada codon
  if name is None:
    frame_info = inspect.currentframe().f_back
    name = [name for name, value in frame_info.f_locals.items() if value is sequence][0]
  y_values_cod=frecCod(sequence)
  x_values_cod = list(codones.keys())
  fig = go.Figure(data=go.Bar(x=x_values_cod, y=y_values_cod,marker_color='#576EF4'))
  fig.update_layout(
    title_text=f"Porcentaje de Cada Codon en la Secuencia: {name}",
    xaxis_title='Codones',
    yaxis_title='Porcentaje (%)',
    yaxis=dict(range=[0, 10]),  # Para que el eje y vaya de 0 a 100%
    template='plotly')
  fig.show()

# 1c) iii)
def translate(sequence): # Traduce una secuencia de ntd a una de aa segun el codigo establecido arriba
  global genetic_code
  amino_acid_sequence_t = []
  for i in range(0, len(sequence) - 2, 3):
    codon = sequence[i:i+3]
    amino_acid = genetic_code.get(codon)
    amino_acid_sequence_t.append(amino_acid) # es la lista con aa de la secuencia dada
    aa_seq_t=''.join(amino_acid_sequence_t)
  return(aa_seq_t)

def graficoFrecAa(aa_sequence,name=None):
  # Calcula la frecuencia de aa en una secuencia dada.
  # Toma en cuenta si la secuencia ingresada es de ntd (en cuyo caso la traduce) o si son aa.
  if name is None:
    frame_info = inspect.currentframe().f_back
    name = [name for name, value in frame_info.f_locals.items() if value is aa_sequence][0]
  if all(char in nucleotides for char in aa_sequence):
    print('perá que hay que traducir esta secuencia de nucleótidos')
    aa_sequence_translated=translate(aa_sequence)
    graficoFrecAa(aa_sequence_translated)
  else:
    for a in range(len(aa_sequence)):
      if aa_sequence[a]=='X' or aa_sequence[a]=='Z':
          continue
      else:
        aminoacidos[aa_sequence[a]] += 1 #Cuenta las veces que aparece cada aa
    x_values_aa = list(aminoacidos.keys())
    aminoacidosTot =sum(list(aminoacidos.values()))
    porcentajesAA = [(aa*100)/aminoacidosTot for aa in list(aminoacidos.values())] #Valores para el eje Y.
    # Ordenar los datos en función del porcentaje
    sorted_data = sorted(zip(x_values_aa, porcentajesAA), key=lambda pair: pair[1])
    # Separar las listas ordenadas
    x_values_sorted, y_values_sorted = zip(*sorted_data)
    y_values_aa_sorted = porcentajesAA
    fig = go.Figure(data=go.Bar(x=x_values_sorted, y=y_values_sorted,marker_color='#576EF4'))
    fig.update_layout(
      title_text=f'Porcentaje de Cada Aminoácido en la Secuencia: {name}',
      xaxis_title='Aminoácidos',
      yaxis_title='Porcentaje (%)',
      yaxis=dict(range=[0, 10]),  # Ajusta el rango del eje y si es necesario
      template='plotly')
    fig.show()

def process_sequence(sequence,name=None):
  if all(char in nucleotides for char in sequence):
    print("Processing nucleotide sequence...")
    graficoFrecBases(sequence,name)
    graficoFredCod(sequence,name)
    graficoFrecAa(sequence,name)
  else:
    print("Processing amino acid sequence...")
    graficoFrecAa(sequence,name)

# 1d)
def count_orf_lengths(secuencia,name=None): # output: una lista con las longitudes de los orfs
  stop = ['TAA', 'TAG', 'TGA']
  for frame in range(3):
    contador=  0
    largoORF=[]
    for cod in range(frame,len(secuencia)-2,3):
      if secuencia[cod:cod+3] in stop:
        largoORF.append(contador)
        contador = 0 # si encuentra un codon de stop resetea el contador a 0
      else:
        contador +=1
  return(largoORF)

def graficoOrfs(secuencia,name=None):
  if name is None:
    frame_info = inspect.currentframe().f_back # para obtener el nombre de la variable y pasarla al titulo del histograma
    name = [name for name, value in frame_info.f_locals.items() if value is secuencia][0]
  y_values=count_orf_lengths(secuencia,name=None) # largo de los distintos ORF
  media=np.mean(y_values)
  fig=go.Figure()
  fig.add_trace(go.Histogram(x=y_values,nbinsx=300))
  fig.add_vline(x=media, line_dash="dash", line_color="red", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
  fig.update_layout(title_text = f'Secuencia: {name}')
  fig.update_xaxes(title_text="Longitud ORF (codones)")
  fig.update_yaxes(title_text="Frecuencia")
  fig.show()

def graficoOrfs3ML(seq,name=None): #cuenta los orfs y hace un histograma para cada marco de lectura
    if name is None:
      frame_info = inspect.currentframe().f_back # para obtener el nombre de la variable y pasarla al titulo del histograma
      name = [name for name, value in frame_info.f_locals.items() if value is seq][0]
    stop = ['TAA', 'TAG', 'TGA']
    for frame in range(3):
      contador=  0
      largoORF=[]
      for cod in range(frame,len(seq)-2,3):
          if seq[cod:cod+3] in stop:
              largoORF.append(contador)
              contador = 0 # si encuentra un codon de stop resetea el contador a 0
          else:
              contador +=1 # si no es un codon de stop suma de a 1 al contador.
      y_values=largoORF # largo de los distintos ORF
      media=np.mean(y_values)
      fig=go.Figure()
      fig.add_trace(go.Histogram(x=y_values,nbinsx=300))
      fig.add_vline(x=media, line_dash="dash", line_color="red", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
      fig.update_layout(title_text = f'Marco de lecura - {frame+1}, Secuencia: {name}')
      fig.update_xaxes(title_text="Longitud ORF (codones)")
      fig.update_yaxes(title_text="Frecuencia")
      fig.show()
    return # muestra los tres 3 histogramas

def graficoOrfs6ML(seq,name=None):
  if name is None:
      frame_info = inspect.currentframe().f_back
      name = [name for name, value in frame_info.f_locals.items() if value is seq][0]
  graficoOrfs3ML(seq,name)
  reverse_complement=Seq(seq).reverse_complement()
  graficoOrfs3ML(str(reverse_complement),name)


# 2a)
def graficoFrecBases_comparacion(seq1,seq2,name1=None,name2=None):
  serie1=frecBases(seq1)
  serie2=frecBases(seq2)
  if name1 is None:
    frame_info = inspect.currentframe().f_back
    name1 = [name for name, value in frame_info.f_locals.items() if value is seq1][0]
  if name2 is None:
    frame_info = inspect.currentframe().f_back
    name2 = [name for name, value in frame_info.f_locals.items() if value is seq2][0]
    fig = go.Figure()
  fig.add_trace(go.Bar(
    x=nucleotides,
    y=serie1,
    name=f'{name1}',
    marker_color='blue'))
  fig.add_trace(go.Bar(
    x=nucleotides,
    y=serie2,
    name=f'{name2}',
    marker_color='lightskyblue'))
  fig.update_layout(
    title=f'Frecuencias de las bases en {name1} y {name2}',
    xaxis_title='Bases',
    yaxis_title='Frecuencias (%)',
    barmode='group')

  fig.show()

def graficoFredCod_comparacion(seq1,seq2,name1=None,name2=None):
  y_values_cod1=frecCod(seq1)
  y_values_cod2=frecCod(seq2)
  x_values_cod = list(codones.keys())
  if name1 is None:
    frame_info = inspect.currentframe().f_back
    name1 = [name for name, value in frame_info.f_locals.items() if value is seq1][0]
  if name2 is None:
    frame_info = inspect.currentframe().f_back
    name2 = [name for name, value in frame_info.f_locals.items() if value is seq2][0]
  fig = go.Figure()
  fig.add_trace(go.Bar(
    x=x_values_cod,
    y=y_values_cod1,
    name=f'{name1}',
    marker_color='blue'))
  fig.add_trace(go.Bar(
    x=x_values_cod,
    y=y_values_cod2,
    name=f'{name2}',
    marker_color='lightskyblue'))
  fig.update_layout(
    title=f'Frecuencias de codones de {name1} y {name2}',
    xaxis_title='Codones',
    yaxis_title='Frecuencias',
    barmode='group')

  fig.show()

def graficoFrecAa_comparacion(seq1,seq2,name1=None,name2=None):
  # prot al azar
  if name1 is None:
    frame_info = inspect.currentframe().f_back
    name1 = [name for name, value in frame_info.f_locals.items() if value is seq1][0]
  for a in range(len(seq1)):
    if seq1[a]=='X' or seq1[a]=='Z':
          continue
    else:
      aminoacidos[seq1[a]] += 1 #Cuenta las veces que aparece cada aa
  x_values_aa = list(aminoacidos.keys())
  aminoacidosTot =sum(list(aminoacidos.values()))
  porcentajesAA_azar = [(aa*100)/aminoacidosTot for aa in list(aminoacidos.values())] #Valores para el eje Y.
  # prot real
  if name2 is None:
    frame_info = inspect.currentframe().f_back
    name2 = [name for name, value in frame_info.f_locals.items() if value is seq2][0]
  for a in range(len(seq2)):
    if seq2[a]=='X' or seq2[a]=='Z':
      continue
    else:
        aminoacidos[seq2[a]] += 1 #Cuenta las veces que aparece cada aa
  x_values_aa = list(aminoacidos.keys())
  aminoacidosTot =sum(list(aminoacidos.values()))
  porcentajesAA_real = [(aa*100)/aminoacidosTot for aa in list(aminoacidos.values())] #Valores para el eje Y.
  fig = go.Figure()
  # Agregar barras para los aminoácidos de la secuencia al azar
  fig.add_trace(go.Bar(
    x=x_values_aa,
    y=porcentajesAA_azar,
    name=f'{name1}',
    marker_color='blue'))
  # Agregar barras para los aminoácidos de las secuencias reales
  fig.add_trace(go.Bar(
    x=x_values_aa,
    y=porcentajesAA_real,
    name=f'{name2}',
    marker_color='lightskyblue'))
  # Configurar el layout del gráfico
  fig.update_layout(
    title=f'Comparación de Porcentajes de Aminoácidos en {name1} y {name2}',
    xaxis_title='Aminoácidos',
    yaxis_title='Porcentaje (%)',
    yaxis=dict(range=[0, 10]),  # Ajusta el rango del eje Y si es necesario
    barmode='group',  # Agrupar barras lado a lado
    template='plotly')

  fig.show()

def graficoOrfs_comparacion(seq1,seq2,name1=None,name2=None):
  # seq 1
  if name1 is None:
    frame_info = inspect.currentframe().f_back
    name1 = [name for name, value in frame_info.f_locals.items() if value is seq1][0]
  serie1=count_orf_lengths(seq1)
  media=np.mean(serie1)
  fig = go.Figure()
  fig.add_trace(go.Histogram(
    x=serie1,
    nbinsx=250,
    name=f'{name1}',
    marker_color='blue',
    opacity=0.55))
  fig.add_vline(x=media, line_dash="dash", line_color="red", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
  # seq 2
  if name2 is None:
    frame_info = inspect.currentframe().f_back
    name2 = [name for name, value in frame_info.f_locals.items() if value is seq2][0]
  serie2=count_orf_lengths(seq2)
  media=np.mean(serie2)
  fig.add_trace(go.Histogram(
    x=serie2,
    nbinsx=250,
    name=f'{name2}',
    marker_color='lightskyblue',
    opacity=0.65))
  fig.add_vline(x=media, line_dash="dash", line_color="violet", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
  # juntos
  fig.update_layout(
    title=f'Longitudes de ORFs en las secuencias {name1} y {name2}',
    xaxis_title='ORF Length',
    yaxis_title='Frequency',
    barmode='overlay')
  fig.show()

def comparacionEjercicio2a(lenght):
  secuencia=generate_random_DNA_sequence(lenght)
  reversa_complementaria=Seq(secuencia)
  reversa_complementaria=reversa_complementaria.reverse_complement()
  secuencia=Seq(secuencia)
  secuencia_traducida=translate(secuencia)
  reversa_complementaria_traducida=translate(reversa_complementaria)
  graficoFrecBases_comparacion(secuencia,reversa_complementaria)
  graficoFredCod_comparacion(secuencia,reversa_complementaria)
  graficoFrecAa_comparacion(secuencia_traducida,reversa_complementaria_traducida)
  graficoOrfs_comparacion(secuencia,reversa_complementaria)

# 2b)
# URL del archivo GenBank
#archivo_gb = 'https://raw.githubusercontent.com/Nehuenpg/Fasta-Sequences/main/sequence.gb'

def atributos(url):
  response = requests.get(url) # Descarga el archivo desde la URL
  response.raise_for_status()  # Lanza un error si la descarga falla
  from io import StringIO # Guarda el contenido del archivo en un objeto StringIO
  file_content = StringIO(response.text)
  record = SeqIO.read(file_content, "genbank") # Lee el archivo GenBank desde el contenido descargado
  print(f"ID: {record.id}")
  print(f"Name: {record.name}")
  print(f"Description: {record.description}")
  print(f"Sequence: {record.seq[:100]}...")  # Muestra los primeros 100 nucleótidos

# 2c)
def abrir_archivo(url): #abre el url, genera una lista con cada secuencia y devuelve un str con todas las secuencias seguidas.
  responseNt = requests.get(url)
  secuencia = responseNt.text #lo transforma en archivo de texto
  archivo = secuencia.splitlines()
  lista_sec = []
  secuencia_parcial = ''
  for line in range(len(archivo)):
    if archivo[line].startswith('>'):
      nro_linea = line+1
      while not archivo[nro_linea].startswith('>') and nro_linea+1 != len(archivo):
        seq = archivo[nro_linea].strip()
        secuencia_parcial +=seq
        nro_linea += 1
      lista_sec.append(secuencia_parcial)
  print('se leyeron las',len(lista_sec),'secuencias')
  return(lista_sec)

# 3a) i)
#url_cRNA='https://raw.githubusercontent.com/agusvidaly/secuencias1/main/seq_codificantes_menos.txt'
#sec_cRNA=''.join(abrir_archivo(url_cRNA))
#sec_RNA_azar = generate_random_DNA_sequence(len(sec_cRNA))
#graficoFrecBases_comparacion(sec_cRNA,sec_RNA_azar)

# 3a ii)  frecs aa en prtos reales vs prots azar
#url_prots = 'https://raw.githubusercontent.com/Nehuenpg/Fasta-Sequences/main/SecuenciaProteinas'
#url_prots = 'https://raw.githubusercontent.com/agusvidaly/secuencias1/main/prots_reales.txt'
#sec_aa_real = ''.join(abrir_archivo(url_prots))
#sec_aa_azar=generate_random_protein_sequence(70000)
#graficoFrecAa_comparacion(sec_aa_real,sec_aa_azar, 'secuencias de proteínas reales', 'secuencia de aminoácidos generada al azar')


# 3a iii) frecs aa en prots reales vs ntds traducidos
#url_prots = 'https://raw.githubusercontent.com/Nehuenpg/Fasta-Sequences/main/SecuenciaProteinas'
#url_prots='https://raw.githubusercontent.com/agusvidaly/secuencias1/main/prots_reales.txt'
#sec_aa_real = ''.join(abrir_archivo(url_prots))
#sec_ntd_azar = generate_random_DNA_sequence(40000)
#sec_ntd_azar_traducida = translate(sec_ntd_azar)
#graficoFrecAa_comparacion(sec_aa_real ,sec_ntd_azar_traducida )

# 3a2)
#seq1 = generate_random_DNA_sequence(100000)
#seq2 = generate_random_DNA_sequence(40000)
#graficoFrecBases_comparacion(seq1,seq2)
#graficoFredCod_comparacion(seq1,seq2)
#graficoOrfs_comparacion(seq1,seq2)
#seq1_translated=translate(seq1)
#seq2_translated=translate(seq2)
#graficoFrecAa_comparacion(seq1_translated,seq2_translated)
#prot1 = generate_random_protein_sequence(70000)
#prot2 = generate_random_protein_sequence(100000)
#graficoFrecAa_comparacion(prot1,prot2)

# 3c)
def marcos_lectura(secuencia):
  # input: una Seq
  #output: diccionario con los 6 marcos de lectura {'ML-1': 'secuencia', 'ML-2':'secuencia2'...}
  marcos_lectura = {}
  for i in range(3):
    marcos_lectura['ML-'+str(i+1)] = secuencia[i:]
    marcos_lectura['ML-'+str(i+4)] = secuencia.reverse_complement()[i:]
  return marcos_lectura

def orfs(secuencia):
  #output: diccionario con los ORFs
  stop = ['TAA','TAG','TGA']
  orfs = {}
  for ML in secuencia:
    seq = secuencia[ML]
    orfs[ML] = []
    contador = 0
    for i in range(0, len(seq), 3):
      codon = seq[i:i+3]
      if codon in stop:
        orfs[ML].append(contador)
        contador = 0
      else:
        contador += 1
  return orfs

def graficos_ml(orfs):
  #output: gráficos de frecuencia de longitud de los ORF de cada ML
  for ML in orfs:
    orf_filtrados = [largoORF for largoORF in orfs[ML] if largoORF > 20]  # largo de los distintos ORF filtrando los que son menores a 20 codones.
    print(ML, orf_filtrados)
    '''Gráfico de frecuencia de longitud de los ORF de cada ML'''
    y_values = orf_filtrados  # largo de los distintos ORF filtrando los que son menores a 20 codones.
    media =np.mean(y_values)
    fig = go.Figure()
    fig.add_trace(go.Histogram(x=y_values,nbinsx=300))
    fig.add_vline(x=media, line_dash="dash", line_color="red", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
    fig.update_xaxes(title_text=f"Longitud ORF {ML}")
    fig.update_yaxes(title_text="Frecuencia")
    fig.show()
  return

#url = 'https://raw.githubusercontent.com/agusvidaly/secuencias1/main/1_seq_mRNA.txt' # G.gallus mRNA for mitochondrial creatine kinase
#url = 'https://raw.githubusercontent.com/Nehuenpg/Fasta-Sequences/main/EF136884.1' # Escherichia coli K-12 lac operon, complete sequence
#secuencia = Seq(''.join(abrir_archivo(url)))
#marcos_lectura=marcos_lectura(secuencia)
#orfs=orfs(marcos_lectura)
#graficos_ml(orfs)

# 3c iv)
'''Trabajo con una secuncia que contiene intrones y exones correspondiente al gen JAK2 en humanos'''
# URL de la secuencia FASTA
url_adn_jak2 = "https://raw.githubusercontent.com/Nehuenpg/Fasta-Sequences/main/JAK2intronexones"

# Obtener el contenido desde la URL de una secuencia Real
response_jak2 = requests.get(url_adn_jak2)
fasta_data = response_jak2.text #lo transforma en archivo de texto

secuencia_jak2 = Seq(''.join(fasta_data.split('\n')[1:]).replace('\n', '')) #configura la secuencia para que pueda ser procesada.
secuencia_jak2_rv_complement = secuencia_jak2.reverse_complement()
secuencia_azar = Seq(''.join([random.choice('ATCG') for _ in range(len(secuencia_jak2))])) #creo una secuencia al azar de la misma longitud que la real.
secuencia_azar_rv_complement = secuencia_azar.reverse_complement()

'''Definos variables que voy a usar para comparar.'''

marcos_jak2 = marcos_lectura(secuencia_jak2)
marcos_azar = marcos_lectura(secuencia_azar)
orfs_jak2 = orfs(marcos_jak2)
orfs_azar = orfs(marcos_azar)
graficos_jak2 = graficos_ml(orfs_jak2)
graficos_azar = graficos_ml(orfs_azar)

'''Diccionarios generados a partir de 3000 secuencias de ADNc con los que voy a comparar mi ADN incognita.'''
#Diccionario de referencia
referencia_aa_reales = {'M': 1.9590386549742596, 'W': 2.082172507068479, 'Y': 2.724825645696219, 'D': 2.8375390794129163, 'H': 3.2431163239494203, 'N': 3.5069485761453874, 'C': 3.8157387905989624, 'Q': 4.401466412238577, 'E': 4.470404987967602, 'I': 4.818283147339483, 'F': 4.914296609245914, 'T': 5.312593212121267, 'V': 5.379484107383094, 'K': 5.4955193338576915, 'A': 5.606549119190247, 'G': 6.125658869630327, 'R': 6.3456252558178585, 'P': 6.515173198501917, 'S': 9.69472041994741, 'L': 10.750845748912965}
referencia_codones_reales = {'AGT': 1.7214622167688194, 'CTC': 1.7007578920449207, 'GCC': 1.65721055410916, 'CAA': 1.8526502787006875, 'ATT': 2.0395124973351715, 'CCC': 1.7505847834134243, 'TTA': 1.6015591494117356, 'CAG': 2.548816133537889, 'ATG': 1.9590386549742596, 'CCG': 0.6062954343323663, 'AGA': 2.4033853163563936, 'CCT': 2.014007499515951, 'TGC': 1.715296423362032, 'TAC': 1.1162361266804346, 'AAA': 3.093681153854292, 'CTT': 1.9761026588675827, 'CTG': 2.397378786985944, 'GCA': 1.7252845536409238, 'AAT': 2.0078644581143545, 'TGT': 2.1004423672369303, 'TCT': 2.0082284901974123, 'ATC': 1.3604789024068664, 'GTT': 1.5060689836246992, 'AGC': 1.854333927084829, 'CCA': 2.1442854812401753, 'TCA': 1.9513484772196685, 'CTA': 1.1384193317417548, 'CAC': 1.5155110657790047, 'TTC': 1.8893265110687367, 'ACT': 1.554052962572724, 'GTC': 1.1138016621249873, 'AAG': 2.401838180003399, 'AAC': 1.4990841180310324, 'TCG': 0.4187279035369585, 'GCT': 1.745442830240236, 'GAT': 1.558512355590179, 'AGG': 1.8885529428922394, 'TTT': 3.0249700981771777, 'ACA': 1.9037285303547016, 'GTA': 0.9944673948976808, 'TCC': 1.7406194051397235, 'TAT': 1.6085895190157846, 'GAA': 2.351601752541456, 'GTG': 1.7651460667357266, 'CGA': 0.4680087467808757, 'ACG': 0.41981999978613116, 'CAT': 1.7276052581704158, 'CGC': 0.5400415952158903, 'GGA': 1.91571883709041, 'GAC': 1.2790267238227373, 'TTG': 1.9366279298610285, 'ATA': 1.418291747597445, 'ACC': 1.4349917194077106, 'TGG': 2.082172507068479, 'GAG': 2.1188032354261463, 'GCG': 0.4786111811999271, 'GGC': 1.52142658712869, 'CGT': 0.42168566421180115, 'GGT': 1.1765971964524164, 'GGG': 1.5119162489588114, 'CGG': 0.6239509903606579}
referencia_aa_azar = {'M': 1.6390142240310601, 'W': 1.6404995835280487, 'Y': 3.2639061069984425, 'C': 3.264500250797238, 'Q': 3.2673567113683704, 'E': 3.2731381875643417, 'N': 3.275629021182369, 'K': 3.2770001222565126, 'F': 3.2817989760160144, 'H': 3.285980834292152, 'D': 3.286460719668102, 'I': 4.9134092542467, 'P': 6.5480816582096395, 'G': 6.5486529503238655, 'T': 6.554754350103804, 'A': 6.5617012622127975, 'V': 6.566408709234023, 'S': 9.844597119088126, 'R': 9.848207685250037, 'L': 9.858902273628356}
referencia_codones_azar = {'GTC': 1.630239177156542, 'TCT': 1.6399739947829604, 'ATA': 1.637368902742088, 'AGA': 1.6464638732005727, 'GGC': 1.6368433139969996, 'AGG': 1.6427847519849545, 'AAC': 1.6439044845288382, 'TTA': 1.6379858982254525, 'CAG': 1.638580042024248, 'CCG': 1.6350608826006132, 'TCG': 1.6488861517648927, 'TAT': 1.6321358669757737, 'GGA': 1.6373231993729498, 'GTT': 1.6527480864570634, 'TGC': 1.635883543245099, 'GCA': 1.6319302018146522, 'CCA': 1.6310389861164591, 'AAT': 1.6317245366535307, 'TGG': 1.6404995835280487, 'CTA': 1.6476293091135945, 'CCT': 1.6432189339917664, 'GCG': 1.6460753945628988, 'ATG': 1.6390142240310601, 'CTT': 1.6384200802322646, 'ACA': 1.6368204623124305, 'ACG': 1.6360435050370825, 'GCC': 1.639448406037872, 'TTG': 1.6460753945628988, 'ATT': 1.6321587186603428, 'CAT': 1.6486576349192021, 'CTG': 1.641756426179347, 'ACC': 1.6409566172194299, 'CAC': 1.6373231993729498, 'TCC': 1.6438587811597, 'CGC': 1.6395169610915792, 'GTA': 1.6379173431717453, 'GGG': 1.6356321747148395, 'GTG': 1.6455041024486723, 'GCT': 1.644247259797374, 'TGT': 1.628616707552139, 'GAC': 1.6382372667557121, 'ACT': 1.6409337655348608, 'AGC': 1.6480406394358376, 'GGT': 1.6388542622390767, 'GAT': 1.64822345291239, 'TTT': 1.6374603094803641, 'CAA': 1.6287766693441224, 'GAG': 1.634603848909232, 'TCA': 1.636980424104414, 'CGA': 1.6375060128495023, 'CGG': 1.645001365388153, 'AGT': 1.6268571278403217, 'GAA': 1.63853433865511, 'CTC': 1.647035165314799, 'TTC': 1.6443386665356503, 'AAG': 1.642213459870728, 'CGT': 1.6369347207352758, 'AAA': 1.6347866623857845, 'ATC': 1.6438816328442691, 'TAC': 1.6317702400226688, 'CCC': 1.6387628555008005}

def graficoProbabilidades(marcos_lectura_incognita, ventana, pasos,referencia):
  for ml in range(1): # Elijo cuantos marcos de lectura a procesar cambiando ese valor.
    marco_lectura = marcos_lectura_incognita['ML-'+str(ml+1)]
    ecm_values = []

    for i in range(0, len(marco_lectura) - ventana + 1, pasos):
      secuencia_incognita = marco_lectura[i:i + ventana].translate()
      frecuencias_incognita = Counter(secuencia_incognita)

      total_aa = sum(frecuencias_incognita.values())
      frecuencias_normalizadas = {aa: (frecuencias_incognita.get(aa, 0) / total_aa) * 100 for aa in referencia.keys() if aa != '*'}

      # Lista para guardar los valores de Error Cuadratico Medio
      # Calcular el ECM entre las frecuencias de la ventana y las del diccionario de referencia
      ecm = np.mean([(frecuencias_normalizadas.get(aa, 0) - referencia[aa]) ** 2 for aa in referencia.keys()])
      ecm_values.append(ecm)

    # Graficar los valores de ECM
    plt.figure(figsize=(12, 6))
    plt.bar(range(len(ecm_values)), ecm_values, color='b', alpha=0.6, label='ECM')


    plt.xlabel('Posición de la Ventana')
    plt.ylabel('ECM (Error Cuadrático Medio)')
    plt.title(f'ECM entre la Muestra Incógnita {ml+1} y la Referencia')
    plt.legend()
    plt.grid(True)
    plt.show()
  return ecm_values

#ecm_values_2 = graficoProbabilidades(marcos_jak2, 5000, 3, referencia_aa_azar)
#ecm_values_3 = graficoProbabilidades(marcos_azar, 5000, 3, referencia_aa_reales)
#ecm_values_4 = graficoProbabilidades(marcos_azar, 5000, 3, referencia_aa_azar)

''' En la comparación de una secuencia incógnita contra una secuencia real de ADN (Gráfico 1) y contra una secuencia de ADN generada al azar (Gráfico 2), 
 representados mediante el error cuadrático medio (ECM), se observan regiones de similitud y discrepancias. 
 Estas diferencias deben ser estudiadas con mayor profundidad para identificar posibles regiones codificantes. 
 El análisis se realizó recorriendo la secuencia incógnita en ventanas de 5000 nucleótidos.
 Al comparar la secuencia incógnita con el ADN real, el ECM varía entre 1 y 7, indicando fluctuaciones en la similitud que podrían reflejar regiones funcionales. 
 Por otro lado, al compararla con secuencias de ADN generadas al azar, el ECM presenta valores más amplios, oscilando entre 2 y 11. 
 Estos valores más elevados sugieren una mayor discrepancia con el ADN aleatorio, lo cual podría indicar que ciertas regiones de la secuencia incógnita 
 no son completamente aleatorias y podrían tener algún grado de estructura similar al ADN funcional.
 Por otro lado, cuando la incógnita es una secuencia de ADN completamente generada al azar, se esperaría que el ECM refleje una baja similitud 
 con secuencias funcionales (grafico 3) y una alta similitud con secuencias generadas al azar (grafico 4). 
 Específicamente, el ECM para las secuencias de ADN aleatorias es de 0.3 ± 0.1, indicando una correspondencia cercana entre secuencias no codificantes. 
 En contraste, el ECM con secuencias funcionales es significativamente mayor, alcanzando 1.5 ± 0.5, lo que sugiere una marcada discrepancia 
 y refuerza la diferencia estructural y funcional entre secuencias aleatorias y codificantes. '''