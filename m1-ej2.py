import random
import requests
from collections import Counter
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
from Bio import ExPASy
from Bio import SwissProt
import inspect
print('hola')

# 2a) Obtenga la reversa complementaria y compara las estadísticas de la misma con la original

length=10000
secuencia=generate_random_DNA_sequence(length)
reversa_complementaria=Seq(secuencia)
reversa_complementaria=str(reversa_complementaria.reverse_complement())

#process_sequence(secuencia,'secuencia al azar')
#process_sequence(reversa_complementaria,'reversa complementaria')


# 2b) Traduzca la secuencia y obtenga la distribución de largo de los ORF (compare con los resultados obtenidos anteriormente). 

length=10000
secuencia=generate_random_DNA_sequence(length)
secuencia_traducida=translate(secuencia)
my_seq=Seq(secuencia)
my_seq_trad=str(secuencia.translate())

#process_sequence(secuencia_traducida,'secuencia traducida con nuestro programa')
#process_sequence(my_seq_trad,'secuencia traducida con biopython')
#FrecAaJuntos(secuencia_traducida,my_seq_trad)

#hacer una process_sequence que compare!
def graficoFrecBases_comparacion(seq1,seq2):
  frec_bases1,cuenta=frecBases(seq1)
  frec_bases2,cuenta=frecBases(seq2)
  fig = go.Figure()
  fig.add_trace(go.Bar(
    x_values = list(cuenta.keys()),
    y_values=frec_bases1,
    name='Original Sequence Base Frequencies',
    marker_color='blue'))
  
  fig.add_trace(go.Bar(
    x_values = list(cuenta.keys()),
    y_values=frec_bases2,
    name='Reverse Complement Base Frequencies',
    marker_color='lightblue'))

  fig.update_layout(
    title='Base Frequencies for Original and Reverse Complement DNA Sequences',
    xaxis_title='Base',
    yaxis_title='Frequency',
    barmode='group')
  fig.show()

def graficoFredCod_comparacion(seq1,seq2):
  x_values_cod1,y_values_cod1=frecCod(seq1)
  x_values_cod2,y_values_cod2=frecCod(seq2)
  fig = go.Figure()
  fig.add_trace(go.Bar(
    x=x_values_cod1,
    y=y_values_cod1,
    name1='Original Sequence Codon Frequencies',
    marker_color='blue'))
  fig.add_trace(go.Bar(
    x=x_values_cod2,
    y=y_values_cod2,
    name2='Reverse Complement Codon Frequencies',
    marker_color='lightblue'))

  fig.update_layout( 
    title='Codon Frequencies for Original and Reverse Complement DNA Sequences',
    xaxis_title='Codon',
    yaxis_title='Frequency',
    barmode='group')
  
  fig.show()

def process_sequence_comparacion(seq_original,seq_reversacomplem):
  if all(char in nucleotides for char in seq_original):
    print("Processing nucleotide sequence...")
    graficoFrecBases_comparacion(seq_original,seq_reversacomplem)
    graficoFredCod_comparacion(seq_original,seq_reversacomplem)
    FrecAaJuntos(seq_original,seq_reversacomplem)
    graficoOrfs_comparacion(seq_original,seq_reversacomplem)

#FALTA DISTRIBUCIÓN DE ORFS

def graficoOrfs_comparacion(seq1,seq2):
  # seq 1
  y_values1=count_orf_lengths(seq1)
  media=np.mean(y_values1)
  fig = go.Figure()
  fig.add_trace(go.Histogram(
    y=y_values1,
    name='Original Sequence ORF Lengths',
    marker_color='blue',
    opacity=0.75))
  fig.add_vline(x=media, line_dash="dash", line_color="violet", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
  
  # seq 2 
  y_values2=count_orf_lengths(seq2)
  media=np.mean(y_values2)
  fig.add_trace(go.Histogram(
    y=y_values2,
    name='Reverse Complement ORF Lengths',
    marker_color='violet',
    opacity=0.75))
  fig.add_vline(x=media, line_dash="dash", line_color="red", annotation_text=f"Media = {media:.2f}", annotation_position="top right")

  # juntos
  fig.update_layout(
    title='Distribution of ORF Lengths Between Stop Codons',
    xaxis_title='ORF Length',
    yaxis_title='Frequency',
    barmode='overlay')

  fig.show()


# 2b) elija un registros de geneBank bajelo en formato fasta y gb, levantelos con
# Biopyhton y explore sus atributos (record, record.id, record.name, record.description, etc.)
# es capaz de obtener “solo” la secuencia y convertirla en un “string” para manipularla?

# URL del archivo GenBank
archivo_gb = 'https://raw.githubusercontent.com/Nehuenpg/Fasta-Sequences/main/sequence.gb'

# Descarga el archivo desde la URL
response = requests.get(archivo_gb)
response.raise_for_status()  # Lanza un error si la descarga falla

# Guarda el contenido del archivo en un objeto StringIO
from io import StringIO
file_content = StringIO(response.text)

# Lee el archivo GenBank desde el contenido descargado
record = SeqIO.read(file_content, "genbank")

# Imprime los atributos del registro
print(f"ID: {record.id}")
print(f"Name: {record.name}")
print(f"Description: {record.description}")
print(f"Sequence: {record.seq[:100]}...")  # Muestra los primeros 100 nucleótidos



# 2c) Combine lo aprendido en la clase previas (for / while loops) y/o la información
# disponible en biopython para cargar un conjunto de secuencias de ADN y/o proteínas. (Puede
# bajarlos manualmente a disco y realizar código para cargarlos y/o hacer el código que los baje
# directamente de internet). Obtenga solo las secuencias como cadenas de caracteres y guardelas
# en una lista tal que cada elemento de la lista sea una de las secuencias.

# URL del archivo GenBank
archivo_gb = 'https://raw.githubusercontent.com/Nehuenpg/Fasta-Sequences/main/sequence.gb'

# Descarga el archivo desde la URL
response = requests.get(archivo_gb)
response.raise_for_status()  # Lanza un error si la descarga falla

# Guarda el contenido del archivo en un objeto StringIO
from io import StringIO
file_content = StringIO(response.text)

# Lee el archivo GenBank desde el contenido descargado
record = SeqIO.read(file_content, "genbank")

# Imprime los atributos del registro
print(f"ID: {record.id}")
print(f"Name: {record.name}")
print(f"Description: {record.description}")
print(f"Sequence: {record.seq[:100]}...")  # Muestra los primeros 100 nucleótidos


# Example ussage


def abrir_archivo(url): #abre el archivo, genera una lista con cada secuencia y devuelve un str con todas las secuencias seguidas.
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

print(abrir_archivo(url))

