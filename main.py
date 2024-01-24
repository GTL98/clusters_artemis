# --- Importar as bibliotecas --- #
from Bio import SeqIO
from json import load
import streamlit as st
from obter_clusters import obter_clusters

# --- Configuração base da página --- #
st.set_page_config(page_title='Clusters Artemis', layout='wide')

# --- Carregar os dados dos genes --- #
with open('genes_bifido_pullorum.json', 'r') as doc:
    dados = load(doc)

# --- Obter as proteínas --- #
proteinas = sorted(list(dados.keys()))

# --- Escolher a proteína --- #
proteina = st.selectbox(
    label='Escolha uma opção',
    options=proteinas,
    placeholder='Selecione uma proteína',
    index=None
)

# --- Obter os IDs dos genes --- #
if proteina is not None:
    indices = list(dados[proteina])[1:]

    # Criar uma caixa de seleção com os IDs
    indice = st.selectbox(
        label='Escolha uma opção',
        options=indices,
        placeholder='Selecione um ID',
        index=None
    )

    # Obter as informações dos genes do cluster
    if indice is not None:
        # Complementar
        complementar_temp = dados[proteina][indice]['complementar']
        if complementar_temp is True:
            complementar = 'Complementar'
        else:
            complementar = 'Não complementar'

        # Obter o cluster
        registro = SeqIO.read('Bifidobacterium_pullorum.gb', 'genbank')
        dados_clusters = obter_clusters(registro, proteina, indice, dados)

        # Obter as chaves dos genes anteriores e posteriores
        chaves = list(dados_clusters[indice].keys())
        id_genes = []
        dic_id_ncbi = {}
        localizacoes = []
        if 'antes' in chaves:
            # Obter os IDs dos genes do cluster
            id_genes_antes = list(dados_clusters[indice]['antes'].keys())
            for id_gene in id_genes_antes:
                if id_gene not in id_genes:
                    id_genes.append(id_gene)

            # Obter o ID do NCBI de cada gene
            for id_gene in id_genes:
                try:
                    id_ncbi = dados_clusters[indice]['antes'][id_gene]['id_ncbi']
                except KeyError:
                    pass
                else:
                    dic_id_ncbi[id_gene] = id_ncbi

            # Obter a localização dos genes
            for id_gene in id_genes:
                try:
                    localizacao_1 = int(dados_clusters[indice]['antes'][id_gene]['localizacao'].split('-')[0])
                    localizacao_2 = int(dados_clusters[indice]['antes'][id_gene]['localizacao'].split('-')[1])
                except KeyError:
                    pass
                else:
                    localizacoes.append(localizacao_1)
                    localizacoes.append(localizacao_2)

        if 'depois' in chaves:
            # Obter os IDs dos genes do cluster
            id_genes_antes = list(dados_clusters[indice]['depois'].keys())
            for id_gene in id_genes_antes:
                if id_gene not in id_genes:
                    id_genes.append(id_gene)

            # Obter o ID do NCBI de cada gene
            for id_gene in id_genes:
                try:
                    id_ncbi = dados_clusters[indice]['depois'][id_gene]['id_ncbi']
                except KeyError:
                    pass
                else:
                    dic_id_ncbi[id_gene] = id_ncbi

            # Obter a localização dos genes
            for id_gene in id_genes:
                try:
                    localizacao_1 = int(dados_clusters[indice]['depois'][id_gene]['localizacao'].split('-')[0])
                    localizacao_2 = int(dados_clusters[indice]['depois'][id_gene]['localizacao'].split('-')[1])
                except KeyError:
                    pass
                else:
                    localizacoes.append(localizacao_1)
                    localizacoes.append(localizacao_2)

        # Organizar os índices pelo número
        id_genes.sort()

        # Organizar as localização pelo número
        localizacoes.sort()

        # Obter a ordem do cluster
        ordem = ''
        for id_gene in id_genes:
            if id_gene == indice:
                ordem += f':red[{id_gene}] --> '
            else:
                ordem += f'{id_gene} --> '

        # Verificar se possui genes antes ou depois
        if ordem != '':
            st.subheader(ordem[:-5], divider='grey')
        else:
            st.subheader('Não há genes anteriores e posteriores.', divider='grey')

# --- Escrever as informações importantes --- #
try:
    complementar = complementar_temp
except NameError:
    pass
else:
    st.subheader(f'Complementar: {complementar}', divider='grey')
    if ordem != '':
        st.subheader(f'Localização: {localizacoes[0]} - {localizacoes[-1]}', divider='grey')
    else:
        loc_alternativa_1 = dados[proteina][indice]['posicao'].split('-')[0]
        loc_alternativa_2 = dados[proteina][indice]['posicao'].split('-')[1]
        st.subheader(f'Localização: {loc_alternativa_1} - {loc_alternativa_2}', divider='grey')
    for id_gene, id_ncbi in dic_id_ncbi.items():
        st.subheader(f'{id_gene}: {id_ncbi}')