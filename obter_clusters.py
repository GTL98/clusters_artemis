from obter_features import obter_features


def obter_clusters(genbank: object, funcao: str, id_gene: str, dados: dict) -> dict:
    """
    Função responsável por retornar os clusters com as features de cada gene do cluster.
    :param genbank: Informações do arquivo GenBank.
    :param funcao: Função probiótica.
    :param id_gene: Identificação do gene.
    :param dados: Dicionário com os dados JSON do gene.
    :return: Dicionário com as features de cada gene do cluster.
    """
    # --- Obter a identificação a fita de leitura --- #
    complementar = dados[funcao][id_gene]['complementar']
    dic_clusters = {id_gene: {}}

    # --- Se a complmentar for "False", os genes posteriores se mantêm na mesma ordem numérica --- #
    if complementar is False:
        # Obter quantos genes vêm antes
        antes = int(dados[funcao][id_gene]['antes'])

        # Obter os genes anteriores
        if antes > 0:
            # Adicionar a chave dos genes anteriores
            dic_clusters[id_gene]['antes'] = {}

            # Separar o id do gene para obter somente o número
            id_bruto = int(id_gene.split('_')[1])

            # Colocar o gene probiótico no dicionário
            features = obter_features(genbank, id_gene)
            dic_clusters[id_gene]['antes'][id_gene] = {}
            dic_clusters[id_gene]['antes'][id_gene]['produto'] = features[0][0]
            dic_clusters[id_gene]['antes'][id_gene]['id_ncbi'] = features[1][0]
            dic_clusters[id_gene]['antes'][id_gene]['localizacao'] = features[2]

            # Obter os outros genes anteriores
            for i in range(antes):
                deslocamento = (i+1) * 5
                id_gene_temp = id_bruto - deslocamento

                # Recriar a string do "locus_tag"
                if id_bruto < 1000:
                    id_gene_anterior = f'ESN35_00{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['antes'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['localizacao'] = features[2]

                elif 1000 < id_bruto < 9999:
                    id_gene_anterior = f'ESN35_0{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['antes'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['localizacao'] = features[2]

                else:
                    id_gene_anterior = f'ESN35_{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['antes'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['localizacao'] = features[2]

        # Obter quantos genes vêm depois
        depois = int(dados[funcao][id_gene]['depois'])

        # Obter os genes posteriores
        if depois > 0:
            # Adicionar a chave dos genes anteriores
            dic_clusters[id_gene]['depois'] = {}

            # Separar o id do gene para obter somente o número
            id_bruto = int(id_gene.split('_')[1])

            # Colocar o gene probiótico no dicionário
            features = obter_features(genbank, id_gene)
            dic_clusters[id_gene]['depois'][id_gene] = {}
            dic_clusters[id_gene]['depois'][id_gene]['produto'] = features[0][0]
            dic_clusters[id_gene]['depois'][id_gene]['id_ncbi'] = features[1][0]
            dic_clusters[id_gene]['depois'][id_gene]['localizacao'] = features[2]

            # Obter os outros genes anteriores
            for i in range(depois):
                deslocamento = (i + 1) * 5
                id_gene_temp = id_bruto + deslocamento

                # Recriar a string do "locus_tag"
                if id_bruto < 1000:
                    id_gene_anterior = f'ESN35_00{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['depois'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['localizacao'] = features[2]

                elif 1000 < id_bruto < 9999:
                    id_gene_anterior = f'ESN35_0{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['depois'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['localizacao'] = features[2]

                else:
                    id_gene_anterior = f'ESN35_{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['depois'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['localizacao'] = features[2]

    # --- Se a complmentar for "True", os genes posteriores NÃO se mantêm na ordem numérica --- #
    if complementar is True:
        # Obter quantos genes vêm antes
        antes = int(dados[funcao][id_gene]['antes'])

        # Obter os genes anteriores
        if antes > 0:
            # Adicionar a chave dos genes anteriores
            dic_clusters[id_gene]['antes'] = {}

            # Separar o id do gene para obter somente o número
            id_bruto = int(id_gene.split('_')[1])

            # Colocar o gene probiótico no dicionário
            features = obter_features(genbank, id_gene)
            dic_clusters[id_gene]['antes'][id_gene] = {}
            dic_clusters[id_gene]['antes'][id_gene]['produto'] = features[0][0]
            dic_clusters[id_gene]['antes'][id_gene]['id_ncbi'] = features[1][0]
            dic_clusters[id_gene]['antes'][id_gene]['localizacao'] = features[2]

            # Obter os outros genes anteriores
            for i in range(antes):
                deslocamento = (i + 1) * 5
                id_gene_temp = id_bruto + deslocamento

                # Recriar a string do "locus_tag"
                if id_bruto < 1000:
                    id_gene_anterior = f'ESN35_00{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['antes'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['localizacao'] = features[2]

                elif 1000 < id_bruto < 9999:
                    id_gene_anterior = f'ESN35_0{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['antes'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['localizacao'] = features[2]

                else:
                    id_gene_anterior = f'ESN35_{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['antes'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['antes'][id_gene_anterior]['localizacao'] = features[2]

        # Obter quantos genes vêm depois
        depois = int(dados[funcao][id_gene]['depois'])

        # Obter os genes posteriores
        if depois > 0:
            # Adicionar a chave dos genes anteriores
            dic_clusters[id_gene]['depois'] = {}

            # Separar o id do gene para obter somente o número
            id_bruto = int(id_gene.split('_')[1])

            # Colocar o gene probiótico no dicionário
            features = obter_features(genbank, id_gene)
            dic_clusters[id_gene]['depois'][id_gene] = {}
            dic_clusters[id_gene]['depois'][id_gene]['produto'] = features[0][0]
            dic_clusters[id_gene]['depois'][id_gene]['id_ncbi'] = features[1][0]
            dic_clusters[id_gene]['depois'][id_gene]['localizacao'] = features[2]

            # Obter os outros genes anteriores
            for i in range(depois):
                deslocamento = (i + 1) * 5
                id_gene_temp = id_bruto - deslocamento

                # Recriar a string do "locus_tag"
                if id_bruto < 1000:
                    id_gene_anterior = f'ESN35_00{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['depois'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['localizacao'] = features[2]

                elif 1000 < id_bruto < 9999:
                    id_gene_anterior = f'ESN35_0{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['depois'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['localizacao'] = features[2]

                else:
                    id_gene_anterior = f'ESN35_{id_gene_temp}'

                    # Adicionar as features dos genes
                    features = obter_features(genbank, id_gene_anterior)
                    dic_clusters[id_gene]['depois'][id_gene_anterior] = {}
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['produto'] = features[0][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['id_ncbi'] = features[1][0]
                    dic_clusters[id_gene]['depois'][id_gene_anterior]['localizacao'] = features[2]

    return dic_clusters