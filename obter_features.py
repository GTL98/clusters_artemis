def obter_features(genbank: object, id_gene: str) -> tuple:
    """
    Função responsável por retornar algumas features do arquivo GenBank.
    :param genbank: Informações do arquivo GenBank.
    :param id_gene: Identificação do gene.
    :return: Tupla com algumas das features.
    """
    # --- Iterar sobre cada feature do GenBank --- #
    for feature in genbank.features:
        # --- Verificar se é uma região codificante e qual é o locus --- #
        if feature.type == 'CDS' and id_gene in feature.qualifiers['locus_tag']:
            try:
                produto = feature.qualifiers['product']
                id_proteina = feature.qualifiers['protein_id']
                localizaco = (str(feature.location)
                              .replace('[', '')
                              .replace(']', '')
                              .replace('(+)', '')
                              .replace('(-)', '')
                              .replace(':', '-'))
            except KeyError:
                pass
            else:
                return produto, id_proteina, localizaco
