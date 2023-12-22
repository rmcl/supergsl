from typing import Dict, List
import requests
from supergsl.core.types.part import Part
from supergsl.core.provider import ProviderConfig, SuperGSLProvider
from supergsl.core.types.protein import Protein
from supergsl.core.exception import ProviderError


class UniprotProvider(SuperGSLProvider):
    """A provider for parts from UniProt."""

    def __init__(self, name : str, config : ProviderConfig):
            self._provider_name = name
            self.sequence_store = config.sequence_store


    def get_nested_key(self, data : dict, key : str) -> dict:
        """Get a nested key from a dictionary."""
        keys = key.split('.')
        current = data
        for key in keys:
            try:
                current = current[key]
                if isinstance(current, list):
                    if len(current) == 0:
                        return None
                    current = current[0]

            except KeyError:
                return None
        return current

    def get_details_from_uniprot_result(self, result : dict) -> dict:
        """Extract details from a UniProt result."""
        function_comments = [
            comment
            for comment in result['comments']
            if comment['commentType'] == 'FUNCTION'
        ]

        if len(function_comments) > 0:
            description = function_comments['texts'][0]['value']
        else:
            description = self.get_nested_key(
                result, 'proteinDescription.submissionNames.fullName.value')


        short_names = self.get_nested_key(
            result, 'proteinDescription.recommendedName.shortNames')

        alternative_names = []
        if short_names:
            alternative_names = [
                name['value']
                for name in short_names
            ]

        name = self.get_nested_key(
            result, 'proteinDescription.recommendedName.fullName.value')
        if not name:
            name = self.get_nested_key(
                result, 'genes.geneName.value')

        return {
            'name': name,
            'alternative_names': alternative_names,
            'description': description
        }

    def get_protein_data(self, uniprot_id : str) -> dict:
        """Get the protein sequence and metadata for the given UniProt identifier."""
        base_url = "https://www.uniprot.org/uniprot/"
        response = requests.get(base_url + uniprot_id + ".json")
        if response.status_code != 200:
            raise Exception("Error: Unable to retrieve protein sequence.")

        result = response.json()

        print(result)
        details = self.get_details_from_uniprot_result(result)

        return {
            'identifier': result['primaryAccession'],
            'name': details['name'],
            'sequence': result['sequence']['value'],
            'alternative_names': details['alternative_names'],
            'part_type': 'protein',
            'description': details['description']
        }

    def search_metadata(self, query : str) -> dict:
        """Search for a protein by name or UniProt accession."""
        query_param=f'protein_name:{query} OR accession:{query}'

        base_url = "https://rest.uniprot.org/uniprotkb/search?"
        response = requests.get(
            base_url,
            params={
                'query': query_param
            })

        if response.status_code != 200:
            raise ProviderError("Error: Unable to perform uniprot search.")

        entries = response.json()['results']

        results = []
        for entry in entries:
            details = self.get_details_from_uniprot_result(entry)
            result = {
                'identifier': entry['primaryAccession'],
                'source': 'uniprot',
                'name': details['name'],
                'description': details['description'],
                'part_type': 'protein',
                'source_details': entry
            }
            results.append(result)

        return results

    def search(self, query : str) -> List[Protein]:
        raise NotImplementedError()

    def get(self, identifier : str) -> Protein:
        """Get the protein with the given UniProt identifier."""
        protein_data = self.get_protein_data(identifier)
        alt_names = protein_data['alternative_names'].copy()
        alt_names.insert(0, protein_data['name'])

        seq_entry = self.sequence_store.add_from_reference(protein_data['sequence'])
        return Protein(
            identifier=protein_data['identifier'],
            sequence_entry=seq_entry,
            alternative_names=alt_names,
            description=protein_data['description']
        )
