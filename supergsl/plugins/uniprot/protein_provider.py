from typing import Dict, List
import requests
from supergsl.core.types.part import Part
from supergsl.core.provider import ProviderConfig, SuperGSLProvider
from supergsl.core.types.builtin import AminoAcidSequence
from supergsl.core.exception import ProviderError


class UniprotProvider(SuperGSLProvider):
    """A provider for parts from UniProt."""

    def __init__(self, name : str, config : ProviderConfig):
            self._provider_name = name
            self._cached_parts: Dict[str, Part] = {}
            self.sequence_store = config.sequence_store

            settings = config.settings

    def get_details_from_uniprot_result(self, result : dict) -> dict:
        """Extract details from a UniProt result."""
        function_comment = [
            comment
            for comment in result['comments']
            if comment['commentType'] == 'FUNCTION'
        ][0]

        try:
            short_names = result['proteinDescription']['recommendedName']['shortNames']
        except KeyError:
            short_names = []

        alternative_names = [
            name['value']
            for name in short_names
        ]

        return {
            'name': result['proteinDescription']['recommendedName']['fullName']['value'],
            'alternative_names': alternative_names,
            'description': function_comment['texts'][0]['value']
        }

    def get_protein_data(self, uniprot_id : str) -> dict:
        """Get the protein sequence and metadata for the given UniProt identifier."""
        base_url = "https://www.uniprot.org/uniprot/"
        response = requests.get(base_url + uniprot_id + ".json")
        if response.status_code != 200:
            raise Exception("Error: Unable to retrieve protein sequence.")

        result = response.json()

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

    def search(self, query : str) -> List[AminoAcidSequence]:
        raise NotImplementedError()

    def get_part(self, identifier : str) -> List[AminoAcidSequence]:
        raise NotImplementedError()
