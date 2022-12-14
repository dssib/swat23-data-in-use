# Tutorial for our submission to SWAT4HCLS "data-in-use"

Here we attempt to prepare the relevant datasets enabling to predict essential proteins, starting from the published article [Predicting essential proteins by integrating orthology, gene expressions, and PPI networks](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0195410). In particular, we aim to compute the number of species where proteins in S. cerevisiae are conserved (also the number of orthologs per species) and the average expression breadth of the orthologs of S. cerevisiae. The rationale for this is: proteins that are highly conserved and whose orthologs are broadly expressed are more likely to be essential proteins.

The federated query below should be executed at the SPARQL endpoint of OMA https://sparql.omabrowser.org/sparql. However, the query does not work as expected, due to problems running aggregations within SERVICE blocks. An extended discussion is provided in our [Jupyter notebook](UniProt-OMA-Bgee.ipynb), including functional alternative queries and lessons learned from our attempts to reproduce the methodology in the paper.


```
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX dct: <http://purl.org/dc/terms/>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX ensembl: <http://rdf.ebi.ac.uk/resource/ensembl/>
PREFIX oma: <http://omabrowser.org/ontology/oma#>
PREFIX orth: <http://purl.org/net/orth#>
PREFIX sio: <http://semanticscience.org/resource/>
PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX void: <http://rdfs.org/ns/void#>
PREFIX lscr: <http://purl.org/lscr#>
PREFIX genex: <http://purl.org/genex#>

SELECT ?p ?taxon (count (DISTINCT ?orthologGeneEns) as ?count)  (avg (?expr_breadth) as ?mean_expr_breadth_per_ortho_in_taxon)
where {
      ?p a orth:Protein.
      ?p orth:organism/obo:RO_0002162 taxon:559292.
      ?cluster a orth:OrthologsCluster.
      ?cluster orth:hasHomologousMember ?node1.
      ?cluster orth:hasHomologousMember ?node2. 
      ?node2 orth:hasHomologousMember* ?ortholog. 
      ?node1 orth:hasHomologousMember* ?p.
      ?ortholog sio:SIO_010079 ?gene . #is encoded by
      ?gene lscr:xrefEnsemblGene ?orthologGeneEns .
      ?ortholog orth:organism/obo:RO_0002162 ?taxon.
       filter(?node1 != ?node2) 

# attempt to get expression breadth per ortholog from Bgee
# this will fail given that ?orthologGeneEns is replaced with a URI at query time and cannot be used as a group by criterion

  SERVICE <https://bgee.org/sparql/> {
	SELECT DISTINCT ?orthologGeneEns  (count (distinct ?anat) as ?expr_breadth) {
	   ?geneB a orth:Gene .
	   ?geneB lscr:xrefEnsemblGene ?orthologGeneEns .
                ?expr <http://purl.org/genex#hasSequenceUnit> ?geneB.
                ?expr a <http://purl.org/genex#Expression> .
                ?expr genex:hasConfidenceLevel obo:CIO_0000029 . # high confidence level
                ?expr genex:hasExpressionLevel ?exprLevel .
                FILTER (?exprLevel > 99) # highly expressed      
                ?expr genex:hasExpressionCondition ?cond .
                ?cond genex:hasAnatomicalEntity ?anat .
	} group by ?orthologGeneEns
  }
} group by ?p ?taxon
```