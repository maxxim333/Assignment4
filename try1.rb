=begin

Disclaimer: In this work I colaborated with Carolina Jiang (https://github.com/carol13jz)

I did not try to run this script using the ENTIRE genome of Arabidopsis because it was simply too large for my PC and I would have
to run it for more than overnight. But I tested it with a subset of the genome large enough to find a couple of orthologues. In principle
it should run smoothly for the rest of genome. If it doesnt, please let me know. I would like to have a chance to try to fix it.

The article I based my parameters on pointed to an e-value threshold of 1*10**-6 and overlap of over 50

=end

require 'bio'

#Here will be the file
orthologues = File.open('./orthologues.txt', 'w')
orthologues.puts "Ortholog pairs: \n" 



# Create databases
system "makeblastdb -in 'pep.fa' -dbtype prot -out prot_db"
system "makeblastdb -in 'TAIR10_seq_20110103_representative_gene_model_updated' -dbtype nucl -out nt_db"

at = File.open('TAIR10_seq_20110103_representative_gene_model_updated', "r")
pombae = File.open('pep.fa', "r")

# Using Bio Blast and creating local factory of alignments.

#This is how its gonna work:
# blastx searches protein databases using a translated nucleotide query, so we use it for protein sequence factory
# BLASTN programs search nucleotide databases using a nucleotide query, so we use it for nucleotide fasta factory

prot_factory = Bio::Blast.local('blastx', './prot_db')
nt_factory = Bio::Blast.local('tblastn','./nt_db')



#Here we transform the file in hash. It is needed to access the ID in form of hash key, when iterating over one of the files 
at = Bio::FlatFile.auto(at)
pombae = Bio::FlatFile.new(Bio::FastaFormat,pombae)

nt_hash = Hash.new 
at.each do |entry|
  nt_hash[(entry.entry_id).to_s] = (entry.seq).to_s 
end



#I will
#1.iterate over records of one databasae,
#2 store the ID of record Im currently iterating
#3. Check if there are hits, if there are, also check if they pass the parameters I applied.
#4. If it meets conditions, store the ID of the hit
#5. Retrieve the sequence of the hit and blast it over the other database
#6. Check conditions of hits as in 3
#7. Compare ID of hit with the innitial record. If they match, write them down.

i=0
pombae.each do |record| 
  prot_id = (record.definition.match(/(\w+\.\w+)/)).to_s #In the definition is the ID of the protein and has the format of \w+\.\w+
  blast1 = nt_factory.query(record)
  i+=1
  puts "Blasting protein number #{i} : #{prot_id}"
  
=begin
#This is the info I can retreive

  blast1.each do |hit|
    puts hit.evalue           # E-value
    puts hit.identity         # % identity
    puts hit.overlap          # length of overlapping region
    puts hit.query_id         # identifier of query sequence
    puts hit.query_def        # definition(comment line) of query sequence
    puts hit.query_len        # length of query sequence
    puts hit.query_seq        # sequence of homologous region
    puts hit.target_id        # identifier of hit sequence
    puts hit.target_def       # definition(comment line) of hit sequence
    puts hit.target_len       # length of hit sequence
    puts hit.target_seq       # hit of homologous region of hit sequence
    puts hit.query_start      # start position of homologous
                              # region in query sequence
    puts hit.query_end        # end position of homologous region
                              # in query sequence
    puts hit.target_start     # start position of homologous region
                              # in hit(target) sequence
    puts hit.target_end       # end position of homologous region
                              # in hit(target) sequence
    puts hit.lap_at           # array of above four numbers
  end
=end
  
  if blast1.hits[0] != nil and blast1.hits[0].evalue <= 10**-6 and blast1.hits[0].overlap.to_i >= 50
    nt_id = (blast1.hits[0].definition.match(/(\w+\.\w+)/)).to_s # store ID of the first hit
    puts "FOUND HIT in #{nt_id}! Checking reciprocal hit"
    sequence = nt_hash[nt_id]
    
    #Second Blast
    blast2 = prot_factory.query("#{sequence}") 
    if blast2.hits[0] != nil and blast2.hits[0].evalue <= 10**-6 and blast2.hits[0].overlap >= 50 
      hit = (blast2.hits[0].definition.match(/(\w+\.\w+)/)).to_s
      if prot_id == hit 
        orthologues.puts "#{prot_id}\t#{nt_id}" 
        puts "#{prot_id}  orthologue to #{nt_id}"
      end
    end
  end
end

=begin
Bonus question:
The next step for identifying orthology would be focusing on phylogenic reconstruction. By analyzing the phylogenetic trees via
tree reconciliation for example, we would derive a collection of fine-graied prediction of all orthology relationship. A crucia parameter
for ths should be defined: how far the evolutionary tree one should go.

Another approach would be to check the positional conservation of genes in the genomes (synteny). Gene and gene clusters that conserve
their relative position in the genome, are more likely to be conserved in function

And of course, we could directly query the function of the protein, either by checking databases or perform a wet-lab experiment.
=end
