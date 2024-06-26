require 'rbbt-util'
require 'rbbt/document/annotation'
require 'rbbt/ner/token_trieNER'
require 'rbbt/ner/rnorm'
require 'rbbt/util/python'

module ExTRI2
  def self.is_tf?(gene)
    return false if gene.code.nil?
    @@codes ||= Rbbt.data.tf_entrez_code.list
    (@@codes & gene.code.split(";")).any?
  end
end

Document.define :pubtator => :single do |dir|
  pmid = self.code
  dir[pmid].read.split("\n").collect do |line|
    _id, start, eend, literal, type, code = line.split("\t")
    code.gsub!('|','-') if code
    NamedEntity.setup(literal, offset: start, code: code, entity_type: type)
  end
end

Document.define :genes => :single do |dir|
  pubtator(dir).select{|e| e.entity_type == "Gene" }
end

Document.define :mutations => :single do |dir|
  pubtator(dir).select{|e| e.entity_type == "ProteinMutation" || e.entity_type == "DNAMutation" || e.entity_type == "SNP" }
end

Document.define :tfs => :single do |dir|
  genes(dir).select{|g| ExTRI2.is_tf?(g) }
end

RbbtPython.add_path Rbbt.python.find(:lib)

Document.define_multiple :sentences do |list|
  list = self if Array === self
  list_sentences = RbbtPython.run :spacy_splitter do
    spacy_splitter.split_texts(list)
  end
  list_sentences = RbbtPython.list2ruby(list_sentences)
  list_sentences.collect do |sentences|
    sentences.collect do |text,start,eend|
      Segment.setup(text, offset: start)
    end
  end
end

#Document.persist :sentences, :annotations, annotation_repo: ExTRI2::ANNOTATION_REPO
#Document.persist :genes, :annotations, annotation_repo: ExTRI2::ANNOTATION_REPO
#Document.persist :tfs, :annotations, annotation_repo: ExTRI2::ANNOTATION_REPO
Document.persist :sentences, :annotations
Document.persist :genes, :memory
Document.persist :tfs, :memory
Document.persist :mutations, :memory

if __FILE__ == $0
  require 'rbbt/workflow'
  require 'rbbt/document/corpus'
  require 'rbbt/document/corpus/pubmed'
  Workflow.require_workflow "ExTRI2"
  pmids = "19522013|20861254".split("|")
  documents = ExTRI2::CORPUS.add_pmid(pmids)
  iii documents.tf
  iii documents.first.tf
end
