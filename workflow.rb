require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/document'
require 'rbbt/document/corpus'
#require 'rbbt/util/python'
#require 'rbbt/vector/model/huggingface'


Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/ExTRI2'

module ExTRI2
  extend Workflow

  CORPUS = Document::Corpus.setup(Rbbt.var.ExTRI2.corpus.find)
  ANNOTATIONS = Rbbt.var.ExTRI2.annotations.find
  ANNOTATION_REPO = Rbbt.var.ExTRI2.annotation_repo.find

  CORPUS.close

  helper :corpus do
    CORPUS
  end

  input :pubtator_file, :file
  task :load_documents => :array do |file|
    Open.read(file).split("\n\n").collect do |chunk|
      title_pre, abstract_pre, *entities = chunk.split("\n")
      pmid,_sep, title = title_pre.partition("|t|")
      pmid,_sep, abstract = abstract_pre.partition("|a|")

      abstract.sub!('|'," ")

      text = [title, abstract] * "\n"

      Document.setup(text, "PMID", pmid, :pubtator_title_and_abstract)

      file(pmid).write(entities * "\n")

      corpus.add_document(text)
      text.docid
    end
  end

  dep :load_documents
  task :tri_candidates => :tsv do
    docids = DocID.setup(step(:load_documents).load, corpus: corpus)
    docids.extend AnnotatedArray

    dir = step(:load_documents).files_dir 

    log :documents, "Loading documents"
    documents = docids.document

    log :sentences, "Loading document sentences"
    documents.sentences
    
    log :entities, "Loading document entities"
    documents.genes(dir)
    documents.tfs(dir)
    documents.mutations(dir)

    tsv = TSV::Dumper.new(:key_field => "SentenceID", :fields => ["Text", "TF", "Gene", "TF Id", "IG Id", "TF offset", "Gene offset", "Mutated Genes", "Mutation offsets"], :type => :list)
    tsv.init
    TSV.traverse docids, :bar => self.progress_bar("Procesing documents for TRI candidates"), into: tsv do |docid|
      doc = docid.document
      genes = doc.genes(dir)
      tfs = doc.tfs(dir)
      mutations = doc.mutations(dir)
      sentences = doc.sentences

      res = []
      sentences.each_with_index do |sentence,i|
        sentence_genes = sentence.overlaps(genes)
        sentence_tfs = sentence.overlaps(tfs)
        sentence_mutations = sentence.overlaps(mutations)

        mutated_genes = sentence_mutations.collect do |m| 
          m.code.split(";").
            select{|e| e.include?("CorrespondingGene") }.
            collect{|e| e.split(":").last }
        end.flatten

        next if sentence_genes.length == 0
        next if sentence_tfs.length == 0

        sentence_tfs.each do |tf|
          sentence_genes.each_with_index do |gene|
            next if tf.overlaps?(gene)
            Transformed.with_transform(sentence, "\n", "  ") do
              Transformed.with_transform(sentence, tf, "[TF]") do
                Transformed.with_transform(sentence, gene, "[TG]") do
                  raise if sentence.include?("\n")
                  id = [docid, i, tf, gene, tf.offset, gene.offset] * ":"
                  res << [id, [sentence.dup.gsub("\n", " ").strip, tf, gene, tf.code, gene.code,tf.offset.to_s,gene.offset.to_s, mutated_genes * ";", sentence_mutations.collect{|m| m.offset } * ";"]]
                end
              end
            end
          end
        end
      end
      res.extend MultipleResult
      res
    end
  end

  dep :tri_candidates
  input :tri_model, :string, "TRI model to load", "TRI_model"
  task :tri_sentences => :tsv do |tri_model|

    tsv = step(:tri_candidates).load

    tri_model = Rbbt.models[tri_model].find unless File.exist?(tri_model)

    model = HuggingfaceModel.new 'SequenceClassification', tri_model, nil,
      :tokenizer_args => {:model_max_length => 512, :truncation => true},
      :return_logits => true

    model.extract_features do |_,feature_list|
      feature_list.collect do |text,tf,tg|
        text.sub("[TF]", "<TF>#{tf}</TF>").
          sub("[TG]", "<TG>#{tg}</TG>")
      end
    end

    model.init

    predictions = model.eval_list tsv.slice(["Text", "TF", "Gene"]).values

    tsv.add_field "Valid score" do 
      non_valid, valid = predictions.shift
      begin
        Misc.softmax([valid, non_valid]).first
      rescue
        0
      end
    end

    tsv.add_field "Valid" do |k,values|
      values.last > 0.5 ? "Valid" : "Non valid"
    end

    tsv
  end

  dep :tri_sentences
  input :mor_model, :string, "TRI model to load", "MoR_model"
  task :tri_MoR => :tsv do |mor_model|

    tsv = step(:tri_sentences).load
    mor_model = Rbbt.models[mor_model].find unless File.exist?(mor_model)
    model = HuggingfaceModel.new 'SequenceClassification', mor_model, nil,
      :tokenizer_args => {:model_max_length => 512, :truncation => true},
      :class_labels => %w(UNDEFINED ACTIVATION REPRESION),
      :return_logits => true

    model.extract_features do |_,feature_list|
      feature_list.collect do |text,tf,tg|
        text.sub("[TF]", "<TF>#{tf}</TF>").
          sub("[TG]", "<TG>#{tg}</TG>")
      end
    end

    model.init

    predictions = model.eval_list tsv.slice(["Text", "TF", "Gene"]).values

    tsv.add_field "MoR scores" do 
      preds = predictions.shift
      begin
        Misc.softmax(preds) * ";"
      rescue
        ""
      end
    end

    tsv.add_field "MoR" do |k,values|
      scores = values.last.split(";").collect{|v| v.to_f }
      %w(UNDEFINED ACTIVATION REPRESION)[scores.index(scores.max)]
    end

    tsv
  end

  dep :tri_MoR, pubtator_file: :placeholder, compute: :produce, canfail: true  do |jobname,options|
    Rbbt.data.pubtator.glob("*.pubtator").collect do |file|
      {task: :tri_MoR, inputs: options.merge(:pubtator_file => file)}
    end
  end
  task :ExTRI2 => :tsv do
    TSV.concat_streams(dependencies.select{|dep| dep.done? })
  end

end
require 'ExTRI2/entities'


require 'ExTRI2/tasks/human.rb'
require 'ExTRI2/tasks/collec_tri2.rb'
require 'ExTRI2/tasks/knocktf.rb'
require 'ExTRI2/tasks/tmp.rb'

#require 'rbbt/knowledge_base/ExTRI2'
#require 'rbbt/entity/ExTRI2'

Workflow.main = ExTRI2
