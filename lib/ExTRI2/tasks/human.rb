require 'rbbt/sources/organism'
module ExTRI2

  input :name, :string, "Taxa name", "Human", required: true
  input :taxa, :string, "Taxa codes, separated by commas", "9606", required: true
  task :to_taxa => :tsv do |name, taxa|
    taxa = taxa.split(",").collect{|t| t.strip }

    ids ={}
    TSV.traverse "https://ftp.ncbi.nih.gov/pub/HomoloGene/last-archive/homologene.data", :type => :array do |line|
      hid, tax, id = line.split("\t")
      ids[hid]||=[tax, id]
      ids[hid]<<[tax, id]
    end
    index = TSV.setup({}, :key_field => "Entrez Gene ID", :fields =>["#{name} (Entrez Gene ID)"], :type => :single)
    ids.each do |hid, pairs|
      tax, translated = pairs.select{|p| taxa.include?(p[0]) }.first
      next if translated.nil?
      pairs.each do |tax,id|
        index[id]= translated
      end

    end
    index
  end

  task_alias :to_human, ExTRI2, :to_taxa, name: "Human", taxa: "9606"

  input :original, :file, "ExTRI2 main result", nil, :nofile => true
  input :organism, :string, "Organism code", "Hsa/feb2014"
  dep :to_taxa, name: :placeholder, taxa: :placeholder do |jobname,options|
    organism = options[:organism]
    taxa = Organism.entrez_taxids(organism).list * ","
    name = organism.split("/").first
    {inputs: options.merge(name: name, taxa: taxa)}
  end
  task :ExTRI2_taxa => :tsv do |original,organism|

    to_human = step(:to_taxa).load

    entrez2name = Organism.identifiers(organism).index :target => "Associated Gene Name", :fields =>["Entrez Gene ID"], persist: true
    parser = TSV::Parser.new original
    dumper = TSV::Dumper.new parser.options
    dumper.init
    TSV.traverse parser, into: dumper, bar: self.progress_bar("Translating gene orthologs") do |k,values|
      text, tf, tg, tf_ids, tg_ids, tf_off, tg_off, mutation, mutation_off, score, valid, mor = values
      next unless valid == 'Valid'

      tf_name = tf_ids.split(";").collect do |tf_id|
        tf_human = to_human[tf_id] || tf_id
        entrez2name[tf_human]
      end.compact *  ";"

      tg_name = tg_ids.split(";").collect do |tg_id|
        tg_human = to_human[tg_id] || tg_id
        entrez2name[tg_human]
      end.compact *  ";"

      if tf_name.empty? && ! tf_ids.empty?
        Log.debug "Missing #{tf_ids}"
      end

      if tg_name.empty? && ! tg_ids.empty?
        Log.debug "Missing #{tg_ids}"
      end

      next if tf_name.empty? or tg_name.empty?

      pmid, sentence = k.split(":").values_at 1, 4

      k = [pmid, sentence, tf_name, tg_name] * ":"
      [k, [text, tf, tg, tf_name, tg_name, tf_off, tg_off, mutation, mutation_off, score, valid, mor]]
    end
  end

  task_alias :ExTRI2_human, ExTRI2, :ExTRI2_taxa, organism: "Hsa/feb2014"

  dep :ExTRI2_human
  task :ExTRI2_clean => :tsv do
    dumper = TSV::Dumper.new(:key_field => "TRI", :fields => ["Transcription Factor (Associated Gene Name)", "Target Gene (Associated Gene Name)", "Interaction score", "Sentence"], :type => :list, :namespace => "Hsa/feb2014")
    dumper.init
    TSV.traverse step(:ExTRI2_human), into: dumper do |k,values|
      text, tf, tg, tf_ids, tg_ids, tf_off, tg_off, mutation, mutation_off, score, valid, mor = values

      sentence = text.sub('[TF]', tf).sub('[TG]', tg)

      [k, [tf, tg, score, sentence]]
    end
  end
end
