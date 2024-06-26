require 'rbbt/sources/organism'
module ExTRI2

  task :to_human => :tsv do

    ids ={}
    TSV.traverse "https://ftp.ncbi.nih.gov/pub/HomoloGene/last-archive/homologene.data", :type => :array do |line|
      hid, tax, id = line.split("\t")
      ids[hid]||=[tax, id]
      ids[hid]<<[tax, id]
    end
    index = TSV.setup({}, :key_field => "Entrez Gene ID", :fields =>["Human(Entrez Gene ID)"], :type => :single)
    ids.each do |hid, pairs|
      tax, human = pairs.select{|p| p[0]== "9606"}.first
      next if human.nil?
      pairs.each do |tax,id|
        index[id]= human
      end

    end
    index
  end


  dep :to_human
  input :original, :file, "ExTRI2 main result", nil, :nofile => true
  task :ExTRI2_human => :tsv do |original|

    to_human = step(:to_human).load

    entrez2name = Organism.identifiers("Hsa/feb2014").index :target => "Associated Gene Name", :fields =>["Entrez Gene ID"], persist: true
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
