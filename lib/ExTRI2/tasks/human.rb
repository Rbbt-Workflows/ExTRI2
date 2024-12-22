require 'rbbt/sources/organism'
module ExTRI2

  input :name, :string, "Taxa name", "Human"
  input :taxa, :string, "Taxa codes, separated by commas", "9606"
  task :to_taxa_homologene => :tsv do |name, taxa|
    taxa = taxa.split(",").collect{|t| t.strip }

    ids ={}
    TSV.traverse "https://ftp.ncbi.nih.gov/pub/HomoloGene/last-archive/homologene.data", :type => :array do |line|
      hid, tax, id = line.split("\t")
      ids[hid]||=[]
      ids[hid]<<[tax, id]
    end
    index = TSV.setup({}, :key_field => "Entrez Gene ID", :fields =>["#{name} (Entrez Gene ID)"], :type => :single)
    ids.each do |hid, pairs|
      tax, translated = pairs.select{|p| taxa.include?(p[0]) }.first
      next if translated.nil?
      pairs.each do |tax,id|
        next if id == translated
        index[id]= translated
      end

    end
    index
  end

  input :name, :string, "Taxa name", "Human"
  input :taxa, :string, "Taxa codes, separated by commas", "9606"
  task :to_taxa_hcop => :tsv do |name, taxa|
    taxa = taxa.split(",").collect{|t| t.strip }

    url = "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/genenames/hcop/human_all_hcop_sixteen_column.txt.gz"
    ids ={}
    TSV.traverse url, :type => :array do |line|
      tax, human_entrez, hid, id = line.split("\t").values_at 0, 1, 3, 8
      ids[hid]||=[]
      ids[hid]<<["9606", human_entrez]
      ids[hid]<<[tax, id]
    end
    index = TSV.setup({}, :key_field => "Entrez Gene ID", :fields =>["#{name} (Entrez Gene ID)"], :type => :single)
    ids.each do |hid, pairs|
      tax, translated = pairs.select{|p| taxa.include?(p[0]) }.first
      next if translated.nil?
      pairs.each do |tax,id|
        next if id == translated
        index[id]= translated
      end
    end
    index
  end

  dep :to_taxa_hcop
  dep :to_taxa_homologene
  task :to_taxa => :tsv do 
    homologene = step(:to_taxa_homologene).load
    hcop = step(:to_taxa_hcop).load
    homologene.merge!(hcop)
    homologene
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
      text, tf, tg, tf_ids, tg_ids, tf_off, tg_off, mutation, mutation_off, score, valid, mor_scores, mor = values

      next unless valid == 'Valid'

      tf_name = tf_ids.split(";").collect do |tf_id|
        tf_human = to_human[tf_id] || tf_id
        entrez2name[tf_human]
      end.compact.uniq *  ";"

      tg_name = tg_ids.split(";").collect do |tg_id|
        tg_human = to_human[tg_id] || tg_id
        entrez2name[tg_human]
      end.compact.uniq *  ";"

      if tf_name.empty? && ! tf_ids.empty?
        Log.debug "Missing #{tf_ids}"
      end

      if tg_name.empty? && ! tg_ids.empty?
        Log.debug "Missing #{tg_ids}"
      end

      next if tf_name.empty? or tg_name.empty?

      pmid, sentence = k.split(":").values_at 1, 4

      k = [pmid, sentence, tf_name, tg_name] * ":"
      [k, [text, tf, tg, tf_name, tg_name, tf_off, tg_off, mutation, mutation_off, score, valid, mor_scores, mor]]
    end
  end

  task_alias :ExTRI2_human, ExTRI2, :ExTRI2_taxa, organism: "Hsa/feb2014"

  dep_task :ExTRI2_final, ExTRI2, :ExTRI2_human, original: Rbbt.share.ExTRI2.sentences.find

  dep :ExTRI2_human, original: Rbbt.share.ExTRI2.sentences.find
  input :remove_auto_regulation, :boolean, "Remove auto-regulation", false
  input :only_authoritative_tfs, :boolean, "Use only authoritative TFs", false
  input :no_MoR, :boolean, "Don't use MoR", false
  task :ExTRI2_clean => :tsv do |remove_auto_regulation,only_authoritative_tfs,no_MoR|
    authoritative_tfs = Rbbt.data.authoritative_tfs.tsv.keys
    dumper = TSV::Dumper.new(:key_field => "TRI", :fields => ["Transcription Factor (Associated Gene Name)", "Target Gene (Associated Gene Name)", "Interaction score", (no_MoR ? "DEACTIVATED" : "Sign"), "Sentence"], :type => :list, :namespace => "Hsa/feb2014")
    dumper.init
    traverse step(:ExTRI2_human), into: dumper do |k,values|
      text, tf_text, tg_text, tfs, tgs, tf_off, tg_off, mutation, mutation_off, score, valid, mor_score, mor = values

      next unless valid == "Valid"

      res = []
      tfs.split(";").uniq.each do |tf|
        next if only_authoritative_tfs && ! authoritative_tfs.include?(tf)
        tgs.split(";").uniq.each do |tg|
          next if remove_auto_regulation && tf == tg

          sentence = text.sub('[TF]', tf_text).sub('[TG]', tg_text)

          sign = case mor
                 when "ACTIVATION"
                   "UP"
                 when "REPRESION"
                   "DOWN"
                 else
                   "UNKNOWN"
                 end

          sign = "" if no_MoR

          pmid,num,_,_ = k.split(":")
          key = [pmid, num, tf, tg, tf_off, tg_off] * ":"
          res << [key, [tf, tg, score, sign, sentence]]
        end
      end
      next if res.empty?
      res.extend MultipleResult
      res
    end
  end
end
