require 'rbbt/sources/organism'
module ExTRI2

  def self.organism
    'Hsa/feb2014'
  end

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

  input :original, :file, "ExTRI2 main result",Rbbt.share.ExTRI2['ExTRI2_final_resource.tsv.gz'].find 
  input :organism, :string, "Organism code", ExTRI2.organism
  input :fix_normalization, :boolean, 'Try to fix normalizations in original file', false
  dep :to_taxa, name: :placeholder, taxa: :placeholder do |jobname,options|
    organism = options[:organism]
    taxa = Organism.entrez_taxids(organism).list * ","
    name = organism.split("/").first
    {inputs: options.merge(name: name, taxa: taxa)}
  end
  task :ExTRI2_taxa => :tsv do |original,organism,fix_normalization|

    to_human = step(:to_taxa).load
    entrez2name = Organism.identifiers(organism).index :target => "Associated Gene Name", :fields =>["Entrez Gene ID"], persist: true
    
    parser = TSV::Parser.new original, type: :list
    dumper = TSV::Dumper.new parser.options
    dumper.init

    valid_pos = parser.fields.index 'Valid'
    tf_id_pos = parser.fields.index 'TF Id'
    tg_id_pos = parser.fields.index('TG Id') || parser.fields.index('IG Id')
    tf_name_orig_pos = parser.fields.index 'TF_human_symbol'
    tg_name_orig_pos = parser.fields.index 'TG_human_symbol'
    tf_human_id_pos = parser.fields.index 'TF_human_Id'
    tg_human_id_pos = parser.fields.index 'TG_human_Id'

    TSV.traverse parser, into: dumper, bar: self.progress_bar("Translating gene orthologs") do |k,values|
      pmid, sentence = k.split(":").values_at 1, 4

      valid = values[valid_pos]
      tf_ids = values[tf_id_pos]
      tg_ids = values[tg_id_pos]

      tf_name_orig = values[tf_name_orig_pos] if tf_name_orig_pos
      tg_name_orig = values[tg_name_orig_pos] if tg_name_orig_pos

      tf_human_id = values[tf_name_orig_pos] if tf_human_id_pos
      tg_human_id = values[tg_name_orig_pos] if tg_human_id_pos

      tf_name_orig = 'None' if tf_name_orig.nil?
      tg_name_orig = 'None' if tg_name_orig.nil?

      if not fix_normalization
        tf_name_orig = 'None' if tf_human_id && tf_human_id.include?(';')
        tg_name_orig = 'None' if tf_human_id && tg_human_id.include?(';')
      end

      tf_name = tf_name_orig
      tg_name = tg_name_orig

      next unless valid == 'Valid'

      #if tf_name_orig == 'None'
      #  next unless fix_normalization
      #  tf_name = tf_name_orig
      #  #tf_name = tf_name_orig.split(";").collect do |tf_id|
      #  #  if tf_id == 'None'
      #  #    tf_human = to_human[tf_id] || tf_id
      #  #    entrez2name[tf_human]
      #  #  else
      #  #    tf_id
      #  #  end
      #  #end *  ";" 
      #  
      #  if tf_name.empty?
      #    Log.debug "Missing #{tf_ids}"
      #    next
      #  end
      #else
      #  tf_name = tf_name_orig
      #end

      #if tg_name_orig == 'None'
      #  next unless fix_normalization
      #  tg_name = tg_name_orig
      #  #tg_name = tg_name_orig.split(";").collect do |tg_id|
      #  #  if tg_id == 'None'
      #  #    tg_human = to_human[tg_id] || tg_id
      #  #    entrez2name[tg_human] || tg_id
      #  #  else
      #  #    tg_id
      #  #  end
      #  #end *  ";"

      #  if tg_name.empty?
      #    Log.debug "Missing #{tg_ids}"
      #    next
      #  end
      #else
      #  tg_name = tg_name_orig
      #end

      k = [pmid, sentence, tf_name, tg_name] * ":"
      #[k, [text, tf, tg, tf_name, tg_name, tf_off, tg_off, mutation, mutation_off, score, valid, mor_scores, mor]]
      new_values = values.dup

      new_values[tf_id_pos] = tf_name
      new_values[tg_id_pos] = tg_name

      [k, new_values]
    end
  end

  task_alias :ExTRI2_human, ExTRI2, :ExTRI2_taxa, organism: ExTRI2.organism

  dep_task :ExTRI2_final, ExTRI2, :ExTRI2_human, original: Rbbt.share.ExTRI2['ExTRI2_final_resource.tsv.gz'].find

  dep :ExTRI2_human
  input :remove_auto_regulation, :boolean, "Remove auto-regulation", true
  input :only_authoritative_tfs, :boolean, "Use only authoritative TFs", false
  input :no_MoR, :boolean, "Don't use MoR", false
  task :ExTRI2_clean => :tsv do |remove_auto_regulation,only_authoritative_tfs,no_MoR|
    authoritative_tfs = Rbbt.data.authoritative_tfs.tsv.keys
    tf_type = Rbbt.data.TF_type.tsv

    dumper = TSV::Dumper.new(:key_field => "TRI", :fields => ["Transcription Factor (Associated Gene Name)", "Target Gene (Associated Gene Name)", "Transcription Factor Type", "Interaction Score", (no_MoR ? "DEACTIVATED" : "Sign"), "Sentence"], :type => :list, :namespace => ExTRI2.organism)
    dumper.init

    parser = TSV::Parser.new step(:ExTRI2_human)

    text_pos = parser.fields.index('Sentence') || parser.fields.index('Text')
    tf_text_pos = parser.fields.index 'TF'
    tg_text_pos = parser.fields.index('TG') || parser.fields.index("Gene")
    tf_pos = parser.fields.index('unique_TF_human_symbol') || parser.fields.index('TF_human_symbol')
    tg_pos = parser.fields.index('unique_TG_human_symbol') || parser.fields.index('TG_human_symbol')
    tf_off_pos = parser.fields.index 'TF offset'
    tg_off_pos = parser.fields.index 'Gene offset'
    score_pos = parser.fields.index('TRI score') || parser.fields.index('Valid score')
    mor_pos = parser.fields.index 'MoR'
    type_pos = parser.fields.index 'human_TF_type' || parser.fields.index("TF_yype")

    traverse parser, into: dumper do |k,values|
      #text, tf_text, tg_text, tfs, tgs, tf_off, tg_off, mutation, mutation_off, score, valid, mor_score, mor = values

      text = values[text_pos]
      tfs = values[tf_pos]
      tgs = values[tg_pos]
      tf_off = values[tf_off_pos]
      tg_off = values[tg_off_pos]
      tf_text = values[tf_text_pos]
      tg_text = values[tg_text_pos]
      score = values[score_pos]
      mor = values[mor_pos]
      type = type_pos ? values[type_pos] : ''

      res = []
      tfs.split(";").uniq.each_with_index do |tf,inf|
        type = tf_type[tf] if type == ''
        next if tf == 'None'

        next unless %w(dbTF coTF).include? type
        next if only_authoritative_tfs && ! authoritative_tfs.include?(tf)
        tgs.split(";").uniq.each do |tg|
          next if tg == 'None'
          next if remove_auto_regulation && tf == tg

          #if type && type.include?(";")
          #  type_ = type.split(';')[inf]
          #else
          #  type_ = type
          #end

          sentence = text.sub('[TF]', tf_text).sub('[TG]', tg_text)

          sign = case mor
                 when "ACTIVATION"
                   "UP"
                 when 'REPRESSION'
                   "DOWN"
                 else
                   "UNKNOWN"
                 end

          sign = "" if no_MoR

          pmid,num,_,_ = k.split(":")
          key = [pmid, num, tf, tg, tf_off, tg_off] * ":"
          res << [key, [tf, tg, type, score, sign, sentence]]
        end
      end
      res.uniq!
      next if res.empty?
      res.extend MultipleResult
      res
    end
  end
end
