module ExTRI2
  task tuberculosis_genes: :array do
    genes = Rbbt.data['ncbi_dataset.tsv'].read.split("\n").collect{|line| line.split("\t").first }
    genes
  end
  
  dep :tuberculosis_genes
  task tuberculosis_sentences: :tsv do
    genes = step(:tuberculosis_genes).load
    tsv = Rbbt.share.ExTRI2.sentences.tsv
    tsv.select("TF Id" => genes)
  end

  task texts: :array do
    Workflow.require_workflow "PubTator"
    entity_repo = PubTator.entity_repo
    entity_repo.unnamed = false
    corpus = PubTator.corpus
    current = nil
    entities = []
    sentences = []
    prompts = []
    Open.open_pipe do |sin|
      TSV.traverse Rbbt.share.ExTRI2.sentences, bar: true do |id,values|
        sentence, tf, tg, tfid, tgid, tfoff, tgoff, mut, mutoff, score, valid, mscore, mor = values
        next unless valid == "Valid"
        docid = id.split(':').values_at(0,1,2,3) * ':'
        ids = entity_repo.prefix docid
        next if ids.empty?

        if current and current != docid
          info = {}
          info[:text] = text = corpus[docid]
          info[:entities] = entities.collect{|e| 
            literal, entity_type = entity_repo[e].values_at 'literal', 'entity_type'
            next unless entity_type == 'Gene'
            literal
          }.compact.uniq
          info[:sentences] = sentences

          sin.puts info.to_json unless info[:entities].empty?
          prompts << info
          entities = []
          sentences = []
        end
        current = docid
        entities.concat ids
        sentences << [sentence, tf, tg, mor]
      end
      prompts
    end
  end

  dep :texts
  task choices: :array do
    TSV.traverse step(:texts), into: :stream do |json|
      info = JSON.parse json
      IndiferentHash.setup info

      genes = info['entities']
      sentences = info['sentences']
      good = sentences.collect do |text, tf, tg, mor|
        [tf, tg, mor, text]
      end
       
      decoys = good.collect{|tf, tg, mor| 
        case mor
        when 'ACTIVATION'
          [tf, tg, 'REPRESION']
        when 'REPRESION'
          [tf, tg, 'ACTIVATION']
        else
          [tg, tf, mor]
        end
      }.compact

      genes.each do |tf|
        genes.each do |tg|
          next if tf == tg
          mor = good.shuffle.first[2]
          next if good.select{|gtf,gtg,mor| tf == gtf && tg == gtg  }.any?
          decoys << [tf, tg, mor]
          break if decoys.length > 7 
        end
      end
      decoys.uniq!
      
      {
        text: info[:text],
        good: good,
        decoys: decoys
      }.to_json
    end
  end

  helper :format_choice do |tf,tg,mor,sentence=nil|
    case mor
    when 'ACTIVATION'
      "The transcription factor #{tf} promotes the expression of #{tg}"
    when 'REPRESION'
      "The transcription factor #{tf} inhibits the expression of #{tg}"
    else
      "The transcription factor #{tf} regulates the expression of #{tg}"
    end
  end

  dep :choices
  task prompts: :array do
    TSV.traverse step(:choices), into: :stream do |json|
      info = JSON.parse json
      IndiferentHash.setup info
      answer = info[:good].shuffle.first
      bad = info[:decoys].shuffle[0..5]

      choices = bad  + [answer]
      choices = choices.shuffle

      question = <<-EOF
Consider the following article abstract:

#{info[:text]}

Which of the following statements is implied in the above text?

      EOF

      choices.each_with_index do |p,i|
        question += " -#{i+1}: " + format_choice(*p) + "\n"
      end

      number = choices.index answer
      answer_text = "The correct answer is #{number+1}: #{format_choice(*answer)}"

      {
        question: question,
        answer: answer_text,
        text: answer.last
      }.to_json
    end
  end

end
